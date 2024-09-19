from os import getenv
import redis
from .db.dao import MoleculesDAO
from rdkit import Chem
from fastapi import FastAPI, UploadFile
from fastapi import status
from .logger import setupLogging
from .generators.sub_search import substructure_search
from .caching.cache_handler import set_cache, get_cached_result, remove_cache
from src.celery_worker import celery
from celery.result import AsyncResult
from src.tasks import sub_search_task

logger = setupLogging()
app = FastAPI()
redis_client = redis.Redis(host='redis', port=6379, db=0, password=getenv('REDIS_PASSWORD'))


@app.post("/tasks/substructure_match")
async def create_task(mol: str, limit: int = 100):
    task = sub_search_task.delay(mol, limit)

    return {"task_id": task.id, "status": task.status}


@app.get("/tasks/{task_id}")
async def get_task_result(task_id: str):
    task_result = AsyncResult(task_id, app=celery)
    if task_result.state == 'PENDING':
        return {"task_id": task_id, "status": "Task is still processing"}
    elif task_result.state == 'SUCCESS':
        return {"task_id": task_id, "status": "Task completed", "result": task_result.result}
    else:
        return {"task_id": task_id, "status": task_result.state}


@app.get("/")
def get_server():
    logger.info('[METHOD] /GET - [PATH] /')
    return {"server_id": getenv("SERVER_ID", "1")}


def read_root():
    return "Hello!"


@app.get('/api/v1/molecules', description="Retrieve all the available molecules")
async def get_molecules(limit=100):
    logger.info('[METHOD] /GET - [PATH] /api/v1/molecules')

    return await MoleculesDAO.get_molecules(limit)


@app.post('/api/v1/molecules/{mol_smiles}', description="Add a new molecule")
async def add_molecule(mol_smiles: str):
    logger.info(f'[METHOD] /POST - [PATH] /api/v1/molecules/{mol_smiles}')

    mol_smiles = mol_smiles.strip()
    elements = await MoleculesDAO.get_molecules()

    new_molecule = Chem.MolFromSmiles(mol_smiles)
    # this ensures that whatever I pass as a parameter is a chemical and not some random string
    if new_molecule:
        identifier = f"PUBCHEM{1 if len(elements) == 0 else int(elements[-1].pubchem_id.split('PUBCHEM')[1]) + 1}"
        result = await MoleculesDAO.add_smiles(pubchem_id=identifier, smiles=mol_smiles)
        if result is None:
            logger.error(f"Failed to add molecule - {mol_smiles} already exists")
            return f"{status.HTTP_400_BAD_REQUEST} BAD REQUEST - already exists"

        logger.info(f"[ACTION] ADD - {mol_smiles}")
        return status.HTTP_201_CREATED

    logger.error(f"Failed to add molecule - {mol_smiles} is not a molecule")
    return f'{status.HTTP_400_BAD_REQUEST} BAD REQUEST - not a molecule'


@app.get('/api/v1/molecules/{mol_id}', description="Get a molecule by its 'PUBCHEM' id")
async def get_molecule(mol_id: str):
    logger.info(f'[METHOD] /GET - [PATH] /api/v1/molecules/{mol_id}')

    cache_key = f"search:{mol_id}"
    cached_result = get_cached_result(redis_client, cache_key)
    if cached_result:
        logger.info(f'[REDIS] - Returned cached result for {mol_id}')
        return cached_result

    mol_id = mol_id.strip()
    molecule = await MoleculesDAO.get_molecule_by_pubchem_id(pubchem_id=mol_id)
    if molecule is None:
        logger.error(f'Failed to get the molecule - {mol_id} not found')
        return f'{status.HTTP_404_NOT_FOUND} - NOT FOUND'
    else:
        logger.info(f'[ACTION] CACHE - cached {cache_key}')
        set_cache(redis_client, cache_key, molecule.smiles)
        return molecule.smiles


@app.delete('/api/v1/molecules/{mol_id}', description="Delete a molecule by its 'PUBCHEM' id")
async def delete_molecule(mol_id: str):
    logger.info(f'[METHOD] /DELETE - [PATH] /api/v1/molecules/{mol_id}')

    mol_id = mol_id.strip()
    result = await MoleculesDAO.delete_molecule_by_pubchem_id(pubchem_id=mol_id)
    if result is None:
        logger.error(f'Failed to delete the molecule with ID {mol_id} - not found')
        return f'{status.HTTP_404_NOT_FOUND} - NOT FOUND'
    else:
        # if molecule with given pubchem id is deleted, the cache will be automatically removed
        # if it exists
        remove_cache(redis_client, f'search:{mol_id}')
        logger.info(f'[ACTION] DELETE - {mol_id}')
        return f'{status.HTTP_204_NO_CONTENT} - DELETED'


@app.get('/api/v1/sub_match/{mol_smiles}',
         description="Match the substructure of given smiles molecule with other saved ones")
async def get_sub_match(mol_smiles: str, limit: int = 100):
    logger.info(f'[METHOD] /GET - [PATH] /api/v1/sub_match/{mol_smiles}')

    cache_key = f'sub_match:{mol_smiles}?limit={limit}'
    cached = get_cached_result(redis_client, cache_key)
    if cached:
        logger.info(f'[REDIS] - Returned cached result for {cache_key}')
        return dict(cached).get('matches')

    sub_matches = []
    mol_smiles = mol_smiles.strip()
    if Chem.MolFromSmiles(mol_smiles):
        all_molecules = [el.getSmiles() for el in (await MoleculesDAO.get_all_molecules())]
        generated_matches = 0
        for molecule in substructure_search(all_molecules, mol_smiles):
            if generated_matches >= limit:
                break
            else:
                sub_matches.append(molecule)
                generated_matches += 1

        logger.info(f'[ACTION] CACHE - cached {cache_key}')
        # putting data to cache for 10 minutes, because sub-substructure matching requires significant amount of time for computation
        # and data doesn't change that frequently (till new molecules are added)
        set_cache(redis_client, cache_key, {'matches': sub_matches}, 600)
        return sub_matches

    logger.error(f'Failed to get submatch - {mol_smiles} is not a molecule')
    return f'{status.HTTP_400_BAD_REQUEST} BAD REQUEST - not a molecule'


@app.put('/api/v1/molecules', description="Update a molecule by its 'PUBCHEM' identifier")
async def update_molecule(mol_id: str, new_mol_smiles: str):
    logger.info(f'[METHOD] /PUT - [PATH] /api/v1/molecules/{mol_id}')

    mol_id = mol_id.strip()
    new_mol_smiles = new_mol_smiles.strip()
    new_molecule = Chem.MolFromSmiles(new_mol_smiles)
    if new_molecule:
        request = await MoleculesDAO.update_molecule(pubchem_id=mol_id, new_mol_smiles=new_mol_smiles)
        if request is None:
            logger.error(
                f'Failed to update the molecule - ID {mol_id} was not found or {new_mol_smiles} already exists')
            return f"{status.HTTP_404_NOT_FOUND} - Either corresponding ID was not found or smiles already exists"

        logger.info(f"[ACTION] UPDATE - {mol_id} with new value of {new_mol_smiles}")

        return status.HTTP_200_OK

    logger.error(f'Failed to update the molecule - {new_mol_smiles} is not a molecule')
    return f'{status.HTTP_400_BAD_REQUEST} BAD REQUEST - not a molecule'


@app.post('/api/v1/upload-molecules')  # format is txt, each smiles is on the new line
async def upload_molecules(molecules: UploadFile):
    logger.info('[METHOD] /POST - [PATH] /api/v1/upload-molecules')
    contents = str(await molecules.read()).split('\\r\\n')
    added = 0
    elements = await MoleculesDAO.get_molecules()
    identifier = 1 if len(elements) == 0 else int(elements[-1].pubchem_id.split("PUBCHEM")[1]) + 1

    for mol_smiles in contents:
        mol_smiles = mol_smiles.strip("'")
        if mol_smiles.startswith("b"):
            mol_smiles = mol_smiles[2:]
        new_molecule = Chem.MolFromSmiles(mol_smiles, sanitize=False)
        if new_molecule and len(mol_smiles) > 0:
            result = await MoleculesDAO.add_smiles(pubchem_id=f"PUBCHEM{identifier}", smiles=mol_smiles)
            if result is not None:
                added += 1
                identifier = identifier + 1
        else:
            continue
    logger.info(f"[ACTION] ADD - {added} molecule{'s' if added > 1 else ''}")
    return f'{status.HTTP_201_CREATED} OK - ADDED {added} molecule{"s" if added > 1 else ""}'
