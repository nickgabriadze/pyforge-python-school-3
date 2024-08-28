from os import getenv
from .dao import MoleculesDAO
from rdkit import Chem
from fastapi import FastAPI, UploadFile
from fastapi import status
from .logger import setupLogging
from .generators.sub_search import substructure_search


logger = setupLogging()
app = FastAPI()


@app.get("/")
def get_server():
    logger.info(f'[METHOD] /GET - [PATH] /')
    return {"server_id": getenv("SERVER_ID", "1")}


@app.get('/api/v1/molecules', description="Retrieve all the available molecules")
async def get_molecules(limit=100):
    logger.info(f'[METHOD] /GET - [PATH] /api/v1/molecules')

    return await MoleculesDAO.get_all_molecules(limit)


@app.post('/api/v1/molecules', description="Add a new molecule")
async def add_molecule(mol_smiles: str):
    logger.info(f'[METHOD] /POST - [PATH] /api/v1/molecules/{mol_smiles}')

    mol_smiles = mol_smiles.strip()
    elements = await MoleculesDAO.get_all_molecules()

    new_molecule = Chem.MolFromSmiles(mol_smiles)
    # this ensures that whatever I pass as a parameter is a chemical and not some random string
    if new_molecule:
        identifier = f"PUBCHEM{1 if len(elements) == 0 else int(elements[-1].pubchem_id.split("PUBCHEM")[1]) + 1}"
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

    mol_id = mol_id.strip()
    molecule = await MoleculesDAO.get_molecule_by_pubchem_id(pubchem_id=mol_id)
    if molecule is None:
        logger.error(f'Failed to get the molecule - {mol_id} not found')
        return f'{status.HTTP_404_NOT_FOUND} - NOT FOUND'
    else:
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
        logger.info(f'[ACTION] DELETE - {mol_id}')
        return f'{status.HTTP_204_NO_CONTENT} - DELETED'


@app.get('/api/v1/sub_match/{mol_smiles}',
         description="Match the substructure of given smiles molecule with other saved ones")
async def get_sub_match(mol_smiles: str, limit=100):
    logger.info(f'[METHOD] /GET - [PATH] /api/v1/sub_match/{mol_smiles}')

    sub_matches = []
    mol_smiles = mol_smiles.strip()
    if Chem.MolFromSmiles(mol_smiles):
        all_molecules = [el.getSmiles() for el in (await MoleculesDAO.get_all_molecules(limit))]

        for molecule in substructure_search(all_molecules, mol_smiles):
            sub_matches.append(molecule)

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
            logger.error(f'Failed to update the molecule - {mol_id} was not found or already exists')
            return f"{status.HTTP_404_NOT_FOUND} - Either corresponding ID was not found or smiles already exists"

        logger.info(f"[ACTION] UPDATE - {mol_id} with new value of {new_mol_smiles}")

        return status.HTTP_200_OK

    logger.error(f'Failed to update the molecule - {new_mol_smiles} is not a molecule')
    return f'{status.HTTP_400_BAD_REQUEST} BAD REQUEST - not a molecule'


@app.post('/api/v1/upload-molecules')  # format is txt, each smiles is on the new line
async def upload_molecules(molecules: UploadFile):
    logger.info(f'[METHOD] /POST - [PATH] /api/v1/upload-molecules')
    contents = str(await molecules.read()).split('\\r\\n')
    added = 0
    elements = await MoleculesDAO.get_all_molecules()
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
