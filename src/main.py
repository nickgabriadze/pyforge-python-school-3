from os import getenv

from rdkit import Chem
from fastapi import FastAPI, UploadFile
from fastapi import status


def substructure_search(mols, mol):
    substructure_matches = []
    initial_mol = Chem.MolFromSmiles(mol)

    for molecule in mols:
        list_mol = Chem.MolFromSmiles(molecule, sanitize=False)

        if initial_mol.HasSubstructMatch(list_mol) \
                or list_mol.HasSubstructMatch(initial_mol):
            substructure_matches.append(molecule)

    return substructure_matches


molecules_db = {
    "PUBCHEM1": "NCC",
    "PUBCHEM2": "c1ccccc1",
    "PUBCHEM3": "CCO",
    "PUBCHEM4": "CC(=O)O",
    "PUBCHEM5": "CC(=O)Oc1ccccc1C(=O)O"
}

ID_COUNTER = len(molecules_db)

app = FastAPI()


@app.get("/")
def get_server():
    return {"server_id": getenv("SERVER_ID", "1")}


@app.get('/api/v1/molecules', description="Retrieve all the available molecules")
def get_molecules():
    return molecules_db


@app.post('/api/v1/molecules', description="Add a new molecule")
def add_molecule(mol_smiles: str):
    global ID_COUNTER
    mol_smiles = mol_smiles.strip()
    # I know this is not the most efficient thing, but it will do the job for now
    if mol_smiles in list(molecules_db.values()):
        return f'{status.HTTP_400_BAD_REQUEST} BAD REQUEST - already exists'

    new_molecule = Chem.MolFromSmiles(mol_smiles)
    # this ensures that whatever I pass as a parameter is a chemical and not some random string
    if new_molecule:
        identifier = f'PUBCHEM{ID_COUNTER + 1}'
        molecules_db[identifier] = mol_smiles
        ID_COUNTER += 1
        return status.HTTP_201_CREATED
    return f'{status.HTTP_400_BAD_REQUEST} BAD REQUEST - not a molecule'


@app.get('/api/v1/molecules/{mol_id}', description="Get a molecule by its 'PUBCHEM' id")
def get_molecule(mol_id: str):
    mol_id = mol_id.strip()

    if not molecules_db.get(mol_id):
        return f'{status.HTTP_404_NOT_FOUND} - NOT FOUND'

    return molecules_db[mol_id]


@app.delete('/api/v1/molecules/{mol_id}', description="Delete a molecule by its 'PUBCHEM' id")
def delete_molecule(mol_id: str):
    mol_id = mol_id.strip()
    if molecules_db.get(mol_id) is not None:
        del molecules_db[mol_id]
        return status.HTTP_200_OK
    else:
        return f'{status.HTTP_404_NOT_FOUND} - NOT FOUND'


@app.get('/api/v1/sub_match/{mol_smiles}',
         description="Match the substructure of given smiles molecule with other saved ones")
def get_sub_match(mol_smiles: str):
    mol_smiles = mol_smiles.strip()
    if Chem.MolFromSmiles(mol_smiles):
        all_molecules = list(molecules_db.values())
        return substructure_search(all_molecules, mol_smiles)

    return f'{status.HTTP_400_BAD_REQUEST} BAD REQUEST - not a molecule'


@app.put('/api/v1/molecules', description="Update a molecule by its 'PUBCHEM' identifier")
def update_molecule(mol_id: str, new_mol_smiles: str):
    mol_id = mol_id.strip()
    new_mol_smiles = new_mol_smiles.strip()
    if new_mol_smiles in list(molecules_db.values()):
        return f'{status.HTTP_400_BAD_REQUEST} BAD REQUEST - already exists'

    if molecules_db.get(mol_id) is not None:
        new_molecule = Chem.MolFromSmiles(new_mol_smiles)
        if new_molecule:
            molecules_db[mol_id] = new_mol_smiles
            return status.HTTP_200_OK
        return f'{status.HTTP_409_CONFLICT} - NOT a molecule!'

    return f'{status.HTTP_404_NOT_FOUND} - NOT FOUND'


@app.post('/api/v1/upload-molecules')  # format is txt, each smiles is on the new line
async def upload_molecules(molecules: UploadFile):
    global ID_COUNTER
    contents = str(await molecules.read()).split('\\r\\n')
    added = 0

    for mol_smiles in contents:
        mol_smiles = mol_smiles.strip("'")
        if mol_smiles.startswith("b"):
            mol_smiles = mol_smiles[2:]
        new_molecule = Chem.MolFromSmiles(mol_smiles, sanitize=False)
        if new_molecule and len(mol_smiles) > 0:
            identifier = f'PUBCHEM{ID_COUNTER + 1}'
            molecules_db[identifier] = mol_smiles
            added += 1
            ID_COUNTER += 1
        else:
            continue

    return f'{status.HTTP_201_CREATED} OK - ADDED {added} molecule{'s' if added > 1 else ''}'
