from src.main import substructure_search
from src.main import app
import pytest
from fastapi.testclient import TestClient
from fastapi import status


@pytest.mark.parametrize('mols,mol,res',
                         [
                             (["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "c1ccccc1",
                              ["c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"]),

                             (["NCC", "c1ccccc1", "CCO", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "CO",
                              ["CCO", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]),

                             (["NCC", "c1ccccc1", "CCO", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "O",
                              ["CCO", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]),

                             (["c1ccccc1", "CCO", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "CC",
                              ["CCO", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]),

                             ([], "NCO", []),

                             (["CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O", "CC[NH2+][C@@H](C)Cc1cccc(c1)C(F)(F)F",
                               "C[C@@H]1C[C@@H](CC([NH+]1C)(C)C)OC(=O)[C@@H](c2ccccc2)O"], "F",
                              ["CC[NH2+][C@@H](C)Cc1cccc(c1)C(F)(F)F"]),

                             (["CCO", "CCC", "CCN", "CNC"], "CC", ["CCO", "CCC", "CCN"]),

                             (["COC", "NCCO", "OCC"], "N", ["NCCO"]),
                         ])
def test_substructure_search(mols, mol, res):
    assert substructure_search(mols, mol) == res


@pytest.mark.xfail(reason="H2O is not in mols and H2O is not substructure of CO, should have been []")
def test_substructure_search_fail_not_in_mols():
    assert substructure_search(['CO'], "H2O") == ["CO"]


@pytest.mark.xfail(reason="result is all of them and shouldn't be an empty list")
def test_substructure_search_empty():
    assert substructure_search(["CCO", "CCC", "CCN", "CNC"], "C") == []


def test_substructure_search_same():
    assert substructure_search(['CCO', 'c1ccccc1'], 'c1ccccc1') == ['c1ccccc1']


def test_substructure_search_exception_mols():
    """ won't be able to parse SMILES "21" as molecule"""
    with pytest.raises(TypeError):
        substructure_search(['21', 'c1ccccc1', 'NCC'], 'CO')


def test_substructure_search_exception_mol():
    """ c1cccc1 is not a molecule """
    with pytest.raises(AttributeError):
        substructure_search(['CCO', 'c1ccccc1'], 'c1cccc1')


def test_substructure_search_exception_mol_number():
    """ 5 is not a molecule """
    with pytest.raises(AttributeError):
        substructure_search(['CCO', 'c1ccccc1'], '5')


""" Other FASTAPI routes' tests """

t_client = TestClient(app=app)


def test_server_running():
    response = t_client.get('/')
    assert response.status_code == status.HTTP_200_OK


def test_get_molecules_success():
    response = t_client.get('/api/v1/molecules')
    assert response.status_code == status.HTTP_200_OK


def test_post_molecule_success():
    response = t_client.post('/api/v1/molecules', params={"mol_smiles": "CCO"})
    assert response.status_code == status.HTTP_200_OK


def test_post_already_exists_success():
    response = t_client.post('/api/v1/molecules', params={"mol_smiles": "CCO"})
    assert response.json() == "400 BAD REQUEST - already exists"


def test_delete_molecule_success():
    response = t_client.delete('/api/v1/molecules/{mol_id}', params={"mol_id": "PUBCHEM3"})
    assert response.status_code == status.HTTP_200_OK


def test_get_sub_match_success():
    response = t_client.get(f'/api/v1/sub_match/{"CO"}')
    assert response.json() == ['CCO', 'CC(=O)O', 'CC(=O)Oc1ccccc1C(=O)O']


def test_put_new_mol_success():
    response = t_client.put('/api/v1/molecules/', params={'mol_id': 'PUBCHEM1', 'new_mol_smiles': 'CO'})
    assert response.status_code == status.HTTP_200_OK
