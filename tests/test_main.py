import pytest
from fastapi.testclient import TestClient
from fastapi import status
from src.generators.sub_search import substructure_search


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
    sub_matches = []
    for mol in substructure_search(mols, mol):
        sub_matches.append(mol)

    assert sub_matches == res


@pytest.mark.xfail(reason="H2O is not in mols and H2O is not substructure of CO, should have been []")
def test_substructure_search_fail_not_in_mols():
    sub_matches = []
    for mol in substructure_search(['CO'], "H2O"):
        sub_matches.append(mol)

    assert sub_matches == ['CO']


@pytest.mark.xfail(reason="result is all of them and shouldn't be an empty list")
def test_substructure_search_empty():
    sub_matches = []
    for mol in substructure_search(["CCO", "CCC", "CCN", "CNC"], "C"):
        sub_matches.append(mol)

    assert sub_matches == []


def test_substructure_search_same():
    sub_matches = []
    for mol in substructure_search(['CCO', 'c1ccccc1'], "c1ccccc1"):
        sub_matches.append(mol)

    assert sub_matches == ['c1ccccc1']


def test_substructure_search_exception_mols():
    """ won't be able to parse SMILES "21" as molecule"""
    with pytest.raises(TypeError):
        for _ in substructure_search(['21', 'c1ccccc1', 'NCC'], 'CO'):
            continue


def test_substructure_search_exception_mol():
    """ c1cccc1 is not a molecule """
    with pytest.raises(AttributeError):
        for _ in substructure_search(['CCO', 'c1ccccc1'], 'c1cccc1'):
            continue


def test_substructure_search_exception_mol_number():
    """ 5 is not a molecule """
    with pytest.raises(AttributeError):
        for _ in substructure_search(['CCO', 'c1ccccc1'], '5'):
            continue
