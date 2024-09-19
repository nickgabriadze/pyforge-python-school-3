from rdkit import Chem


def substructure_search(mols, mol):
    initial_mol = Chem.MolFromSmiles(mol)

    for molecule in mols:
        list_mol = Chem.MolFromSmiles(molecule, sanitize=False)

        if initial_mol.HasSubstructMatch(list_mol) \
                or list_mol.HasSubstructMatch(initial_mol):
            yield molecule
