from rdkit import Chem


def substructure_search(mols, mol):
    substructure_matches = []
    initial_mol = Chem.MolFromSmiles(mol)

    for molecule in mols:
        list_mol = Chem.MolFromSmiles(molecule, sanitize=False)

        if initial_mol.HasSubstructMatch(list_mol) \
                or list_mol.HasSubstructMatch(initial_mol):
            substructure_matches.append(molecule)

    return substructure_matches


print(substructure_search(["CCO", "CCC", "CCN", "CNC"], "C"))
