from rdkit import Chem, RDLogger
from rdchiral.chiral import copy_chirality
from rdkit.Chem.AllChem import AssignStereochemistry


RDLogger.DisableLog('rdApp.*')


def canonicalize2(smiles, isomeric=False, canonical=True, kekulize=False):
    # When canonicalizing a SMILES string, we typically want to
    # run Chem.RemoveHs(mol), but this will try to kekulize the mol
    # which is not required for canonical SMILES. Instead, we make a
    # copy of the mol retaining only the information we desire (not explicit Hs).
    # Then, we sanitize the mol without kekulization. copy_atom and copy_edit_mol
    # are used to create this clean copy of the mol.
    
    def copy_atom(atom):
        new_atom = Chem.Atom(atom.GetSymbol())
        new_atom.SetFormalCharge(atom.GetFormalCharge())
        
        # Preserve atom mapping number
        atom_map_num = atom.GetAtomMapNum()
        if atom_map_num:
            new_atom.SetAtomMapNum(atom_map_num)
        
        if atom.GetIsAromatic() and atom.GetNoImplicit():
            new_atom.SetNumExplicitHs(atom.GetNumExplicitHs())

        return new_atom

    def copy_edit_mol(mol):
        new_mol = Chem.RWMol(Chem.MolFromSmiles(''))
        atom_mapping = {}

        for atom in mol.GetAtoms():
            new_atom = copy_atom(atom)
            new_idx = new_mol.AddAtom(new_atom)
            atom_mapping[atom.GetIdx()] = new_idx

        for bond in mol.GetBonds():
            a1 = atom_mapping[bond.GetBeginAtom().GetIdx()]
            a2 = atom_mapping[bond.GetEndAtom().GetIdx()]
            bt = bond.GetBondType()
            new_mol.AddBond(a1, a2, bt)
            new_bond = new_mol.GetBondBetweenAtoms(a1, a2)
            new_bond.SetBondDir(bond.GetBondDir())
            new_bond.SetStereo(bond.GetStereo())
        
        for new_atom in new_mol.GetAtoms():
            atom = mol.GetAtomWithIdx(new_atom.GetIdx())
            copy_chirality(atom, new_atom)

        for atom in mol.GetAtoms():
            new_atom = new_mol.GetAtomWithIdx(atom_mapping[atom.GetIdx()])
            copy_chirality(atom, new_atom)

        return new_mol

    smiles = smiles.replace(" ", "")  
    tmp = Chem.MolFromSmiles(smiles, sanitize=False)
    tmp.UpdatePropertyCache()
    new_mol = copy_edit_mol(tmp)

    if not kekulize:
        Chem.SanitizeMol(new_mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_SETAROMATICITY | 
                         Chem.SanitizeFlags.SANITIZE_PROPERTIES | 
                         Chem.SanitizeFlags.SANITIZE_ADJUSTHS, catchErrors=True)
    else:
        Chem.SanitizeMol(new_mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_KEKULIZE | 
                         Chem.SanitizeFlags.SANITIZE_PROPERTIES | 
                         Chem.SanitizeFlags.SANITIZE_ADJUSTHS, catchErrors=True)
   
    AssignStereochemistry(new_mol, cleanIt=False, force=True, flagPossibleStereoCenters=True)
    
    new_smiles = Chem.MolToSmiles(new_mol, isomericSmiles=isomeric, canonical=canonical)
    return new_smiles


def canonicalize_molecule_smiles(smiles, return_none_for_error=True, skip_mol=False, sort_things=True, isomeric=True, kekulization=True, allow_empty_part=False, keep_atom_map=True):
    things = smiles.split('.')
    if skip_mol:
        new_things = things
    else:
        new_things = []
        for thing in things:
            try:
                if thing == '' and not allow_empty_part:
                    raise ValueError('SMILES contains empty part.')

                mol = Chem.MolFromSmiles(thing)
                assert mol is not None
                if not keep_atom_map:
                    for atom in mol.GetAtoms():
                        atom.SetAtomMapNum(0)
                thing_smiles = Chem.MolToSmiles(mol, kekuleSmiles=False, isomericSmiles=isomeric)
                thing_smiles = Chem.MolFromSmiles(thing_smiles)
                thing_smiles = Chem.MolToSmiles(thing_smiles, kekuleSmiles=False, isomericSmiles=isomeric)
                thing_smiles = Chem.MolFromSmiles(thing_smiles)
                thing_smiles = Chem.MolToSmiles(thing_smiles, kekuleSmiles=False, isomericSmiles=isomeric)
                assert thing_smiles is not None
                can_in = thing_smiles
                can_out = canonicalize2(thing_smiles, isomeric=isomeric)
                assert can_out is not None, can_in
                thing_smiles = can_out
                if kekulization:
                    thing_smiles = keku_mid = Chem.MolFromSmiles(thing_smiles)
                    assert keku_mid is not None, 'Before can: %s\nAfter can: %s' % (can_in, can_out)
                    thing_smiles = Chem.MolToSmiles(thing_smiles, kekuleSmiles=True, isomericSmiles=isomeric)
            except KeyboardInterrupt:
                raise
            except:
                if return_none_for_error:
                    return None
                else:
                    raise
            new_things.append(thing_smiles)
    if sort_things:
        new_things = sorted(new_things)
    new_things = '.'.join(new_things)
    return new_things
