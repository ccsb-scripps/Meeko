from rdkit import Chem
from rdkit.Chem import AllChem


def find_smarts(mol, smarts, smarts_indices):
    """ find occurrences of the SMARTS indices atoms in the requested SMARTS"""
    indices = []
    patt = Chem.MolFromSmarts(smarts)
    smarts_size = patt.GetNumAtoms()
    if smarts_indices[0] >= smarts_size or smarts_indices[1] >= smarts_size:
        raise ValueError("SMARTS index exceeds number of atoms in SMARTS (%d)" % (smarts_size))
    found = mol.GetSubstructMatches(patt)
    if len(found)>1:
        print(f"WARNING: the specified pattern ({smarts}) returned {len(found)} matches. " )
    for f in found:
        indices.append([f[x] for x in smarts_indices])
    return indices


def transform(ligand, index_pair, attractors_p3d):
    """ generate translatead and aligned molecules for each of the indices requested
        and for all the residues defined in the class constructor
        SOURCE: https://sourceforge.net/p/rdkit/mailman/message/36750909/
    """
    # make a copy of the ligand
    # TODO: maybe define new conformers?
    mol = Chem.Mol(ligand)
    target = Chem.MolFromSmiles("CC")
    # add hydrogens
    target = Chem.AddHs(target)
    # generate 3D coords
    AllChem.EmbedMolecule(target)
    # get the first conformer
    conf = target.GetConformer()
    # set coordinates to the ones of the actual target residue
    c1, c2 = attractors_p3d
    conf.SetAtomPosition(0, c1)
    conf.SetAtomPosition(1, c2)
    # perform alignment
    Chem.rdMolAlign.AlignMol(mol, target, -1, -1,[(index_pair[0],0), (index_pair[1], 1)] )
    return mol

def get_fragments_by_atom_indices(mol: Chem.Mol,
                                atom_idx1: int, atom_idx2: int) -> tuple[tuple[int], tuple[int]]:

    """
    Splits mol into two frags by bond between atoms w/ atom_idx1 and atom_idx2
    """

    bond = mol.GetBondBetweenAtoms(atom_idx1, atom_idx2)
    if bond is None:
        raise ValueError(f"No bond exists between the specified atom indices {atom_idx1, atom_idx2}.")

    emol = Chem.EditableMol(mol)
    emol.RemoveBond(atom_idx1, atom_idx2)
    mol = emol.GetMol()

    index_fragments = Chem.GetMolFrags(mol, asMols=False)
    if len(index_fragments) != 2:
        raise ValueError(f"Expected 2 fragments, but got {len(index_fragments)}.")

    frag1_indices, frag2_indices = index_fragments
    if atom_idx1 in frag1_indices:
        return (frag1_indices, frag2_indices)
    else:
        return (frag2_indices, frag1_indices)