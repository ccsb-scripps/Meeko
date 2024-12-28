from rdkit import Chem
from rdkit.Chem import rdChemReactions
from .utils import mini_periodic_table
from .pdbutils import PDBAtomInfo
from rdkit.Geometry import Point3D
from rdkit.Chem import rdDetermineBonds

periodic_table = Chem.GetPeriodicTable()


"""
create new RDKIT residue

mi  =  Chem.AtomPDBResidueInfo()
mi.SetResidueName('MOL')
mi.SetResidueNumber(1)
mi.SetOccupancy(0.0)
mi.SetTempFactor(0.0)

source: https://sourceforge.net/p/rdkit/mailman/message/36404394/
"""

def set_h_isotope_atom_coords(mol: Chem.Mol, conf: Chem.Conformer, return_mol: bool = False) -> dict[int, Point3D]: 

    """
    use AddHs() and RemoveHs() to generate H isotopes coordinates
    in a mol copy as if they were regular Hs
    returns a dict of H isotope atom idx and assigned position as Point3D
    """

    # check if molecule has H isotopes
    def is_h_isotope(atom: Chem.Atom) -> bool:
        return atom.GetAtomicNum() == 1 and atom.GetIsotope() > 0
    def has_h_isotopes(mol: Chem.Mol) -> bool:
        for atom in mol.GetAtoms():
            if is_h_isotope(atom):
                return True
        return False
    
    if not has_h_isotopes(mol): 
        return {}

    """
    create a nested dictionary for index mapping, isotope type and chirality flag
    
    isotope_data = {
        'H_isotope_idx': {  
            1: {'parent': 0, 'Isotope': 2}, 
            2: {'parent': 0, 'Isotope': 3}, 
            3: {'parent': 4, 'Isotope': 3}, 
        },
        'parent_idx': {  
            0: {'kids': [1, 2], 'CIPCode': 'R'}, 
            4: {'kids': [3], 'CIPCode': None}, 
        }
    }
    """

    # initialize 
    isotope_data = {
        'H_isotope_idx': {}, 
        'parent_idx': {}
    }
    H_isotope_idxs = [
        atom.GetIdx() for atom in mol.GetAtoms() 
        if is_h_isotope(atom)
    ]
    parent_idxs = []
    
    # populate H_isotope section
    for idx in H_isotope_idxs: 
        atom = mol.GetAtomWithIdx(idx)
        parent_idx = [nei.GetIdx() for nei in atom.GetNeighbors()][0]
        isotope_data['H_isotope_idx'][idx] = {
            'parent': parent_idx, 'Isotope': atom.GetIsotope()
        }
        if parent_idx not in parent_idxs:
            parent_idxs.append(parent_idx)
    
    # populate parent section
    for idx in parent_idxs: 
        atom = mol.GetAtomWithIdx(idx)
        kids = [nei.GetIdx() for nei in atom.GetNeighbors() if is_h_isotope(nei)]
        cip_code = atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None
        isotope_data['parent_idx'][idx] = {
            'kids': kids, 'CIPCode': cip_code
        }

    # in an editable copy, remove isotopes and add back as regular Hs
    editable_mol = Chem.RWMol(mol)
    indices_to_remove = sorted(H_isotope_idxs, reverse=True)
    parent_idxs_shifted = parent_idxs.copy()
    atom_idxs_original = [i for i in range(mol.GetNumAtoms()) if i not in indices_to_remove]
    atom_idxs_shifted = atom_idxs_original.copy()
    for idx in indices_to_remove:
        parent_atom = editable_mol.GetAtomWithIdx(
            isotope_data['H_isotope_idx'][idx]['parent']
        )
        editable_mol.RemoveAtom(idx)
        parent_atom.SetNumExplicitHs(parent_atom.GetNumExplicitHs() + 1)

        # update parent idx
        for i, parent_idx in enumerate(parent_idxs_shifted):
            if parent_idx > idx: 
                parent_idxs_shifted[i] += -1
        for i, atom_idx_shifted in enumerate(atom_idxs_shifted): 
            if atom_idx_shifted > idx:
                atom_idxs_shifted[i] += -1
        
    copy_mol = editable_mol.GetMol()
    shifted_conf = Chem.Conformer(copy_mol.GetNumAtoms())
    for i, atom_idx_shifted in enumerate(atom_idxs_shifted): 
        atom_idx_original = atom_idxs_original[i]
        shifted_conf.SetAtomPosition(atom_idx_shifted, conf.GetAtomPosition(atom_idx_original))
    copy_mol.RemoveAllConformers()
    copy_mol.AddConformer(shifted_conf)
    copy_mol = Chem.AddHs(copy_mol, addCoords=True)
    copy_conf = copy_mol.GetConformer()
    # assign H isotope coordinates from the regularized equivalent H
    assigned_isotope_pos = {}
    for ip, parent_idx in enumerate(parent_idxs):
        parent_atom = copy_mol.GetAtomWithIdx(parent_idxs_shifted[ip])
        kids_all = [nei for nei in parent_atom.GetNeighbors() if nei.GetAtomicNum() == 1]
        
        # have just enough regular Hs
        num_kids_needed = len(isotope_data['parent_idx'][parent_idx]['kids'])
        kids_all = kids_all[:num_kids_needed]
        kids_idxs_all = [nei.GetIdx() for nei in kids_all]
        for ik, kid in enumerate(isotope_data['parent_idx'][parent_idx]['kids']):
            # initial assignment, in original order
            copy_mol.GetAtomWithIdx(kids_idxs_all[ik]).SetIsotope(
                isotope_data["H_isotope_idx"][kid]['Isotope']
            )
            assigned_isotope_pos[kid] = copy_conf.GetAtomPosition(kids_idxs_all[ik])
        
        expected_cip_code = isotope_data['parent_idx'][parent_idx]['CIPCode']
        if len(kids_all) >1 and expected_cip_code is not None: 
            # check CIPCode
            Chem.rdmolops.AssignStereochemistry(copy_mol, force=True, cleanIt=True)
            if parent_atom.GetProp('_CIPCode')==expected_cip_code: 
                continue
            # catch some unsolvable problems
            kids_types = [nei.GetIsotope() for nei in kids_all]
            if any(len(lst) <= 1 for lst in (
                kids_all, 
                isotope_data['parent_idx'][parent_idx]['kids'],
                set(kids_types),
            )): 
                raise RuntimeError(
                    f"Unable to recover original chirality by manipulating H sotope positions: \n"
                    f"Atom # ({parent_idx}) \n"
                    f"Current CIPCode ({parent_atom.GetProp('_CIPCode')}) \n" 
                    f"Expected CIPCode: ({expected_cip_code})\n"
                    "Its chirality might have changed due to re-arrangements of heavy atoms, "
                    "or become ambiguous to Chem.rdmolops.AssignStereochemistry due to strained geometry. "
                )
            # swap assignment of max and min isotope, re-evaluate CIPCode
            max_idx = kids_types.index(max(kids_types))
            min_idx = kids_types.index(min(kids_types))
            copy_mol.GetAtomWithIdx(kids_idxs_all[max_idx]).SetIsotope(kids_types[min_idx])
            copy_mol.GetAtomWithIdx(kids_idxs_all[min_idx]).SetIsotope(kids_types[max_idx])
            Chem.rdmolops.AssignStereochemistry(copy_mol, force=True, cleanIt=True)
            if parent_atom.GetProp('_CIPCode')!=expected_cip_code: 
                raise RuntimeError(
                    "Failed to recover original chirality after attempts to manipulate H sotope positions: \n"
                    f"Atom # ({parent_idx}) \n"
                    f"Current CIPCode ({parent_atom.GetProp('_CIPCode')}) \n" 
                    f"Expected CIPCode: ({expected_cip_code})\n"
                    "Its chirality might have changed due to re-arrangements of heavy atoms, "
                    "or become ambiguous to Chem.rdmolops.AssignStereochemistry due to strained geometry. "
                )
            max_pos = copy_conf.GetAtomPosition(kids_idxs_all[max_pos])
            min_pos = copy_conf.GetAtomPosition(kids_idxs_all[min_pos])
            H_isotope_idx_with_maxpos = isotope_data['parent_idx'][parent_idx]['kids'][max_pos]
            assigned_isotope_pos[H_isotope_idx_with_maxpos] = min_pos
            H_isotope_idx_with_minpos = isotope_data['parent_idx'][parent_idx]['kids'][min_pos]
            assigned_isotope_pos[H_isotope_idx_with_minpos] = max_pos

    if return_mol is False: 
        return assigned_isotope_pos
    
    return copy_mol
        
def getPdbInfoNoNull(atom):
    """extract information for populating an ATOM/HETATM line
    in the PDB"""
    minfo = atom.GetMonomerInfo()  # same as GetPDBResidueInfo
    if minfo is None:
        atomic_number = atom.GetAtomicNum()
        if atomic_number == 0:
            name = "%-2s" % "*"
        else:
            name = "%-2s" % mini_periodic_table[atomic_number]
        chain = " "
        resNum = 1
        icode = ""
        resName = "UNL"
    else:
        name = minfo.GetName()
        chain = minfo.GetChainId()
        resNum = minfo.GetResidueNumber()
        icode = minfo.GetInsertionCode()
        resName = minfo.GetResidueName()
    return PDBAtomInfo(
        name=name, resName=resName, resNum=resNum, icode=icode, chain=chain
    )


class Mol2MolSupplier:
    """RDKit Mol2 molecule supplier.
    Parameters
        sanitize: perform RDKit sanitization of Mol2 molecule"""

    def __init__(
        self, filename, sanitize=True, removeHs=False, cleanupSubstructures=True
    ):
        self.fp = open(filename, "r")
        self._opts = {
            "sanitize": sanitize,
            "removeHs": removeHs,
            "cleanupSubstructures": cleanupSubstructures,
        }
        self.buff = []

    def __iter__(self):
        return self

    def __next__(self):
        """iterator step"""
        while True:
            line = self.fp.readline()
            # empty line
            if not line:
                if len(self.buff):
                    # buffer full, returning last molecule
                    mol = Chem.MolFromMol2Block("".join(self.buff), **self._opts)
                    self.buff = []
                    return mol
                # buffer empty, stopping the iteration
                self.fp.close()
                raise StopIteration
            if "@<TRIPOS>MOLECULE" in line:
                # first molecule parsed
                if len(self.buff) == 0:
                    self.buff.append(line)
                else:
                    # found the next molecule, breaking to return the complete one
                    break
            else:
                # adding another line in the current molecule
                self.buff.append(line)
        # found a complete molecule, returning it
        mol = Chem.MolFromMol2Block("".join(self.buff), **self._opts)
        self.buff = [line]
        return mol

class AtomField:
    """Stores data parsed from PDB or mmCIF"""

    def __init__(
        self,
        atomname: str,
        altloc: str,
        resname: str,
        chain: str,
        resnum: int,
        icode: str,
        x: float,
        y: float,
        z: float,
        element: str,
    ):
        self.atomname = atomname
        self.altloc = altloc
        self.resname = resname
        self.chain = chain
        self.resnum = resnum
        self.icode = icode
        self.x = x
        self.y = y
        self.z = z
        if len(element) > 1:
            element = f"{element[0].upper()}{element[1].lower()}"
        else:
            element = f"{element.upper()}"
        self.atomic_nr = periodic_table.GetAtomicNumber(element)


def _build_rdkit_mol_for_altloc(atom_fields_list, wanted_altloc:str=None):
    mol = Chem.EditableMol(Chem.Mol())
    mol.BeginBatchEdit() 
    positions = []
    idx_to_rdkit = {}
    for index_list, atom in enumerate(atom_fields_list):
        if wanted_altloc is not None:
            if atom.altloc and atom.altloc != wanted_altloc:
                # if atom.altloc is "" we still want to consider this atom
                continue
        rdkit_atom = Chem.Atom(atom.atomic_nr)
        positions.append(Point3D(atom.x, atom.y, atom.z))
        res_info = Chem.AtomPDBResidueInfo()
        res_info.SetName(atom.atomname)
        res_info.SetResidueName(atom.resname)
        res_info.SetResidueNumber(atom.resnum)
        res_info.SetChainId(atom.chain)
        res_info.SetInsertionCode(atom.icode)
        rdkit_atom.SetPDBResidueInfo(res_info)
        index_rdkit = mol.AddAtom(rdkit_atom)
        idx_to_rdkit[index_list] = index_rdkit
    mol.CommitBatchEdit()
    mol = mol.GetMol()
    conformer = Chem.Conformer(mol.GetNumAtoms())
    for index, position in enumerate(positions):
        conformer.SetAtomPosition(index, position)
    mol.AddConformer(conformer, assignId=True)
    return mol, idx_to_rdkit
        

def build_one_rdkit_mol_per_altloc(atom_fields_list):
    """ if no altlocs, the only key in the output dict is None
        if altlocs exist, None is not a key: the keys are the altloc IDs
    """
    altlocs = set([atom.altloc for atom in atom_fields_list if atom.altloc])
    rdkit_mol_dict = {}
    if not altlocs:
        altlocs = {None}
    for altloc in altlocs:
        mol, idx_to_rdkit = _build_rdkit_mol_for_altloc(atom_fields_list, altloc)
        rdkit_mol_dict[altloc] = (mol, idx_to_rdkit)
    return rdkit_mol_dict


def _aux_altloc_mol_build(atom_field_list, requested_altloc, default_altloc):
    missed_altloc = False
    needed_altloc = False
    mols_dict = build_one_rdkit_mol_per_altloc(atom_field_list) 
    has_altloc = None not in mols_dict
    if has_altloc and requested_altloc is None and default_altloc is None:
        pdbmol = None
        missed_altloc = False 
        needed_altloc = True
    elif requested_altloc and requested_altloc in mols_dict:
        pdbmol, idx_to_rdkit = mols_dict[requested_altloc]
    elif requested_altloc and requested_altloc not in mols_dict:
        pdbmol = None
        missed_altloc = True
        needed_altloc = False
    elif default_altloc and default_altloc in mols_dict:
        pdbmol, idx_to_rdkit = mols_dict[default_altloc]
    elif has_altloc and default_altloc not in mols_dict:
        pdbmol = None
        missed_altloc = True
        needed_altloc = False
    elif not has_altloc and requested_altloc is None:
        pdbmol, idx_to_rdkit = mols_dict[None]
    else:
        raise RuntimeError("programming bug, please post full error on github")
    if pdbmol is None: 
        idx_to_rdkit = None
        return pdbmol, idx_to_rdkit, missed_altloc, needed_altloc
    else:
        rdDetermineBonds.DetermineConnectivity(pdbmol)
        for atom in pdbmol.GetAtoms():
            if atom.GetAtomicNum() == 7 and len(atom.GetNeighbors()) == 4:
                atom.SetFormalCharge(1)
        _ = Chem.SanitizeMol(pdbmol)

    return pdbmol, idx_to_rdkit, missed_altloc, needed_altloc

def react_and_map(reactants: tuple[Chem.Mol], rxn: rdChemReactions.ChemicalReaction):
    """
    Run a reaction and keep track of atom indices from reactants to products.
    
    Parameters
    ----------
    reactants : tuple[Chem.Mol]
        A tuple of RDKit molecule objects representing the reactants.
    rxn : rdChemReactions.ChemicalReaction
        The RDKit reaction object.
        
    Returns
    -------
    list[tuple[Chem.Mol, dict[str, list[Optional[int]]]]]
        A list of tuples where each tuple contains a product molecule and a dictionary.
        The dictionary has keys 'atom_idx' and 'new_atom_label', which are ordered lists for product atoms:
        - 'atom_idx' holds the corresponding atom indices in reactant. None for newly added atoms. 
        - 'new_atom_label' holds the reaction mapping number, only for newly added atoms. 
    """

    # Prepare for multiple possible outcomes resulted from multiple matched reactive sites in reactant
    outcomes = []
    for products in rxn.RunReactants(reactants): 
        # Assumes single product 
        product = products[0]
        # For each atom, get react_atom_idx if they were in reactant
        atom_idxmap = [
            atom.GetIntProp("react_atom_idx") if atom.HasProp("react_atom_idx")
            else None
            for atom in product.GetAtoms()
        ]
        # For each atom, get the rxn mapping number if the were added in the rxn
        new_atom_label = [
            atom.GetIntProp("old_mapno") if atom.HasProp("old_mapno") and not atom.HasProp("react_atom_idx")
            else None
            for atom in product.GetAtoms()
        ]
        # Collect product and index_map
        index_map = {"atom_idx": atom_idxmap, "new_atom_label": new_atom_label}
        outcomes.append((product, index_map))

    return outcomes


covalent_radius = {  # from wikipedia
    1: 0.31,
    5: 0.84,
    6: 0.76,
    7: 0.71,
    8: 0.66,
    9: 0.57,
    12: 0.00,  # hack to avoid bonds with metals
    14: 1.11,
    15: 1.07,
    16: 1.05,
    17: 1.02,
    # 19: 2.03,
    20: 0.00,
    # 24: 1.39,
    25: 0.00,  # hack to avoid bonds with metals
    26: 0.00,
    30: 0.00,  # hack to avoid bonds with metals
    # 34: 1.20,
    35: 1.20,
    53: 1.39,
}
