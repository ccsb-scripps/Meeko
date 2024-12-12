#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#

import argparse
import gzip
import pathlib
import sys
eol="\n"
import warnings
import copy
import numpy as np

from rdkit import Chem

from meeko import PDBQTMolecule
from meeko import RDKitMolCreate
from meeko import Polymer
from meeko.export_flexres import export_pdb_updated_flexres, pdb_updated_flexres_from_rdkit
from meeko.utils.utils import parse_begin_res
from meeko.covalentbuilder import get_fragments_by_atom_indices
from meeko.polymer import Monomer


def cmd_lineparser():
    parser = argparse.ArgumentParser(
        description='Export docked ligand to SDF, and receptor to PDB',
    )
    parser.add_argument(dest='docking_results_filename', nargs = "+",
                        help='Docking output file(s), either PDBQT \
                        file from Vina or DLG file from AD-GPU.')
    parser.add_argument('-s', '--write_sdf', metavar='output_SDF_filename',
                        help="defaults to input filename with suffix from --sufix")
    parser.add_argument(
        '-p',
        '--write_pdb',
        metavar='output_PDB_filename',
        help="defaults to input filename with suffix from --suffix",
    )
    parser.add_argument("--suffix", default="_docked",
                        help="suffix for output filesnames that are not explicitly set")
    parser.add_argument('-j', '--read_json', metavar='input_JSON_filename',
                        help="receptor written by mk_prepare_receptor -j/--write_json")
    parser.add_argument('--all_dlg_poses', action='store_true',
                        help="write all AutoDock-GPU poses, not just cluster leads.")
    parser.add_argument('-k', '--keep_flexres_sdf', action='store_true',
                        help="add flexres, if any, to SDF ouput")
    parser.add_argument('-', '--',  dest='redirect_stdout', action='store_true',
                        help="do not write SDF file, just print it to STDOUT")
    return parser.parse_args()

def main():
    args = cmd_lineparser()
    
    docking_results_filenames = args.docking_results_filename
    write_sdf = args.write_sdf
    write_pdb = args.write_pdb
    read_json = args.read_json
    suffix = args.suffix
    all_dlg_poses = args.all_dlg_poses
    keep_flexres_sdf = args.keep_flexres_sdf
    redirect_stdout = args.redirect_stdout
    
    if (
        (write_sdf is not None or write_pdb is not None) and
        len(docking_results_filenames) > 1
    ):
        msg = "With multiple input files, the output filenames are based on the"
        msg += "input filename. The suffix can be controlled with option --suffix."
        msg += "--write options are incompatible with multiple input files"
        print("--write options incompatible with multiple input files", file=sys.stderr)
        sys.exit(2)
    
    if redirect_stdout and len(docking_results_filenames) > 1:
        print("option -/-- incompatible with multiple input files", file=sys.stderr)
        sys.exit(2)
    
    if read_json is not None:
        with open(read_json) as f:
            json_string = f.read()
        polymer = Polymer.from_json(json_string)
    else:
        polymer = None
        if write_pdb is not None:
            print("option -p (write pdb) requires -j (receptor receptor file)", file=sys.stderr)
            sys.exit(2)
    
    
    for filename in docking_results_filenames:
        is_dlg = filename.endswith('.dlg') or filename.endswith(".dlg.gz")
        if filename.endswith(".gz"):
            with gzip.open(filename) as f:
                string = f.read().decode()
        else:
            with open(filename) as f:
                string = f.read()
        mol_name = pathlib.Path(filename).with_suffix("").name
        pdbqt_mol = PDBQTMolecule(
            string,
            name=mol_name,
            is_dlg=is_dlg,
            skip_typing=True,
        )
        only_cluster_leads = is_dlg and not all_dlg_poses
        sdf_string, failures = RDKitMolCreate.write_sd_string(
                pdbqt_mol,
                only_cluster_leads=only_cluster_leads,
                keep_flexres=keep_flexres_sdf,
        )
        for i in failures:
            warnings.warn("molecule %d not converted to RDKit/SD File" % i)
        if sdf_string == "": 
            warnings.warn("sdf_string does not contain molecular data.")
            if redirect_stdout or write_pdb: 
                pass
            else:
                print("Output SDF will not be created because there is no pose data for ligand. \n"
                        + "Maybe the input poses only contain flexible sidechains and \n"
                        + "keep_flexres_sdf is set to False. \n"
                        + "Use -k with mk_export.py to retain the flexres and write to output SDF File. ")
                sys.exit(2)
        if len(failures) == len(pdbqt_mol._atom_annotations["mol_index"]):
            msg = "\nCould not convert to RDKit. Maybe meeko was not used for preparing\n"
            msg += "the input PDBQT for docking, and the SMILES string is missing?\n"
            msg += "Except for standard protein sidechains, all ligands and flexible residues\n"
            msg += "require a REMARK SMILES line in the PDBQT, which is added automatically by meeko."
            raise RuntimeError(msg)
        if not redirect_stdout:
            if write_sdf is None:
                fn = pathlib.Path(filename).with_suffix("").name + f"{suffix}.sdf"
            else:
                fn = write_sdf
            with open(fn, "w") as f:
                f.write(sdf_string)
        else:
            print(sdf_string)
    
        # write receptor with updated flexres
        if read_json is not None:
            has_covlig = False
            
            # region modifies polymer if there's a covalent ligand
            flexres_id = pdbqt_mol._pose_data["mol_index_to_flexible_residue"]
            for mol_idx, _ in pdbqt_mol._atom_annotations["mol_index"].items(): 
                has_covlig = True
                pdbqt_mol_copy = copy.deepcopy(pdbqt_mol)
                # if it's a covalent ligand as flexres, it will have smiles and index_map in pose_data
                if flexres_id[mol_idx] is not None and pdbqt_mol._pose_data["smiles"][mol_idx] is not None: 
                    covligmol = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol_copy,
                                                              only_cluster_leads=False, 
                                                              keep_flexres=True)[mol_idx]
                    index_map = pdbqt_mol_copy._pose_data['smiles_index_map'][mol_idx]
                    
                    # locate tethered atoms in ligmol and get xyz
                    tethered_mol_index = [int(index_map[i * 2]) - 1 for i in (0,1)]
                    tethered_xyz = [covligmol.GetConformer().GetPositions()[i] for i in tethered_mol_index]

                    # match with ignored monomer atoms
                    res_id = parse_begin_res(flexres_id[mol_idx])
                    tolerance = 1e-2
                    monomer_rdkit_mol = polymer.monomers[res_id].rdkit_mol
                    monomer_conformer = monomer_rdkit_mol.GetConformer()

                    attractor_mol_index = []
                    for target_xyz in tethered_xyz: 
                        target_x, target_y, target_z = target_xyz
                        for atom in monomer_rdkit_mol.GetAtoms():
                            idx = atom.GetIdx()
                            pos = monomer_conformer.GetAtomPosition(idx)
                            if (abs(pos.x - target_x) <= tolerance and
                                abs(pos.y - target_y) <= tolerance and
                                abs(pos.z - target_z) <= tolerance):
                                attractor_mol_index.append(idx)
                                break
                    
                    # breaks monomer's rdkit_mol by the attractors bond and keep only the receptor-side fragment
                    keep_mol, _ = get_fragments_by_atom_indices(monomer_rdkit_mol, attractor_mol_index[0],
                                                                attractor_mol_index[1], get_as_mols=True)
                    polymer.monomers[res_id].rdkit_mol = keep_mol

                    # creates a new monomer from pdbqt_mol and keeps only the ligand-side fragment
                    _, keep_lig = get_fragments_by_atom_indices(covligmol, tethered_mol_index[0],
                                                                tethered_mol_index[1], get_as_mols=True)
                    # _, keep_indices = get_fragments_by_atom_indices(covligmol, tethered_mol_index[0],
                    #                                            tethered_mol_index[1], get_as_mols=False)
                    
                    # only if we need to export the updated polymer...
                    # new_chid = res_id.split(":")[0]
                    # new_resnum = max([int(key.split(":")[1]) for key in polymer.monomers if key.split(":")[0] == new_chid]) + 1
                    # new_res_id = f"{new_chid}:{new_resnum}"
                    # polymer.monomers[new_res_id] = Monomer(covligmol, keep_lig, 
                    #                                        {i:idx for i,idx in enumerate(keep_indices)},
                    #                                         atom_names=["x"]*len(keep_indices))
                    # polymer.monomers[new_res_id].input_resname = "LIG"
            # endregion 
            
            pdb_string = ""
            pose_id_to_iter = [pose.pose_id for pose in pdbqt_mol]
            for pose_id in pose_id_to_iter:
                model_nr = pose_id + 1
                pdb_string += "MODEL " + f"{model_nr:8}" + eol
                pdb_string += export_pdb_updated_flexres(polymer, pdbqt_mol)
                if has_covlig: 
                    pdb_string += Chem.MolToPDBBlock(keep_lig, confId=pose_id)
                pdb_string += "ENDMDL" + eol
            if write_pdb is None:
                fn = pathlib.Path(filename).with_suffix("").name + f"{suffix}.pdb"
            else:
                fn = write_pdb
            with open(fn, "w") as f:
                f.write(pdb_string)
    return

if __name__ == "__main__":
    sys.exit(main())
