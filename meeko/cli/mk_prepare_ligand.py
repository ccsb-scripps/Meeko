#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#

import argparse
from datetime import datetime
import io
import os
eol="\n"
import sys
import json
import tarfile
import warnings

from rdkit import Chem

from meeko.utils.utils import parse_cmdline_res_assign
from meeko import MoleculePreparation
from meeko import rdkitutils
from meeko import PDBQTWriterLegacy
from meeko import Polymer
from meeko import PolymerCreationError
from meeko import ResidueChemTemplates

try:
    import prody
    from meeko import CovalentBuilder

    _prody_parsers = {"pdb": prody.parsePDB, "mmcif": prody.parseMMCIF}
except ImportError as err:
    _has_prody = False
    _prody_parsers = {}
else:
    _has_prody = True
supported_recfile_formats = set(["json", "pdb"]) | set(_prody_parsers.keys())
need_prody_msg = ""
if not _has_prody:
    need_prody_msg = ". Parsing of receptor mmCIF files needs Prody which can be installed from PyPI or conda-forge"


def cmd_lineparser():

    conf_parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False,
    )
    conf_parser.add_argument(
        "-c",
        "--config_file",
        help="configure MoleculePreparation from JSON file. Overriden by command line args.",
    )
    confargs, remaining_argv = conf_parser.parse_known_args()

    parser = (
        argparse.ArgumentParser()
    )  # parents=[conf_parser]) # parents shows --config_file in help msg
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="print information about molecule setup",
    )

    io_group = parser.add_argument_group("Input/Output")
    io_group.add_argument(
        "-i",
        "--mol",
        dest="input_molecule_filename",
        required=True,
        action="store",
        help="molecule file (MOL2, SDF,...)",
    )
    io_group.add_argument(
        "--name_from_prop",
        help="set molecule name from RDKit/SDF property",
    )
    io_group.add_argument(
        "-o",
        "--out",
        dest="output_pdbqt_filename",
        action="store",
        help="output pdbqt filename. Single molecule input only.",
    )
    io_group.add_argument(
        "--multimol_outdir",
        dest="multimol_output_dir",
        action="store",
        help="folder to write output pdbqt for multi-mol inputs. Incompatible with -o/--out and -/--.",
    )
    io_group.add_argument(
        "--multimol_prefix",
        dest="multimol_prefix",
        action="store",
        help="replace internal molecule name in multi-molecule input by specified prefix. Incompatible with -o/--out and -/--.",
    )
    io_group.add_argument(
        "-z",
        "--multimol_targz",
        action="store_true",
        help="compress output files in .tar.gz",
    )
    io_group.add_argument(
        "--multimol_targz_size",
        type=int,
        default=10000,
        help="number of PDBQT files per .tar.gz",
    )
    io_group.add_argument(
        "-",
        "--",
        dest="redirect_stdout",
        action="store_true",
        help="do not write file, redirect output to STDOUT. Argument -o/--out is ignored. Single molecule input only.",
    )

    config_group = parser.add_argument_group("Molecule preparation")
    config_group.add_argument(
        "-c",
        "--config_file",
        help="configure MoleculePreparation from JSON file. Overriden by command line args.",
    )  # parsed above by conf_parser, here for help msg
    config_group.add_argument(
        "--rigid_macrocycles",
        dest="rigid_macrocycles",
        action="store_true",
        help="keep macrocycles rigid in input conformation",
    )
    config_group.add_argument(
        "--macrocycle_allow_A",
        action="store_true",
        help="allow bond break with atom type A, retyped as C",
    )
    config_group.add_argument(
        "--keep_chorded_rings",
        dest="keep_chorded_rings",
        action="store_true",
        help="return all rings from exhaustive perception",
    )
    config_group.add_argument(
        "--keep_equivalent_rings",
        dest="keep_equivalent_rings",
        action="store_true",
        help="equivalent rings have the same size and neighbors",
    )
    config_group.add_argument(
        "--min_ring_size",
        dest="min_ring_size",
        type=int,
        help="min nr of atoms in ring for opening",
    )
    config_group.add_argument(
        "-w",
        "--hydrate",
        dest="hydrate",
        action="store_true",
        help="add water molecules for hydrated docking",
    )
    config_group.add_argument(
        "--merge_these_atom_types",
        dest="merge_these_atom_types",
        nargs="*",
        help='list of atom types to merge, default is "H"',
        default=["H"],
    )
    config_group.add_argument(
        "-r",
        "--rigidify_bonds_smarts",
        dest="rigidify_bonds_smarts",
        action="append",
        help="SMARTS patterns to rigidify bonds",
        metavar="SMARTS",
    )
    config_group.add_argument(
        "-b",
        "--rigidify_bonds_indices",
        dest="rigidify_bonds_indices",
        action="append",
        help="indices of two atoms (in the SMARTS) that define a bond (start at 1)",
        nargs="+",
        type=int,
        metavar="i j",
    )
    config_group.add_argument(
        "-a",
        "--flexible_amides",
        dest="flexible_amides",
        action="store_true",
        help="allow amide bonds to rotate and be non-planar, which is bad",
    )
    config_group.add_argument(
        "-p",
        "--load_atom_params",
        nargs="+",
        action="store",
        default=["ad4_types"],
        help="filename with SMARTS defined atom types (JSON format)",
    )
    config_group.add_argument(
        "-aa",
        "--add_atom_types",
        dest="add_atom_types_json",
        action="append",
        help="Additional atom types to assign (JSON formated)",
        metavar="[{'smarts': '<smarts pattern>', 'atype': ',atomtype name>'}, {'smarts': '<smarts pattern>', 'atype': ',atomtype name>'}]",
    )
    config_group.add_argument(
        "--double_bond_penalty",
        help="penalty > 100 prevents breaking double bonds",
        type=int,
    )
    config_group.add_argument(
        "--charge_model",
        choices=("gasteiger", "espaloma", "zero", "read"),
        help="default is 'gasteiger', 'zero' sets all zeros",
        default="gasteiger",
    )
    config_group.add_argument(
        "--charge_atom_prop",
        help="set atom partial charges from an RDKit atom property based on the input file. The default is 'PartialCharge' for SDF and '_TriposPartialCharge' for MOL2 unless overriden by a user defined property name. ",
    )
    config_group.add_argument(
        "--bad_charge_ok",
        help="NaN and Inf charges allowed in PDBQT",
        action="store_true",
    )
    config_group.add_argument(
        "--add_index_map",
        dest="add_index_map",
        action="store_true",
        help="write map of atom indices from input to pdbqt",
    )
    config_group.add_argument(
        "--remove_smiles",
        dest="remove_smiles",
        action="store_true",
        help="do not write smiles as remark to pdbqt",
    )
    config_group.add_argument(
        "--rename_atoms",
        dest="rename_atoms",
        action="store_true",
        help="rename atoms: the new name will be the original name and its (1-based) index in MoleculeSetup",
    )
    reactive_group = parser.add_argument_group("Reactive docking")
    reactive_group.add_argument(
        "--reactive_smarts", help="SMARTS pattern for reactive group"
    )
    reactive_group.add_argument(
        "--reactive_smarts_idx",
        help="index (1-based) of the reactive atom in --reactive_smarts",
        type=int,
    )

    covalent_group = parser.add_argument_group(
        "Covalent docking (tethered)"
    )
    covalent_group.add_argument(
        "--receptor",
        help="receptor filename. Supported formats: [%s]%s"
        % ("/".join(supported_recfile_formats), need_prody_msg),
    )
    covalent_group.add_argument("--add_templates", help="include additional residue chemical template files [.json]", metavar="JSON_FILENAME", nargs="+", default=[])
    covalent_group.add_argument("--default_altloc", help="default alternate location (overridden by --wanted_altloc)")
    covalent_group.add_argument("--wanted_altloc", help="require altloc for specific residues, e.g. :5=B,B:17=A")
    covalent_group.add_argument(
        "--rec_residue", help='examples: "A:LYS:204", "A:HIS:" (all HIS in chain A), ":LYS:" (all LYS)'
    )
    covalent_group.add_argument(
        "--rec_attractor_names", 
        type=str,
        nargs=2,
        metavar="ATOM_NAME",
        default=["CA", "CB"],
        help=("Atom names of the two \"attractor\" atoms within the receptor. "
        " (default: CA CB) " 
        "  The closer atom to the receptor core needs to be put first, "
        "  followed by the other atom that is closer to the ligand. "
        )
    )
    covalent_group.add_argument(
        "--rec_attractor_smarts", 
        help=("SMARTS pattern to define receptor atoms for ligand attachment. "
        "Overrides the default atom names (CA CB) if --rec_attractor_names is not explicitly provided. "
        )
    )
    covalent_group.add_argument(
        "--rec_smarts_indices",
        type=int,
        nargs=2,
        default=[1, 2],
        help=("indices (1-based) of the receptor SMARTS atoms that will be attached. "
        " (default: 1 2) " 
        "  The closer atom to the receptor core needs to be put first, "
        "  followed by the other atom that is closer to the ligand. "
        )
    )
    covalent_group.add_argument(
        "--tether_smarts",
        help="SMARTS pattern to define ligand atoms for receptor attachment. "
    )
    covalent_group.add_argument(
        "--tether_smarts_indices",
        type=int,
        nargs=2,
        metavar="IDX",
        default=[1, 2],
        help=("indices (1-based) of the ligand SMARTS atoms that will be attached. "
        " (default: 1 2) " 
        "  The closer atom to the receptor core needs to be put first, "
        "  followed by the other atom that is closer to the ligand. "
        )
    )

    config = MoleculePreparation.get_defaults_dict()

    if confargs.config_file is not None:
        with open(confargs.config_file) as f:
            c = json.load(f)
            config.update(c)

    # Command line arguments should override the config file only for options
    # set explicitly by the user. The config file still has priority over the
    # defaults set in argparse. To achieve this, we reset the argparse
    # defaults to the values from the config file. Then, we can just update
    # variable `config` with the values parsed with argparse
    parser.set_defaults(**config)
    args = parser.parse_args(remaining_argv)

    # check reactive arguments
    if (args.reactive_smarts is None) != (args.reactive_smarts_idx is None):
        print(
            "Arguments --reactive_smarts and --reactive_smarts_idx require each other",
            file=sys.stderr,
        )
        sys.exit(2)
    elif args.reactive_smarts_idx is not None:
        if args.reactive_smarts_idx < 1:
            print(
                "--reactive_smarts_idx is 1-indexed, but got %d"
                % args.reactive_smarts_idx,
                file=sys.stderr,
            )
            sys.exit(2)
        args.reactive_smarts_idx -= 1  # convert from 1- to 0-index

    # This is where command line arguments override config file.
    # Relies on key/parameter names being equal.
    # Deliberate mismatch for add_atom_types/add_atom_types_json, as these
    # are extended below instead of being replaced
    for key in config:
        if key in args.__dict__:
            config[key] = args.__dict__[key]

    if args.add_atom_types_json is not None:
        additional_ats = []
        for at in args.add_atom_types_json:
            at = json.loads(at)
            if type(at) == dict:
                additional_ats.append(at)
            elif type(at) == list:
                additional_ats.extend(at)
        config["add_atom_types"] += additional_ats

    if args.multimol_output_dir is not None or args.multimol_prefix is not None:
        if args.output_pdbqt_filename is not None:
            print(
                "Warning: -o/--out ignored with --multimol_outdir or --multimol_prefix",
                file=sys.stderr,
            )
        if args.redirect_stdout:
            print(
                "Warning: -/-- ignored with --multimol_outdir or --multimol_prefix",
                file=sys.stderr,
            )

    # verify sanity of covalent docking input
    num_required_covalent_args = 0
    num_required_covalent_args += int(args.receptor is not None)
    num_required_covalent_args += int(args.rec_residue is not None)
    num_required_covalent_args += int(args.tether_smarts is not None)
    if num_required_covalent_args not in [0, 3]:
        print(
            "Error: --receptor, --rec_residue, and --tether_smarts are all required for covalent docking.",
            file=sys.stderr,
        )
        sys.exit(2)
    is_covalent = num_required_covalent_args == 3

    if is_covalent: 
        # verify ways to specify receptor attractor atoms
        provided_names = '--rec_attractor_names' in remaining_argv
        provided_smarts = '--rec_attractor_smarts' in remaining_argv

        if provided_names and provided_smarts: 
        # when both are explicitly provided
            print(
                "Error: --rec_attractor_names and --rec_attractor_smarts are conflicting arguments. ",
                file=sys.stderr,
            )
            sys.exit(2)

        if not provided_names and not provided_smarts: 
        # when both are not explicitly provided
        # use the default of --rec_attractor_names 
            provided_names = True

        if provided_smarts: 
            if min(args.rec_smarts_indices) < 1:
                print(
                "--rec_smarts_indices is 1-indexed, all values must be greater than zero",
                file=sys.stderr,
                )
                sys.exit(2)
            args.rec_smarts_indices = [
                i - 1 for i in args.rec_smarts_indices
            ]  # convert to 0-index

        if min(args.tether_smarts_indices) < 1:
            print(
                "--tether_smarts_indices is 1-indexed, all values must be greater than zero",
                file=sys.stderr,
            )
            sys.exit(2)
        args.tether_smarts_indices = [
            i - 1 for i in args.tether_smarts_indices
        ]  # convert to 0-index

    # verify sanity of SMARTS patterns to make bonds rigid and convert to 0-based indices
    rigidify_bonds_smarts = config["rigidify_bonds_smarts"]
    rigidify_bonds_indices = config["rigidify_bonds_indices"]
    if len(rigidify_bonds_indices) != len(rigidify_bonds_smarts):
        raise RuntimeError(
            "length of --rigidify_bonds_indices differs from length of --rigidify_bonds_smarts"
        )
    for indices in rigidify_bonds_indices:
        if len(indices) != 2:
            raise RuntimeError(
                "--rigidify_bonds_indices must specify pairs, e.g. -b 1 2 -b 3 4"
            )
        indices[0] = indices[0] - 1  # convert from 1- to 0-index
        indices[1] = indices[1] - 1

    return args, config, is_covalent, provided_names


class Output:
    def __init__(self, multimol_output_dir, multimol_targz, multimol_targz_size, multimol_prefix, redirect_stdout, output_filename, name_from_prop):
        is_multimol = (multimol_prefix is not None) or (multimol_output_dir is not None)
        self._mkdir(multimol_output_dir)

        if multimol_output_dir is None:
            multimol_output_dir = "."
        self.multimol_output_dir = multimol_output_dir
        self.multimol_targz = multimol_targz
        self.multimol_targz_size = max(int(multimol_targz_size), 1)
        self.tarf = None
        self.tar_pdbqt_count = 0
        self.tarf_index = 0
        self.multimol_prefix = multimol_prefix
        self.redirect_stdout = redirect_stdout
        self.output_filename = output_filename
        self.is_multimol = is_multimol
        self.visited_filenames = set()
        self.duplicate_filenames = set()
        self.visited_names = set()
        self.duplicate_names = set()
        self.num_files_written = 0
        self.counter = 0
        self.name_from_prop = name_from_prop

    def _open_new_tar(self):
        self.tarf_index += 1
        prefix = "" if self.multimol_prefix is None else self.multimol_prefix
        tgz_path = os.path.join(
            self.multimol_output_dir,
            f"{prefix}{self.tarf_index:07d}.tar.gz"
        )
        if self.tarf is not None:
            self.tarf.close()
        tarf = tarfile.open(tgz_path, "w:gz")
        return tarf

    def _add_to_tar(self, pdbqt_string, name):
        if self.tarf is None or self.tar_pdbqt_count >= self.multimol_targz_size:
            self.tarf = self._open_new_tar()
            self.tar_pdbqt_count = 0
        tarinfo = tarfile.TarInfo(name=f"{name}.pdbqt")
        tarinfo.size = len(pdbqt_string)
        tarinfo.mtime = datetime.timestamp(datetime.now())
        self.tarf.addfile(tarinfo, io.BytesIO(pdbqt_string.encode()))
        self.tar_pdbqt_count += 1
        return

    def __call__(self, pdbqt_string, name, suffixes=()):
        self.counter += 1
        if name in self.visited_names:
            self.duplicate_names.add(name)
        self.visited_names.add(name)
        if self.multimol_prefix is not None:
            name = "%s-%d" % (self.multimol_prefix, self.counter)
        for suffix in suffixes:
            if suffix is not None and len(suffix) > 0:
                name += "_" + suffix

        if self.is_multimol:
            if name in self.visited_filenames:
                self.duplicate_filenames.add(name)
                repeat_id = 1
                newname = name + "-again%d" % repeat_id
                while newname in self.visited_filenames:
                    repeat_id += 1
                    newname = name + "-again%d" % repeat_id
                print(
                    "Renaming %s to %s to disambiguate filename" % (name, newname),
                    file=sys.stderr,
                )
                name = newname

            self.visited_filenames.add(name)

            if self.multimol_targz:
                self._add_to_tar(pdbqt_string, name)
            else:
                fpath = os.path.join(self.multimol_output_dir, name + ".pdbqt")
                with open(fpath, "w") as f:
                    print(pdbqt_string, end="", file=f)
            self.num_files_written += 1

        elif self.redirect_stdout:
            print(pdbqt_string, end="")
        else:
            if self.output_filename is None:
                filename = "%s.pdbqt" % name
            else:
                filename = self.output_filename
            with open(filename, "w") as f:
                print(pdbqt_string, end="", file=f)
            self.num_files_written += 1

    def _mkdir(self, multimol_output_dir):
        """make directory if it doesn't exist yet"""
        if multimol_output_dir is not None:
            if not os.path.exists(multimol_output_dir):
                os.mkdir(multimol_output_dir)

    def get_duplicates_info_string(self):
        if not self.is_multimol:
            return None
        if len(self.duplicate_names):
            # with multimol_prefix with can have duplicate molecule names,
            # but not duplicate filenames. This warning is for such cases.
            string = "Warning: %d molecules have duplicated names" % len(
                self.duplicate_names
            )
        else:
            string = "No duplicate molecule molecule names were found"
        # if we have duplicate_filenames, we also have duplicate molecule names,
        # but it suffices to return the string below
        if len(self.duplicate_filenames):
            string = (
                "Warning: %d molecules with repeated names were suffixed with -again<n>"
                % len(self.duplicate_filenames)
            )
        return string

    @staticmethod
    def get_suffixes(molsetups):
        if len(molsetups) == 1:
            return ("",)  # no suffix needed
        else:
            return tuple("mk%d" % (i + 1) for i in range(len(molsetups)))

def main():
    args, config, is_covalent, provided_names = cmd_lineparser()
    input_molecule_filename = args.input_molecule_filename

    # read input
    input_fname, ext = os.path.splitext(input_molecule_filename)
    ext = ext[1:].lower()

    parsers = {
        "sdf": Chem.SDMolSupplier,
        "mol2": rdkitutils.Mol2MolSupplier,
        "mol": Chem.SDMolSupplier,
    }
    if not ext in parsers:
        print(
            "Error: Format [%s] not in supported formats [%s]"
            % (ext, "/".join(list(parsers.keys())))
        )
        sys.exit(1)
    mol_supplier = parsers[ext](
        input_molecule_filename, removeHs=False
    )  # input must have explicit H
    
    # configure output writer
    if args.output_pdbqt_filename is None:
        output_filename = input_fname + ".pdbqt"
    else:
        output_filename = args.output_pdbqt_filename


    output = Output(
        args.multimol_output_dir,
        args.multimol_targz,
        args.multimol_targz_size,
        args.multimol_prefix,
        args.redirect_stdout,
        output_filename,
        args.name_from_prop,
    )

    # initialize covalent object for receptor
    if is_covalent:

        # region parsing receptor file and args into target_monomers
        rec_filename = args.receptor
        _, rec_extension = os.path.splitext(rec_filename)
        rec_extension = rec_extension[1:].lower()

        if rec_extension not in supported_recfile_formats: 
            print(
            "Error: Format [%s] not in supported receptor formats [%s]%s"
            % (rec_extension, "/".join(supported_recfile_formats), need_prody_msg)
            )
            sys.exit(1)

        if rec_extension=="json": 
            try: 
                with open(rec_filename, "r") as file:
                    json_string = file.read()
                polymer = Polymer.from_json(json_string)

            except Exception as e: 
                print(e)
                sys.exit(1)

        else: 
            templates = ResidueChemTemplates.create_from_defaults()
            for fn in args.add_templates:
                templates.add_json_file(fn)
            mk_prep = MoleculePreparation()

            if args.wanted_altloc is None:
                wanted_altloc = None
            else:
                wanted_altloc = parse_cmdline_res_assign(args.wanted_altloc)
                # Ensure meaningful wanted_altloc
                for key, value in wanted_altloc.items():
                    if isinstance(value, str) and value.strip() == "":
                        msg = "Wanted atloc cannot be an empty string or a string with just space"
                        print("Command line error: " + msg, file=sys.stderr)
                        sys.exit(2)
        
            # Ensure meaningful default_altloc
            if args.default_altloc is not None and args.default_altloc.strip()=="":
                msg = "Allowed atloc cannot be an empty string or a string with just space"
                print("Command line error: " + msg, file=sys.stderr)
                sys.exit(2)


            if _has_prody: 
                prody_parser = _prody_parsers[rec_extension]
                rec_prody_mol = prody_parser(rec_filename, altloc="all")
                # covalent_builder = CovalentBuilder(rec_prody_mol, args.rec_residue)
                try:
                    polymer = Polymer.from_prody(
                        rec_prody_mol,
                        templates,
                        mk_prep,
                        wanted_altloc=wanted_altloc,
                        default_altloc=args.default_altloc,
                    )
                except PolymerCreationError as e:
                    print(e)
                    sys.exit(1)
            
            else: 
                try:
                    with open(rec_filename, "r") as f:
                        pdb_string = f.read()
                    polymer = Polymer.from_pdb_string(
                        pdb_string,
                        templates,
                        mk_prep,
                        wanted_altloc=wanted_altloc,
                        default_altloc=args.default_altloc,
                    )
                except PolymerCreationError as e:
                    print(e)
                    sys.exit(1)
                pass

        residue_string = f"{args.rec_residue}"
        chid, res_type, res_num = residue_string.split(":")

        def get_target_monomers(polymer, chid, res_type, res_num): 
            keys = polymer.monomers.keys()
            if chid: 
                keys = [key for key in keys if key.split()[0]==chid]
            if res_type: 
                keys = [key for key in keys if polymer.monomers[key].input_resname==res_type]
            if res_num: 
                keys = [key for key in keys if key.split()[1]==res_num]
            return keys

        target_monomers = get_target_monomers(polymer, chid, res_type, res_num)
                
        if not target_monomers: 
            print("ERROR: no residue found with the following specification: \n"
                  "chain[%s] residue[%s] number[%s]\n" % (chid, res_type, res_num))
            sys.exit(1)
        # endregion

    input_mol_skipped = 0
    input_mol_with_failure = (
        0  # if reactive or covalent, each mol can yield multiple PDBQT
    )
    nr_failures = 0
    is_after_first = False

    if  config["charge_atom_prop"] is not None: 
        if config["charge_model"] != "read": 
            print(f'Error: --charge_atom_prop must be used with --charge_model "read", but the current charge_model is "{config["charge_model"]}". ',
                  file=sys.stderr,)
            sys.exit(1)
    elif config["charge_model"] == "read":  
        if ext=="sdf":
            config["charge_atom_prop"] = "PartialCharge"
        elif ext=="mol2": 
            config["charge_atom_prop"] = "_TriposPartialCharge"

    preparator = MoleculePreparation.from_config(config)
    for mol in mol_supplier:
        if is_after_first and output.is_multimol == False:
            print("Processed only the first molecule of multiple molecule input.")
            print(
                "Use --multimol_prefix and/or --multimol_outdir to process all molecules in %s."
                % (input_molecule_filename)
            )
            break

        # check that molecule was successfully loaded
        is_valid = mol is not None
        input_mol_skipped += int(is_valid == False)
        if not is_valid:
            continue

        this_mol_had_failure = False

        if args.name_from_prop is not None:
            if mol.HasProp(args.name_from_prop):
                name = mol.GetProp(args.name_from_prop)
            else:
                continue  # TODO log this event
        else:
            name = mol.GetProp("_Name")
        is_after_first = True

        if is_covalent:
            connect_pattern = 1
            for cov_lig in covalent_builder.process(
                mol, args.tether_smarts, args.tether_smarts_indices
            ):
                root_atom_index = cov_lig.indices[0]
                molsetups = preparator.prepare(
                    cov_lig.mol,
                    root_atom_index=root_atom_index,
                    not_terminal_atoms=[root_atom_index],
                    rename_atoms=args.rename_atoms,
                )
                chain, res, num = cov_lig.res_id
                suffixes = output.get_suffixes(molsetups)
                for molsetup, suffix in zip(molsetups, suffixes):
                    pdbqt_string, success, error_msg = PDBQTWriterLegacy.write_string(
                        molsetup,
                        bad_charge_ok=args.bad_charge_ok,
                        add_index_map=args.add_index_map,
                    )
                    if success:
                        pdbqt_string = (
                            PDBQTWriterLegacy.adapt_pdbqt_for_autodock4_flexres(
                                pdbqt_string, res, chain, num
                            )
                        )
                        name = molsetup.name
                        output.output_filename = molsetup.name + f"_{connect_pattern}_{chain}_{res}_{num}.pdbqt"
                        output(pdbqt_string, name, (cov_lig.label, suffix))
                    else:
                        nr_failures += 1
                        this_mol_had_failure = True
                        print(error_msg, file=sys.stderr)
                    connect_pattern += 1

        else:
            try: 
                molsetups = preparator.prepare(mol, rename_atoms=args.rename_atoms)
            except Exception as error_msg: 
                nr_failures += 1
                this_mol_had_failure = True
                print(error_msg, file=sys.stderr)
                input_mol_with_failure += int(this_mol_had_failure)
                continue

            if len(molsetups) > 1:
                output.is_multimol = True
            suffixes = output.get_suffixes(molsetups)
            for molsetup, suffix in zip(molsetups, suffixes):
                pdbqt_string, success, error_msg = PDBQTWriterLegacy.write_string(
                    molsetup,
                    bad_charge_ok=args.bad_charge_ok,
                    add_index_map=args.add_index_map,
                )
                if success:
                    output(pdbqt_string, name, (suffix,))
                    if args.verbose:
                        molsetup.show()
                else:
                    nr_failures += 1
                    this_mol_had_failure = True
                    print(error_msg, file=sys.stderr)

        input_mol_with_failure += int(this_mol_had_failure)

    # Close the last tarf opened
    if output.tarf is not None:
        output.tarf.close()

    # Print final status
    print(
        "Input molecules processed: %d, skipped: %d"
        % (output.counter, input_mol_skipped)
    )
    print("PDBQT files written: %d" % (output.num_files_written))
    print("PDBQT files not written due to error: %d" % (nr_failures))
    print("Input molecules with errors: %d" % (input_mol_with_failure))
    if output.is_multimol:
        # would be None if not is_multimol
        print(output.get_duplicates_info_string())

    # Determine if exit code should be non-zero based on processing results
    if output.num_files_written == 0 and not output.redirect_stdout:
        print("No PDBQT files were written due to errors!")
        sys.exit(3)  # full failure
    elif input_mol_with_failure > 0:
        print("Some molecules encountered errors!")
        sys.exit(4)  # partial failure
    return

if __name__ == '__main__':
    sys.exit(main())
