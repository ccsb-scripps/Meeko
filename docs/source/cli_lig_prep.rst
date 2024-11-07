mk_prepare_ligand.py
====================

A command-line script for ligand preparation to generate the ligand PDBQT file. Currently support SD files (.sdf), Mol2 files (.mol2) and Mol files (.mol), but SDF is strongly preferred. 

Basic usage
-----------

.. code-block:: bash

    # write a single ligand PDBQT file from a single-molecule input file
    mk_prepare_ligand.py -i molecule.sdf -o molecule.pdbqt

    # prepare ligand PDBQT files in batch from a multiple-molecule input file
    mk_prepare_ligand.py -i multi_mol.sdf --multimol_outdir folder_for_pdbqt_files

Options
-------

Input/Output Options
~~~~~~~~~~~~~~~~~~~~

.. option:: -i, --mol <input_molecule_filename>

   The input molecule file, in formats such as MOL2, SDF, etc. This option is required.

.. option:: --name_from_prop <property_name>

   Set the molecule name using a specified RDKit or SDF property.

.. option:: -o, --out <output_pdbqt_filename>

   Specify the output PDBQT filename. Only compatible with single-molecule input.

.. option:: --multimol_outdir <output_directory>

   Specify the directory to write PDBQT output files for multi-molecule inputs. Incompatible with `-o/--out` and `-`/`--`.

.. option:: --multimol_prefix <prefix>

   Replace the internal molecule name in multi-molecule input with the specified prefix. Incompatible with `-o/--out` and `-`/`--`.

.. option:: -z, --multimol_targz

   (Flag) Compress output files into a `.tar.gz` archive.

.. option:: --multimol_targz_size <size>

   Define the number of PDBQT files per `.tar.gz` archive. Default is 10000.

.. option:: -, --

   (Flag) Redirect output to standard output (STDOUT) instead of writing a file. Ignored if `-o/--out` is specified. Only compatible with single-molecule input.

Molecule Preparation Options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. option:: -c, --config_file <config_file>

   Configure `MoleculePreparation` from a JSON file. Command-line arguments will override settings in the file.

.. option:: --rigid_macrocycles

   (Flag) Keep macrocycles rigid in their input conformation.

.. option:: --macrocycle_allow_A

   (Flag) Allow bond break with atom type A, which will be retyped as carbon (C).

.. option:: --keep_chorded_rings

   (Flag) Retain all rings from exhaustive ring perception.

.. option:: --keep_equivalent_rings

   (Flag) Retain rings with equivalent sizes and neighboring atoms.

.. option:: --min_ring_size <size>

   Define the minimum number of atoms required in a ring for it to be considered for opening.

.. option:: -w, --hydrate

   (Flag) Add water molecules to the structure for hydrated docking.

.. option:: --merge_these_atom_types <types> [*]

   Specify a list of atom types to merge. The default is `"H"`.

.. option:: -r, --rigidify_bonds_smarts <SMARTS>

   Provide SMARTS patterns to rigidify specific bonds in the molecule.

.. option:: -b, --rigidify_bonds_indices <i j>

   Specify the indices of two atoms that define a bond in the SMARTS pattern (starting from 1).

.. option:: -a, --flexible_amides

   (Flag) Allow amide bonds to rotate, making them non-planar (not recommended).

.. option:: -p, --atom_type_smarts <JSON_FILENAME>

   Specify SMARTS-based atom typing in JSON format.

.. option:: -v, --verbose

   (Flag) Print detailed information about the molecule setup process.
