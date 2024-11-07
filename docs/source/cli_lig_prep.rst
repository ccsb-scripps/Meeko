mk_prepare_ligand.py
====================

Command line tool to prepare small organic molecules.

About
-----

Write PDBQT files
~~~~~~~~~~~~~~~~~

AutoDock-GPU and Vina read molecules in the PDBQT format. These can be prepared
by Meeko from SD files, or from Mol2 files, but SDF is strongly preferred.

.. code-block:: bash

    mk_prepare_ligand.py -i molecule.sdf -o molecule.pdbqt
    mk_prepare_ligand.py -i multi_mol.sdf --multimol_outdir folder_for_pdbqt_files

Usage
-----

.. code-block:: bash

   python mk_prepare_ligand.py [OPTIONS]

Positional Argument
~~~~~~~~~~~~~~~~~~~

.. option:: -i, --mol <input_molecule_filename>

   The input molecule file, in formats such as MOL2, SDF, etc.

Options
~~~~~~~

.. option:: -c, --config_file <config_file>

   Configure `MoleculePreparation` from a JSON file. Command-line arguments will override settings in the file.

   **Example:**

   .. code-block:: bash

      python mk_prepare_ligand.py -c config.json -i ligand.sdf

.. option:: -v, --verbose

   (Flag) Print detailed information about the molecule setup process.

.. option:: --name_from_prop <property_name>

   Set the molecule name using a specified RDKit or SDF property.

   **Example:**

   .. code-block:: bash

      python mk_prepare_ligand.py -i ligand.sdf --name_from_prop compound_name

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
----------------------------

.. option:: --rigid_macrocycles

   (Flag) Keep macrocycles rigid in their input conformation.

.. option:: --macrocycle_allow_A

   (Flag) Allow bond break with atom type A, which will be retyped as carbon (C).

.. option:: --keep_chorded_rings

   (Flag) Retain all rings from exhaustive ring perception.

