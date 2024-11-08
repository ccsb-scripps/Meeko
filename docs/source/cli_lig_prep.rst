mk_prepare_ligand.py
====================

A command-line script for ligand preparation that generates the ligand PDBQT file(s). Currently supports SD files (.sdf), Mol2 files (.mol2) and Mol files (.mol), but SDF is strongly preferred as input files. 

Basic usage
-----------

.. code-block:: bash

    # write a single ligand PDBQT file from a single-molecule input file
    mk_prepare_ligand.py -i molecule.sdf -o molecule.pdbqt

    # prepare ligand PDBQT files in batch from a multiple-molecule input file
    mk_prepare_ligand.py -i multi_mol.sdf --multimol_outdir folder_for_pdbqt_files

Example: Prepare ligand from a multi-molecule SDF, output PDBQT in (multiple) tar.gz of a certain size
~~~~~~~~

.. code-block:: bash

   lig_sdf="Meeko/example/tutorial1/mols.sdf"
   mk_prepare_ligand.py -i $lig_sdf --multimol_prefix mol_batch -z --multimol_targz_size 100

Example: Prepare ligand with the macrocycle-typing option turned off
~~~~~~~~
.. code-block:: bash

   lig_sdf="Meeko/example/cli_lig_prep/Rifamycin_PubChem.sdf"
   mk_prepare_ligand.py -i $lig_sdf -o Rifamycin_rigidmacro.pdbqt --rigid_macrocycles

   # current default allows flexible macrocycles
   mk_prepare_ligand.py -i $lig_sdf -o Rifamycin_flexmacro.pdbqt 

Example: Prepare ligand with the espaloma charge model 
~~~~~~~~

Python module espaloma and its dependencies are required. Visit `the official documentation for installation guide <https://espaloma.wangyq.net/install.html>`_. 

.. code-block:: bash

   mk_prepare_ligand.py -i $lig_sdf -o AMP_espaloma.pdbqt --charge_model espaloma

   # current default charge model is gasteiger
   mk_prepare_ligand.py -i $lig_sdf -o AMP_gasteiger.pdbqt 


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

   Define the number of PDBQT files per `.tar.gz` archive. Default is 10000. Only effective when used with `-z, --multimol_targz`. 

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

.. option:: -aa, --add_atom_types <JSON>

   Specify additional atom types to assign in JSON format, with SMARTS patterns and atom type names.

.. option:: --double_bond_penalty <penalty>

   Set a penalty value; values greater than 100 prevent breaking double bonds.

.. option:: --charge_model <model>

   Choose the charge model: `gasteiger`, `espaloma`, or `zero`. Default is `gasteiger`; `zero` sets all charges to zero.

.. option:: --bad_charge_ok

   (Flag) Allow NaN and Inf charges in the PDBQT output.

.. option:: --add_index_map

   (Flag) Include a map of atom indices from the input to the PDBQT file.

.. option:: --remove_smiles

   (Flag) Exclude SMILES from being written as a remark in the PDBQT output.

Reactive Docking Options
~~~~~~~~~~~~~~~~~~~~~~~~

.. option:: --reactive_smarts <SMARTS>

   Provide a SMARTS pattern for defining the reactive group.

.. option:: --reactive_smarts_idx <index>

   Specify the 1-based index of the reactive atom within the SMARTS pattern provided by `--reactive_smarts`.

Covalent Docking (Tethered) Options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. option:: --receptor <filename>

   Specify the receptor file. Supported formats depend on ProDy availability, such as `.pdb` and `.mmcif`.

.. option:: --rec_residue <residue>

   Specify the residue in the receptor for attachment, e.g., `A:LYS:204`.

.. option:: --tether_smarts <SMARTS>

   Provide a SMARTS pattern defining the ligand atoms used for attachment to the receptor.

.. option:: --tether_smarts_indices <IDX IDX>

   Specify the 1-based indices of the two atoms in the SMARTS pattern that will be attached (default: `1 2`).