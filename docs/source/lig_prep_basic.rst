Basic ligand preparation
========================

Command line script
-------------------

Upon installing Meeko, ``mk_prepare_ligand.py`` is available in the terminal
and running it should print a help message. This script writes one or more
PDBQT files for AutoDock-GPU or AutoDock-Vina. It reads SD files, often referred
to as SDF, and these can contain multiple molecules. The input molecules must
have all hydrogens as real atoms, and 3D positions. Mol2 files are supported,
but the SDF is strongly preferred.

To print the help message pass the ``-h`` option:

.. code-block:: bash

   mk_prepare_ligand.py -h


Writing a single PDBQT file:

.. code-block:: bash

    mk_prepare_ligand.py -i molecule.sdf -o molecule.pdbqt

If the ``-o`` option is omitted, the output filename will be the same as the
input but with ``.pdbqt`` extension. Option ``-`` prints the PDBQT
string to standard output instead of writing to a file.

Writing multiple PDBQT files from an SD file with multiple molecules:

.. code-block:: bash

    mk_prepare_ligand.py -i multi_mol.sdf --multimol_outdir folder_for_pdbqt_files

Optionally, passing ``-z`` or ``--multimol_targz`` will compress the output as
``.tar.gz`` files with 10000 PDBQT files each.

Python API
----------

The following script reads a file named ``molecules.sdf`` using RDKit,
configures an instance of ``MoleculePreparation`` with default options,
and creates PDBQT strings for all molecules in the input file:

.. code-block:: python

    from meeko import MoleculePreparation
    from meeko import PDBQTWriterLegacy
    from rdkit import Chem
    
    input_filename = "molecules.sdf"
    
    # iterate over molecules in SD file
    for mol in Chem.SDMolSupplier(input_filename, removeHs=False):
        mk_prep = MoleculePreparation()
        molsetup_list = mk_prep(mol)

        molsetup = molsetup_list[0]

        pdbqt_string = PDBQTWriterLegacy.write_string(molsetup)


The ``pdbqt_string`` can be written to a file for docking with AutoDock-GPU or
AutoDock-Vina, or passed directly to Vina within Python using Vina's Python API,
and avoiding writing PDBQT files to the filesystem.

Note that calling ``mk_prep`` returns a list of molecule setups.
As of v0.6.0, this list contains only one element unless ``mk_prep`` is
configured for reactive docking, which is not the case in this example. This is
why we are considering the first (and only) molecule setup in the list.
