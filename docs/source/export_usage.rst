Exporting docking results
=========================

AutoDock writes results in PDBQT format, which lacks bond orders and hydrogens
bonded to carbon. When other software reads PDBQT files, it needs to infer
bond orders, which can be impossible to do accurately for some molecules.
As an example, consider the following conversion, in which OpenBabel adds an
extra double bond:

.. code-block:: bash

    obabel -:"C1C=CCO1" -o pdbqt --gen3d | obabel -i pdbqt -o smi
    [C]1=[C][C]=[C]O1

Meeko can export docking results to SDF while preserving bond orders because
when it writes the PDBQT file for a ligand (as input for docking) it writes
a SMILES string and index mapping information in remark lines.
Both AutoDock-GPU and AutoDock-Vina preserve
the remarks in the output files, and Meeko then uses them to create an accurate
RDKit molecule and write an SD file.

When docking with flexible sidechains,
Meeko can write a PDB file for the entire receptor with the modified
sidechain positions, or include the sidechains as additional fragments in
the SDF for the ligand. 

Since hydrogens bonded to carbon are excluded by default in AutoDock,
exported positions of these hydrogens are calculated by RDKit. This can
be annoying if a careful forcefield minimization is employed before
docking, as probably rigorous Hs positions will be replaced by the
RDKit geometry rules, which are empirical and much simpler than most
force fields.


Command line usage
------------------

The script reads vina output files (``.pdbqt`` extension):

.. code-block:: bash

    mk_export.py vina_results.pdbqt -s vina_results.sdf

And AutoDock-GPU output files (``.dlg`` extension):

.. code-block:: bash

    mk_export.py autodock-gpu_results.dlg -s autodock-gpu_results.sdf

It can also read files zipped with gzip (``.gz`` extension), and multiple
filenames can be passed, for example using wildcards. In this case, the output
file names will be the input filenames with an added suffix. The suffix is
``_docked`` by default but can be modified with the ``--suffix`` flag:

.. code-block:: bash

   mk_export.py docking_results/*.dlg.gz --suffix _docking_result

To write a PDB file with the entire receptor including modified sidechain
positions, the prepared receptor written by ``mk_prepare_receptor.py`` or by
the equivalent Python code is required using the ``-j/--read_json`` option,
and the ``-p/--write_pdb`` option sets the filename of the output PDB file:

.. code-block:: bash

   mk_export.py results.pdbqt -s results.sdf -p results.pdb -j rec.json

To write flexible sidechains as additional fragments to the ligand SDF,
use the ``-k/--keep_flexres_sdf`` flag.


In Python
---------

This example converts an output file from AutoDock-GPU into a list of
RDKit molecules, each with multiple conformers corresponding to a different
docked pose:

.. code-block:: python
   
   from meeko import PDBQTMolecule
   from meeko import RDKitMolCreate
   
   fn = "autodock-gpu_results.dlg"
   pdbqt_mol = PDBQTMolecule.from_file(fn, is_dlg=True, skip_typing=True)
   rdkitmol_list = RDKitMolCreate.from_pdbqt_mol(
       pdbqt_mol,
       only_cluster_leads=True,
       keep_flexres=False,
   )

The length of `rdkitmol_list` is one if there are no sidechains and only one
ligand was docked.
If multiple ligands and/or sidechains are docked simultaneously, each will be
an individual RDKit molecule in `rdkitmol_list`.
In this example, ``keep_flexres=False`` would exclude sidechains from the
RDKit molecules even if such sidechains existed.
Note that docking multiple
ligands simultaneously is only available in Vina, and it differs from docking
multiple ligands one after the other. Each failed creation of an RDKit molecule
for a ligand or sidechain results in a `None` in `rdkitmol_list`.

For Vina's output PDBQT files, omit `is_dlg=True`:

.. code-block:: python

   pdbqt_mol = PDBQTMolecule.from_file("vina_results.pdbqt", skip_typing=True)

When using Vina from Python, the output string can be passed directly.
See [the docs](https://autodock-vina.readthedocs.io/en/latest/docking_python.html)
for context on the `v` object.

.. code-block:: python

    vina_output_string = v.poses()
    pdbqt_mol = PDBQTMolecule(vina_output_string, is_dlg=True, skip_typing=True)
    rdkitmol_list = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
