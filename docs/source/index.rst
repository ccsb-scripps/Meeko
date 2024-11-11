Meeko: interface for AutoDock
=============================

Parameterization of molecules
-----------------------------

Meeko assigns parameters to small organic molecules, referred to as ligands,
and to proteins and nucleic acids, referred to as receptors.
Parameterization includes assigning atom types, partial charges, setting
bonds as rotatable, and making receptor sidechains flexible.

Prepare input and convert output
--------------------------------

Meeko writes the input PDBQT files (or strings in Python) for AutoDock-Vina
and AutoDock-GPU, and exports docking results in SDF format for ligands and
in PDB format for receptor.

Python API
----------

Meeko is written in Python and exposes functions and classes that operate on
RDKit molecules for the ligands, leveraging RDKit's popularity to facilitate
integration with external software that also interfaces with RDKit.

Command line scripts
--------------------

There are two scripts to write the input files for docking,
``mk_prepare_ligand.py`` and ``mk_prepare_receptor.py``, and one script to
convert docking output files  ``mk_export.py``.

Running a docking
-----------------

To run a docking, more packages are required besides Meeko:

 * AutoDock-Vina
 * AutoDock-GPU
 * Ringtail

Check the tutorials page to learn about using meeko with these other packages
to run molecular docking and virtual screening.



.. toctree::
   :maxdepth: 3
   :hidden:
   :caption: Getting started

   installation
   colab_examples
   tutorials

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Ligand preparation

   Overview <lig_overview>
   Basic usage <lig_prep_basic>
   Advanced usage <lig_prep_advanced>
   mk_prepare_ligand.py options <lig_cli_options>

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Receptor preparation

   Overview <rec_overview>
   Info on templates  <py_build_temp>
   Command line usage <cli_rec_prep>
   mk_prepare_receptor.py options <rec_cli_options>

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Exporting results

   Usage <export_usage>
   mk_export.py options <export_cli_options>
