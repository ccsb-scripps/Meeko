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
integration with external software. Command line scripts are also available.

AutoDock ecosystem
------------------

To run a docking, more packages are required besides Meeko:

 * AutoDock-Vina
 * AutoDock-GPU
 * Ringtail



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
   cli_lig_prep
   In Python <py_lig_prep>

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Receptor preparation

   cli_rec_prep

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Exporting results

   cli_export_result
   In Python <py_export_result>

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Building residue templates

   In Python <py_build_temp>



