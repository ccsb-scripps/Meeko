Meeko
=====

Parameterization of molecules for AutoDock
------------------------------------------

Meeko assigns parameters to small organic molecules, often called ligands,
and to proteins and to nucleic acids, often called receptors.
This includes assigning atom types, partial charges, setting
bonds as rotatable or fixed, and making receptor sidechains flexible.

Write input and process output
------------------------------

Meeko writes the input PDBQT files for AutoDock-Vina and AutoDock-GPU, and it
also converts the output files from docking, which are PDBQT for Vina and
DLG for AutoDock-GPU, into SDF for ligands and PDB for receptor.

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



