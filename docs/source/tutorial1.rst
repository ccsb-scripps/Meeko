.. _tutorial1:

Basic Docking 
-------------------------------------

This tutorial provides practice examples and a step-by-step guide for the four basic procedures with Meeko for molecular docking and virtual screening with `AutoDock Vina <https://github.com/ccsb-scripps/AutoDock-Vina>`_ and `AutoDock-GPU <https://github.com/ccsb-scripps/AutoDock-GPU>`_: 

- Ligand Preparation 
- Receptor Preparation 
- Molecular Docking (Single Ligand) 
- Virtual Screening (Batch Docking) 

It is based on, but not a full version of the tutorial materials in `Forlilab tutorials <https://github.com/forlilab/tutorials>`_. 

.. contents::
   :local:
   :depth: 2

Prerequisites and Environment Setup
===================================

Create a new virtual environment (recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   micromamba create -c conda-forge -n meeko_tutorial_py39 python=3.9 -y
   micromamba activate meeko_tutorial_py39         

In this tutorial, we will use ``micromamba`` as the example package manager. Visit `this official guide  <https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html>`_ for a quick install and setup of micromamba. There are many equivalent ways to manage Python packages, such as ``conda`` and ``mamba``. You can easily adapt the commands to your preferred tool, as the syntax is largely compatible across these package managers. 

Install the required Python packages through ``conda-forge``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   micromamba install -c conda-forge numpy scipy rdkit gemmi vina -y

Install the additional packages and data from GitHub repositories
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- (Python package) Meeko 

.. code-block:: bash

   git clone --single-branch --branch develop https://github.com/forlilab/Meeko.git
   cd Meeko; pip install --use-pep517 -e .; cd ..

- (Python package) scrubber 

.. code-block:: bash

   git clone --single-branch --branch develop https://github.com/forlilab/scrubber.git
   cd scrubber; pip install --use-pep517 -e .; cd ..

- (Python package) Ringtail

.. code-block:: bash

   git clone --single-branch --branch develop https://github.com/forlilab/Ringtail.git
   cd Ringtail; pip install --use-pep517 -e .; cd ..

- (Example files for this tutorial) Forlilab Tutorials

.. code-block:: bash

   git clone https://github.com/forlilab/tutorials.git

Ligand Preparation
==================

Ligand Preparation is the process that generates ligand input files for docking calculation and virtual screening. At present, AutoDock Vina and AutoDock-GPU need the ligand input files in the PDBQT format. In this example, we will use ``mk_prepare_ligand.py``, a command-line script in Meeko, to prepare such ligand PDBQT files. 

Prepare a Single Ligand from a Smiles string
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`Imatinib <https://pubchem.ncbi.nlm.nih.gov/compound/Imatinib>_` is a small-molecule drug. You can find the SMILES string for Imatinib from various reliable chemical databases and resources, including but not limited to `PubChem <https://pubchem.ncbi.nlm.nih.gov/>_` and `DrugBank <https://go.drugbank.com/>`_. 

``scrub.py`` is a command-line script in Scrubber that generates 3D conformers of protomers and tautomers for given small molecules at a specified (range of) pH. Given a pH range of 5 to 9, the output protomers will include those which make up no less than 1% of the total population at pH = 7. Based on the reference pKa values, the amine nitrogens and the pyridine nitrogen will be considered for acid/base enumeration. With the ``meeko_tutorial_py39`` micromamba environment active, run ``scrub.py`` to generate 3D conformers of Imatinib from the SMILES string. 

.. code-block:: bash

    smiles_string="CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5"
    scrub.py $smiles_string -o imatinib.sdf --skip_tautomers --ph_low 5 --ph_high 9

The output file ``imatinib.sdf`` will contain two protomers of Imatinib, one with a neutral pyridine group and the other with a (+1) pyridinium group. All of the aliphatic amininium nitrogens will be protonated. 

.. code-block:: bash

    Scrub completed.
    Summary of what happened:
    Input molecules supplied: 1
    mols processed: 1, skipped by rdkit: 0, failed: 0
    nr isomers (tautomers and acid/base conjugates): 2 (avg. 2.000 per mol)
    nr conformers:  2 (avg. 1.000 per isomer, 2.000 per mol)

In case there are multiple molecules in the SDF file, ``mk_prepare_ligand.py`` needs to know the prefix of filenames (by ``--multimol_prefix``) or alternatively where to output (by ``--multimol_outdir``) the multiple PDBQT files. Here, we will give the PDBQT files a prefix ``imatinib_protomer`` in the names. The output PDBQT files will be ``imatinib_protomer-1.pdbqt`` and ``imatinib_protomer-2.pdbqt``. 

.. code-block:: bash

    mk_prepare_ligand.py -i imatinib.sdf --multimol_prefix imatinib_protomer


Receptor Preparation
====================

Receptor Preparation is 



Molecular Docking (Single Ligand)
=================================

Docking with AutoDock-Vina
~~~~~~~~~~~~~~~~~~~~~~~~~~

Docking with AutoDock-GPU
~~~~~~~~~~~~~~~~~~~~~~~~~

Processing the Docking Results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Virtual Screening (Batch Docking)
=================================

Batch Docking with AutoDock-Vina
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Batch Docking with AutoDock-GPU
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Processing the Screening Results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


What's Next?
------------

Now that you've completed this tutorial, you're ready to move on to :ref:`tutorial2` and :ref:`tutorial3` where we dive deeper into more advanced docking methods: reactive docking and tethered docking.
