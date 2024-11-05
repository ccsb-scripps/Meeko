.. _tutorial1:

Basic Docking 
-------------------------------------

This tutorial provides a step-by-step guide for the five basic procedures with Meeko for molecular docking and virtual screening with `AutoDock Vina <https://github.com/ccsb-scripps/AutoDock-Vina>`_ and `AutoDock-GPU <https://github.com/ccsb-scripps/AutoDock-GPU>`_: 

- Ligand Preparation
- Receptor Preparation
- Molecular Docking (Single Ligand) and Analysis
- Virtual Screening (Batch Docking) and Analysis

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

Ligand Preparation
==================

Ligand Preparation is 


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
