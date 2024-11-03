.. _tutorial2:

=========================
Reactive Docking
=========================

This is a reactive docking example that uses the AutoDock-GPU executable to generate the near-attack conformation of a small molecule and a protein receptor. 

Follow the instructions to set up the environment and run this example on your own device (Linux, MacOS or WSL). To run this example in a Colab notebook, see :ref:`colab_examples`. 

.. contents::
   :local:
   :depth: 2

Introduction
------------

The reactive docking example is based on reactive docking method that has been developed for high-throughput virtual screenings of reactive species. This method is currently only implemented in AutoDock-GPU. In this example, a small molecule substrate (Adenosine monophosphate, PDB token AMP) is targeting at the catalytic histidine residue of a hollow protein structure of bacteria RNA 3' cyclase (PDB token 3KGD) to generate the near-attack conformation for the formation of the phosphoamide bond. A docked pose that closely resembles the original position of the ligand is expected among the top-ranked poses. 

This tutorial is intended to showcase the Meeko usage in the preparation of receptor and ligand for reactive docking. 

Prerequisites and Environment Setup
-----------------------------------

1. **Create a new virtual environment (recommended)**

.. code-block:: bash

   micromamba create -c conda-forge -n meeko_tutorial python=3.10 -y
   micromamba activate meeko_tutorial         

In this tutorial, we will use ``micromamba`` as the example package manager. Visit `this official guide  <https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html>`_ for a quick install and setup of micromamba. There are many equivalent ways to manage Python packages, such as ``conda`` and ``mamba``. You can easily adapt the commands to your preferred tool, as the syntax is largely compatible across these package managers. 

2. Install the required Python packages through `conda-forge`

.. code-block:: bash

   micromamba install -c conda-forge cctbx-base numpy scipy rdkit gemmi -y

Expose ``reduce2.py`` to system ``PATH``

3. Install the additional packages and data from GitHub repositories

- (Python package) Meeko 

.. code-block:: bash

   git clone --single-branch --branch develop https://github.com/forlilab/Meeko.git
   cd Meeko; pip install --use-pep517 -e .; cd ..

- (Python package) scrubber 

.. code-block:: bash

   git clone --single-branch --branch develop https://github.com/forlilab/scrubber.git
   cd scrubber; pip install --use-pep517 -e .; cd ..

- (Python package) ProDy 

.. code-block:: bash

   pip install prody

- (Required data package for reduce2) Phenix Project geostd (restraint) Library 

.. code-block:: bash

   git clone https://github.com/phenix-project/geostd.git


Ligand Peparation
-----------------

.. code-block:: bash

    ligand_smiles="c1nc(c2c(n1)n(cn2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)([O-])[O-])O)O)N"
    scrub.py $ligand_smiles -o AMP.sdf --ph 6.5 --skip_tautomer --skip_acidbase

.. code-block:: bash

    reactive_smarts="COP(=O)([O-])[O-]"
    reactive_smarts_idx=3
    mk_prepare_ligand.py -i AMP.sdf -o AMP.pdbqt \
    --reactive_smarts $reactive_smarts \
    --reactive_smarts_idx $reactive_smarts_idx

Receptor Peparation
-----------------

.. code-block:: bash

    pdb_token="3kgd"
    curl "http://files.rcsb.org/view/${pdb_token}.pdb" -o "${pdb_token}.pdb"

.. code-block:: bash

    python - <<EOF
    from prody import parsePDB, writePDB

    pdb_token = "3kgd"
    atoms_from_pdb = parsePDB(pdb_token)
    receptor_selection = "chain A and not water and not hetero and not resname AMP"
    receptor_atoms = atoms_from_pdb.select(receptor_selection)
    prody_receptorPDB = f"{pdb_token}_receptor_atoms.pdb"
    writePDB(prody_receptorPDB, receptor_atoms)
    EOF

    # Add CRYST1 card (temporarily required for reduce2)
    cat <(grep "CRYST1" "${pdb_token}.pdb") "${pdb_token}_receptor_atoms.pdb" > "${pdb_token}_receptor.pdb"

.. code-block:: bash

   reduce2="$(python -c "import site; print(site.getsitepackages()[0])")/mmtbx/command_line/reduce2.py"
   chmod +x $reduce2
   geostd="$(realpath geostd)"
   export MMTBX_CCP4_MONOMER_LIB=$geostd
   reduce_opts="approach=add add_flip_movers=True"
   python $reduce2 "${pdb_token}_receptor.pdb" $reduce_opts

.. code-block:: bash

    python - <<EOF
    from prody import parsePDB, writePDB, calcCenter

    pdb_token = "3kgd"
    atoms_from_pdb = parsePDB(pdb_token)
    ligand_selection = "chain A and resname AMP"
    ligand_atoms = atoms_from_pdb.select(ligand_selection)
    prody_ligandPDB = "LIG.pdb"
    writePDB(prody_ligandPDB, ligand_atoms)
    EOF

    reactive_name_specific="A:309=NE2"
    mk_prepare_receptor.py -i "${pdb_token}_receptorH.pdb" -o "${pdb_token}_receptorH" -p -g \
    --default_altloc A --reactive_name_specific $reactive_name_specific \
    --box_enveloping "LIG.pdb" --padding 8.0 

.. code-block:: bash

    @> 2510 atoms and 1 coordinate set(s) were parsed in 0.01s.

    Flexible residues:
    chain resnum is_reactive reactive_atom
        A    309        True           NE2
    reactive_flexres={'A:309'}

    For reactive docking, pass the configuration file to AutoDock-GPU:
        autodock_gpu -C 1 --import_dpf 3kgd_receptorH.reactive_config --flexres 3kgd_receptorH_flex.pdbqt -L <ligand_filename>


    Files written:
        3kgd_receptorH_flex.pdbqt <-- flexible receptor input file
        3kgd_receptorH_rigid.pdbqt <-- static (i.e., rigid) receptor input file
        boron-silicon-atom_par.dat <-- atomic parameters for B and Si (for autogrid)
        3kgd_receptorH_rigid.gpf <-- autogrid input file
            3kgd_receptorH.box.pdb <-- PDB file to visualize the grid box
    3kgd_receptorH.reactive_config <-- reactive parameters for AutoDock-GPU