Overview of ligand preparation
==============================

Meeko provides both command line scripts and a Python API. Here, we briefly
describe how data flows under the hood, and some of the key methods and data
structures. While this information is primarily for those using Meeko from
Python, it provides context for users of the command line script
``mk_prepare_ligand.py``. When using this script, the input are SD files
(.sdf extension), and the output are PDBQT files.

Use of RDKit
------------

Meeko takes as input an RDKit molecule that has 3D positions and
all hydrogens as real atoms, and creates an object called MoleculeSetup
that stores all the parameters, such as atom types, partial charges, or
rotatable bonds. Adding hydrogens and 3D positions is not performed by Meeko.

Parameterization using SMARTS patterns
--------------------------------------

Many of the parameters are assigned using SMARTS strings, which are a compact
and versatile query language for identigying chemical substructures. This makes
it easier to define custom atom types and other parameters. The SMARTS
that define the AutoDock parameters are included in Meeko and set as the
default, This goal of this feature is to facilitate the implementation of
custom docking workflows.

The MoleculePreparation class stores all the configuration required to
parameterize a ligand, and contains the methods that take an RDKit molecule
as input, and return a parameterized MoleculeSetup. The MoleculePreparation
offers several option to control how molecules are parameterized, and many of
which are exposed as command line options in ``mk_prepare_ligand.py``.

Preparing the input for docking
-------------------------------

AutoDock-GPU and Vina use the PDBQT format for input molecules.
PDBQT files, or strings in Python, are produced from a parameterized instance
of MoleculeSetup. The methods for that live in the PDBQTWriterLegacy class.
PDBQT strings can be passed directly to Vina using its Python API.
