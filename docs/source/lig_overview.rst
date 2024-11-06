Overview of ligand preparation
==============================

Meeko takes as input an RDKit molecule that has 3D positions and
all hydrogens as real atoms, and creates an object called MoleculeSetup
that stores all the parameters, such as atom types, partial charges, or
rotatable bonds. Adding hydrogens and 3D positions is not performed by Meeko.

The MoleculePreparation class stores all the configuration required to
parameterize a ligand, and containes the methods that take an RDKit molecule
as input, and return a parameterized MoleculeSetup. The MoleculePreparation
offers several option to control how molecules are parameterized, and many of
which are exposed as command line options in ``mk_prepare_ligand.py``.

PDBQT files, or strings in Python, are produced using a parameterized instance
of MoleculeSetup. The methods for that live in the PDBQTWriterLegacy class.
PDBQT strings can be passed directly to Vina using its Python API.
