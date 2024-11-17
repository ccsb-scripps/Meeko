Overview of receptor preparation
================================

Polymers and monomers
---------------------

Receptors are represented as a collection of subunits where each subunit is
treated as an individual molecule. The Python class for receptors is called
``Polymer`` and the class for subunits is ``Monomer``.
The code for ligand parameterization is also used for the receptor.

During most of the v0.6.0 development, the name of the ``Polymer`` class was
``LinkedRDKitChorizo``. Many of the commit messages, as well as GitHub issues and pull
requests mention chorizos rather then polymers. The analogy is that each of the
links in the chorizo is an RDKit molecule, and that these classes often represent protein.

The input structure is matched against templates to
guarantee chemical correctness and identify problems with the input structures.
This allows the user to identify and fix problems, resulting in a molecular
model that is correct with respect to heavy atoms, protonation state,
connectivity, bond orders, and formal charges.

The matching algorithm uses the connectivity and elements, but not bond orders
or atom names. Hydrogens are optional. This makes it compatible with input
files from various sources.

Templates are matched on a per residue basis. Each residue is represented
as an instance of a PolymerResidue object, which contains:

 * an RDKit molecule that represents the actual state
 * a padded RDKit molecule containing a few atoms from the adjacent residues
 * parameters such as partial charges

The positions are set by the input, and the connectivity and formal charges
are defined by the templates. Heavy atoms must match exactly. If heavy atoms
are missing or in excess, the templates will fail to match.

Missing hydrogens are added by RDKit, but are not subjected to minimization
with a force field. Thus, their bond lengths are not super accurate.

Different states of the same residue are stored as different templates,
for example different protonation states of HIS, N-term, LYN/LYS, etc.
Residue name is primary key unless user overrides.

Currently not supported: capped residues from charmm-gui.


Templates
---------

The templates contain SMILES strings that are used to create the RDKit
molecules that constitute every residue in the processed model. In this way,
the chemistry of the processed model is fully defined by the templates,
and the only thing that is preserved from the input are the atom positions
and the connectivity between residues.

The SMILES strings contain all atoms that exist in the final model,
and none that do not exist. This also applies to hydrogens,
meaning that the SMILES are expected to have real hydrogens. Note that
real hydrogens are different from explicit hydrogens. Real hydrogens will be
represented as an actual atom in an RDKit molecule, while explicit hydrogens
are a just property of heavy atoms. In the SMILES, real hydrogens are defined
with square brackets "[H]" and explicit hydrogens without, e.g. "[nH]" to set
the number of explicit hydrogens on an aromatic nitrogen to one.

Residues that are part of a polymer, which is often all of them, will have
bonds to adjacent residues. The heavy atoms involved in the bonds will miss
a real hydrogen and have an implicit (or explicit) one instead. As an
example, consider modeling an alkyl chain as a polymer, in which the monomer
is a single carbon atom. Our template SMILES would be "[H]C[H]". The RDKit
molecule will have three atoms and the carbon will have two implicit hydrogens.
The implicit hydrogens correspond to bonds to adjacent residues in the
processed polymer. 
