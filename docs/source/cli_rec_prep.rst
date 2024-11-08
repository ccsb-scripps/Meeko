mk_prepare_receptor.py
======================

A command-line script for receptor preparation from a PDB/CIF file, setting up various properties (flexible/reactive residues, box specifications, etc.), and writing useful input files (PDBQT, GPF, box configuration file for Vina) for docking. 

Basic usage
-----------

.. code-block:: bash

    mk_prepare_receptor.py -i examples/system.pdb --write_pdbqt prepared.pdbqt

This is equivalent to: 

.. code-block:: bash

    mk_prepare_receptor.py -i examples/system.pdb -o prepared -p

Read more about the syntax in :ref:`Write flags` and :ref:`Options`. 

About
-----

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

.. _templates:

Templates
~~~~~~~~~

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

Write flags
~~~~~~~~~~~

The option flags starting with ``--write`` in  ``mk_prepare_receptor`` can
be used both with an argument to specify the outpuf filename: 

.. code-block:: bash

    --write_pdbqt myenzyme.pdbqt --write_json myenzyme.json

and without the filename argument as long as a default basename is provided:

.. code-block:: bash

    --output_basename myenzyme --write_pdbqt --write_json

It is also possible to combine the two types of usage:

.. code-block:: bash

    --output_basename myenzyme --write_pdbqt --write_json --write_vina_box box_for_myenzyme.txt

in which case the specified filenames have priority over the default basename. 

Residue selection and assignment language
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Meeko uses the **chain ID** and **residue number** to identify a residue. The arguments involving selection of residues: 

.. code-block:: bash

    -d, --delete_residues <residues>

    -f, --flexres <residues>

    -r, --reactive_flexres <residues>

use the compact selection language that specify residues efficiently. The chain ID and the residue number(s) are separated by a colon (``:``) delimiter. Each residue number is combined with the most recent chain ID that precedes it, resulting in an expanded list of chain-residue pairs. 

For an input like ``A:5,7,BB:12C``, this selection language represents: ``residues (number) 5 and 7 in Chain A`` and ``residue (number) 12C in Chain BB``. 

The arguments involving assignment of residues to properties: 

.. code-block:: bash

    -n, --set_template <template>

    -b, --blunt_ends <positions>

    --wanted_altloc <location>

    -s, --reactive_name_specific <residue:atom>

use the residue selection lanaguge described above, followed by an equal sign (``=``) as the delimiter and the assigned value, which could be the name of a residue template, the atom index for the blunt end, the wanted altloc ID, or the atom name of the reactive atom. Each residue selection is comibned with the most recent assignment that precedes it, resulting in a further expanded list of residue-assignment pairs. 

For an input like ``"A:5,7=CYX,A:19A,B:17=HID``, this assignment language represents: ``residues (number) 5 in Chain A are set to (template name) CYX`` and ``residue (number) 19 A in Chain A, and residue (number) 17 in Chain B are set to (template name) HID``. 

Options
-------

Input/Output Options
~~~~~~~~~~~~~~~~~~~~

.. option:: --read_pdb <PDB_FILENAME>

   Read a PDB file using the PDB parser in RDKit.

.. option:: -i, --read_with_prody <MACROMOL_FILENAME>

   Read a PDB or mmCIF file using ProDy (if installed). ProDy can be installed from PyPI or conda-forge.

.. option:: -o, --output_basename <basename>

   Specify a default basename for output files created by `--write` options when no filename is specified.

.. option:: -p, --write_pdbqt <PDBQT_FILENAME> [*]

   Output PDBQT files with `_rigid` or `_flex` suffixes for flexible residues. Defaults to `--output_basename` if no filename is provided.

.. option:: -j, --write_json <JSON_FILENAME> [*]

   Save the receptor's parameterized configuration to JSON format. Defaults to `--output_basename` if unspecified.

.. option:: -g, --write_gpf <GPF_FILENAME> [*]

   Output an AutoGrid input file (GPF). Defaults to `--output_basename` if not specified.

.. option:: -v, --write_vina_box <VINA_BOX_FILENAME> [*]

   Generate a configuration file for Vina with grid box dimensions. Defaults to `--output_basename` if not specified.

.. option:: --write_pdb <PDB_FILENAME> [*]

   Save the prepared receptor in PDB format. Must specify the filename.

Receptor Perception Options
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. option:: -n, --set_template <template>

   Assign templates to residues, e.g., `A:5,7=CYX,B:17=HID`.

.. option:: -d, --delete_residues <residues>

   Specify residues to delete, e.g., `A:350,B:15,16,17`.

.. option:: -b, --blunt_ends <positions>

   Blunt end definitions, e.g., `A:123,200=2,A:1=0`.

.. option:: --add_templates <JSON_FILENAME> [*]

   Load additional templates from one or more JSON files.

.. option:: -a, --allow_bad_res

   (Flag) Ignore residues with missing atoms instead of raising an error.

.. option:: --default_altloc <location>

   Define a default alternate location for residues, overridden by `--wanted_altloc`.

.. option:: --wanted_altloc <location>

   Specify alternate locations for particular residues, e.g., `:5=B,B:17=A`.

.. option:: --mk_config <JSON_FILENAME>

   Specify a JSON configuration file for receptor preparation.

Grid Box Options
~~~~~~~~~~~~~~~~

.. option:: --box_size <X Y Z>

   Set the size of the grid box in Angstroms (x, y, z).

.. option:: --box_center <X Y Z>

   Define the center of the grid box in Angstroms (x, y, z).

.. option:: --box_center_off_reactive_res

   (Flag) Shift the grid box center 5 Å along the CA-CB bond from CB. Applicable only when there is one reactive flexible residue.

.. option:: --box_enveloping <FILENAME>

   Adjust the grid box to enclose atoms in the specified file. Supported formats: `.sdf`, `.mol`, `.mol2`, `.pdb`, `.pdbqt`.

.. option:: --padding <value>

   Set padding around atoms specified in `--box_enveloping` (in Å).

Flexible and/or Reactive Options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. option:: -f, --flexres <residues>

   Define flexible residues by chain ID and residue number, e.g., `-f ":42,B:23"`. 

.. option:: -r, --reactive_flexres <residues>

   Define reactive flexible residues by chain ID and residue number, e.g., `-r ":42,B:23"`. Maximum of 8 reactive residues.

.. option:: --reactive_name <residue:atom>

   Specify the reactive atom name for a residue type, e.g., `--reactive_name "TRP:NE1"`. Can be repeated for multiple assignments.

.. option:: -s, --reactive_name_specific <residue:atom>

   Specify the reactive atom for an individual residue by residue ID, e.g., `-s "A:42=NE2"`. The residue becomes reactive.

.. option:: --r_eq_12 <value>

   Set the equilibrium distance (r_eq) for reactive atoms in 1-2 interactions. Default is 1.8 Å.

.. option:: --eps_12 <value>

   Set the epsilon value for reactive atoms in 1-2 interactions. Default is 2.5.

.. option:: --r_eq_13_scaling <factor>

   Scale r_eq for 1-3 interactions across reactive atoms. Default is 0.5.

.. option:: --r_eq_14_scaling <factor>

   Scale r_eq for 1-4 interactions across reactive atoms. Default is 0.5.
