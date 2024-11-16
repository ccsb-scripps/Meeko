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


Write flags
~~~~~~~~~~~

The option flags starting with ``--write`` in  ``mk_prepare_receptor`` can
be used both with an argument to specify the output filename: 

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

use the residue selection lanaguge described above, followed by an equal sign (``=``) as the delimiter and the assigned value, which could be the name of a residue template, the atom index for the blunt end, the wanted altloc ID, or the atom name of the reactive atom. Each residue selection is combined with the most recent assignment that precedes it, resulting in a further expanded list of residue-assignment pairs. 

For an input like ``"A:5,7=CYX,A:19A,B:17=HID``, this assignment language represents: ``residues (number) 5 in Chain A are set to (template name) CYX`` and ``residue (number) 19 A in Chain A, and residue (number) 17 in Chain B are set to (template name) HID``. 

At present, Meeko always requires **chain ID** and **residue number** to identify a residue (cofactor, ligand, or ion). The only exception occurs when the chain ID in the input file is empty; in this case, the chain ID should be omitted from the selection language, i.e. ``"A:19,:17"`` represents: ``residue (number) 19 in Chain A`` and ``residue (number) 17 with an empty chain ID``. 

