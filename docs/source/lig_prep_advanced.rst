Advanced ligand preparation
===========================

Parameterization of RDKit molecules is controlled by several parameters
that can be used when creating an instance of ``MoleculePreparation``.
Many of such parameters are exposed to the command line script, 
``mk_prepare_ligand.py``. Here, we cover the most relevant options, and
how to use them both from command line and Python scripts.

Configuring MoleculePreparation in Python and command line
----------------------------------------------------------

In Python, an instance of ``MoleculePreparation`` is configured by passing
any number of optional parameters during initialization:

.. code-block:: python

   from meeko import MoleculePreparation

   mk_prep = MoleculePreparation(
       merge_these_atom_types=("H"),
       charge_model="gasteiger",
   )

Here we are passing two parameters, one that will merge atoms typed "H" into
their parent atom, and another to use Gasteiger charges. Both are the default
parameter, so we would have gotten the same configuration by passing zero
parameters: ``mk_prep = MoleculePreparation()``.

Parameters can be in a dictionary and use the ``from_config`` constructor:

.. code-block:: python
   config_dict = {"merge_these_atom_types": (), "charge_model": "gasteiger"}
   mk_prep = MoleculePreapration.form_config(config_dict)

The configuration dictionary can be written to a JSON file:

.. code-block:: python
   import json
   with open("my_config.json", "w") as f:
       json.dump(config_dict, f)

And the command line script ``mk_prepare_ligand.py`` can read the JSON file
using the ``-c`` or ``--config_file`` option:

.. code-block:: bash
   mk_prepare_ligand.py -i mol.sdf -c my_config.json

Here, ``mol.sdf`` is some file in the current working directory.
Alternatively, each parameter can be passed directly:

.. code-block:: bash
   mk_prepare_ligand.py -i mol.sdf --charge_model gasteiger --merge_these_atom_types H



Merging hydrogens (or not)
--------------------------

By default, atoms of type ``H`` are merged. Merging consists of adding
their partial charge to the parent atom, and setting the ``is_ignore``
attribute to ``True`` in the instance of ``MoleculeSetup``, which in turn
prevents the atom from being written to PDBQT. Since the atoms are merged
based on their atom type, changing the atom typing can change the set of
atoms that is merged.

Option name is ``merge_these_atom_types``. To prevent merging of any atom
types set the variable to an empty tuple in Python:

.. code-block:: python
   mk_prep = MoleculePreparation(merge_these_atom_types=())

or pass no parameters in command line
.. code-block:: bash
   mk_prepare_ligand.py -i mol.sdf --merge_these_atom_types



Modifying atom types
--------------------

Atom typing relies on


No rigid macrocycles
--------------------

Hydrated docking
----------------
