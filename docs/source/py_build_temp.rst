Building residue templates in Python
===================================

.. code-block:: python

    from meeko.chemtempgen import *
    from rdkit import Chem
    from rdkit.Chem import Draw
    from rdkit import RDLogger
    from PIL import Image
    import io
    import copy
    import logging
    import sys
    
    rdkit_logger = RDLogger.logger()
    rdkit_logger.setLevel(RDLogger.CRITICAL)

    # Create a chemical component from the definition CIF file
    basename = "CRO" # must be a recognized name in CCD (chemical component library)
    CRO_from_cif = ChemicalComponent.from_cif(fetch_from_pdb(basename), basename) # needs network connection to download the definition CIF file

    def draw_cc_mol(cc_mol: Chem.Mol): 
        # Label atoms by atom name
        for atom in cc_mol.GetAtoms():
            atom.SetProp("atomNote", atom.GetProp("atom_id"))

        # Draw the molecule
        drawer = Draw.MolDraw2DCairo(600, 600) 
        drawer.DrawMolecule(cc_mol)
        drawer.FinishDrawing()

        # Get the image as PNG
        png_data = drawer.GetDrawingText()
        img = Image.open(io.BytesIO(png_data))
        img.show()

    draw_cc_mol(CRO_from_cif.rdkit_mol)

.. code-block:: python

    cc = (
        cc
        .make_pretty_smiles()
        .make_link_labels_from_patterns(pattern_to_label_mapping = AA_recipe.pattern_to_label_mapping_standard)
        )
    cc.ResidueTemplate_check()
    export_chem_templates_to_json([cc])

.. code-block:: bash

    ******************** New Template Built ********************
    {
        "ambiguous": {
            "CRO": ["CRO"]
        },
        "residue_templates": {
            "CRO": {
                "smiles": "[H]NC([H])(C1=NC(=C([H])C2=C([H])C([H])=C(O[H])C([H])=C2[H])C(=O)N1C([H])([H])C=O)C([H])(O[H])C([H])([H])[H]",
                "atom_name": ["H", "N1", "CA1", "HA1", "C1", "N2", "CA2", "CB2", "HB2", "CG2", "CD1", "HD1", "CE1", "HE1", "CZ", "OH", "HOH", "CE2", "HE2", "CD2", "HD2", "C2", "O2", "N3", "CA3", "HA31", "HA32", "C3", "O3", "CB1", "HB1", "OG1", "HOG1", "CG1", "HG11", "HG12", "HG13"],
                "link_labels": {"1": "N-term", "27": "C-term"}
            }
        }
    }
    ************************************************************

