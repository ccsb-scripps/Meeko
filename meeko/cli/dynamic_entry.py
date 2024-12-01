#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# wrapper functions for each console script

def prepare_receptor():
    import importlib
    main = importlib.import_module('meeko.cli.mk_prepare_receptor').main
    main()

def prepare_ligand():
    import importlib
    main = importlib.import_module('meeko.cli.mk_prepare_ligand').main
    main()

def export():
    import importlib
    main = importlib.import_module('meeko.cli.mk_export').main
    main()