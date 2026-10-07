#
#  Copyright (C) 2000-2026 Greg Landrum and other RDKit contributors
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
# rdoxide: PROVISIONAL version of upstream rdkit/Chem/__init__.py.
# It only imports the native modules ported so far; it is replaced by the
# verbatim upstream file once rdmolfiles/rdmolops/... are available
# (see docs/STATUS.md).
""" A module for molecules and stuff

 see Chem/index.html in the doc tree for documentation

"""
from rdkit import rdBase
from rdkit.Chem import rdchem
from rdkit.Chem.rdchem import *
from rdkit.Chem.rdmolfiles import *
