# rdoxide: this RDKit extension module is implemented in Rust (rdkit._rdoxide).
# The shim replaces itself in sys.modules with the native submodule, so
# `import rdkit.Chem.rdSubstructLibrary` behaves exactly like importing the upstream extension.
import sys

from rdkit import _rdoxide

sys.modules[__name__] = _rdoxide.rdSubstructLibrary
