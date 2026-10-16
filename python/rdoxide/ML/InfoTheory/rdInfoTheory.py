# rdoxide: this RDKit extension module is implemented in Rust (rdoxide._rdoxide).
# The shim replaces itself in sys.modules with the native submodule, so
# `import rdoxide.ML.InfoTheory.rdInfoTheory` behaves exactly like importing the upstream extension.
import sys

from rdoxide import _rdoxide

sys.modules[__name__] = _rdoxide.rdInfoTheory
