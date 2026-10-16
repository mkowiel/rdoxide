"""Makes `import rdkit...` resolve to the installed `rdoxide` package.

rdoxide does not install a package called `rdkit` (so both can be installed
together), but upstream's Python tests import `rdkit`. The upstream test
runner and the differential tools call `install()` to alias the names.
"""

import importlib
import importlib.abc
import importlib.util
import sys


class _Alias(importlib.abc.Loader):
    def __init__(self, target):
        self.target = target

    def create_module(self, spec):
        return importlib.import_module(self.target)

    def exec_module(self, module):
        pass


class _Finder(importlib.abc.MetaPathFinder):
    def find_spec(self, name, path=None, target=None):
        if name == "rdkit" or name.startswith("rdkit."):
            real = "rdoxide" + name[len("rdkit"):]
            try:
                mod = importlib.import_module(real)
            except ImportError:
                return None
            spec = importlib.util.spec_from_loader(
                name, _Alias(real), is_package=hasattr(mod, "__path__"))
            return spec
        return None


def install() -> None:
    if not any(isinstance(f, _Finder) for f in sys.meta_path):
        sys.meta_path.insert(0, _Finder())
