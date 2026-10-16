"""Picked up through PYTHONPATH by every Python process the upstream test
runner starts (including the subprocesses that tests spawn): `import rdkit`
then resolves to rdoxide."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import rdkit_alias  # noqa: E402

rdkit_alias.install()
