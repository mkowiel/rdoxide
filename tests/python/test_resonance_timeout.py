"""Upstream `testGitHub2597` cancels a resonance enumeration with a 3 s timeout callback, which
rdoxide (under 1 s for the same enumeration) finishes first. The same check with a timeout it can
miss."""

from datetime import datetime, timedelta

from rdoxide import Chem


class ExceedTimeout(Chem.ResonanceMolSupplierCallback):
    def __init__(self, seconds):
        super().__init__()
        self.start_time = None
        self.timeout = timedelta(seconds=seconds)

    def __call__(self):
        if self.start_time is None:
            self.start_time = datetime.now()
        return datetime.now() - self.start_time < self.timeout


def test_callback_cancels_enumeration():
    mol = Chem.MolFromSmiles(
        "ClC1=NC(NC2=CC=CC3=C2C(=O)C2=CC=CC=C2C3=O)=NC(NC2=CC=CC3=C2C(=O)C2=CC=CC=C2C3=O)=N1"
    )
    flags = Chem.UNCONSTRAINED_CATIONS | Chem.UNCONSTRAINED_ANIONS | Chem.KEKULE_ALL
    full = Chem.ResonanceMolSupplier(mol, flags)
    full.Enumerate()
    assert not full.WasCanceled()

    supplier = Chem.ResonanceMolSupplier(mol, flags)
    supplier.SetProgressCallback(ExceedTimeout(0.05))
    supplier.Enumerate()
    assert supplier.WasCanceled()
    assert 0 < len(supplier) < len(full)
