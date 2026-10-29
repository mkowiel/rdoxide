"""PNG drawing against the official wheel's MolDraw2DCairo (tests/data/png_ref, tools/make_png_reference.py).

Anti-aliased edge pixels cannot be bit-identical to Cairo's, so the checks are: same size, the ink bounding box within two
pixels, the amount of ink within 10%, and a mean absolute grey difference below 1 (of 255; measured 0.13-0.35).
"""
import io
import os
import sys

import pytest

Image = pytest.importorskip("PIL.Image")
from PIL import ImageChops, ImageStat  # noqa: E402

from rdoxide import Chem  # noqa: E402
from rdoxide.Chem.Draw import rdMolDraw2D  # noqa: E402

sys.path.insert(0, os.path.dirname(__file__))
from png_cases import CASES, draw  # noqa: E402

REF = os.path.join(os.path.dirname(__file__), "..", "data", "png_ref")


def _ink(img):
    return img.point(lambda p: 255 if p < 250 else 0)


@pytest.mark.parametrize("name", sorted(CASES))
def test_png_close_to_cairo(name):
    got = Image.open(io.BytesIO(draw(Chem, rdMolDraw2D, name))).convert("L")
    want = Image.open(os.path.join(REF, name + ".png")).convert("L")
    assert got.size == want.size
    gi, wi = _ink(got), _ink(want)
    gb, wb = gi.getbbox(), wi.getbbox()
    assert all(abs(a - b) <= 2 for a, b in zip(gb, wb)), (gb, wb)
    ratio = ImageStat.Stat(gi).mean[0] / ImageStat.Stat(wi).mean[0]
    assert 0.9 < ratio < 1.1, ratio
    mad = ImageStat.Stat(ImageChops.difference(got, want)).mean[0]
    assert mad < 1.0, mad
