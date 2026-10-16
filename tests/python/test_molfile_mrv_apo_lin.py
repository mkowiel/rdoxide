"""V2000 `M  APO`, `M  LIN` and `M  MRV SMA` lines; expectations validated against RDKit 2026.09.1."""

import pytest

from rdoxide import Chem


def block(extra, n=5):
    def atom(x):
        return "%10.4f%10.4f%10.4f C   0  0  0  0  0  0  0  0  0  0  0  0" % (x, 0.0, 0.0)

    lines = ["", "  probe", "", "%3d%3d  0  0  0  0            999 V2000" % (n, n - 1)]
    lines += [atom(float(i)) for i in range(n)]
    lines += ["%3d%3d  1  0" % (i + 1, i + 2) for i in range(n - 1)]
    return "\n".join(lines + extra + ["M  END", ""])


def parse(extra, strict=True):
    return Chem.MolFromMolBlock(block(extra), sanitize=False, strictParsing=strict)


def test_attach_points():
    m = parse(["M  APO  2   1   1   5   3"])
    assert m.GetAtomWithIdx(0).GetIntProp("molAttchpt") == 1
    assert m.GetAtomWithIdx(4).GetIntProp("molAttchpt") == -1
    assert not m.GetAtomWithIdx(2).HasProp("molAttchpt")


@pytest.mark.parametrize(
    "line",
    ["M  APO  1   1   4", "M  APO  1   9   1", "M  APO  1   0   1", "M  APO  1   x   1", "M  APO  x"],
)
def test_bad_attach_points_are_rejected(line):
    assert parse([line]) is None


def test_duplicate_attach_point_is_strict_only():
    lines = ["M  APO  2   1   1   1   2"]
    assert parse(lines, strict=True) is None
    assert parse(lines, strict=False).GetAtomWithIdx(0).GetIntProp("molAttchpt") == 1


def test_link_nodes_become_v3000_text():
    m = parse(["M  LIN  1   2   3   2   1"])
    assert m.GetProp("_molLinkNodes") == "1 3 2 2 2 2 1"
    m = parse(["M  LIN  2   2   3   2   1   4   3   2   2"])
    assert m.GetProp("_molLinkNodes") == "1 3 2 2 2 2 1|1 3 2 4 2 4 2"


@pytest.mark.parametrize(
    "line",
    [
        "M  LIN  1   2   3   2",  # RDKit rejects a missing substituent C
        "M  LIN  1   2    ",
        "M  LIN  1   2   1   1",
        "M  LIN  1   2   3   9",
        "M  LIN  1   9   3   1",
        "M  LIN  q",
    ],
)
def test_bad_link_nodes_are_rejected(line):
    assert parse([line]) is None


def test_marvin_smarts_adds_a_recursive_query():
    m = parse(["M  MRV SMA   1 [*;A]"])
    atom = m.GetAtomWithIdx(0)
    assert atom.GetProp("MRV SMA") == "[*;A]"
    assert atom.HasQuery()
    assert "RecursiveStructure" in atom.DescribeQuery()
    assert atom.GetIntProp("_MolFileAtomQuery") == 1


@pytest.mark.parametrize(
    "line", ["M  MRV SMA   1 [[[", "M  MRV SMA   9 [C]", "M  MRV SMA   x [C]"]
)
def test_bad_marvin_smarts_are_rejected(line):
    assert parse([line]) is None


def test_other_mrv_lines_are_ignored():
    assert parse(["M  MRV FOO   1 [C]"]).GetNumAtoms() == 5
