#!/usr/bin/env python3
"""Keep rdoxide in sync with an upstream RDKit release.

Usage (Python >= 3.11, stdlib only):

    python3 tools/sync_upstream.py fetch    # clone/checkout UPSTREAM.toml tag into upstream/rdkit
    python3 tools/sync_upstream.py regen    # re-extract data tables and copy pure-python files/test data
    python3 tools/sync_upstream.py diff     # upstream files changed since they were last ported
    python3 tools/sync_upstream.py stamp [MODULE...]  # record current upstream blob hashes as ported

`porting/manifest.toml` maps every Rust module to the upstream files it was
ported from, together with the git blob hash of each file at port time. When a
new RDKit release is checked out, `diff` lists exactly which ported upstream
files changed so the corresponding Rust code can be updated.
"""

from __future__ import annotations

import argparse
import math
import re
import shutil
import subprocess
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
UPSTREAM_DIR = ROOT / "upstream" / "rdkit"
MANIFEST = ROOT / "porting" / "manifest.toml"
DATA = ROOT / "crates" / "rdoxide" / "data"
COPIES = ROOT / "porting" / "copies.toml"


def load_upstream_cfg() -> dict:
    with open(ROOT / "UPSTREAM.toml", "rb") as fh:
        return tomllib.load(fh)["rdkit"]


def load_externals() -> dict:
    with open(ROOT / "UPSTREAM.toml", "rb") as fh:
        return tomllib.load(fh).get("external", {})


def fetch_externals() -> None:
    """Download pinned third-party sources and verify their MD5."""
    import hashlib
    import tarfile
    import urllib.request

    base = ROOT / "upstream" / "external"
    base.mkdir(parents=True, exist_ok=True)
    for name, spec in load_externals().items():
        dest = base / name
        if dest.exists():
            print(f"external {name}: present")
            continue
        archive = base / f"{name}.tar.gz"
        urllib.request.urlretrieve(spec["url"], archive)
        digest = hashlib.md5(archive.read_bytes()).hexdigest()
        if digest != spec["md5"]:
            archive.unlink()
            sys.exit(f"external {name}: md5 mismatch ({digest} != {spec['md5']})")
        with tarfile.open(archive) as tf:
            top = tf.getnames()[0].split("/")[0]
            tf.extractall(base, filter="data")
        (base / top).rename(dest)
        archive.unlink()
        print(f"external {name}: fetched")


def git(*args: str, cwd: Path = UPSTREAM_DIR) -> str:
    return subprocess.run(
        ["git", *args], cwd=cwd, check=True, capture_output=True, text=True
    ).stdout.strip()


# --------------------------------------------------------------------------- fetch
def cmd_fetch(_args) -> None:
    cfg = load_upstream_cfg()
    if not UPSTREAM_DIR.exists():
        UPSTREAM_DIR.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            ["git", "clone", "--depth", "1", "--branch", cfg["tag"], cfg["repository"],
             str(UPSTREAM_DIR)],
            check=True,
        )
    else:
        git("fetch", "--depth", "1", "origin", "tag", cfg["tag"])
        git("checkout", "--quiet", cfg["tag"])
    head = git("rev-parse", "HEAD")
    if head != cfg["commit"]:
        sys.exit(f"upstream HEAD {head} does not match UPSTREAM.toml commit {cfg['commit']}")
    print(f"upstream/rdkit at {cfg['tag']} ({head})")
    fetch_externals()


# --------------------------------------------------------------------------- regen
_RAW_RE = re.compile(r'R"DAT\((.*?)\)DAT"', re.S)


def _cxx_raw_concat(src: str, start_marker: str) -> str:
    """Concatenate the R"DAT(...)DAT" literals of the initializer after start_marker."""
    # drop whole-line // comments (they may contain ';')
    src = re.sub(r"(?m)^\s*//.*$", "", src)
    start = pos = src.index(start_marker)
    while True:
        # the initializer ends at the first ';' that is not inside a raw literal
        lit = src.find('R"DAT(', pos)
        semi = src.find(";", pos)
        if lit == -1 or semi < lit:
            return "".join(_RAW_RE.findall(src[start:semi]))
        pos = src.index(')DAT"', lit) + 5


def regen_atomic_data() -> None:
    src = (UPSTREAM_DIR / "Code/GraphMol/atomic_data.cpp").read_text()
    header = (
        "# Extracted from RDKit Code/GraphMol/atomic_data.cpp by tools/sync_upstream.py.\n"
        "# Copyright (C) 2001-2024 Greg Landrum and other RDKit contributors (BSD-3-Clause).\n"
        "# DO NOT EDIT - regenerate with `python3 tools/sync_upstream.py regen`.\n"
    )
    # periodic table: one element per line
    table = _cxx_raw_concat(src, "periodicTableAtomData =")
    lines = [ln.strip() for ln in table.splitlines() if ln.strip()]
    (DATA / "atomic_data.tsv").write_text(
        header
        + "# anum symbol row rCov rB0 rVdw mass nOuterElecs commonIsotope commonIsotopeMass valences...\n"
        + "\n".join("\t".join(ln.split()) for ln in lines)
        + "\n"
    )
    # element names
    m = re.search(r"elementNames\s*=\s*\{(.*?)\};", src, re.S)
    names = re.findall(r'"([^"]*)"', m.group(1))
    (DATA / "element_names.txt").write_text(header + "\n".join(names) + "\n")
    # isotopes (lines in upstream end with a literal "\n \" sequence)
    iso = _cxx_raw_concat(src, "isotopesAtomData[] =")
    out = []
    for ln in iso.splitlines():
        ln = ln.replace("\\n", " ").replace("\\", " ").strip()
        if ln:
            out.append("\t".join(ln.split()[:5]))
    (DATA / "isotopes.tsv").write_text(
        header + "# anum symbol isotope mass abundance\n" + "\n".join(out) + "\n"
    )
    print(f"atomic_data.tsv: {len(lines)} rows; data/isotopes.tsv: {len(out)} rows; "
          f"data/element_names.txt: {len(names)} names")


_MOLSTANDARDIZE_TABLES = [
    # (upstream .in file, initializer name, output file, fields per entry)
    ("AcidBaseCatalog/acid_base_pairs.in", "defaultAcidBasePairs", "acid_base_pairs.tsv", 3),
    ("FragmentCatalog/fragmentPatterns.in", "defaultFragments", "fragment_patterns.tsv", 2),
    ("TransformCatalog/normalizations.in", "defaultNormalizations", "normalizations.tsv", 2),
    ("TautomerCatalog/tautomerTransforms.in", "defaultTautomerTransforms",
     "tautomer_transforms.tsv", 4),
    ("TautomerCatalog/tautomerTransforms.v1.in", "defaultTautomerTransformsv1",
     "tautomer_transforms_v1.tsv", 4),
]


def _cxx_initializer_entries(src: str, name: str) -> list[list[str]]:
    """The string literals of each top-level element of the brace initializer
    of `name` (adjacent literals concatenated, as the C++ compiler does)."""
    toks = re.findall(r'//[^\n]*|"(?:[^"\\]|\\.)*"|[{}(),]|[^\s{}(),"/]+|/', src)
    i = next(k for k, t in enumerate(toks) if t == name)
    while toks[i] != "{":
        i += 1
    entries, cur, depth, prev_lit = [], [], 0, False
    for t in toks[i:]:
        if t.startswith("//"):
            continue
        if t.startswith('"'):
            lit = t[1:-1]
            if "\\" in lit:
                raise ValueError(f"{name}: escape sequences are not handled: {t}")
            if prev_lit:
                cur[-1] += lit
            else:
                cur.append(lit)
            prev_lit = True
            continue
        prev_lit = False
        if t in "{(":
            depth += 1
        elif t in "})":
            depth -= 1
            if depth == 0:
                break
        elif t == "," and depth == 1:
            entries.append(cur)
            cur = []
    if cur:
        entries.append(cur)
    return entries


def regen_molstandardize() -> None:
    """The default MolStandardize tables (compiled-in C++ initializers upstream)
    as tab-separated files for `graphmol::mol_standardize::catalogs`."""
    srcdir = UPSTREAM_DIR / "Code/GraphMol/MolStandardize"
    outdir = DATA / "molstandardize"
    outdir.mkdir(exist_ok=True)
    for path, name, out, nfields in _MOLSTANDARDIZE_TABLES:
        entries = _cxx_initializer_entries((srcdir / path).read_text(), name)
        for e in entries:
            if len(e) != nfields or any("\t" in f or "\n" in f for f in e):
                raise ValueError(f"{path}: unexpected entry {e}")
        (outdir / out).write_text(
            f"# Extracted from RDKit Code/GraphMol/MolStandardize/{path} ({name})\n"
            "# by tools/sync_upstream.py.\n"
            "# Copyright (C) RDKit contributors (BSD-3-Clause).\n"
            "# DO NOT EDIT - regenerate with `python3 tools/sync_upstream.py regen`.\n"
            "# One entry per line, fields separated by tabs (empty fields allowed).\n"
            + "".join("\t".join(e) + "\n" for e in entries)
        )
        print(f"molstandardize/{out}: {len(entries)} entries")


def regen_uff_params() -> None:
    """The default UFF parameter table (a concatenated C string literal
    upstream) as `crates/rdoxide/data/uff_params.tsv`."""
    src = (UPSTREAM_DIR / "Code/ForceField/UFF/Params.cpp").read_text()
    body = src[src.index("defaultParamData =") :]
    body = body[: body.index(";")]
    literals = re.findall(r'"((?:[^"\\]|\\.)*)"', body)
    text = "".join(literals).replace("\\t", "\t").replace("\\n", "\n")
    (DATA / "uff_params.tsv").write_text(
        "# Extracted from RDKit Code/ForceField/UFF/Params.cpp (defaultParamData)\n"
        "# by tools/sync_upstream.py.\n"
        "# Copyright (C) RDKit contributors (BSD-3-Clause).\n"
        "# DO NOT EDIT - regenerate with `python3 tools/sync_upstream.py regen`.\n"
        + text
    )
    print(f"uff_params.tsv: {text.count(chr(10))} lines")


_MMFF_TABLES = [
    "Def", "Prop", "PBCI", "Chg", "Bond", "Bndk", "HerschbachLaurie",
    "CovRadPauEle", "Stbn", "Dfsb", "Oop", "sOop", "Tor", "sTor", "VdW",
]


def _c_string_literals(text: str) -> str:
    """The decoded concatenation of the C string literals in `text`."""
    literals = re.findall(r'"((?:[^"\\]|\\.)*)"', text)
    return "".join(literals).replace("\\t", "\t").replace("\\n", "\n")


def _literal_run(text: str, start: int, separators: str = r"\s*") -> str:
    """The run of adjacent C string literals beginning at `start`; a `;` or
    `}` inside a literal does not end it."""
    literal = r'"(?:[^"\\]|\\.)*"'
    m = re.compile(rf"{separators}(?:{literal}{separators})+").match(text, start)
    assert m, "no string literal at the given position"
    return m.group(0)


def regen_mmff_params() -> None:
    """The default MMFF parameter tables (concatenated C string literals
    upstream) as `crates/rdoxide/data/mmff/*.tsv`."""
    src = (UPSTREAM_DIR / "Code/ForceField/MMFF/Params.cpp").read_text()
    outdir = DATA / "mmff"
    outdir.mkdir(exist_ok=True)
    header = (
        "# Extracted from RDKit Code/ForceField/MMFF/Params.cpp (defaultMMFF{name})\n"
        "# by tools/sync_upstream.py.\n"
        "# Copyright (C) RDKit contributors (BSD-3-Clause).\n"
        "# DO NOT EDIT - regenerate with `python3 tools/sync_upstream.py regen`.\n"
    )
    tables = {}
    for name in _MMFF_TABLES:
        start = src.index(f"const std::string defaultMMFF{name} =")
        tables[name] = _c_string_literals(_literal_run(src, src.index("=", start) + 1))
    # the angle table is an array of strings ending with "EOS"
    start = src.index("const std::string defaultMMFFAngleData[] = {")
    run = _literal_run(src, src.index("{", start) + 1, separators=r"[\s,]*")
    tables["Angle"] = _c_string_literals(run).removesuffix("EOS")
    for name, text in tables.items():
        (outdir / f"{name}.tsv").write_text(header.format(name=name) + text)
        print(f"mmff/{name}.tsv: {text.count(chr(10))} lines")


def _strip_cpp_comments(text: str) -> str:
    """`text` without `//` line comments (outside of string literals)."""
    out = []
    for line in text.splitlines():
        in_str = False
        i = 0
        while i < len(line):
            c = line[i]
            if c == "\\" and in_str:
                i += 2
                continue
            if c == '"':
                in_str = not in_str
            elif not in_str and line.startswith("//", i):
                line = line[:i]
                break
            i += 1
        out.append(line)
    return "\n".join(out)


def _cpp_tokens(text: str):
    """Tokens of C++ source (comments and preprocessor lines dropped):
    ('s', decoded string) | ('n', int) | ('i', identifier) | ('p', punctuation)."""
    i, n = 0, len(text)
    toks = []
    while i < n:
        c = text[i]
        if c in " \t\r\n":
            i += 1
        elif text.startswith("//", i):
            i = text.find("\n", i)
            i = n if i < 0 else i
        elif text.startswith("/*", i):
            i = text.index("*/", i) + 2
        elif c == "#":
            i = text.find("\n", i)
            i = n if i < 0 else i
        elif c == '"':
            j, out = i + 1, []
            while text[j] != '"':
                if text[j] == "\\":
                    j += 1
                    out.append({"n": "\n", "t": "\t"}.get(text[j], text[j]))
                else:
                    out.append(text[j])
                j += 1
            toks.append(("s", "".join(out)))
            i = j + 1
        elif c.isdigit():
            j = i
            while j < n and text[j].isalnum():
                j += 1
            toks.append(("n", int(text[i:j].rstrip("uU"))))
            i = j
        elif c.isalpha() or c == "_":
            j = i
            while j < n and (text[j].isalnum() or text[j] == "_"):
                j += 1
            toks.append(("i", text[i:j]))
            i = j
        else:
            toks.append(("p", c))
            i += 1
    return toks


def regen_filtercatalog() -> None:
    """The built-in filter catalogs (Brenk, NIH, PAINS, ZINC, ChEMBL) from
    `Filters.cpp` and its `*.in` files as `crates/rdoxide/data/filtercatalog.txt`:
    `@NAME` starts a catalog, `P<TAB>key<TAB>value` is a property and
    `E<TAB>name<TAB>smarts<TAB>max<TAB>comment` an entry."""
    srcdir = UPSTREAM_DIR / "Code/GraphMol/FilterCatalog"
    text = (srcdir / "Filters.cpp").read_text()
    text = re.sub(
        r'^#include "([\w.]+\.in)"',
        lambda m: (srcdir / m.group(1)).read_text(),
        text,
        flags=re.M,
    )
    toks = _cpp_tokens(text)
    entries, props = {}, {}
    i = 0
    while i < len(toks):
        if toks[i] == ("i", "const") and toks[i + 1][1] in ("FilterData_t", "FilterProperty_t"):
            kind, name = toks[i + 1][1], toks[i + 2][1]
            if toks[i + 3] != ("p", "[") or toks[i + 5] != ("p", "="):
                i += 1
                continue
            j = i + 7  # after the opening brace
            rows = []
            while toks[j] != ("p", "}") or toks[j - 1] == ("p", "}") and False:
                assert toks[j] == ("p", "{"), (name, toks[j])
                j += 1
                fields = [[]]
                while toks[j] != ("p", "}"):
                    if toks[j] == ("p", ","):
                        fields.append([])
                    else:
                        fields[-1].append(toks[j])
                    j += 1
                j += 1
                if toks[j] == ("p", ","):
                    j += 1
                rows.append(["".join(str(t[1]) for t in f) for f in fields])
            (entries if kind == "FilterData_t" else props)[name] = rows
            i = j
        else:
            i += 1
    out = [
        "# Extracted from RDKit Code/GraphMol/FilterCatalog/Filters.cpp and *.in",
        "# by tools/sync_upstream.py.",
        "# Copyright (c) 2015, Novartis Institutes for BioMedical Research Inc. (BSD-3-Clause).",
        "# DO NOT EDIT - regenerate with `python3 tools/sync_upstream.py regen`.",
    ]
    for name, rows in entries.items():
        out.append(f"@{name}")
        for key, value in props[f"{name}_PROPS"]:
            out.append(f"P\t{key}\t{value}")
        for ename, smarts, mx, comment in rows:
            assert "\t" not in ename + smarts + comment and "\n" not in ename + smarts + comment
            out.append(f"E\t{ename}\t{smarts}\t{mx}\t{comment}")
    (DATA / "filtercatalog.txt").write_text("\n".join(out) + "\n")
    print(f"filtercatalog.txt: {len(entries)} catalogs, {sum(map(len, entries.values()))} entries")


def regen_crystalff() -> None:
    """The experimental torsion-angle preferences of CrystalFF/ETKDG as
    `crates/rdoxide/data/crystalff/*.txt`: SMARTS and cosine-series constants
    (`torsion_*.txt`, one pattern per line) or Gaussian fits
    (`gauss_*.tsv`: SMARTS, heights, positions and widths)."""
    srcdir = UPSTREAM_DIR / "Code/GraphMol/ForceFieldHelpers/CrystalFF/torsionData"
    outdir = DATA / "crystalff"
    outdir.mkdir(exist_ok=True)
    header = (
        "# Extracted from RDKit Code/GraphMol/ForceFieldHelpers/CrystalFF/torsionData/{src}\n"
        "# by tools/sync_upstream.py.\n"
        "# Copyright (C) RDKit contributors (BSD-3-Clause).\n"
        "# DO NOT EDIT - regenerate with `python3 tools/sync_upstream.py regen`.\n"
    )
    cosine = {
        "torsion_v1": "torsionPreferences_v1.in",
        "torsion_v2": "torsionPreferences_v2.in",
        "torsion_smallrings": "torsionPreferences_smallrings.in",
        "torsion_macrocycles": "torsionPreferences_macrocycles.in",
    }
    for name, src in cosine.items():
        text = (srcdir / src).read_text()
        # skip the block comment, then take the literals after the `=`
        body = _strip_cpp_comments(text[text.index("const std::string") :])
        body = body[body.index("=") + 1 :]
        data = _c_string_literals(_literal_run(body, 0))
        (outdir / f"{name}.txt").write_text(header.format(src=src) + data)
        print(f"crystalff/{name}.txt: {data.count(chr(10))} lines")
    gauss = {
        "gauss_v4": "gaussTorsionPreferences_v4.h",
        "gauss_smallrings": "gaussTorsionPreferences_v4_smallrings.h",
        "gauss_macrocycles": "gaussTorsionPreferences_v4_macrocycles.h",
    }
    entry = re.compile(
        r'^\{("(?:[^"\\]|\\.)*"),\s*\{([^}]*)\},\s*\{([^}]*)\},\s*\{([^}]*)\}\},?$'
    )
    for name, src in gauss.items():
        lines = _strip_cpp_comments((srcdir / src).read_text()).splitlines()
        rows = []
        for line in lines:
            line = line.strip()
            if not line.startswith('{"'):
                continue
            m = entry.match(line)
            assert m, f"unparsable torsion entry in {src}: {line}"
            smarts = m.group(1)[1:-1]
            cols = [smarts]
            for group in m.groups()[1:]:
                values = [eval(v, {"PI": math.pi}) for v in group.split(",")]
                cols.append(",".join(repr(float(v)) for v in values))
            rows.append("\t".join(cols))
        (outdir / f"{name}.tsv").write_text(header.format(src=src) + "\n".join(rows) + "\n")
        print(f"crystalff/{name}.tsv: {len(rows)} entries")


def regen_common_properties() -> None:
    """Generate rdgeneral/common_properties.rs from Code/RDGeneral/types.h."""
    src = (UPSTREAM_DIR / "Code/RDGeneral/types.h").read_text()
    block = src[src.index("namespace common_properties {"):]
    block = block[: block.index("}  // namespace common_properties")]
    out = [
        "//! Property names shared across RDKit (`RDKit::common_properties`).",
        "//!",
        "//! GENERATED from RDKit `Code/RDGeneral/types.h` by tools/sync_upstream.py regen.",
        "//! DO NOT EDIT. Copyright (C) RDKit contributors, BSD-3-Clause.",
        "#![allow(non_upper_case_globals)]",
        "",
    ]
    # declarations may wrap onto a second line
    block = re.sub(r"//[^\n]*", "", block)
    tokens = re.findall(
        r'inline constexpr std::string_view (\w+)\s*=\s*("(?:[^"\\]|\\.)*");'
        r"|namespace (\w+) \{|(\})",
        block,
    )
    n, depth = 0, 0
    for name, value, ns, close in tokens:
        if name:
            out.append(f"{'    ' * depth}pub const {name}: &str = {value};")
            n += 1
        elif ns == "common_properties":
            continue
        elif ns:
            out.append(f"{'    ' * depth}pub mod {ns} {{")
            depth += 1
        elif close and depth:
            depth -= 1
            out.append(f"{'    ' * depth}}}")
    path = ROOT / "crates/rdoxide/src/rdgeneral/common_properties.rs"
    path.write_text("\n".join(out) + "\n")
    print(f"common_properties.rs: {n} names")


def regen_copies() -> None:
    """Copy verbatim files (pure-python modules, test data) listed in porting/copies.toml."""
    if not COPIES.exists():
        return
    with open(COPIES, "rb") as fh:
        cfg = tomllib.load(fh)
    n = 0
    for entry in cfg.get("copy", []):
        src = UPSTREAM_DIR / entry["from"]
        dst = ROOT / entry["to"]
        if src.is_dir():
            if dst.exists():
                shutil.rmtree(dst)
            shutil.copytree(src, dst, ignore=shutil.ignore_patterns(*entry.get("exclude", [])))
        else:
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
        n += 1
    print(f"copied {n} upstream paths")


def cmd_regen(_args) -> None:
    if not UPSTREAM_DIR.exists():
        sys.exit("upstream/rdkit missing - run `sync_upstream.py fetch` first")
    regen_atomic_data()
    regen_molstandardize()
    regen_uff_params()
    regen_mmff_params()
    regen_crystalff()
    regen_filtercatalog()
    regen_common_properties()
    regen_copies()


# --------------------------------------------------------------------------- diff / stamp
def _load_manifest() -> dict:
    with open(MANIFEST, "rb") as fh:
        return tomllib.load(fh)


def _blob(path: str) -> str:
    return git("rev-parse", f"HEAD:{path}")


def cmd_diff(_args) -> None:
    manifest = _load_manifest()
    stale = 0
    for module, spec in manifest.get("module", {}).items():
        for path, blob in spec.get("sources", {}).items():
            try:
                cur = _blob(path)
            except subprocess.CalledProcessError:
                cur = "<deleted>"
            if cur != blob:
                stale += 1
                print(f"{module}: {path} changed ({blob[:10]} -> {cur[:10]})")
    print(f"{stale} upstream file(s) changed since last port" if stale else "all ported files up to date")
    sys.exit(1 if stale else 0)


def cmd_stamp(args) -> None:
    """Rewrite blob hashes in manifest.toml (simple line-based rewrite keeps comments)."""
    text = MANIFEST.read_text()
    out, module = [], None
    for line in text.splitlines():
        mh = re.match(r'\[module\."([^"]+)"', line) or re.match(r"\[module\.([\w.]+)", line)
        if mh:
            module = mh.group(1)
        ms = re.match(r'^(\s*)"([^"]+)"\s*=\s*"[^"]*"(.*)$', line)
        if ms and module and (not args.modules or module in args.modules):
            indent, path, rest = ms.groups()
            line = f'{indent}"{path}" = "{_blob(path)}"{rest}'
        out.append(line)
    MANIFEST.write_text("\n".join(out) + "\n")
    print("manifest stamped")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("fetch").set_defaults(fn=cmd_fetch)
    sub.add_parser("regen").set_defaults(fn=cmd_regen)
    sub.add_parser("diff").set_defaults(fn=cmd_diff)
    s = sub.add_parser("stamp")
    s.add_argument("modules", nargs="*")
    s.set_defaults(fn=cmd_stamp)
    args = p.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
