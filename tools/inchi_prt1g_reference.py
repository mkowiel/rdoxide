"""Generates tests/data/inchi_prt1g.json from the InChI C code: the original
structure text and the reversibility part of the AuxInfo of ichiprt1.c
(`OrigStruct_FillOut`, `OutputAUXINFO_ReversibilityInfo`). The harness
includes ichiprt1.c itself (the latter is static). Links the upstream InChI
base library (see inchi_markalt_reference.py).

    python tools/inchi_prt1g_reference.py
"""
import importlib.util
import json
import random
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load(name, rel):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


p1e = load("prt1e", "tools/inchi_prt1e_reference.py")
p1f = load("prt1f", "tools/inchi_prt1f_reference.py")
p3c = p1e.p3c
ma = p3c.ma
SRC = ma.SRC

HARNESS = p3c.HARNESS.split("int str_AuxEqu(")[0] + r"""
#include "ichiprt1.c"
static void rd_coord(MOL_COORD c) {
    char t[64]; scanf("%s", t);
    memset(c, 'X', sizeof(MOL_COORD));
    int n = (int) strlen(t);
    for (int i = 0; i < n && i < (int) sizeof(MOL_COORD); i++) c[i] = t[i] == '_' ? ' ' : t[i] == '~' ? '\0' : t[i];
}
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        CANON_GLOBALS cg; memset(&cg, 0, sizeof(cg));
        STRUCT_DATA sd; memset(&sd, 0, sizeof(sd));
        ORIG_ATOM_DATA od; memset(&od, 0, sizeof(od));
        ORIG_STRUCT os; memset(&os, 0, sizeof(os));
        INCHI_IOS_STRING sb; INCHI_IOSTREAM out; INCHI_OUT_CTL io; memset(&io, 0, sizeof(io));
        printf("#CASE %d\n", c);
        int n = rd(), has_coord = rd(), alloc = rd(), plain = rd(), second = rd(), ntag = rd(), fh = rd(), lfopt = rd();
        sd.bChiralFlag = rd(); od.n_zy = rd();
        inp_ATOM *at = inchi_calloc(n + 1, sizeof(inp_ATOM));
        od.num_inp_atoms = n; od.at = at;
        if (has_coord) od.szCoord = inchi_calloc(n + 1, sizeof(MOL_COORD));
        for (int k = 0; k < n; k++) {
            char el[16]; scanf("%s", el);
            strcpy(at[k].elname, el);
            at[k].el_number = rd();
            at[k].valence = rd();
            for (int m = 0; m < at[k].valence; m++) at[k].neighbor[m] = rd();
            for (int m = 0; m < at[k].valence; m++) at[k].bond_type[m] = rd();
            for (int m = 0; m < at[k].valence; m++) at[k].bond_stereo[m] = rd();
            at[k].chem_bonds_valence = rd();
            for (int m = 0; m < 3; m++) at[k].sb_parity[m] = rd();
            for (int m = 0; m < 3; m++) at[k].sb_ord[m] = rd();
            for (int m = 0; m < 3; m++) at[k].sn_ord[m] = rd();
            at[k].orig_at_number = rd();
            at[k].num_H = rd();
            at[k].charge = rd();
            at[k].radical = rd();
            for (int m = 0; m < 3; m++) at[k].num_iso_H[m] = rd();
            at[k].iso_atw_diff = rd();
            for (int m = 0; m < 3; m++) { MOL_COORD mc; rd_coord(mc); if (has_coord && m == 0) memcpy(od.szCoord[k], mc, sizeof(MOL_COORD)); }
        }
        int ret = OrigStruct_FillOut(&cg, &od, &os, &sd);
        printf("%d %d %d [%s] [%s] [%s] [%s]\n", ret, os.num_atoms, os.n_zy, os.szAtoms ? os.szAtoms : "-", os.szBonds ? os.szBonds : "-", os.szCoord ? os.szCoord : "-", sd.pStrErrStruct);
        inchi_strbuf_init(&sb, alloc, 16);
        inchi_ios_init(&out, INCHI_IOS_TYPE_STRING, NULL);
        io.bPlainTextTags = plain; io.bSecondNonTautPass = second; io.nTag = ntag; io.bFhTag = fh;
        char *pLF, *pTAB;
        set_line_separators(lfopt ? 0x80 : 0, &pLF, &pTAB);
        zero_stack();
        ret = OutputAUXINFO_ReversibilityInfo(&cg, &out, &sb, ret ? NULL : &os, &io, pLF, pTAB);
        printf("%d [%s] [%s] %d\n", ret, io.szTag2, out.s.pStr ? out.s.pStr : "", sb.nUsedLength);
    }
    return 0;
}
"""

R = random.Random(2024)
SYM = {1: "H", 6: "C", 7: "N", 8: "O", 11: "Na", 15: "P", 16: "S", 17: "Cl", 26: "Fe", 30: "Zn"}


def main():
    cases = []
    for _ in range(6000):
        atoms = p1f.mk_atoms(R)
        n = len(atoms)
        parts = []
        has_coord = R.choice([1, 1, 1, 0])
        for idx, (el, val, nb, types, st, cbv, sb, so, sn) in enumerate(atoms):
            extra = [idx + 1 if R.random() < 0.97 else R.randrange(0, 5), R.choice([0, 0, 1, 2]), R.choice([0, 0, 0, 1, -1]),
                     R.choice([0, 0, 0, 0, 2]), *[R.choice([0, 0, 0, 1]) for _ in range(3)], R.choice([0, 0, 0, 1, -1])]
            sym = SYM.get(el, "C")
            if R.random() < 0.05:
                sym = "Zy"
            parts.append(" ".join(map(str, [sym, el, val, *nb, *types, *st, cbv, *sb, *so, *sn, *extra])) + " "
                         + " ".join(p1e.mk_coord(R) if m == 0 else "0" for m in range(3)))
        cases.append(" ".join(map(str, [n, has_coord, R.choice([100, 100, 120, 200, 1000]), R.choice([0, 1]), R.choice([0, 0, 0, 1]),
                                        R.choice([0, 1, 2, 3]), R.choice([0, 1]), R.choice([0, 1]), R.choice([0, 1, 2, 3]), R.choice([0, 0, 1, 2])]))
                     + " " + " ".join(parts))
    objs = ma.build_lib()
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "h.c").write_text(HARNESS)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), *objs, "-lm"], check=True)
        r = subprocess.run([str(exe)], input=(f"{len(cases)}\n" + "\n".join(cases) + "\n").encode(), capture_output=True)
        if r.returncode == 0:
            out = r.stdout.decode(errors="replace").split("\n")
        else:
            kept, out = [], []
            for c in cases:
                r1 = subprocess.run([str(exe)], input=f"1\n{c}\n".encode(), capture_output=True)
                if r1.returncode == 0:
                    kept.append(c)
                    out += [f"#CASE {len(kept) - 1}"] + r1.stdout.decode(errors="replace").split("\n")[1:]
            print(len(cases) - len(kept), "cases dropped")
            cases = kept
    starts = [i for i, l in enumerate(out) if l.startswith("#CASE ")]
    starts.append(len(out))
    assert len(starts) - 1 == len(cases), (len(starts), len(cases))
    data = [{"in": c, "out": [x for x in out[a + 1:b] if x != ""]} for c, a, b in zip(cases, starts, starts[1:])]
    data = [x for x in data if all(ch.isprintable() and ord(ch) < 127 for l in x["out"] for ch in l)]
    path = ROOT / "tests/data/inchi_prt1g.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
