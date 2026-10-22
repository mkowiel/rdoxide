"""Generates tests/data/inchi_rvr_cmp2.json from the InChI C code:
`FillOutCMP2FHINCHI` and `FillOutCMP2MHINCHI` (ichirvr4.c) on random
`StrFromINChI` data.

    python tools/inchi_rvr_cmp2_reference.py
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


r1 = load("rvr1", "tools/inchi_rvr1_reference.py")
r2, ma, SRC = r1.r2, r1.ma, r1.SRC

MAIN = r1.PRE + r"""
static int rd(void) { int x; if (scanf("%d", &x) != 1) exit(3); return x; }
static void read_ti(T_GROUP_INFO *ti, int n) {
    memset(ti, 0, sizeof(*ti));
    int nt = rd();
    ti->num_t_groups = nt;
    ti->t_group = calloc(nt + 1, sizeof(T_GROUP));
    for (int i = 0; i < nt; i++) {
        ti->t_group[i].num[0] = rd(); ti->t_group[i].num[1] = rd();
        ti->t_group[i].nNumEndpoints = rd(); ti->t_group[i].nFirstEndpointAtNoPos = rd();
    }
    int ne = rd();
    ti->nNumEndpoints = ne;
    ti->nEndpointAtomNumber = calloc(ne + 1, sizeof(AT_NUMB));
    for (int i = 0; i < ne; i++) ti->nEndpointAtomNumber[i] = rd();
}
static INChI *read_inchi(int with_fixed, int *present) {
    *present = rd();
    if (!*present) return NULL;
    INChI *c = calloc(1, sizeof(INChI));
    c->nNumberOfAtoms = rd(); c->bDeleted = rd();
    if (rd()) { int n = rd(); c->nNum_H = calloc(n + 1, 1); for (int i = 0; i < n; i++) c->nNum_H[i] = rd(); }
    if (with_fixed && rd()) { int n = rd(); c->nNum_H_fixed = calloc(n + 1, 1); for (int i = 0; i < n; i++) c->nNum_H_fixed[i] = rd(); }
    c->nTotalCharge = rd();
    return c;
}
int main(void) {
    int ncases = rd();
    for (int cs = 0; cs < ncases; cs++) {
        printf("#CASE %d\n", cs);
        int mode = rd(), n = rd();
        StrFromINChI s; memset(&s, 0, sizeof(s));
        s.num_atoms = n;
        read_ti(&s.ti, n);
        read_ti(&s.One_ti, n);
        s.One_ti.tni.nNumRemovedProtons = rd();
        s.nNumRemovedProtonsMobHInChI = rd();
        s.endpoint = calloc(n + 1, sizeof(AT_NUMB)); s.fixed_H = calloc(n + 1, 1);
        for (int i = 0; i < n; i++) s.endpoint[i] = rd();
        for (int i = 0; i < n; i++) s.fixed_H[i] = rd();
        VAL_AT *va = calloc(n + 1, sizeof(VAL_AT));
        for (int i = 0; i < n; i++) { va[i].cNumValenceElectrons = rd(); va[i].cPeriodicRowNumber = rd(); va[i].cMetal = rd(); }
        inp_ATOM *at2 = calloc(n + 1, sizeof(inp_ATOM));
        for (int i = 0; i < n; i++) { at2[i].num_H = rd(); at2[i].charge = rd(); at2[i].valence = rd(); at2[i].chem_bonds_valence = rd(); at2[i].endpoint = rd(); }
        ALL_TC_GROUPS g; memset(&g, 0, sizeof(g));
        g.tgroup_charge = rd();
        INChI *pin[2]; int pr[2];
        pin[0] = read_inchi(0, &pr[0]); pin[1] = read_inchi(0, &pr[1]);
        for (int k = 0; k < 2; k++) s.pOneINChI[k] = read_inchi(1, &pr[k]);
        for (int k = 0; k < 2; k++) {
            if (rd()) {
                s.pOne_norm_data[k] = calloc(1, sizeof(INP_ATOM_DATA));
                s.pOne_norm_data[k]->at = calloc(n + 1, sizeof(inp_ATOM));
                for (int i = 0; i < n; i++) { s.pOne_norm_data[k]->at[i].endpoint = rd(); s.pOne_norm_data[k]->at[i].charge = rd(); }
            }
        }
        s.nAtno2Canon[0] = calloc(n + 1, sizeof(AT_NUMB)); s.nCanon2Atno[0] = calloc(n + 1, sizeof(AT_NUMB));
        for (int i = 0; i < n; i++) { s.nAtno2Canon[0][i] = rd(); s.nCanon2Atno[0][s.nAtno2Canon[0][i]] = i; }
        if (mode == 0) {
            CMP2FHINCHI c; int r = FillOutCMP2FHINCHI(&s, at2, va, pin, &c);
            printf("%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", r, c.len_c2at, c.nNumRemHInChI, c.nNumRemHRevrs, c.nNumTgInChI, c.nNumTgRevrs,
                   c.nNumEndpInChI, c.nNumEndpRevrs, c.nNumTgDiffMinus, c.nNumTgDiffH, c.nNumTgMInChI, c.nNumTgHInChI, c.nNumTgMRevrs, c.nNumTgHRevrs,
                   c.nChargeFixHInChI, c.nChargeMobHInChI, c.nChargeFixHRevrs, c.nChargeMobHRevrs, c.nChargeFixHRevrsNonMetal, c.nChargeMobHRevrsNonMetal,
                   c.bFixedHLayerExistsRevrs, c.bHasDifference, c.nNumDiffMobH);
            for (int i = 0; i < c.len_c2at; i++) {
                CMP2FHATOMS *a = &c.c2at[i];
                printf("a %d %d %d %d %d %d %d %d %d %d %d\n", a->endptInChI, a->endptRevrs, a->atomNumber, a->nValElectr, a->nPeriodNum, a->nFixHInChI, a->nFixHRevrs, a->nMobHInChI, a->nMobHRevrs, a->nNumHRevrs, a->nAtChargeRevrs);
            }
        } else {
            CMP2MHINCHI c; int r = FillOutCMP2MHINCHI(&s, &g, at2, va, pin, &c);
            printf("%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d", r, c.len_c2at, c.nNumRemHInChI, c.nNumRemHRevrs, c.nNumTgInChI, c.nNumTgRevrs,
                   c.nNumEndpInChI, c.nNumEndpRevrs, c.nNumTgDiffMinus, c.nNumTgDiffH, c.nNumTgMInChI, c.nNumTgHInChI, c.nNumTgOInChI, c.nNumTgNInChI,
                   c.nNumTgMRevrs, c.nNumTgHRevrs, c.nNumTgORevrs, c.nNumTgNRevrs);
            printf(" %d %d %d %d %d %d %d %d %d %d", c.nNumTgOMinusRevrs, c.nNumTgOHRevrs, c.nNumTgDBORevrs, c.nNumTgNMinusRevrs, c.nNumTgNHMinusRevrs, c.nNumTgNHRevrs, c.nNumTgNH2Revrs, c.nNumTgDBNHRevrs, c.nNumTgDBNMinusRevrs, c.nNumTgDBNRevrs);
            printf(" %d %d %d %d %d %d %d %d %d %d", c.nNumTgOMinusInChI, c.nNumTgOHInChI, c.nNumTgDBOInChI, c.nNumTgNMinusInChI, c.nNumTgNHMinusInChI, c.nNumTgNHInChI, c.nNumTgNH2InChI, c.nNumTgDBNHInChI, c.nNumTgDBNMinusInChI, c.nNumTgDBNInChI);
            printf(" %d %d %d %d %d %d\n", c.nChargeMobHInChI, c.nChargeMobHRevrs, c.nChargeMobHRevrsNonMetal, c.bFixedHLayerExistsRevrs, c.bHasDifference, c.nNumDiffMobH);
            for (int i = 0; i < c.len_c2at; i++) {
                CMP2MHATOMS *a = &c.c2at[i];
                printf("a %d %d %d %d %d %d %d %d %d\n", a->endptInChI, a->endptRevrs, a->atomNumber, a->nValElectr, a->nPeriodNum, a->nMobHInChI, a->nMobHRevrs, a->nNumHRevrs, a->nAtChargeRevrs);
            }
        }
    }
    return 0;
}
"""

HARNESS = r2.HARNESS.split("int main(void)")[0] + MAIN

R = random.Random(20261022)


def gen_ti(n, mh_valid_atoms):
    nt = R.choice([0, 1, 1, 2, 3])
    tgs, eps = [], []
    for _ in range(nt):
        ne = R.randrange(0, 4)
        h = R.randrange(0, 4)
        m = R.randrange(0, 3)
        tgs.append((h + m, m, ne, len(eps)))
        eps += [R.choice(mh_valid_atoms) for _ in range(ne)]
    return [nt] + [x for t in tgs for x in t] + [len(eps)] + eps


def inchi_part(n, with_fixed):
    if R.random() < 0.15:
        return [0]
    out = [1, R.choice([0, n, n]), R.choice([0, 0, 0, 1])]
    if R.random() < 0.8:
        out += [1, n] + [R.randrange(0, 3) for _ in range(n)]
    else:
        out += [0]
    if with_fixed:
        if R.random() < 0.7:
            out += [1, n] + [R.randrange(0, 3) for _ in range(n)]
        else:
            out += [0]
    out += [R.randrange(-2, 3)]
    return out


def case():
    n = R.randrange(1, 9)
    mode = R.randrange(0, 2)
    ve = [R.choice([5, 6, 6, 5, 4, 7]) for _ in range(n)]
    row = [1 if v == 5 else R.choice([1, 2, 3, 4]) for v in ve]
    metal = [R.choice([0, 0, 0, 1]) for _ in range(n)]
    good = [i for i in range(n) if (ve[i] == 5 and row[i] == 1) or ve[i] == 6] or list(range(n))
    v = [mode, n]
    v += gen_ti(n, good if R.random() < 0.9 else list(range(n)))
    v += gen_ti(n, good if R.random() < 0.9 else list(range(n)))
    v += [R.randrange(0, 3), R.randrange(0, 3)]
    v += [R.choice([0, 0, 1, 2]) for _ in range(n)]
    v += [R.choice([0, 0, 1]) for _ in range(n)]
    for i in range(n):
        v += [ve[i], row[i], metal[i]]
    for i in range(n):
        val = R.randrange(1, 4)
        v += [R.randrange(0, 3), R.choice([-1, 0, 0, 1]), val, val + R.choice([0, 0, 1, 1, 2]), R.choice([0, 0, 1, 2])]
    v += [R.randrange(-2, 3)]
    v += inchi_part(n, False) + inchi_part(n, False)
    v += inchi_part(n, True) + inchi_part(n, True)
    for _ in range(2):
        if R.random() < 0.7:
            v += [1] + [x for _ in range(n) for x in (R.choice([0, 0, 1, 2]), R.choice([-1, 0, 0, 1]))]
        else:
            v += [0]
    perm = list(range(n))
    R.shuffle(perm)
    v += perm
    return " ".join(map(str, v))


def main():
    cases = [case() for _ in range(6000)]
    objs = ma.build_lib()
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "h.c").write_text(HARNESS)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), *objs, "-lm"], check=True)
        data = []
        for c in cases:
            try:
                r = subprocess.run([str(exe)], input=f"1\n{c}\n".encode(), capture_output=True, timeout=5)
            except subprocess.TimeoutExpired:
                continue
            if r.returncode != 0:
                continue
            lines = [x for x in r.stdout.decode(errors="replace").split("\n")[1:] if x]
            data.append({"in": c, "out": lines})
    print(len(cases) - len(data), "cases dropped")
    path = ROOT / "tests/data/inchi_rvr_cmp2.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
