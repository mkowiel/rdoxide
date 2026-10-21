"""Generates tests/data/inchi_onecomp.json from the InChI C code:
`CreateOneComponentINChI` of runichi.c (allocation of the results via
`Alloc_INChI`/`Alloc_INChI_Aux`, `Create_INChI`, error detection, warnings and
component bookkeeping) on random small molecules. Links the upstream InChI base
library (see inchi_markalt_reference.py).

    python tools/inchi_onecomp_reference.py
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


cr = load("create", "tools/inchi_create_reference.py")
ma = cr.ma
SRC = ma.SRC

HARNESS = cr.HARNESS.split("int main(void)")[0] + r"""
#include "ichimain.h"
int CreateOneComponentINChI( CANON_GLOBALS *pCG, INCHI_CLOCK *ic, STRUCT_DATA *sd, INPUT_PARMS *ip, INP_ATOM_DATA *inp_cur_data,
                             ORIG_ATOM_DATA *orig_inp_data, PINChI2 *pINChI, PINChI_Aux2 *pINChI_Aux, int iINChI, int i, long num_inp,
                             INP_ATOM_DATA **inp_norm_data, NORM_CANON_FLAGS *pncFlags, INCHI_IOSTREAM *log_file );
int main(void) {
    int ncases;
    setvbuf(stdout, NULL, _IONBF, 0);
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        int n, which, sddone, nowarn, ii, outopt;
        unsigned long long mode, flags, done;
        int loose;
        scanf("%d %d %d %d %d", &sddone, &nowarn, &ii, &outopt, &which);
        scanf("%d %d %llu %llu %llu %d", &n, &which, &mode, &flags, &done, &loose);
        inp_ATOM *at = calloc(n + 8, sizeof(inp_ATOM));
        for (int i = 0; i < n; i++) {
            int v, cbv, el, rad, chg, nh, i0, i1, i2, aw;
            double x, y, z;
            scanf("%d %d %d %d %d %d %d %d %d %d %lf %lf %lf", &v, &cbv, &el, &rad, &chg, &nh, &i0, &i1, &i2, &aw, &x, &y, &z);
            at[i].valence = v; at[i].chem_bonds_valence = cbv; at[i].el_number = el; strcpy(at[i].elname, elsym(el));
            at[i].radical = rad; at[i].charge = chg; at[i].num_H = nh;
            at[i].num_iso_H[0] = i0; at[i].num_iso_H[1] = i1; at[i].num_iso_H[2] = i2; at[i].iso_atw_diff = aw;
            at[i].orig_at_number = i + 1; at[i].x = x; at[i].y = y; at[i].z = z;
            for (int j = 0; j < v; j++) { int q, b; scanf("%d %d", &q, &b); at[i].neighbor[j] = q; at[i].bond_type[j] = b; }
        }
        CANON_GLOBALS cg;
        memset(&cg, 0, sizeof(cg));
        cg.m_num_bit = 16;
        cg.m_bBit = calloc(16, sizeof(bitWord));
        for (int i = 0; i < 16; i++) cg.m_bBit[i] = 1 << i;
        INCHI_CLOCK ic; memset(&ic, 0, sizeof(ic));
        INPUT_PARMS ip; memset(&ip, 0, sizeof(ip));
        ip.bLooseTSACheck = loose; ip.nMode = mode; ip.bTautFlags = flags; ip.bTautFlagsDone = done;
        ip.bNoWarnings = nowarn; ip.bINChIOutputOptions = outopt;
        STRUCT_DATA sd; memset(&sd, 0, sizeof(sd));
        sd.bTautFlagsDone[0] = sddone;
        INP_ATOM_DATA cur; memset(&cur, 0, sizeof(cur));
        cur.at = at; cur.num_at = n;
        INP_ATOM_DATA nd[2]; memset(nd, 0, sizeof(nd));
        INP_ATOM_DATA *pnd[2] = { &nd[0], &nd[1] };
        PINChI2 pi[3]; PINChI_Aux2 pa[3];
        memset(pi, 0, sizeof(pi)); memset(pa, 0, sizeof(pa));
        NORM_CANON_FLAGS nc; memset(&nc, 0, sizeof(nc));
        printf("#CASE %d\n", c);
        int ret = CreateOneComponentINChI(&cg, &ic, &sd, &ip, &cur, NULL, pi, pa, ii, 0, 1, pnd, &nc, NULL);
        printf("R %d %d %d %d %d %d\n", ret, sd.nErrorCode, sd.num_non_taut[ii], sd.num_taut[ii], cur.num_bonds, cur.num_isotopic);
        printf("E %s\n", sd.pStrErrStruct);
        printf("C %d\n", n ? at[0].component : -1);
        for (int k = 0; k < 2; k++) {
            printf("B%d %d %d %d %d\n", k, nd[k].bExists, nd[k].bHasIsotopicLayer, nd[k].num_at, nd[k].at != NULL);
            printf("F%d %llu %llu %llu\n", k, (unsigned long long)nc.bTautFlags[ii][k], (unsigned long long)nc.bTautFlagsDone[ii][k],
                   (unsigned long long)nc.bNormalizationFlags[ii][k]);
            INChI *p = pi[0][k]; INChI_Aux *a = pa[0][k];
            if (!p) { printf("I%d -\n", k); continue; }
            printf("I%d %d %d %d %s %d %d %d %d\n", k, p->nErrorCode, p->nNumberOfAtoms, p->nTotalCharge, p->szHillFormula ? p->szHillFormula : "NULL",
                   p->lenConnTable, p->lenTautomer, p->nNumberOfIsotopicAtoms, p->nNumberOfIsotopicTGroups);
            printf("A%d %d %d %d %d\n", k, a ? a->nNumberOfAtoms : -1, a ? a->nNumRemovedProtons : -1,
                   a ? a->nNumRemovedIsotopicH[0] + a->nNumRemovedIsotopicH[1] + a->nNumRemovedIsotopicH[2] : -1, a ? a->nCanonFlags : -1);
        }
    }
    return 0;
}
"""


def main():
    rng = random.Random(777)
    cases = []
    for _ in range(3000):
        atoms = cr.gen(rng)
        flags = 0
        for f in cr.FLAGS:
            if rng.random() < 0.5:
                flags |= f
        flags |= 0x1
        which = 3
        c = {"atoms": atoms, "which": which, "mode": rng.choice(cr.MODES), "flags": flags,
             "done": rng.choice([0, 0, 0, 0x400]), "loose": rng.randrange(2)}
        head = f"{rng.choice([0, 0, 0x4000, 0x8000, 0xC000])} {rng.randrange(2)} {rng.randrange(3)} {rng.choice([0, 0, 1, 2])} 0\n"
        cases.append((head, c))
    objs = ma.build_lib()
    fmt = lambda hc: hc[0] + cr.fmt(hc[1])
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "h.c").write_text(HARNESS)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), *objs, "-lm"], check=True)
        r = subprocess.run([str(exe)], input=(f"{len(cases)}\n" + "".join(fmt(c) for c in cases)).encode(), capture_output=True)
        if r.returncode == 0:
            out = r.stdout.decode().strip().split("\n")
        else:
            kept, lines = [], []
            for c in cases:
                r1 = subprocess.run([str(exe)], input=("1\n" + fmt(c)).encode(), capture_output=True)
                if r1.returncode == 0:
                    kept.append(c)
                    lines += [f"#CASE {len(kept) - 1}"] + r1.stdout.decode().strip().split("\n")[1:]
            print(len(cases) - len(kept), "cases dropped")
            cases = kept
            out = lines
    starts = [i for i, l in enumerate(out) if l.startswith("#CASE ")]
    starts.append(len(out))
    assert len(starts) - 1 == len(cases)
    data = [{"in": fmt(c), "out": out[a + 1:b]} for c, a, b in zip(cases, starts, starts[1:])]
    path = ROOT / "tests/data/inchi_onecomp.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
