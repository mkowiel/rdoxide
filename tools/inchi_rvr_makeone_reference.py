"""Generates tests/data/inchi_rvr_makeone.json from the InChI C code:
`MakeOneInChIOutOfStrFromINChI` (ichirvr1.c) on the random small molecules of
`inchi_create_reference.py`, as restored structures (`StrFromINChI.at2`).

    python tools/inchi_rvr_makeone_reference.py
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
r1 = load("rvr1", "tools/inchi_rvr1_reference.py")
r2, ma, SRC = r1.r2, r1.ma, r1.SRC

MAIN = r1.PRE + r"""
static const char *elsym(int el) {
    switch (el) { case 6: return "C"; case 7: return "N"; case 8: return "O"; case 15: return "P"; case 16: return "S"; case 17: return "Cl"; default: return "Xx"; }
}
static void dump_one_inchi(int k, INChI *p) {
    if (!p) { printf("I%d -\n", k); return; }
    printf("I%d %d %lu %d %d %s\n", k, p->nErrorCode, (unsigned long)p->nFlags, p->nTotalCharge, p->nNumberOfAtoms, p->szHillFormula ? p->szHillFormula : "NULL");
    printf("IA%d", k); for (int i = 0; i < p->nNumberOfAtoms; i++) printf(" %d", p->nAtom[i]); printf("\n");
    printf("IC%d %d", k, p->lenConnTable); for (int i = 0; i < p->lenConnTable; i++) printf(" %d", p->nConnTable[i]); printf("\n");
    printf("IT%d %d", k, p->lenTautomer); for (int i = 0; i < p->lenTautomer; i++) printf(" %d", p->nTautomer[i]); printf("\n");
    printf("IH%d", k); for (int i = 0; i < p->nNumberOfAtoms; i++) printf(" %d/%d", p->nNum_H[i], p->nNum_H_fixed ? p->nNum_H_fixed[i] : -1); printf("\n");
}
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        int n, bmob, imob;
        unsigned long long mode, flags;
        int loose;
        scanf("%d %d %d %llu %llu %d", &n, &bmob, &imob, &mode, &flags, &loose);
        inp_ATOM *at2 = calloc(n + 8, sizeof(inp_ATOM));
        for (int i = 0; i < n; i++) {
            int v, cbv, el, rad, chg, nh, i0, i1, i2, aw;
            double x, y, z;
            scanf("%d %d %d %d %d %d %d %d %d %d %lf %lf %lf", &v, &cbv, &el, &rad, &chg, &nh, &i0, &i1, &i2, &aw, &x, &y, &z);
            at2[i].valence = v; at2[i].chem_bonds_valence = cbv; at2[i].el_number = el; strcpy(at2[i].elname, elsym(el));
            at2[i].radical = rad; at2[i].charge = chg; at2[i].num_H = nh;
            at2[i].num_iso_H[0] = i0; at2[i].num_iso_H[1] = i1; at2[i].num_iso_H[2] = i2; at2[i].iso_atw_diff = aw;
            at2[i].orig_at_number = i + 1; at2[i].x = x; at2[i].y = y; at2[i].z = z;
            for (int j = 0; j < v; j++) { int q, b; scanf("%d %d", &q, &b); at2[i].neighbor[j] = q; at2[i].bond_type[j] = b; }
        }
        inp_ATOM *at3 = calloc(n + 8, sizeof(inp_ATOM));
        CANON_GLOBALS cg;
        memset(&cg, 0, sizeof(cg));
        cg.m_num_bit = 16;
        cg.m_bBit = calloc(16, sizeof(bitWord));
        for (int i = 0; i < 16; i++) cg.m_bBit[i] = 1 << i;
        INCHI_CLOCK ic; memset(&ic, 0, sizeof(ic));
        INPUT_PARMS ip; memset(&ip, 0, sizeof(ip));
        ip.bLooseTSACheck = loose; ip.nMode = mode; ip.bTautFlags = flags;
        STRUCT_DATA sd; memset(&sd, 0, sizeof(sd));
        StrFromINChI s; memset(&s, 0, sizeof(s));
        s.num_atoms = n; s.bMobileH = bmob; s.iMobileH = imob;
        ALL_TC_GROUPS g; memset(&g, 0, sizeof(g));
        printf("#CASE %d\n", c);
        int ret = MakeOneInChIOutOfStrFromINChI(&cg, &ic, &ip, &sd, &s, at2, at3, &g);
        printf("R %d %d %d %d %d\n", ret, s.nChargeRevrs, s.nOneINChI_bMobileH, s.nNumRemovedProtons, s.One_ti.num_t_groups);
        printf("TI %d", s.One_ti.nNumEndpoints);
        for (int i = 0; i < s.One_ti.num_t_groups && s.One_ti.t_group; i++)
            printf(" (%d %d %d %d %d)", s.One_ti.t_group[i].num[0], s.One_ti.t_group[i].num[1], s.One_ti.t_group[i].nGroupNumber, s.One_ti.t_group[i].nNumEndpoints, s.One_ti.t_group[i].nFirstEndpointAtNoPos);
        printf(" |");
        for (int i = 0; i < s.One_ti.nNumEndpoints && s.One_ti.nEndpointAtomNumber; i++) printf(" %d", s.One_ti.nEndpointAtomNumber[i]);
        printf("\n");
        for (int k = 0; k < 2; k++) {
            dump_one_inchi(k, s.pOneINChI[k]);
            INP_ATOM_DATA *d = s.pOne_norm_data[k];
            if (!d) { printf("D%d -\n", k); continue; }
            printf("D%d %d %d %d %d", k, d->num_at, d->num_removed_H, d->bTautomeric, d->bTautPreprocessed);
            for (int i = 0; d->at && i < d->num_at; i++) printf(" %d,%d,%d,%d,%d", d->at[i].charge, d->at[i].num_H, d->at[i].chem_bonds_valence, d->at[i].endpoint, d->at[i].valence);
            printf("\n");
        }
    }
    return 0;
}
"""

HARNESS = r2.HARNESS.split("int main(void)")[0] + MAIN


def main():
    rng = random.Random(777)
    cases = []
    for _ in range(3000):
        atoms = cr.gen(rng)
        flags = 0
        for f in cr.FLAGS:
            if rng.random() < 0.5:
                flags |= f
        cases.append((atoms, rng.randrange(2), rng.randrange(2), rng.choice(cr.MODES), flags, rng.randrange(2)))

    def fmt(c):
        atoms, bmob, imob, mode, flags, loose = c
        t = f"{len(atoms)} {bmob} {imob} {mode} {flags} {loose}\n"
        for a in atoms:
            t += f"{a['v']} {a['cbv']} {a['el']} {a['rad']} {a['chg']} {a['nh']} {a['iso'][0]} {a['iso'][1]} {a['iso'][2]} {a['aw']} " \
                 f"{a['xyz'][0]!r} {a['xyz'][1]!r} {a['xyz'][2]!r}\n"
            t += " ".join(f"{x} {t2}" for x, t2 in a["nb"]) + "\n"
        return t

    objs = ma.build_lib()
    data = []
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "h.c").write_text(HARNESS)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), *objs, "-lm"], check=True)
        for c in cases:
            r = subprocess.run([str(exe)], input=("1\n" + fmt(c)).encode(), capture_output=True)
            if r.returncode != 0:
                continue
            data.append({"in": fmt(c), "out": r.stdout.decode().strip().split("\n")[1:]})
    print(len(cases) - len(data), "cases dropped")
    path = ROOT / "tests/data/inchi_rvr_makeone.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
