"""Generates tests/data/inchi_cgroups.json from the InChI C code: the charge
group and salt charge type helpers of ichitaut.c (`RegisterCPoints`,
`CmpCCandidates`, `GetSaltChargeType`, `GetOtherSaltChargeType`,
`GetOtherSaltType`, `bDoNotMergeNonTautAtom`, `comp_candidates`) on random
atom arrays, t-group lists and registration sequences. The function texts
are cut out of ichitaut.c. Needs a C compiler and the fetched InChI sources."""

import json
import random
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "upstream/external/InChI/INCHI_BASE/src"

HARNESS = r"""
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mode.h"
#include "incomdef.h"
#include "ichidrp.h"
#include "inpdef.h"
#include "ichi.h"
#include "strutil.h"
#include "util.h"
#include "ichierr.h"
#include "extr_ct.h"
#include "ichitaut.h"
@FUNCS@
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        int n;
        scanf("%d", &n);
        inp_ATOM *at = calloc(n, sizeof(inp_ATOM));
        for (int i = 0; i < n; i++) {
            int v, cbv, el, chg, rad, nh, cp, ep, ring; char en[8];
            scanf("%d %d %d %d %d %d %d %d %d %s", &v, &cbv, &el, &chg, &rad, &nh, &cp, &ep, &ring, en);
            at[i].valence = v; at[i].chem_bonds_valence = cbv; at[i].el_number = el; at[i].charge = chg;
            at[i].radical = rad; at[i].num_H = nh; at[i].c_point = cp; at[i].endpoint = ep; at[i].nNumAtInRingSystem = ring;
            strcpy(at[i].elname, en);
            for (int j = 0; j < v; j++) { int q, b; scanf("%d %d", &q, &b); at[i].neighbor[j] = q; at[i].bond_type[j] = b; }
        }
        int ntg, hastg;
        scanf("%d %d", &hastg, &ntg);
        T_GROUP tg[8]; memset(tg, 0, sizeof(tg));
        for (int g = 0; g < ntg; g++) { int a, b, num; scanf("%d %d %d", &num, &a, &b); tg[g].nGroupNumber = num; tg[g].num[0] = a; tg[g].num[1] = b; }
        T_GROUP_INFO tgi; memset(&tgi, 0, sizeof(tgi));
        tgi.t_group = hastg ? tg : NULL; tgi.num_t_groups = ntg;
        for (int i = 0; i < n; i++) {
            int s1, s2, s3, s4;
            int r1 = GetSaltChargeType(at, i, hastg == 2 ? NULL : &tgi, &s1);
            int r2 = GetOtherSaltChargeType(at, i, hastg == 2 ? NULL : &tgi, &s2, 0);
            int r3 = GetOtherSaltChargeType(at, i, hastg == 2 ? NULL : &tgi, &s3, 1);
            s4 = 0;
            int r4 = GetOtherSaltType(at, i, &s4);
            printf("S %d %d %d %d %d %d %d %d %d\n", r1, s1, r2, s2, r3, s3, r4, s4, bDoNotMergeNonTautAtom(at, i));
        }
        /* RegisterCPoints sequence */
        int nreg, maxc;
        scanf("%d %d", &maxc, &nreg);
        C_GROUP cg[8]; memset(cg, 0, sizeof(cg));
        int numc = 0;
        for (int r = 0; r < nreg; r++) {
            int p1, p2, ct; scanf("%d %d %d", &p1, &p2, &ct);
            int ret = RegisterCPoints(cg, &numc, maxc, NULL, p1, p2, ct, at, n);
            printf("C %d %d |", ret, numc);
            for (int g = 0; g < numc; g++) printf(" %d,%d,%d,%d,%d", cg[g].num[0], cg[g].num[1], cg[g].num_CPoints, cg[g].nGroupNumber, cg[g].cGroupType);
            printf(" |");
            for (int i = 0; i < n; i++) printf(" %d", at[i].c_point);
            printf("\n");
        }
        /* sorting */
        int nc; scanf("%d", &nc);
        C_CANDIDATE cc[16];
        for (int i = 0; i < nc; i++) { int a, t, s; scanf("%d %d %d", &a, &t, &s); cc[i].atnumber = a; cc[i].type = t; cc[i].subtype = s; }
        qsort(cc, nc, sizeof(cc[0]), CmpCCandidates);
        printf("Q");
        for (int i = 0; i < nc; i++) printf(" %d,%d,%d", cc[i].atnumber, cc[i].type, cc[i].subtype);
        printf("\n");
        int ns; scanf("%d", &ns);
        S_CANDIDATE sc[16];
        for (int i = 0; i < ns; i++) { int a, t, s, e; scanf("%d %d %d %d", &a, &t, &s, &e); sc[i].atnumber = a; sc[i].type = t; sc[i].subtype = s; sc[i].endpoint = e; }
        qsort(sc, ns, sizeof(sc[0]), comp_candidates);
        printf("T");
        for (int i = 0; i < ns; i++) printf(" %d,%d,%d,%d", sc[i].atnumber, sc[i].type, sc[i].subtype, sc[i].endpoint);
        printf("\n");
        free(at);
    }
    return 0;
}
"""

SYM = {6: "C", 7: "N", 8: "O", 15: "P", 16: "S", 34: "Se", 52: "Te", 17: "Cl", 1: "H"}
ELS = [6, 6, 7, 7, 8, 8, 8, 8, 16, 34, 52, 15, 17, 1]


def gen(rng):
    n = rng.randrange(2, 7)
    nb = [[] for _ in range(n)]
    for i in range(1, n):
        if rng.random() < 0.92:
            j = rng.randrange(i)
            t = rng.choice([1, 1, 2, 2])
            nb[i].append((j, t))
            nb[j].append((i, t))
    atoms = []
    for i in range(n):
        el = rng.choice(ELS)
        v = len(nb[i])
        cbv = sum(t for _, t in nb[i])
        if el == 6 and rng.random() < 0.7:
            cbv = 4 if v else cbv
        atoms.append({"el": el, "sym": SYM[el], "v": v, "cbv": cbv, "chg": rng.choice([0, 0, 0, 1, -1, -1]),
                      "rad": rng.choice([0] * 14 + [1, 2]), "nh": rng.choice([0, 0, 1, 1, 2]),
                      "cp": rng.choice([0, 0, 0, 1]), "ep": rng.choice([0, 0, 0, 1, 2, 3]),
                      "ring": rng.choice([0, 5, 6]), "nb": nb[i]})
        a = atoms[-1]
        if el in (8, 16, 34, 52) and v == 1 and rng.random() < 0.6:
            a["cbv"] = 1 if a["chg"] == -1 else rng.choice([1, 2])
            a["nh"] = 1 if (a["chg"] == 0 and a["cbv"] == 1) else 0
    return atoms


def main():
    rng = random.Random(77)
    cases = []
    for _ in range(12000):
        atoms = gen(rng)
        n = len(atoms)
        hastg = rng.choice([0, 1, 1, 2])
        tgs = [[g + 1, rng.randrange(0, 5), rng.randrange(0, 3)] for g in range(rng.randrange(0, 4))]
        maxc = rng.randrange(1, 4)
        regs = []
        for _ in range(rng.randrange(0, 7)):
            a, b = rng.randrange(n), rng.randrange(n)
            if a != b:
                regs.append([a, b, rng.randrange(0, 3)])
        nc = rng.randrange(0, 8)
        cands = [[i, rng.choice([0, 1, 2]), rng.choice([0, 1, 4, 5, 16, 20])] for i in rng.sample(range(40), nc)]
        ns = rng.randrange(0, 8)
        scands = [[i, rng.choice([-1, 0, 1, 2]), rng.choice([0, 1, 2, 4]), rng.choice([0, 0, 1, 2, 3])] for i in rng.sample(range(40), ns)]
        cases.append({"atoms": atoms, "hastg": hastg, "tgs": tgs, "maxc": maxc, "regs": regs, "cands": cands, "scands": scands})
    stdin = f"{len(cases)}\n"
    for c in cases:
        stdin += f"{len(c['atoms'])}\n"
        for a in c["atoms"]:
            stdin += (f"{a['v']} {a['cbv']} {a['el']} {a['chg']} {a['rad']} {a['nh']} {a['cp']} {a['ep']} {a['ring']} {a['sym']} "
                      + " ".join(f"{x} {t}" for x, t in a["nb"]) + "\n")
        stdin += f"{c['hastg']} {len(c['tgs'])}\n" + "".join(f"{g[0]} {g[1]} {g[2]}\n" for g in c["tgs"])
        stdin += f"{c['maxc']} {len(c['regs'])}\n" + "".join(f"{r[0]} {r[1]} {r[2]}\n" for r in c["regs"])
        stdin += f"{len(c['cands'])}\n" + "".join(f"{x[0]} {x[1]} {x[2]}\n" for x in c["cands"])
        stdin += f"{len(c['scands'])}\n" + "".join(f"{x[0]} {x[1]} {x[2]} {x[3]}\n" for x in c["scands"])
    text = (SRC / "ichitaut.c").read_text(errors="replace")

    def cut(start):
        a = text.rindex(start)
        return text[a:text.index("\n}\n", a) + 3]

    ctype_start = text.index("typedef struct tagChargeType")
    ctype_end = text.index("#define NUM_C_TYPES")
    funcs = [text[ctype_start:ctype_end], "#define NUM_C_TYPES  (int)(sizeof( CType )/sizeof(CType[0]))\n",
             cut("int is_centerpoint_elem( U_CHAR el_number )\n{"),
             cut("int bCanBeACPoint( inp_ATOM *at,\n"),
             cut("int GetChargeType( inp_ATOM *atom, int iat, S_CHAR *cChargeSubtype )\n{"),
             cut("int nGetEndpointInfo( inp_ATOM *atom, int iat, ENDPOINT_INFO *eif )\n{"),
             cut("int CmpCCandidates( const void *a1, const void *a2 )\n{"),
             cut("int RegisterCPoints( C_GROUP *c_group,\n"),
             cut("int GetSaltChargeType( inp_ATOM *at,\n"),
             cut("int bDoNotMergeNonTautAtom( inp_ATOM *at, int at_no )\n{"),
             cut("int GetOtherSaltChargeType( inp_ATOM *at,\n"),
             cut("int GetOtherSaltType( inp_ATOM *at, int at_no, int *s_subtype )\n{"),
             cut("int comp_candidates( const void *a1, const void *a2 )\n{")]
    protos = """int is_centerpoint_elem( U_CHAR el_number );
int GetChargeType( inp_ATOM *atom, int iat, S_CHAR *cChargeSubtype );
int bCanBeACPoint( inp_ATOM *at, S_CHAR cCharge, S_CHAR cChangeValence, S_CHAR neutral_bonds_valence, S_CHAR neutral_valence, S_CHAR nEndpointValence, S_CHAR *cChargeSubtype );
"""
    with tempfile.TemporaryDirectory() as d:
        exe = Path(d) / "h"
        (Path(d) / "h.c").write_text(HARNESS.replace("@FUNCS@", protos + "\n".join(funcs)))
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), str(SRC / "util.c"), str(SRC / "ichierr.c")], check=True)
        out = subprocess.run([str(exe)], input=stdin.encode(), capture_output=True, check=True).stdout.decode().strip().split("\n")
    k = 0
    data = []
    for c in cases:
        cnt = len(c["atoms"]) + len(c["regs"]) + 2
        data.append({"in": c, "out": out[k:k + cnt]})
        k += cnt
    path = ROOT / "tests/data/inchi_cgroups.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


main()
