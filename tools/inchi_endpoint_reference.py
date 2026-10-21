"""Generates tests/data/inchi_endpoint.json from the InChI C code: the
tautomer endpoint helpers of ichitaut.c (`is_centerpoint_elem*`,
`AddAtom2num`, `AddAtom2DA`, `AddEndPoint`, `GetChargeType`,
`bCanBeACPoint`, `nGetEndpointInfo` and its PT_22_00/16_00/06_00/39_00/13_00/
18_00/KET variants) on random small atom arrays. The function texts are cut
out of ichitaut.c. Needs a C compiler and the fetched InChI sources."""

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
static void pe(ENDPOINT_INFO *e, int r) {
    if (r > 0) printf(" %d:%d,%d,%d,%d,%d,%d,%d", r, e->cMoveableCharge, e->cNeutralBondsValence, e->cMobile, e->cDonor, e->cAcceptor, e->cKetoEnolCode, 0);
    else printf(" 0");
}
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        int n;
        scanf("%d", &n);
        inp_ATOM *at = calloc(n, sizeof(inp_ATOM));
        for (int i = 0; i < n; i++) {
            int v, cbv, el, chg, rad, nh, ih0, ih1, ih2, cp, ep, att, ring; char en[8];
            scanf("%d %d %d %d %d %d %d %d %d %d %d %d %d %s", &v, &cbv, &el, &chg, &rad, &nh, &ih0, &ih1, &ih2, &cp, &ep, &att, &ring, en);
            at[i].valence = v; at[i].chem_bonds_valence = cbv; at[i].el_number = el; at[i].charge = chg;
            at[i].radical = rad; at[i].num_H = nh; at[i].c_point = cp; at[i].endpoint = ep; at[i].at_type = att;
            at[i].num_iso_H[0] = ih0; at[i].num_iso_H[1] = ih1; at[i].num_iso_H[2] = ih2; at[i].nNumAtInRingSystem = ring;
            strcpy(at[i].elname, en);
            for (int j = 0; j < v; j++) { int q, b; scanf("%d %d", &q, &b); at[i].neighbor[j] = q; at[i].bond_type[j] = b; }
        }
        for (int i = 0; i < n; i++) {
            S_CHAR sub = 99;
            int ct = GetChargeType(at, i, &sub);
            printf("G %d %d %d %d %d %d", ct, sub, is_centerpoint_elem(at[i].el_number), is_centerpoint_elem_KET(at[i].el_number), is_centerpoint_elem_strict(at[i].el_number), get_endpoint_valence_KET(at[i].el_number));
            ENDPOINT_INFO e;
            int r;
            memset(&e, 0, sizeof(e)); r = nGetEndpointInfo(at, i, &e); pe(&e, r);
            memset(&e, 0, sizeof(e)); r = nGetEndpointInfo_PT_22_00(at, i, &e); pe(&e, r);
            memset(&e, 0, sizeof(e)); r = nGetEndpointInfo_PT_16_00(at, i, &e); pe(&e, r);
            memset(&e, 0, sizeof(e)); r = nGetEndpointInfo_PT_06_00(at, i, &e); pe(&e, r);
            memset(&e, 0, sizeof(e)); r = nGetEndpointInfo_PT_39_00(at, i, &e); pe(&e, r);
            memset(&e, 0, sizeof(e)); r = nGetEndpointInfo_PT_13_00(at, i, &e); pe(&e, r);
            memset(&e, 0, sizeof(e)); r = nGetEndpointInfo_PT_18_00(at, i, &e); pe(&e, r);
            memset(&e, 0, sizeof(e)); r = nGetEndpointInfo_KET(at, i, &e); pe(&e, r);
            printf("\n");
            T_ENDPOINT ep; memset(&ep, 0, sizeof(ep));
            AddEndPoint(&ep, at, i);
            printf("E %d %d %d", ep.nAtomNumber, ep.nEquNumber, ep.nGroupNumber);
            for (int k = 0; k < T_NUM_NO_ISOTOPIC + T_NUM_ISOTOPIC; k++) printf(" %d", ep.num[k]);
            for (int k = 0; k < TG_NUM_DA; k++) printf(" %d", ep.num_DA[k]);
            AT_RANK num[T_NUM_NO_ISOTOPIC + T_NUM_ISOTOPIC] = {9, 8, 7, 6, 5};
            AT_RANK da[TG_NUM_DA] = {9, 8, 7, 6, 5, 4};
            for (int mode = 0; mode < 3; mode++) {
                int m = AddAtom2num(num, at, i, mode);
                AddAtom2DA(da, at, i, mode);
                printf(" |%d", m);
                for (int k = 0; k < T_NUM_NO_ISOTOPIC + T_NUM_ISOTOPIC; k++) printf(" %d", num[k]);
                for (int k = 0; k < TG_NUM_DA; k++) printf(" %d", da[k]);
                if (mode == 1) { for (int k = 0; k < 5; k++) num[k] = 9 - k; for (int k = 0; k < 6; k++) da[k] = 9 - k; }
            }
            printf("\n");
        }
        free(at);
    }
    return 0;
}
"""

SYM = {6: "C", 7: "N", 8: "O", 15: "P", 16: "S", 34: "Se", 52: "Te", 33: "As", 51: "Sb", 17: "Cl", 35: "Br", 53: "I", 1: "H", 9: "F"}
ELS = [6, 6, 7, 7, 8, 8, 8, 15, 16, 34, 52, 33, 51, 17, 35, 53, 1, 9]


def gen(rng):
    n = rng.randrange(1, 6)
    nb = [[] for _ in range(n)]
    for i in range(1, n):
        if rng.random() < 0.9:
            j = rng.randrange(i)
            t = rng.choice([1, 1, 2, 2, 3])
            nb[i].append((j, t))
            nb[j].append((i, t))
    atoms = []
    for i in range(n):
        el = rng.choice(ELS)
        v = len(nb[i])
        cbv = sum(t for _, t in nb[i])
        if rng.random() < 0.3:
            cbv += rng.choice([0, 1])
        chg = rng.choice([0, 0, 0, 1, 1, -1, -1, 2])
        nh = rng.choice([0, 0, 1, 1, 2, 3])
        atoms.append({"el": el, "sym": SYM[el], "v": v, "cbv": cbv, "chg": chg, "rad": rng.choice([0] * 12 + [1, 2]),
                      "nh": nh, "ih": [rng.choice([0, 0, 0, 1]) for _ in range(3)],
                      "cp": rng.choice([0, 0, 1, 1]), "ep": rng.choice([0, 0, 0, 1, 2]),
                      "att": rng.choice([0, 0, 1, 1, 2]), "ring": rng.choice([0, 1, 5, 6]), "nb": nb[i]})
    return atoms


def main():
    rng = random.Random(61)
    cases = [{"atoms": gen(rng)} for _ in range(12000)]
    stdin = f"{len(cases)}\n"
    for c in cases:
        stdin += f"{len(c['atoms'])}\n"
        for a in c["atoms"]:
            stdin += (f"{a['v']} {a['cbv']} {a['el']} {a['chg']} {a['rad']} {a['nh']} {a['ih'][0]} {a['ih'][1]} {a['ih'][2]} "
                      f"{a['cp']} {a['ep']} {a['att']} {a['ring']} {a['sym']} " + " ".join(f"{x} {t}" for x, t in a["nb"]) + "\n")
    text = (SRC / "ichitaut.c").read_text(errors="replace")
    util = (SRC / "util.c").read_text(errors="replace")

    def cut(start):
        a = text.rindex(start)
        return text[a:text.index("\n}\n", a) + 3]

    ctype_start = text.index("typedef struct tagChargeType")
    ctype_end = text.index("#define NUM_C_TYPES")
    funcs = [text[ctype_start:ctype_end],
             cut("int is_centerpoint_elem( U_CHAR el_number )\n{"),
             cut("int is_centerpoint_elem_KET( U_CHAR el_number )\n{"),
             cut("int is_centerpoint_elem_strict( U_CHAR el_number )\n{"),
             "#define NUM_C_TYPES  (int)(sizeof( CType )/sizeof(CType[0]))\n",
             cut("int AddAtom2num( AT_RANK num[], inp_ATOM *atom, int at_no, int bSubtract )\n{"),
             cut("void AddAtom2DA( AT_RANK num_DA[], inp_ATOM *atom, int at_no, int bSubtract )\n{"),
             cut("int AddEndPoint( T_ENDPOINT *pEndPoint, inp_ATOM *at, int iat )\n{"),
             cut("int bCanBeACPoint( inp_ATOM *at,\n"),
             cut("int GetChargeType( inp_ATOM *atom, int iat, S_CHAR *cChargeSubtype )\n{")]
    for v in ["", "_PT_22_00", "_PT_16_00", "_PT_06_00", "_PT_39_00", "_PT_13_00", "_PT_18_00"]:
        sig = "int nGetEndpointInfo%s(inp_ATOM *atom, int iat, ENDPOINT_INFO *eif)" % v if v else "int nGetEndpointInfo( inp_ATOM *atom, int iat, ENDPOINT_INFO *eif )"
        funcs.append(cut(sig))
    funcs.append(cut("int nGetEndpointInfo_KET( inp_ATOM *atom,\n"))
    protos = """int is_centerpoint_elem( U_CHAR el_number );
int is_centerpoint_elem_KET( U_CHAR el_number );
int is_centerpoint_elem_strict( U_CHAR el_number );
int GetChargeType( inp_ATOM *atom, int iat, S_CHAR *cChargeSubtype );
int bCanBeACPoint( inp_ATOM *at, S_CHAR cCharge, S_CHAR cChangeValence, S_CHAR neutral_bonds_valence, S_CHAR neutral_valence, S_CHAR nEndpointValence, S_CHAR *cChargeSubtype );
int nGetEndpointInfo_PT_22_00(inp_ATOM *atom, int iat, ENDPOINT_INFO *eif);
int nGetEndpointInfo_PT_16_00(inp_ATOM *atom, int iat, ENDPOINT_INFO *eif);
int nGetEndpointInfo_PT_06_00(inp_ATOM *atom, int iat, ENDPOINT_INFO *eif);
int nGetEndpointInfo_PT_39_00(inp_ATOM *atom, int iat, ENDPOINT_INFO *eif);
int nGetEndpointInfo_PT_13_00(inp_ATOM *atom, int iat, ENDPOINT_INFO *eif);
int nGetEndpointInfo_PT_18_00(inp_ATOM *atom, int iat, ENDPOINT_INFO *eif);
int nGetEndpointInfo_KET( inp_ATOM *atom, int iat, ENDPOINT_INFO *eif );
"""
    with tempfile.TemporaryDirectory() as d:
        exe = Path(d) / "h"
        (Path(d) / "h.c").write_text(HARNESS.replace("@FUNCS@", protos + "\n".join(funcs)))
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), str(SRC / "util.c"), str(SRC / "ichierr.c")], check=True)
        out = subprocess.run([str(exe)], input=stdin.encode(), capture_output=True, check=True).stdout.decode().strip().split("\n")
    k = 0
    data = []
    for c in cases:
        n = len(c["atoms"])
        data.append({"in": c, "out": out[k:k + 2 * n]})
        k += 2 * n
    path = ROOT / "tests/data/inchi_endpoint.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


main()
