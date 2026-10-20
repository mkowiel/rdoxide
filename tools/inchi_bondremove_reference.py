"""Generates tests/data/inchi_bondremove.json from the InChI C code:
`RemoveInpAtBond` and `DisconnectInpAtBond` of strutil.c on random atom arrays
with CML tetrahedral and stereobond parities. The function texts are cut out
of strutil.c and ichister.c. Needs a C compiler and the fetched InChI sources."""

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
@FUNCS@
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        int n;
        scanf("%d", &n);
        inp_ATOM *at = calloc(n, sizeof(inp_ATOM));
        AT_NUMB comp[16];
        for (int q = 0; q < 16; q++) comp[q] = 1;
        for (int i = 0; i < n; i++) {
            int v, cbv, pp, cmp;
            scanf("%d %d %d %d", &v, &cbv, &pp, &cmp);
            at[i].valence = v; at[i].chem_bonds_valence = cbv; at[i].p_parity = pp; at[i].component = cmp;
            at[i].orig_at_number = i + 1;
            for (int m = 0; m < 4; m++) { int x; scanf("%d", &x); at[i].p_orig_at_num[m] = x; }
            for (int j = 0; j < v; j++) { int x, b, s; scanf("%d %d %d", &x, &b, &s); at[i].neighbor[j] = x; at[i].bond_type[j] = b; at[i].bond_stereo[j] = s; }
            for (int m = 0; m < MAX_NUM_STEREO_BONDS; m++) {
                int a, b, p, o; scanf("%d %d %d %d", &a, &b, &p, &o);
                at[i].sb_ord[m] = a; at[i].sn_ord[m] = b; at[i].sb_parity[m] = p; at[i].sn_orig_at_num[m] = o;
            }
        }
        int nr; scanf("%d", &nr);
        for (int r = 0; r < nr; r++) {
            int mode, iat, k;
            scanf("%d %d %d", &mode, &iat, &k);
            int ret = mode ? DisconnectInpAtBond(at, comp, iat, k) : RemoveInpAtBond(at, iat, k);
            printf("R %d", ret);
            for (int q = 0; q < 3; q++) printf(" %d", comp[q]);
            printf("\n");
        }
        for (int i = 0; i < n; i++) {
            printf("A %d %d %d", at[i].valence, at[i].chem_bonds_valence, at[i].p_parity);
            for (int m = 0; m < 4; m++) printf(" %d", at[i].p_orig_at_num[m]);
            for (int k = 0; k < at[i].valence; k++) printf(" %d:%d:%d", at[i].neighbor[k], at[i].bond_type[k], at[i].bond_stereo[k]);
            for (int m = 0; m < MAX_NUM_STEREO_BONDS; m++) printf(" S%d:%d:%d:%d", at[i].sb_ord[m], at[i].sn_ord[m], at[i].sb_parity[m], at[i].sn_orig_at_num[m]);
            printf("\n");
        }
        free(at);
    }
    return 0;
}
"""


def gen(rng):
    n = rng.randrange(2, 9)
    nb = [[] for _ in range(n)]
    for _ in range(rng.randrange(1, 12)):
        a, b = rng.randrange(n), rng.randrange(n)
        if a != b and all(x != b for x, _, _ in nb[a]) and len(nb[a]) < 5 and len(nb[b]) < 5:
            t = rng.choice([1, 1, 2, 2, 3, 4])
            st = rng.choice([0, 0, 1, -1])
            nb[a].append((b, t, st))
            nb[b].append((a, t, -st))
    atoms = []
    for i in range(n):
        atoms.append({"nb": nb[i], "pp": 0, "po": [0, 0, 0, 0], "cmp": rng.choice([0, 1, 2, 3]), "sb": []})
    for i in range(n):
        v = len(nb[i])
        if 3 <= v <= 4 and rng.random() < 0.5:
            atoms[i]["pp"] = rng.choice([1, 2, 3, 4])
            ids = [nb[i][k][0] + 1 for k in range(v)]
            if v == 3:
                ids.insert(rng.randrange(4), i + 1 if rng.random() < 0.6 else rng.randrange(1, n + 1))
            if rng.random() < 0.1:
                ids[0] = rng.randrange(1, n + 1)
            atoms[i]["po"] = ids[:4] + [0] * (4 - len(ids[:4]))
    # stereobonds: pairs of atoms with a parity each
    for i in range(n):
        for k, (j, t, _) in enumerate(nb[i]):
            if i < j and t == 2 and rng.random() < 0.8 and len(atoms[i]["sb"]) < 3 and len(atoms[j]["sb"]) < 3:
                kj = [x for x, _, _ in nb[j]].index(i)
                for (a, ko, kk) in ((i, k, kj), (j, kj, k)):
                    others = [q for q in range(len(nb[a])) if q != ko]
                    if others and rng.random() < 0.9:
                        sn = rng.choice(others)
                        so = nb[nb[a][sn][0]][0][0] + 1 if False else nb[a][sn][0] + 1
                    else:
                        sn, so = -1, 0
                    par = rng.choice([1, 2, 1, 2, 3, 4])
                    atoms[a]["sb"].append([ko, sn, par, so])
    for a in atoms:
        a["cbv"] = sum(t if t <= 3 else 1 for _, t, _ in a["nb"])
    return atoms


def main():
    rng = random.Random(5)
    cases = []
    for _ in range(8000):
        atoms = gen(rng)
        rem = []
        for _ in range(rng.randrange(1, 4)):
            i = rng.randrange(len(atoms))
            if atoms[i]["nb"]:
                rem.append([rng.choice([0, 1, 1]), i, rng.randrange(len(atoms[i]["nb"]))])
        cases.append({"atoms": atoms, "rem": rem})
    # the removal sequence is applied on a changing graph: keep only valid positions
    # (positions are re-checked in the harness input by simulating nothing; invalid ones are dropped by the
    # generator through the rule that an atom's valence only shrinks, so indexes are clamped here)
    for c in cases:
        val = [len(a["nb"]) for a in c["atoms"]]
        ok = []
        for mode, i, k in c["rem"]:
            if k < val[i]:
                ok.append([mode, i, k])
                val[i] -= 1
                if mode:
                    val_j = None
                    # the opposite atom shrinks too; make later removals conservative
                    j = c["atoms"][i]["nb"][k][0]
                    val[j] = max(0, val[j] - 1)
        # only the first removal is guaranteed valid on the original indexes
        c["rem"] = ok[:1]
    stdin = f"{len(cases)}\n"
    for c in cases:
        atoms = c["atoms"]
        stdin += f"{len(atoms)}\n"
        for a in atoms:
            stdin += f"{len(a['nb'])} {a['cbv']} {a['pp']} {a['cmp']} " + " ".join(map(str, a["po"])) + " " + " ".join(f"{x} {t} {s}" for x, t, s in a["nb"]) + " "
            sb = a["sb"][:3] + [[0, 0, 0, 0]] * (3 - len(a["sb"][:3]))
            stdin += " ".join(" ".join(map(str, e)) for e in sb) + "\n"
        stdin += f"{len(c['rem'])}\n" + "".join(f"{m} {i} {k}\n" for m, i, k in c["rem"])
    text = (SRC / "strutil.c").read_text(errors="replace")
    ister = (SRC / "ichister.c").read_text(errors="replace")

    def cut(t, start):
        a = t.rindex(start)
        return t[a:t.index("\n}\n", a) + 3]

    funcs = [cut(ister, "int get_opposite_sb_atom( inp_ATOM *at,\n"),
             cut(text, "int RemoveInpAtBond( inp_ATOM *atom, int iat, int k )\n{"),
             cut(text, "int DisconnectInpAtBond( inp_ATOM *at,\n")]
    protos = """int get_opposite_sb_atom( inp_ATOM *at, int cur_atom, int icur2nxt, int *pnxt_atom, int *pinxt2cur, int *pinxt_sb_parity_ord );
int RemoveInpAtBond( inp_ATOM *atom, int iat, int k );
int DisconnectInpAtBond( inp_ATOM *at, AT_NUMB *nOldCompNumber, int iat, int neigh_ord );
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
        cnt = len(c["rem"]) + n
        data.append({"in": c, "out": out[k:k + cnt]})
        k += cnt
    path = ROOT / "tests/data/inchi_bondremove.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


main()
