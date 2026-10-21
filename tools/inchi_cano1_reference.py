"""Generates tests/data/inchi_cano1.json from the InChI C code: the linear
connection table helpers of ichicano.c (`GetCanonLengths`,
`FillIsotopicAtLinearCT`, `FillTautLinearCT2`, `UpdateFullLinearCT`),
compiled with the upstream source included.

    python tools/inchi_cano1_reference.py
"""
import importlib.util
import json
import random
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / f"tools/{name}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


ma = load("inchi_markalt_reference")
d2 = load("inchi_can2d_reference")
SRC = ma.SRC

HARNESS = r"""
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ichicano.c"
int CreateCompositeNormAtom() { return 0; }
static void pi_(const char *t, int *a, int n) { printf("%s", t); for (int i = 0; i < n; i++) printf(" %d", a[i]); printf("\n"); }
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        CANON_GLOBALS cg;
        memset(&cg, 0, sizeof(cg));
        int n, m, ntg, nE;
        scanf("%d %d %d %d", &n, &m, &ntg, &nE);
        int ntot = n + ntg;
        sp_ATOM *at = calloc(n + 1, sizeof(sp_ATOM));
        int edges[200];
        for (int i = 0; i < 2 * m; i++) scanf("%d", &edges[i]);
        for (int i = 0; i < n; i++) {
            long isk; int par, sb0, sb1, sb2, ep, h0, h1, h2, aw, cf;
            scanf("%ld %d %d %d %d %d %d %d %d %d %d", &isk, &par, &sb0, &sb1, &sb2, &ep, &h0, &h1, &h2, &aw, &cf);
            at[i].iso_sort_key = isk; at[i].parity = par;
            at[i].stereo_bond_neighbor[0] = sb0; at[i].stereo_bond_neighbor[1] = sb1; at[i].stereo_bond_neighbor[2] = sb2;
            at[i].endpoint = ep; at[i].num_iso_H[0] = h0; at[i].num_iso_H[1] = h1; at[i].num_iso_H[2] = h2;
            at[i].iso_atw_diff = aw; at[i].cFlags = cf;
        }
        for (int i = 0; i < m; i++) {
            int a = edges[2 * i], b = edges[2 * i + 1];
            at[a].neighbor[at[a].valence++] = b;
            at[b].neighbor[at[b].valence++] = a;
        }
        T_GROUP *tg = calloc(ntg + 2, sizeof(T_GROUP));
        AT_NUMB *epl = calloc(nE + 2, sizeof(AT_NUMB));
        for (int k = 0; k < ntg; k++) {
            int ne, fp; long iw;
            scanf("%d %d", &ne, &fp);
            tg[k].nNumEndpoints = ne; tg[k].nFirstEndpointAtNoPos = fp;
            for (int j = 0; j < 5; j++) { int x; scanf("%d", &x); tg[k].num[j] = x; }
            scanf("%ld", &iw); tg[k].iWeight = iw;
        }
        for (int i = 0; i < nE; i++) { int x; scanf("%d", &x); epl[i] = x; }
        int nIsoEnd, sv[8];
        scanf("%d", &nIsoEnd);
        for (int i = 0; i < 8; i++) scanf("%d", &sv[i]);
        T_GROUP_INFO tgi; memset(&tgi, 0, sizeof(tgi));
        tgi.t_group = ntg ? tg : NULL; tgi.nEndpointAtomNumber = epl; tgi.num_t_groups = ntg; tgi.nNumEndpoints = nE;
        tgi.nNumIsotopicEndpoints = nIsoEnd;
        tgi.tGroupNumber = calloc(4 * ntg + 4, sizeof(AT_NUMB));
        int hasTgi; scanf("%d", &hasTgi);
        printf("N %d\n", c);
        /* GetCanonLengths */
        {
            ATOM_SIZES s; memset(&s, 0, sizeof(s));
            s.nLenBonds = sv[0]; s.nLenCTAtOnly = sv[1]; s.nLenCT = sv[2]; s.nLenIsotopic = sv[3];
            s.nLenLinearCTStereoDble = sv[4]; s.nLenLinearCTStereoCarb = sv[5]; s.nLenIsotopicEndpoints = sv[6];
            s.nLenLinearCTTautomer = sv[7];
            int r = GetCanonLengths(n, at, &s, hasTgi ? &tgi : NULL);
            printf("L %d %d %d %d %d %d %d %d\n", r, s.nLenBonds, s.nLenCTAtOnly, s.nLenCT, s.nLenIsotopic, s.nLenLinearCTStereoDble, s.nLenLinearCTStereoCarb, s.nLenIsotopicEndpoints);
        }
        /* rank vectors */
        AT_RANK *anum[3], *rnk[3], *symm[3];
        for (int v = 0; v < 3; v++) {
            anum[v] = calloc(ntot + 2, 2); rnk[v] = calloc(ntot + 2, 2); symm[v] = calloc(ntot + 2, 2);
            for (int i = 0; i < ntot; i++) { int x; scanf("%d", &x); anum[v][i] = x; }
            for (int i = 0; i < ntot; i++) { int x; scanf("%d", &x); symm[v][i] = x; }
            for (int i = 0; i < ntot; i++) rnk[v][anum[v][i]] = i + 1;
        }
        int nmax, plen, nmaxIsoAt, plenIsoAt, nmaxT, plenT, nmaxIT, plenIT, bIso;
        scanf("%d %d %d %d %d %d %d %d %d", &nmax, &plen, &nmaxIsoAt, &plenIsoAt, &nmaxT, &plenT, &nmaxIT, &plenIT, &bIso);
        /* FillIsotopicAtLinearCT */
        {
            AT_ISOTOPIC *iso = calloc(nmaxIsoAt + 2, sizeof(AT_ISOTOPIC));
            int len = plenIsoAt;
            int r = FillIsotopicAtLinearCT(n, at, anum[0], nmaxIsoAt ? iso : NULL, nmaxIsoAt, &len);
            printf("I %d %d |", r, len);
            for (int i = 0; i < nmaxIsoAt; i++) printf(" %d,%d,%d,%d,%d", iso[i].at_num, iso[i].num_1H, iso[i].num_D, iso[i].num_T, iso[i].iso_atw_diff);
            printf("\n");
        }
        /* FillTautLinearCT2 */
        if (hasTgi) {
            AT_TAUTOMER *lt = calloc(nmaxT + 2, sizeof(AT_TAUTOMER));
            AT_ISO_TGROUP *lit = calloc(nmaxIT + 2, sizeof(AT_ISO_TGROUP));
            int lenT = plenT, lenIT = plenIT;
            int r = FillTautLinearCT2(&cg, n, ntot, bIso, rnk[0], anum[0], symm[0], rnk[1], anum[1], symm[1],
                                      lt, nmaxT, &lenT, lit, nmaxIT, &lenIT, &tgi);
            printf("T %d %d %d |", r, lenT, lenIT);
            for (int i = 0; i < nmaxT; i++) printf(" %d", lt[i]);
            printf(" |");
            for (int i = 0; i < lenIT && i < nmaxIT; i++) printf(" %d:%d,%d,%d", lit[i].tgroup_num, lit[i].num[0], lit[i].num[1], lit[i].num[2]);
            printf(" |");
            for (int i = 0; i < 4 * ntg; i++) printf(" %d", tgi.tGroupNumber[i]);
            printf(" |");
            for (int i = 0; i < nE; i++) printf(" %d", epl[i]);
            printf("\n");
        }
        /* UpdateFullLinearCT */
        {
            CANON_STAT cs; memset(&cs, 0, sizeof(cs));
            int maxlen = nmax;
            cs.nMaxLenLinearCT = maxlen;
            cs.LinearCT = calloc(maxlen + 2, sizeof(AT_NUMB));
            cs.t_group_info = hasTgi ? &tgi : NULL;
            cs.nLenLinearCT = plen;
            int tot = ntot;
            int nn = hasTgi ? tot : n;
            for (int v = 0; v < 3; v++) {
                int r = UpdateFullLinearCT(n, nn, at, rnk[v], anum[v], &cs, &cg, v == 0);
                printf("U%d %d %d %d |", v, r, cs.nLenLinearCT, cs.nLenLinearCTAtOnly);
                for (int i = 0; i < maxlen; i++) printf(" %d", cs.LinearCT[i]);
                printf("\n");
            }
            printf("UE");
            for (int i = 0; i < nE; i++) printf(" %d", epl[i]);
            printf("\n");
        }
    }
    return 0;
}
"""


def perm_ranks(rng, n, ntg):
    a = list(range(n))
    rng.shuffle(a)
    t = list(range(n, n + ntg))
    rng.shuffle(t)
    anum = a + t
    # symm ranks: class representative style ranks (arbitrary but in range)
    symm = []
    tot = n + ntg
    for i in range(tot):
        symm.append(rng.randrange(1, tot + 1) if rng.random() < 0.5 else (anum.index(i) + 1))
    return anum, symm


def gen(rng):
    name = rng.choice(d2.NAMES)
    nn, e = d2.edges_of(name, rng)
    if nn > 14:
        nn, e = 6, [(i, (i + 1) % 6) for i in range(6)]
    n = max(nn, 1)
    perm = list(range(n))
    rng.shuffle(perm)
    e = [(perm[a], perm[b]) for a, b in e]
    m = len(e)
    ntg = rng.choice([0, 0, 1, 2, 3])
    avail = list(range(n))
    rng.shuffle(avail)
    groups = []
    for _ in range(ntg):
        k = rng.randrange(1, 4)
        if len(avail) < k:
            break
        groups.append(sorted(avail[:k]))
        avail = avail[k:]
    ntg = len(groups)
    ep = [0] * n
    for g, mem in enumerate(groups):
        for a in mem:
            ep[a] = g + 1
    epl = [a for mem in groups for a in mem]
    nE = len(epl)
    first, pos = [], 0
    for mem in groups:
        first.append(pos)
        pos += len(mem)
    atoms = []
    for i in range(n):
        par = rng.choice([0, 0, 1, 2, 3, 4])
        sb = [rng.choice([0, 0] + list(range(1, n + 1))) if par else 0 for _ in range(3)]
        sb = [x for x in sb if x] + [0] * 3
        sb = sb[:3]
        atoms.append([rng.choice([0, 0, 0, 5, 100, 4097]), par] + sb + [ep[i], rng.choice([0, 0, 1, 2]),
                      rng.choice([0, 0, 1]), rng.choice([0, 0, 1]), rng.choice([0, 0, 0, 1, -1, 2]),
                      1 if rng.random() < 0.3 else 0])
    tgs = []
    for g, mem in enumerate(groups):
        tgs.append([len(mem), first[g]] + [rng.randrange(0, 4), rng.randrange(0, 2)] + [rng.choice([0, 0, 1, 2]) for _ in range(3)]
                   + [rng.choice([0, 0, 5, 100])])
    sv = [rng.randrange(0, 12) for _ in range(7)] + [rng.choice([0, 1])]
    ranks = []
    for _ in range(3):
        anum, symm = perm_ranks(rng, n, ntg)
        ranks.append((anum, symm))
    # keep the 2nd/3rd ranking identical to the first sometimes (equal CT)
    if rng.random() < 0.4:
        ranks[1] = ranks[0]
    ctlen = n + m + (ntg + nE if ntg else 0)
    nmax = rng.choice([ctlen + 1, ctlen + 3, ctlen + 3, max(ctlen - 1, 0), ctlen + 8])
    plen = rng.choice([0, 0, 0, ctlen, ctlen + 1])
    tmax = 3 * ntg + nE + 1
    return {
        "n": n, "m": m, "ntg": ntg, "nE": nE, "edges": e, "atoms": atoms, "tgs": tgs, "epl": epl,
        "nIsoEnd": rng.choice([0, 0, 1, 2]), "sv": sv, "hasTgi": 1 if (ntg or rng.random() < 0.3) and True else 0,
        "ranks": ranks,
        "lens": [nmax, plen, rng.choice([0, n, n + 2, 1]), rng.choice([0, 0, 1, 3]),
                 rng.choice([tmax, tmax, tmax + 2, tmax - 1, 0]), rng.choice([0, 0, tmax]),
                 rng.choice([0, ntg, ntg + 1, max(ntg - 1, 0)]), rng.choice([0, 0, ntg, 1]), rng.randrange(2)],
    }


def fmt(c):
    o = [f"{c['n']} {c['m']} {c['ntg']} {c['nE']}"]
    o.append(" ".join(f"{a} {b}" for a, b in c["edges"]))
    o += [" ".join(map(str, a)) for a in c["atoms"]]
    o += [" ".join(map(str, t)) for t in c["tgs"]]
    o.append(" ".join(map(str, c["epl"])))
    o.append(str(c["nIsoEnd"]))
    o.append(" ".join(map(str, c["sv"])))
    o.append(str(c["hasTgi"]))
    for anum, symm in c["ranks"]:
        o.append(" ".join(map(str, anum)))
        o.append(" ".join(map(str, symm)))
    o.append(" ".join(map(str, c["lens"])))
    return "\n".join(o) + "\n"


def main():
    rng = random.Random(2468)
    cases = [gen(rng) for _ in range(4000)]
    stdin = f"{len(cases)}\n" + "".join(fmt(c) for c in cases)
    objs = ma.build_lib()
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "h.c").write_text(HARNESS)
        exe = Path(d) / "h"
        subprocess.run(["cc", "-w", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"), *objs, "-lm"], check=True)
        r = subprocess.run([str(exe)], input=stdin.encode(), capture_output=True, timeout=900)
        assert r.returncode == 0, (r.returncode, r.stderr.decode()[:2000])
        out = r.stdout.decode().strip().split("\n")
    starts = [i for i, l in enumerate(out) if l.startswith("N ")]
    starts.append(len(out))
    assert len(starts) - 1 == len(cases)
    data = [{"in": c, "out": out[a + 1:b]} for c, a, b in zip(cases, starts, starts[1:])]
    path = ROOT / "tests/data/inchi_cano1.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases ->", path)


if __name__ == "__main__":
    main()
