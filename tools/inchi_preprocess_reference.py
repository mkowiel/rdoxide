"""Generates tests/data/inchi_preprocess.json from the InChI C code:
`PreprocessOneStructure` of runichi3.c with every function it calls
(strutil.c, ichister.c; `OrigAtData_Duplicate`,
`OrigAtData_bCheckUnusualValences`) on random salts, metal complexes,
ion pairs and plain molecules. The function texts are cut out of the C
sources. Needs a C compiler and the fetched InChI sources."""

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
#include <math.h>
#include <stdint.h>
#include "mode.h"
#include "incomdef.h"
#include "ichidrp.h"
#include "inpdef.h"
#include "ichi.h"
#include "strutil.h"
#include "util.h"
#include "ichierr.h"
#include "extr_ct.h"
#include "ichimain.h"
#include "ichicomn.h"
OAD_PolymerUnit * OAD_PolymerUnit_CreateCopy( OAD_PolymerUnit *u ) { return NULL; }
void inchi_swap( char *a, char *b, size_t width ) { char tmp; if (a != b) { while (width--) { tmp = *a; *a++ = *b; *b++ = tmp; } } }
@FUNCS@
static uint64_t bits(double d) { uint64_t u; memcpy(&u, &d, 8); return u; }
static void dump(const char *tag, ORIG_ATOM_DATA *od) {
    printf("%s %d %d %d %d\n", tag, od->num_inp_atoms, od->bDisconnectSalts, od->bDisconnectCoord, od->num_components);
    if (!od->at) { printf("-\n"); return; }
    printf("L");
    if (od->nCurAtLen) for (int i = 0; i < od->num_components; i++) printf(" %d", od->nCurAtLen[i]);
    printf("\nO");
    if (od->nOldCompNumber) for (int i = 0; i < od->num_components; i++) printf(" %d", od->nOldCompNumber[i]);
    printf("\n");
    inp_ATOM *at = od->at;
    for (int i = 0; i < od->num_inp_atoms; i++) {
        printf("A %d %d %d %d %d %d %d %d %d %d %d %d %llu %llu %llu", at[i].el_number, at[i].valence, at[i].chem_bonds_valence, at[i].charge, at[i].num_H,
               at[i].num_iso_H[0], at[i].num_iso_H[1], at[i].num_iso_H[2], at[i].iso_atw_diff, at[i].radical, at[i].component, at[i].orig_at_number,
               (unsigned long long)bits(at[i].x), (unsigned long long)bits(at[i].y), (unsigned long long)bits(at[i].z));
        for (int k = 0; k < at[i].valence; k++) printf(" %d:%d:%d", at[i].neighbor[k], at[i].bond_type[k], at[i].bond_stereo[k]);
        for (int m = 0; m < MAX_NUM_STEREO_BONDS; m++) printf(" S%d:%d:%d:%d", at[i].sb_ord[m], at[i].sn_ord[m], at[i].sb_parity[m], at[i].sn_orig_at_num[m]);
        printf("\n");
    }
}
int main(void) {
    int ncases;
    scanf("%d", &ncases);
    for (int c = 0; c < ncases; c++) {
        int n, nbonds;
        unsigned long long tf; int fnu, nowarn;
        scanf("%d %d %llu %d %d", &n, &nbonds, &tf, &fnu, &nowarn);
        inp_ATOM *at = calloc(n + 1, sizeof(inp_ATOM));
        for (int i = 0; i < n; i++) {
            int v, cbv, el, chg, rad, iso, nh, ih0, ih1, ih2;
            double x, y, z; char en[8];
            scanf("%d %d %d %d %d %d %d %d %d %d %lf %lf %lf %s", &v, &cbv, &el, &chg, &rad, &iso, &nh, &ih0, &ih1, &ih2, &x, &y, &z, en);
            at[i].valence = v; at[i].chem_bonds_valence = cbv; at[i].el_number = el; at[i].charge = chg;
            at[i].radical = rad; at[i].iso_atw_diff = iso; at[i].num_H = nh;
            at[i].num_iso_H[0] = ih0; at[i].num_iso_H[1] = ih1; at[i].num_iso_H[2] = ih2;
            at[i].x = x; at[i].y = y; at[i].z = z; at[i].orig_at_number = i + 1;
            strcpy(at[i].elname, en);
            for (int j = 0; j < v; j++) { int q, b, s; scanf("%d %d %d", &q, &b, &s); at[i].neighbor[j] = q; at[i].bond_type[j] = b; at[i].bond_stereo[j] = s; }
            for (int m = 0; m < MAX_NUM_STEREO_BONDS; m++) {
                int a, b2, p, o; scanf("%d %d %d %d", &a, &b2, &p, &o);
                at[i].sb_ord[m] = a; at[i].sn_ord[m] = b2; at[i].sb_parity[m] = p; at[i].sn_orig_at_num[m] = o;
            }
        }
        ORIG_ATOM_DATA orig; memset(&orig, 0, sizeof(orig));
        ORIG_ATOM_DATA prep[2]; memset(prep, 0, sizeof(prep));
        orig.at = at; orig.num_inp_atoms = n; orig.num_inp_bonds = nbonds;
        STRUCT_DATA sd; memset(&sd, 0, sizeof(sd));
        INPUT_PARMS ip; memset(&ip, 0, sizeof(ip));
        ip.bTautFlags = tf; ip.bFixNonUniformDraw = fnu; ip.bNoWarnings = nowarn;
        int r = PreprocessOneStructure(NULL, &sd, &ip, &orig, prep);
        printf("R %d %d %d [%s] %llu %llu %llu %llu\n", r, sd.nErrorType, sd.nStructReadError, sd.pStrErrStruct,
               (unsigned long long)sd.bTautFlags[0], (unsigned long long)sd.bTautFlags[1],
               (unsigned long long)sd.bTautFlagsDone[0], (unsigned long long)sd.bTautFlagsDone[1]);
        dump("O", &orig);
        dump("P", &prep[0]);
        dump("Q", &prep[1]);
    }
    return 0;
}
"""

METALS = [(3, "Li"), (11, "Na"), (19, "K"), (12, "Mg"), (20, "Ca"), (30, "Zn"), (26, "Fe"), (29, "Cu"), (13, "Al")]
SYM = {1: "H", 6: "C", 7: "N", 8: "O", 9: "F", 15: "P", 16: "S", 17: "Cl", 35: "Br", 53: "I"}
RNG = random.Random(5)


def new(atoms, el, chg=0, nh=0, rad=0, sym=None):
    atoms.append({"el": el, "sym": sym or SYM.get(el, "X"), "chg": chg, "rad": rad, "iso": 0, "nh": nh, "ih": [0, 0, 0],
                  "xyz": [round(RNG.uniform(-5, 5), 3), round(RNG.uniform(-5, 5), 3), RNG.choice([0.0, 0.0, 0.5])],
                  "nb": [], "sb": []})
    return len(atoms) - 1


def link(atoms, a, b, t=1):
    atoms[a]["nb"].append((b, t, 0))
    atoms[b]["nb"].append((a, t, 0))


def chain(atoms, rng):
    n = rng.randrange(1, 5)
    prev = None
    for _ in range(n):
        a = new(atoms, rng.choice([6, 6, 6, 7, 8]))
        if prev is not None:
            link(atoms, prev, a, rng.choice([1, 1, 2]))
        prev = a
    return prev


def fragment(atoms, rng):
    k = rng.choice(["chain", "nitro", "amm", "metalsalt", "complex", "pair", "iso", "ion3", "hcat"])
    if k == "chain":
        chain(atoms, rng)
    elif k == "nitro":
        c = new(atoms, 6)
        nn = new(atoms, 7, 1)
        link(atoms, c, nn)
        link(atoms, nn, new(atoms, 8), 2)
        link(atoms, nn, new(atoms, 8, -1))
    elif k == "amm":
        nitr = new(atoms, 7, rng.choice([0, 1]), rng.choice([2, 3, 4]))
        o = new(atoms, 8, rng.choice([0, -1]))
        link(atoms, nitr, o)
        if rng.random() < 0.8:
            c = new(atoms, 6, 0, 3)
            link(atoms, o, c)
    elif k == "metalsalt":
        el, sy = rng.choice(METALS[:5])
        m = new(atoms, el, rng.choice([0, 1]), 0, sym=sy)
        for _ in range(rng.randrange(1, 3)):
            o = new(atoms, 8)
            link(atoms, m, o)
            c = new(atoms, 6)
            link(atoms, o, c)
            link(atoms, c, new(atoms, 8), 2)
            link(atoms, c, new(atoms, 6, 0, 3))
    elif k == "complex":
        el, sy = rng.choice(METALS)
        m = new(atoms, el, rng.choice([0, 0, 1, 2]), rng.choice([0, 0, 1]), sym=sy)
        for _ in range(rng.randrange(1, 4)):
            lig = new(atoms, rng.choice([7, 8, 8, 17, 16]), rng.choice([0, 0, -1]), rng.choice([0, 1]))
            link(atoms, m, lig)
            if rng.random() < 0.6:
                link(atoms, lig, new(atoms, 6, 0, 3))
    elif k == "pair":
        a = new(atoms, rng.choice([7, 8, 6]), 1)
        b = new(atoms, rng.choice([8, 6, 7]), -1)
        link(atoms, a, b, rng.choice([1, 1, 2]))
        for x in (a, b):
            for _ in range(rng.randrange(0, 3)):
                link(atoms, x, new(atoms, rng.choice([6, 8, 7])), rng.choice([1, 1, 2]))
    elif k == "iso":
        o = new(atoms, 8, 0, 1)
        atoms[o]["ih"][rng.randrange(3)] = 1
        atoms[o]["nh"] = 0
        link(atoms, o, new(atoms, 6, 0, 3))
        if rng.random() < 0.4:
            atoms[o]["iso"] = 1
    elif k == "ion3":
        c1 = new(atoms, 6, 1)
        o = new(atoms, 8)
        n1 = new(atoms, 7, -1)
        link(atoms, c1, o)
        link(atoms, o, n1)
        link(atoms, c1, new(atoms, 6, 0, 3))
    else:
        h = new(atoms, 1, 1)


def gen(rng):
    atoms = []
    for _ in range(rng.randrange(1, 4)):
        fragment(atoms, rng)
    # a few cross bonds
    for _ in range(rng.randrange(0, 2)):
        a, b = rng.randrange(len(atoms)), rng.randrange(len(atoms))
        if a != b and all(x != b for x, _, _ in atoms[a]["nb"]) and len(atoms[a]["nb"]) < 5 and len(atoms[b]["nb"]) < 5 \
                and atoms[a]["el"] not in (1,) and atoms[b]["el"] not in (1,):
            link(atoms, a, b)
    # explicit hydrogens
    for i in range(len(atoms)):
        if rng.random() < 0.1 and atoms[i]["el"] in (6, 7, 8):
            link(atoms, i, new(atoms, 1))
    for a in atoms:
        a["cbv"] = sum(t for _, t, _ in a["nb"])
    return atoms


def main():
    rng = random.Random(99)
    cases = []
    for _ in range(4000):
        flags = 0
        for f in (0x2, 0x80, 0x100, 0x200, 0x40000):
            if rng.random() < 0.6:
                flags |= f
        cases.append({"atoms": gen(rng), "flags": flags, "fnu": rng.choice([0, 1]), "nowarn": rng.choice([0, 0, 1])})
    stdin = f"{len(cases)}\n"
    for c in cases:
        atoms = c["atoms"]
        nbonds = sum(len(a["nb"]) for a in atoms) // 2
        stdin += f"{len(atoms)} {nbonds} {c['flags']} {c['fnu']} {c['nowarn']}\n"
        for a in atoms:
            stdin += (f"{len(a['nb'])} {a['cbv']} {a['el']} {a['chg']} {a['rad']} {a['iso']} {a['nh']} "
                      f"{a['ih'][0]} {a['ih'][1]} {a['ih'][2]} {a['xyz'][0]} {a['xyz'][1]} {a['xyz'][2]} {a['sym']} "
                      + " ".join(f"{x} {t} {s}" for x, t, s in a["nb"]) + " 0 0 0 0 0 0 0 0 0 0 0 0\n")
        c["nbonds"] = nbonds
    text = (SRC / "strutil.c").read_text(errors="replace")
    ister = (SRC / "ichister.c").read_text(errors="replace")
    util = (SRC / "util.c").read_text(errors="replace")
    run3 = (SRC / "runichi3.c").read_text(errors="replace")

    def cut(t, start):
        a = t.rindex(start)
        return t[a:t.index("\n}\n", a) + 3]

    st = lambda s: cut(text, s)
    funcs = [
        cut(ister, "int get_opposite_sb_atom( inp_ATOM *at,\n"),
        cut(ister, "int ReconcileCmlIncidentBondParities( inp_ATOM *at,\n"),
        cut(ister, "int ReconcileAllCmlBondParities( inp_ATOM *at,\n"),
        st("int RemoveInpAtBond( inp_ATOM *atom, int iat, int k )\n{"),
        st("int DisconnectInpAtBond( inp_ATOM *at,\n"),
        st("int the_only_doublet_neigh( inp_ATOM *at,\n"),
        st("int fix_non_uniform_drawn_oxoanions( int num_atoms,\n"),
        st("int fix_non_uniform_drawn_amidiniums( int num_atoms,\n"),
        st("int nFindOneOM( inp_ATOM *at, int at_no, int ord_OM[], int num_OM )\n{"),
        st("int remove_ion_pairs( int num_atoms, inp_ATOM *at )\n{"),
        st("int fix_odd_things( int num_atoms,\n"),
        st("int post_fix_odd_things( int num_atoms, inp_ATOM *at )\n{"),
        st("int bIsAmmoniumSalt( inp_ATOM *at,\n"),
        st("int DisconnectAmmoniumSalt( inp_ATOM *at,\n"),
        st("int bIsMetalSalt( inp_ATOM *at, int i )\n{"),
        st("int DisconnectMetalSalt( inp_ATOM *at, int i )\n{"),
        st("int DisconnectSalts( ORIG_ATOM_DATA *orig_inp_data, int bDisconnect )\n{"),
        st("int bIsMetalToDisconnect( inp_ATOM *at, int i, int bCheckMetalValence )\n{"),
        st("int bMayDisconnectMetals( ORIG_ATOM_DATA *orig_inp_data,\n"),
        st("int DisconnectMetals( ORIG_ATOM_DATA *orig_inp_data,\n"),
        st("int DisconnectOneLigand( inp_ATOM *at,\n"),
        st("double dist3D( inp_ATOM *at1, inp_ATOM *at2 )\n{"),
        st("double GetMinDistDistribution( inp_ATOM *at,\n"),
        st("int move_explicit_Hcation( inp_ATOM *at,\n"),
        st("int get_iat_number( int el_number )\n{"),
        st("int bNumHeterAtomHasIsotopicH( inp_ATOM *atom, int num_atoms )\n{"),
        st("int cmp_components( const void *a1, const void *a2 )\n{"),
        st("int MarkDisconnectedComponents( ORIG_ATOM_DATA *orig_at_data,\n"),
        cut(run3, "int OrigAtData_bCheckUnusualValences( ORIG_ATOM_DATA *orig_at_data,\n"),
        cut(run3, "int OrigAtData_Duplicate( ORIG_ATOM_DATA *new_orig_atom,\n"),
        cut(run3, "int PreprocessOneStructure( struct tagINCHI_CLOCK *ic,\n"),
    ]
    protos = """int get_opposite_sb_atom( inp_ATOM *at, int cur_atom, int icur2nxt, int *pnxt_atom, int *pinxt2cur, int *pinxt_sb_parity_ord );
int ReconcileCmlIncidentBondParities( inp_ATOM *at, int cur_atom, int prev_atom, S_CHAR *visited, int bDisconnected );
int RemoveInpAtBond( inp_ATOM *atom, int iat, int k );
int DisconnectInpAtBond( inp_ATOM *at, AT_NUMB *nOldCompNumber, int iat, int neigh_ord );
int remove_ion_pairs( int num_atoms, inp_ATOM *at );
int nFindOneOM( inp_ATOM *at, int at_no, int ord_OM[], int num_OM );
int post_fix_odd_things( int num_atoms, inp_ATOM *at );
int bIsAmmoniumSalt( inp_ATOM *at, int i, int *piO, int *pk, S_CHAR *num_explicit_H );
int bIsMetalSalt( inp_ATOM *at, int i );
int bIsMetalToDisconnect( inp_ATOM *at, int i, int bCheckMetalValence );
int DisconnectOneLigand( inp_ATOM *at, AT_NUMB *nOldCompNumber, S_CHAR *bMetal, char *elnumber_Heteroat, int num_halogens, int num_atoms, int iMetal, int jLigand, INCHI_MODE *bTautFlagsDone );
double dist3D( inp_ATOM *at1, inp_ATOM *at2 );
double GetMinDistDistribution( inp_ATOM *at, int num_at, int iat, int iat_H, int bInAllComponents, double min_dist[], int num_segm );
int move_explicit_Hcation( inp_ATOM *at, int num_at, int iat, int iat_H, int bInAllComponents );
int get_iat_number( int el_number );
int cmp_components( const void *a1, const void *a2 );
int OrigAtData_bCheckUnusualValences( ORIG_ATOM_DATA *orig_at_data, int bAddIsoH, char *pStrErrStruct, int bNoWarnings);
enum { IAT_H = 0, IAT_C, IAT_N, IAT_P, IAT_O, IAT_S, IAT_Se, IAT_Te, IAT_F, IAT_Cl, IAT_Br, IAT_I, IAT_MAX = 12 };
#define MIN_BOND_LENGTH   (1.0e-6)
#define MIN_COS           (1.0e-6)
#define MIN_BOND_LENGTH2  (MIN_BOND_LENGTH*MIN_BOND_LENGTH)
#define MAX_BOND_LENGTH   (1.0e30)
"""
    with tempfile.TemporaryDirectory() as d:
        exe = Path(d) / "h"
        (Path(d) / "h.c").write_text(HARNESS.replace("@FUNCS@", protos + "\n".join(funcs)))
        subprocess.run(["cc", "-w", "-ffp-contract=off", "-DTARGET_API_LIB", f"-I{SRC}", "-o", str(exe), str(Path(d) / "h.c"),
                        str(SRC / "util.c"), str(SRC / "ichierr.c"), "-lm"], check=True)
        out = subprocess.run([str(exe)], input=stdin.encode(), capture_output=True, check=True).stdout.decode().strip().split("\n")
    # record boundaries: lines starting with "R " begin a case
    starts = [i for i, l in enumerate(out) if l.startswith("R ")]
    starts.append(len(out))
    assert len(starts) - 1 == len(cases), (len(starts), len(cases))
    data = []
    changed = 0
    for c, a, b in zip(cases, starts, starts[1:]):
        lines = out[a:b]
        changed += lines[0].split()[3] != "0" or "disconnected" in lines[0]
        data.append({"in": c, "out": lines})
    path = ROOT / "tests/data/inchi_preprocess.json"
    path.write_text(json.dumps(data, separators=(",", ":")))
    print(len(data), "cases,", changed, "with messages ->", path)


main()
