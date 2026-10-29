/* Driver for tests/data/freesasa_sasa.txt: freesasa 2.0.3 Lee & Richards / Shrake & Rupley on coordinates and radii.
     cc -O1 -ffp-contract=off -I SRC -DNDEBUG sasa.c SRC/{coord,nb,sasa_sr,sasa_lr,util}.c -lm -o sasa
   Input per case: "CASE alg n probe" then n lines "x y z r" as hex doubles; output the total and the areas as hex doubles. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "freesasa.h"
#include "freesasa_internal.h"

freesasa_verbosity freesasa_get_verbosity(void) { return FREESASA_V_NORMAL; }
const freesasa_parameters freesasa_default_parameters = {FREESASA_LEE_RICHARDS, 1.4, 100, 20, 1};

static double fromhex(const char *s) { uint64_t u = strtoull(s, NULL, 16); double d; memcpy(&d, &u, 8); return d; }
static void tohex(double d, char *out) { uint64_t u; memcpy(&u, &d, 8); sprintf(out, "%016llx", (unsigned long long)u); }

int main(void) {
    char line[512];
    while (fgets(line, sizeof line, stdin)) {
        if (strncmp(line, "CASE", 4)) continue;
        int alg, n; char pb[64];
        sscanf(line, "CASE %d %d %63s", &alg, &n, pb);
        double *xyz = malloc(sizeof(double) * 3 * n), *r = malloc(sizeof(double) * n), *sasa = malloc(sizeof(double) * n);
        for (int i = 0; i < n; i++) {
            char a[40], b[40], c[40], d[40];
            fgets(line, sizeof line, stdin);
            sscanf(line, "%39s %39s %39s %39s", a, b, c, d);
            xyz[3 * i] = fromhex(a); xyz[3 * i + 1] = fromhex(b); xyz[3 * i + 2] = fromhex(c); r[i] = fromhex(d);
        }
        freesasa_parameters p = freesasa_default_parameters;
        p.alg = alg ? FREESASA_SHRAKE_RUPLEY : FREESASA_LEE_RICHARDS;
        p.probe_radius = fromhex(pb);
        p.n_threads = 1;
        coord_t *c = freesasa_coord_new_linked(xyz, n);
        int rc = alg ? freesasa_shrake_rupley(sasa, c, r, &p) : freesasa_lee_richards(sasa, c, r, &p);
        double total = 0;
        for (int i = 0; i < n; i++) total += sasa[i];
        char h[32];
        tohex(total, h);
        printf("rc %d total %s\n", rc, h);
        for (int i = 0; i < n; i++) { tohex(sasa[i], h); printf("s %d %s\n", i, h); }
        printf("$$\n");
        freesasa_coord_free(c);
        free(xyz); free(r); free(sasa);
    }
    return 0;
}
