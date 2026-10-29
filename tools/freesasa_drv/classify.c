/* Driver for tests/data/freesasa_classify.txt: freesasa 2.0.3 classifier lookups.
     cc -O1 -I SRC classify.c SRC/{classifier,classifier_protor,classifier_naccess,classifier_oons,util}.c -o classify
   Input lines "P|N|O residue atom"; output "radius class" (radius as %.17g). */
#include <stdio.h>
#include <string.h>
#include "freesasa.h"
#include "freesasa_internal.h"
#include "classifier.h"

freesasa_verbosity freesasa_get_verbosity(void) { return FREESASA_V_NORMAL; }

int main(void) {
    char k[8], r[64], a[64];
    while (scanf("%7s %63s %63s", k, r, a) == 3) {
        const freesasa_classifier *c = k[0] == 'P' ? &freesasa_protor_classifier : k[0] == 'N' ? &freesasa_naccess_classifier : &freesasa_oons_classifier;
        printf("%.17g %d\n", freesasa_classifier_radius(c, r, a), (int)freesasa_classifier_class(c, r, a));
    }
    return 0;
}

/* the structure functions classifier.c refers to for per-structure summaries are not used here */
void freesasa_range_nodearea(freesasa_nodearea *a, const freesasa_structure *s, const freesasa_result *r, int f, int l) {}
int freesasa_structure_n(const freesasa_structure *s) { return 0; }
