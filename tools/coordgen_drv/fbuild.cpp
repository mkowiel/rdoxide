// Driver for tests/data/coordgen_fbuild.txt: CoordgenFragmentBuilder::initializeCoordinates of every fragment of a molecule
// (after forceUpdateStruct and CoordgenFragmenter::splitIntoFragments); prints coordinates, ring flags and DOFs.
//   c++ -std=c++17 -O1 -ffp-contract=off -I SRC fbuild.cpp OBJDIR/*.o -o fbuild
// Input per molecule: "MOL n m", "A Z charge xbits ybits [f][c]", "B a b order", "END". Molecules with macrocycles or
// systems the builder would hand to the macrocycle builder print "skip".
#define private public
#define protected public
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>
#include "CoordgenFragmenter.h"
#include "CoordgenFragmentBuilder.h"
#include "sketcherMinimizerFragment.h"
#include "sketcherMinimizerAtom.h"
#include "sketcherMinimizerBond.h"
#include "sketcherMinimizerMolecule.h"
#include "sketcherMinimizerRing.h"

static float fromBits(uint32_t u)
{
    float f;
    memcpy(&f, &u, 4);
    return f;
}
static uint32_t bits(float f)
{
    uint32_t u;
    memcpy(&u, &f, 4);
    return u;
}

int main()
{
    std::string line;
    while (std::getline(std::cin, line)) {
        if (line.rfind("MOL", 0) != 0) continue;
        int n, m;
        sscanf(line.c_str(), "MOL %d %d", &n, &m);
        sketcherMinimizerMolecule mol;
        for (int i = 0; i < n; i++) {
            std::getline(std::cin, line);
            int z, ch;
            unsigned xb, yb;
            char fl[8] = "";
            sscanf(line.c_str(), "A %d %d %x %x %7s", &z, &ch, &xb, &yb, fl);
            auto* a = mol.addNewAtom();
            a->setAtomicNumber(z);
            a->charge = ch;
            a->fixed = strchr(fl, 'f') != nullptr;
            a->constrained = strchr(fl, 'c') != nullptr;
            a->setCoordinates(sketcherMinimizerPointF(fromBits(xb), fromBits(yb)));
        }
        for (int i = 0; i < m; i++) {
            std::getline(std::cin, line);
            int a, b, o;
            sscanf(line.c_str(), "B %d %d %d", &a, &b, &o);
            auto* bd = mol.addNewBond(mol._atoms[a], mol._atoms[b]);
            bd->setBondOrder(o);
        }
        std::getline(std::cin, line); // END
        sketcherMinimizerMolecule::forceUpdateStruct(mol._atoms, mol._bonds, mol._rings);
        auto& A = mol._atoms;
        auto& R = mol._rings;
        CoordgenFragmentBuilder fb;
        fb.m_evenAngles = false;
        bool skip = false;
        for (auto* r : R)
            if (r->isMacrocycle()) skip = true;
        CoordgenFragmenter::splitIntoFragments(&mol);
        auto& F = mol._fragments;
        if (!skip) {
            for (auto* f : F)
                if (fb.newScorePlanarity(f->getRings()) >= 1000.f) skip = true;
        }
        if (skip) {
            printf("skip\n$$\n");
            continue;
        }
        auto ai = [&](sketcherMinimizerAtom* a) -> long { return std::find(A.begin(), A.end(), a) - A.begin(); };
        auto ri = [&](sketcherMinimizerRing* r) -> long { return std::find(R.begin(), R.end(), r) - R.begin(); };
        for (auto* f : F) fb.initializeCoordinates(f);
        for (size_t i = 0; i < A.size(); i++) {
            auto* a = A[i];
            printf("a %zu %08x %08x rigid %d inner %d set %d\n", i, bits(a->coordinates.x()), bits(a->coordinates.y()), a->rigid,
                   a->isSharedAndInner, a->coordinatesSet);
        }
        for (size_t i = 0; i < R.size(); i++) {
            auto* r = R[i];
            printf("r %zu gen %d side %d fused", i, r->coordinatesGenerated, r->side);
            for (size_t k = 0; k < r->fusedWith.size(); k++) {
                printf(" %ld:", ri(r->fusedWith[k]));
                for (auto* a : r->fusionAtoms[k]) printf("%ld,", ai(a));
            }
            printf("\n");
        }
        for (size_t i = 0; i < F.size(); i++) {
            auto* f = F[i];
            printf("f %zu templated %d dofs", i, f->isTemplated);
            for (auto* d : f->m_dofs) {
                printf(" [%d", d->tier());
                for (auto* a : d->m_atoms) printf(" %ld", ai(a));
                printf("]");
            }
            printf(" coords");
            std::vector<std::pair<long, sketcherMinimizerPointF>> cs; // the map is ordered by pointer: sort by atom index
            for (auto& kv : f->_coordinates) cs.push_back({ai(kv.first), kv.second});
            std::sort(cs.begin(), cs.end(), [](auto& a, auto& b) { return a.first < b.first; });
            for (auto& kv : cs) printf(" %ld:%08x:%08x", kv.first, bits(kv.second.x()), bits(kv.second.y()));
            printf("\n");
        }
        printf("$$\n");
    }
    return 0;
}
