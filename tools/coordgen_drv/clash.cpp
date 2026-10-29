// Driver for tests/data/coordgen_clash.txt: fragment builder, buildFromFragments and CoordgenMinimizer::avoidClashes.
//   c++ -std=c++17 -O1 -ffp-contract=off -I SRC clash.cpp OBJDIR/*.o -o clash
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
#include "CoordgenMinimizer.h"
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
        for (auto* f : F) fb.initializeCoordinates(f);
        CoordgenMinimizer min;
        min.m_evenAngles = false;
        min.m_atoms = mol._atoms;
        min.m_bonds = mol._bonds;
        min.m_molecules.push_back(&mol);
        min.buildFromFragments(true);
        bool clean = min.avoidClashes();
        float e = min.scoreClashes(&mol);
        printf("clean %d energy %08x\n", clean, bits(e));
        for (size_t i = 0; i < A.size(); i++)
            printf("a %zu %08x %08x\n", i, bits(A[i]->coordinates.x()), bits(A[i]->coordinates.y()));
        printf("$$\n");
    }
    return 0;
}
