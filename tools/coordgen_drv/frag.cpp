// Driver for tests/data/coordgen_frag.txt: builds molecules from a description on stdin, runs forceUpdateStruct and
// CoordgenFragmenter::splitIntoFragments and prints the fragments of coordgenlibs 3.0.2.
//   c++ -std=c++17 -O1 -ffp-contract=off -I SRC model.cpp OBJDIR/*.o -o model
// Input per molecule: "MOL n m", n lines "A Z charge xbits ybits [f][c]", m lines "B a b order", "END".
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include "CoordgenFragmenter.h"
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
template <class T> static std::string idx(const std::vector<T*>& v, const std::vector<T*>& all)
{
    std::string s;
    for (auto* p : v) {
        s += " " + std::to_string(std::find(all.begin(), all.end(), p) - all.begin());
    }
    return s;
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
        auto& B = mol._bonds;
        auto& R = mol._rings;
        CoordgenFragmenter::splitIntoFragments(&mol);
        auto& F = mol._fragments;
        auto fi = [&](sketcherMinimizerFragment* f) -> long { return f ? std::find(F.begin(), F.end(), f) - F.begin() : -1; };
        printf("main %ld hasFixed %d hasConstrained %d nfrag %zu\n", fi(mol.getMainFragment()), mol.hasFixedFragments,
               mol.hasConstrainedFragments, F.size());
        for (size_t i = 0; i < F.size(); i++) {
            auto* f = F[i];
            printf("f %zu atoms%s bonds%s rings%s inter%s children", i, idx(f->getAtoms(), A).c_str(), idx(f->getBonds(), B).c_str(),
                   idx(f->getRings(), R).c_str(), idx(f->_interFragmentBonds, B).c_str());
            for (auto* c : f->_children) printf(" %ld", fi(c));
            printf(" parent %ld btp %ld chain %d fixed %d constrained %d cflip %d weight %u dbl %u heavy %u nconstr %u nfixed %u\n",
                   fi(f->getParent()), f->_bondToParent ? std::find(B.begin(), B.end(), f->_bondToParent) - B.begin() : -1L, f->isChain,
                   f->fixed, f->constrained, f->constrainedFlip, f->totalWeight(), f->countDoubleBonds(), f->countHeavyAtoms(),
                   f->countConstrainedAtoms(), f->countFixedAtoms());
        }
        for (size_t i = 0; i < B.size(); i++) {
            printf("b %zu %ld %ld rev %d\n", i, std::find(A.begin(), A.end(), B[i]->startAtom) - A.begin(),
                   std::find(A.begin(), A.end(), B[i]->endAtom) - A.begin(), B[i]->isReversed);
        }
        printf("$$\n");
    }
    return 0;
}
