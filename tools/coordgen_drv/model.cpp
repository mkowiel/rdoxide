// Driver for tests/data/coordgen_model.txt: builds molecules from a description on stdin, runs forceUpdateStruct and
// prints what the model classes of coordgenlibs 3.0.2 answer.
//   c++ -std=c++17 -O1 -ffp-contract=off -I SRC model.cpp OBJDIR/*.o -o model
// Input per molecule: "MOL n m", n lines "A Z charge xbits ybits", m lines "B a b order", "END".
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
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
            sscanf(line.c_str(), "A %d %d %x %x", &z, &ch, &xb, &yb);
            auto* a = mol.addNewAtom();
            a->setAtomicNumber(z);
            a->charge = ch;
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
        for (size_t i = 0; i < A.size(); i++) {
            auto* a = A[i];
            printf("a %zu hs %d nb%s bd%s rg%s", i, a->_implicitHs, idx(a->neighbors, A).c_str(), idx(a->bonds, B).c_str(),
                   idx(a->rings, R).c_str());
            printf(" chiral %d stereoActive %d", a->canBeChiral(), !a->hasNoStereoActiveBonds());
            printf(" cw%s", idx(a->clockwiseOrderedNeighbors(), A).c_str());
            if (!a->neighbors.empty()) {
                auto sub = a->getSubmolecule(a->neighbors[0]);
                printf(" sub%s", idx(sub, A).c_str());
            }
            auto v = a->getSingleAdditionVector();
            printf(" sav %08x %08x\n", bits(v.x()), bits(v.y()));
        }
        for (size_t i = 0; i < B.size(); i++) {
            auto* b = B[i];
            printf("b %zu small %d macro %d term %d inter %d stereo %d rings%s\n", i, b->isInSmallRing(), b->isInMacrocycle(),
                   b->isTerminal(), b->isInterFragment(), b->isStereo(), idx(b->rings, R).c_str());
        }
        for (size_t i = 0; i < R.size(); i++) {
            auto* r = R[i];
            auto c = r->findCenter();
            printf("r %zu size %d benz %d arom %d mac %d atoms%s bonds%s center %08x %08x\n", i, r->size(), r->isBenzene(),
                   r->isAromatic(), r->isMacrocycle(), idx(r->_atoms, A).c_str(), idx(r->_bonds, B).c_str(), bits(c.x()),
                   bits(c.y()));
        }
        auto tc = mol.totalCharge();
        sketcherMinimizerPointF mn, mx;
        mol.boundingBox(mn, mx);
        auto ce = mol.center();
        printf("m charge %d bbox %08x %08x %08x %08x center %08x %08x\n", tc, bits(mn.x()), bits(mn.y()), bits(mx.x()), bits(mx.y()),
               bits(ce.x()), bits(ce.y()));
        // flip every bond that is not in a ring, on a copy of the coordinates
        for (size_t i = 0; i < B.size(); i++) {
            if (!B[i]->rings.empty() || B[i]->isTerminal()) continue;
            std::vector<sketcherMinimizerPointF> saved;
            for (auto* a : A) saved.push_back(a->coordinates);
            B[i]->flip();
            printf("flip %zu", i);
            for (auto* a : A) printf(" %08x %08x", bits(a->coordinates.x()), bits(a->coordinates.y()));
            printf(" wedge");
            for (auto* b : B) printf(" %d", b->isWedge);
            printf("\n");
            for (size_t k = 0; k < A.size(); k++) A[k]->coordinates = saved[k];
            for (auto* b : B) b->isWedge = false;
        }
        printf("$$\n");
    }
    return 0;
}
