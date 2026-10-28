// Driver for tests/data/coordgen_templates.txt: the normalized templates ("dump" argument) and, for molecules read from
// stdin (same blocks as model.cpp), the Morgan scores, the template matching of the ring atoms (as findTemplate does it)
// and an alignment matrix of coordgenlibs 3.0.2.
//   c++ -std=c++17 -O1 -ffp-contract=off -I SRC tmpl.cpp OBJDIR/*.o -o tmpl
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <map>
#include <string>
#include "sketcherMinimizer.h"
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

int main(int argc, char** argv)
{
    sketcherMinimizer::loadTemplates();
    auto& T = sketcherMinimizer::m_templates.getTemplates();
    if (argc > 1) {
        for (size_t i = 0; i < T.size(); i++) {
            auto* t = T[i];
            std::vector<int> scores;
            int it = sketcherMinimizer::morganScores(t->_atoms, t->_bonds, scores);
            printf("T %zu %zu %zu iter %d scores", i, t->_atoms.size(), t->_bonds.size(), it);
            for (int s : scores) printf(" %d", s);
            printf(" coords");
            for (auto* a : t->_atoms) printf(" %08x %08x", bits(a->coordinates.x()), bits(a->coordinates.y()));
            printf("\n");
        }
        return 0;
    }
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
            mol.addNewBond(mol._atoms[a], mol._atoms[b])->setBondOrder(o);
        }
        std::getline(std::cin, line); // END
        sketcherMinimizerMolecule::forceUpdateStruct(mol._atoms, mol._bonds, mol._rings);
        // as findTemplate: ring atoms and ring bonds in ring order
        std::vector<sketcherMinimizerAtom*> atoms;
        std::vector<sketcherMinimizerBond*> bonds;
        std::map<sketcherMinimizerAtom*, bool> seen;
        for (auto* r : mol._rings)
            for (auto* a : r->_atoms)
                if (!seen[a]) {
                    seen[a] = true;
                    a->_generalUseN = static_cast<int>(atoms.size());
                    atoms.push_back(a);
                }
        for (auto* r : mol._rings)
            for (auto* b : r->_bonds)
                if (std::find(bonds.begin(), bonds.end(), b) == bonds.end()) bonds.push_back(b);
        std::vector<int> scores;
        int it = sketcherMinimizer::morganScores(atoms, bonds, scores);
        printf("morgan %d", it);
        for (int s : scores) printf(" %d", s);
        printf("\n");
        std::vector<unsigned> mapping;
        bool found = false;
        for (size_t t = 0; t < T.size(); t++) {
            if (sketcherMinimizer::compare(atoms, bonds, T[t], mapping)) {
                printf("match %zu", t);
                for (auto v : mapping) printf(" %u", v);
                printf("\n");
                found = true;
                break;
            }
        }
        if (!found) printf("match none\n");
        // alignment matrix between the coordinates and a rotated, shifted copy
        std::vector<sketcherMinimizerPointF> ref, pts;
        float ang = 0.7f, s = std::sin(ang), c = std::cos(ang);
        for (size_t i = 0; i < mol._atoms.size() && i < 12; i++) {
            auto p = mol._atoms[i]->coordinates;
            ref.push_back(p);
            pts.push_back(sketcherMinimizerPointF(p.x() * c + p.y() * s + (i % 3) * 0.37f, -p.x() * s + p.y() * c - (i % 4) * 0.21f));
        }
        float am[4] = {0, 0, 0, 0};
        if (!ref.empty()) sketcherMinimizer::alignmentMatrix(ref, pts, am);
        printf("align %08x %08x %08x %08x\n", bits(am[0]), bits(am[1]), bits(am[2]), bits(am[3]));
        printf("$$\n");
    }
    return 0;
}
