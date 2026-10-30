#!/usr/bin/env python3
"""Times similarity and bit operations on large collections of ExplicitBitVects.

    python tools/bench_bitvect.py rdoxide
    python tools/bench_bitvect.py rdkit      # original RDKit, other venv
Prints one JSON object {operation: milliseconds}.
"""
import json
import random
import sys
import time

lib = sys.argv[1] if len(sys.argv) > 1 else "rdoxide"
if lib=='rdkit':
    from rdkit import DataStructs
else:
    from rdoxide import DataStructs
random.seed(1)
def mk(n, dens, N):
    out=[]
    for _ in range(N):
        bv=DataStructs.ExplicitBitVect(n)
        bv.SetBitsFromList(random.sample(range(n), int(n*dens)))
        out.append(bv)
    return out
res={}
for nbits in (2048, 16384):
    N=100000 if nbits==2048 else 20000
    fps=mk(nbits,0.1,N); q=fps[0]
    def t(f,rep=3):
        b=1e9
        for _ in range(rep):
            s=time.perf_counter(); f(); b=min(b,time.perf_counter()-s)
        return b*1e3
    res[f'BulkTanimoto {nbits} bits, 1 x {N}']=t(lambda: DataStructs.BulkTanimotoSimilarity(q,fps))
    res[f'BulkDice {nbits} bits, 1 x {N}']=t(lambda: DataStructs.BulkDiceSimilarity(q,fps))
    sub=fps[:5000]
    res[f'bitwise AND + OnBitCount {nbits} bits, 5000 pairs']=t(lambda: [(q&f).GetNumOnBits() for f in sub])
    res[f'NumBitsInCommon {nbits} bits, {N}']=t(lambda: [DataStructs.AllProbeBitsMatch(q,f) for f in fps])
print(json.dumps(res))
