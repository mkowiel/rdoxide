#!/usr/bin/env python3
"""Avalon Toolkit reachability: which C functions of the toolkit are reachable from the entry
points RDKit's AvalonTools wrapper uses. Run from `src/main/C` of the unpacked
`ava-formake-AvalonToolkit_2.0.5-pre.3` (clang needed). Three steps in one script; see
docs/AVALON_PORT.md. Writes defs.json and reachable.json to the current directory.
"""

# ---- step 1: function definitions and their callees (clang AST)
import json,subprocess,sys,glob,collections
files=sorted(glob.glob('common/*.c'))+['programs/struchk.c']
files=[f for f in files if f.split('/')[-1] in {l.split('/')[-1].strip().rstrip(')') for l in """layout.c symboltable.c patclean.c utilities.c symbol_lists.c stereo.c set.c perceive.c local.c graph.c geometry.c forio.c depictutil.c denormal.c casutils.c ssmatch.c rtutils.c smi2mol.c didepict.c pattern.c canonizer.c aacheck.c fixcharges.c struchk.c reaccsio.c hashcode.c""".split()}]
defs={}  # name -> (file, lines, callees)
def walk(n, cur, file):
    if isinstance(n,dict):
        if n.get('kind')=='DeclRefExpr':
            r=n.get('referencedDecl',{})
            if r.get('kind')=='FunctionDecl' and cur is not None: cur.add(r['name'])
        for v in n.values(): walk(v,cur,file)
    elif isinstance(n,list):
        for v in n: walk(v,cur,file)
for f in files:
    out=subprocess.run(['clang','-fsyntax-only','-w','-Xclang','-ast-dump=json','-Iinclude','-Wno-everything',f],capture_output=True,text=True)
    if not out.stdout: print('FAIL',f,out.stderr[:200],file=sys.stderr); continue
    ast=json.loads(out.stdout)
    for n in ast['inner']:
        if n.get('kind')=='FunctionDecl' and any(i.get('kind')=='CompoundStmt' for i in n.get('inner',[])):
            loc=n.get('loc',{}); rng=n.get('range',{})
            cs=set(); walk([i for i in n['inner'] if i.get('kind')=='CompoundStmt'],cs,f)
            b=rng.get('begin',{}).get('line'); e=rng.get('end',{}).get('line')
            # lines may be omitted when equal to previous; recompute below
            defs.setdefault(n['name'],[]).append((f,b,e,cs))
json.dump({k:[(a,b,c,sorted(d)) for a,b,c,d in v] for k,v in defs.items()},open('defs.json','w'))
print(len(defs),'function definitions')

# ---- step 2: reachability
import json,collections
defs=json.load(open('defs.json'))
entry="""SMIToMOL MolStr2Mol MolToMolStr FreeMolecule MyFree SetFingerprintBits SetFingerprintCountsWithFocus MOLToSMI CanSmiles RecolorMolecule LayoutMolecule RunStruchk ClearParameters InitCheckMol GetMsgList FreeMsgList CloseOpenFiles""".split()
missing=[e for e in entry if e not in defs]; print('entry not defined:',missing)
seen=set(); st=[e for e in entry if e in defs]
while st:
    f=st.pop()
    if f in seen: continue
    seen.add(f)
    for (file,b,e,cs) in defs[f]:
        for c in cs:
            if c in defs and c not in seen: st.append(c)
print(len(seen),'reachable of',len(defs))
byfile=collections.Counter()
for f in seen:
    for (file,*_) in defs[f]: byfile[file]+=1
tot=collections.Counter()
for f in defs:
    for (file,*_) in defs[f]: tot[file]+=1
for file in sorted(tot): print(f"{file:28s} {byfile[file]:4d}/{tot[file]:4d}")
json.dump(sorted(seen),open('reachable.json','w'))

# ---- step 3: line counts
import json,re,collections
root='./'
defs=json.load(open('defs.json')); seen=set(json.load(open('reachable.json')))
tot=collections.Counter(); alll=collections.Counter(); per={}
cache={}
def lines(f):
    if f not in cache: cache[f]=open(root+f,errors='replace').read().split('\n')
    return cache[f]
def size(name,f):
    L=lines(f); pat=re.compile(r'^[A-Za-z_].*\b'+re.escape(name)+r'\s*\(')
    for i,l in enumerate(L):
        if pat.match(l) and not l.rstrip().endswith(';'):
            # find opening brace then closing at column 0
            j=i
            while j<len(L) and not L[j].startswith('{') and '{' not in L[j]: j+=1
            k=j
            while k<len(L) and not L[k].startswith('}'): k+=1
            return k-i+1
    return 0
for n,v in defs.items():
    for (f,*_) in v:
        s=size(n,f); alll[f]+=s
        if n in seen: tot[f]+=s
print(sum(tot.values()),'lines in reachable functions of',sum(alll.values()))
for f in sorted(alll): print(f"{f:28s} {tot[f]:5d}/{alll[f]:5d}")
