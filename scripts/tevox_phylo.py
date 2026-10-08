#!/usr/bin/env python3
"""Fitch-parsimony event candidates from a TEvoX locus-state matrix."""
from __future__ import annotations

import gzip
from pathlib import Path
import argparse,csv,dataclasses,pathlib,re,sys
@dataclasses.dataclass(eq=False)
class Node:
    name:str
    children:list['Node']=dataclasses.field(default_factory=list)
    parent:'Node|None'=None
    @property
    def tip(self): return not self.children
class Parser:
    def __init__(self,text): self.t=[x for x in re.findall(r'\s*([(),:;]|[^(),:;\s]+)',text) if x];self.i=0;self.k=0
    def peek(self): return self.t[self.i] if self.i<len(self.t) else None
    def take(self,x=None):
        y=self.peek()
        if y is None or x is not None and y!=x: raise ValueError(f'expected {x!r}, found {y!r}')
        self.i+=1;return y
    def sub(self):
        if self.peek()=='(':
            self.take('(');ch=[self.sub()]
            while self.peek()==',':self.take(',');ch.append(self.sub())
            self.take(')');name=self.take() if self.peek() not in {None,',',')',':',';'} else '';n=Node(name,ch)
            for c in ch:c.parent=n
        else:
            if self.peek() in {None,',',')',':',';'}:raise ValueError('missing tip name')
            n=Node(self.take())
        if self.peek()==':':self.take(':');float(self.take())
        return n
    def parse(self):
        root=self.sub()
        if self.peek()==';':self.take(';')
        if self.i!=len(self.t):raise ValueError('unexpected Newick token')
        def names(n):
            for c in n.children:names(c)
            if not n.name:self.k+=1;n.name=f'internal_{self.k}'
        names(root);return root
def table_path(path):
    path = Path(path)
    compressed = Path(str(path) + ".gz")
    if path.suffix == ".tsv" and compressed.is_file():
        if path.is_file():
            raise ValueError("ambiguous plain/compressed table: " + str(path))
        return compressed
    return path


def open_table(path):
    path = table_path(path)
    return gzip.open(path, "rt", newline="", encoding="utf-8") if path.suffix == ".gz" else path.open(newline="", encoding="utf-8")


def post(n):
    for c in n.children:yield from post(c)
    yield n
def pre(n):
    yield n
    for c in n.children:yield from pre(c)
def load(path):
    loci=[];m={}
    with open_table(path) as h:
        r=csv.DictReader(h,delimiter='\t')
        if not r.fieldnames or not {'locus_id','genome_id','state'}<=set(r.fieldnames):raise ValueError('states TSV requires locus_id, genome_id and state')
        for row in r:
            l,g,s=row['locus_id'],row['genome_id'],row['state']
            if l not in m:m[l]={};loci.append(l)
            if g in m[l]:raise ValueError(f'duplicate {l}/{g}')
            m[l][g]={1} if s in {'PRESENT_ANNOTATED','PRESENT_UNANNOTATED'} else {0} if s=='EMPTY_SITE_CONFIRMED' else {0,1}
    if not loci:raise ValueError('states TSV has no rows')
    return loci,m
def fitch(root,tips):
    poss={};score=0
    for n in post(root):
        if n.tip:poss[n]=set(tips.get(n.name,{0,1}));continue
        cur=set(poss[n.children[0]])
        for c in n.children[1:]:
            inter=cur&poss[c]
            if inter:cur=inter
            else:cur|=poss[c];score+=1
        poss[n]=cur
    assign={root:min(poss[root])}
    for n in pre(root):
        for c in n.children:assign[c]=assign[n] if assign[n] in poss[c] else min(poss[c])
    return score,poss,assign
def main():
    ap=argparse.ArgumentParser(description='Infer candidate TE gains/losses with Fitch parsimony')
    ap.add_argument('--states',type=pathlib.Path,required=True);ap.add_argument('--tree',type=pathlib.Path,required=True);ap.add_argument('-o','--output',default='tevox_phylo');a=ap.parse_args()
    try:
        root=Parser(a.tree.read_text()).parse();loci,m=load(a.states);tips=[n.name for n in pre(root) if n.tip]
        if len(tips)!=len(set(tips)):raise ValueError('duplicate tree tips')
        extra={g for x in m.values() for g in x}-set(tips)
        if extra:raise ValueError('genomes absent from tree: '+','.join(sorted(extra)))
        ep=pathlib.Path(a.output+'.events.tsv');sp=pathlib.Path(a.output+'.ancestral_states.tsv');count=0
        with ep.open('w',newline='') as eh,sp.open('w',newline='') as sh:
            ew=csv.writer(eh,delimiter='\t',lineterminator='\n');sw=csv.writer(sh,delimiter='\t',lineterminator='\n');ew.writerow(['locus_id','event','parent','child','parsimony_score','ambiguous']);sw.writerow(['locus_id','node','assigned_state','possible_states'])
            for l in loci:
                score,p,x=fitch(root,m[l])
                for n in pre(root):
                    sw.writerow([l,n.name,'present' if x[n] else 'absent',','.join('present' if z else 'absent' for z in sorted(p[n]))])
                    if n.parent and x[n]!=x[n.parent]:count+=1;ew.writerow([l,'gain' if x[n] else 'loss',n.parent.name,n.name,score,str(len(p[n])>1 or len(p[n.parent])>1).lower()])
        print(f'TEvoX phylo inferred {count} candidate event(s) across {len(loci)} locus/loci.')
        return 0
    except (OSError,ValueError) as e:print(f'tevox_phylo: {e}',file=sys.stderr);return 1
if __name__=='__main__':raise SystemExit(main())
