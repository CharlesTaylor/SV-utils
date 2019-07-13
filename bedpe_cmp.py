#!/usr/bin/python3
import sys
import numpy as np
#Assumes same types
def reciprocal(a,b):
    q = min(a[1],b[1])-max(a[0],b[0])
    r1 = q/(a[1]-a[0])
    r2 = q/(b[1]-b[0])
    return min(r1,r2) if  min(r1,r2) > 0 else 0
def acceptable(a1, a2, b1, b2, rang):
    if a1 < b1:
        return a2 > b1 - rang
    return a1 < b2 + rang
def middle(a,b):
    return (a-b)/2+b
class bedpel:
    R=0.75
    def __init__(self,line):
        parts = line.split("\t")
        self.chr1 = parts[0]
        self.start1 = min(int(parts[1]),int(parts[2]))
        self.end1 = max(int(parts[1]),int(parts[2]))
        self.chr2 = parts[3]
        self.start2 = min(int(parts[4]),int(parts[5]))
        self.end2 = max(int(parts[4]),int(parts[5]))
        self.type = parts[6].strip()
        self.hit = False
        self.rest = "\t".join(parts[7:]).rstrip() 
    def __repr__(self):
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(self.chr1,
            self.start1,self.end1,self.chr2,self.start2,self.end2,self.type,self.rest)
    def __str__(self):
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(self.chr1,
            self.start1,self.end1,self.chr2,self.start2,self.end2,self.type,self.rest)
    def later(self,other,rang):
        return self.start1 > other.start1 + rang
    

    def overlaps(self, other):
        if self.type != other.type:
            return False

        if "tandem" not in self.type and ( "dup" in self.type or "tra" in self.type):

            return self.chr1 == other.chr1 and \
                self.chr2 == other.chr2 and \
                reciprocal((self.start1,self.end1),(other.start1,other.end1)) > bedpel.R and \
                acceptable(self.start2/2+self.end2/2,self.start2/2+self.end2/2,other.end2/2+other.start2/2,other.end2/2+other.start2/2,10000)
           #     acceptable(self.start2,self.end2,other.start2,other.end2,10000)
        if "delet" in self.type or "tandem" in self.type:
            return self.chr1 == other.chr1 and \
                self.chr2 == other.chr2 and \
                reciprocal((middle(self.start1,self.end1),middle(self.start2,self.end2)),(middle(other.start1,other.end1),middle(other.start2,other.end2))) > bedpel.R
        return self.chr1 == other.chr1 and \
                self.chr2 == other.chr2 and \
                reciprocal((self.end1,self.start2),(other.end1,other.start2))> bedpel.R


# argv[1] is ground truth, argv[2] is the predictions, argv[3] is the confidence-range
if __name__ == "__main__":

    if sys.argv[1] == "-":
        stream1 = sys.stdin
    else:
        stream1 = open(sys.argv[1],"r");    

    if sys.argv[2] == "-":
        stream2 = sys.stdin
    else:
        stream2 = open(sys.argv[2],"r");
    if len(sys.argv) > 3:
        reci = float(sys.argv[3])
        bedpel.R = reci
    if len(sys.argv) > 4:
        print_policy = sys.argv[4]
    else:
        print_policy = "NO"
    gt = []
    pr = []

    types = {}
    for line in stream1:
        bp = bedpel(line)
        gt.append(bp)
        if bp.type not in types:
            types[bp.type] = 1
    for line in stream2:
        bp = bedpel(line)
        pr.append(bp)

        if bp.type not in types:
            types[bp.type] = 1
    
    for real in gt:
        for pred in pr:
            #if( pred.later(real,MEAN)):
            #    break
            hit = real.overlaps(pred)
            real.hit = real.hit or hit
            pred.hit = pred.hit or hit


    rhitcount = 0
    phitcount = 0

    if print_policy == "ec":
        chrs = {}
        for real in gt:
            if real.chr1 not in chrs:
                chrs[real.chr1] = []
            chrs[real.chr1].append(real)
        chri = {}
        for pred in pr:
            if pred.chr1 not in chri:
                chri[pred.chr1] = []
            chri[pred.chr1].append(pred)
        for ch in sorted(chrs.keys()):
            rhitcount = 0 
            for q in chrs[ch]:
                if q.hit:
                    rhitcount+=1
            print("Truth chr: {} -> {}/{}".format(chrs[ch][0].chr1,rhitcount,len(chrs[ch])))
        for ch in sorted(chri.keys()):
            phitcount = 0 
            for q in chri[ch]:
                if q.hit:
                    phitcount+=1
            print("Prediction chr: {} -> {}/{}".format(chri[ch][0].chr1,phitcount,len(chri[ch])))
    elif print_policy == "et":
        

        for t in types.keys():
        
            rrhitcount = 0
            rtotal = 0
            for real in gt:
                if real.type != t:
                    continue
                rtotal+=1
                if real.hit:
                    rrhitcount+=1

           
            ptotal = 0
            pphitcount = 0
            for pred in pr:

                if pred.type != t:
                    continue
                ptotal+=1 
                if pred.hit:
                    pphitcount+=1

            print("{}:".format(t))
            if(rtotal>0):
                print("GT:/ {}/{}={}".format(rrhitcount,rtotal,rrhitcount/rtotal),file=sys.stderr)
            
            if(ptotal>0):
                print("PR: {}/{}={}".format(pphitcount,ptotal,pphitcount/ptotal),file=sys.stderr)
            phitcount+=pphitcount
            rhitcount+=rrhitcount

    
    else:

        rhitcount = 0
        for real in gt:
            if real.hit:
                rhitcount+=1
                if print_policy == "gt":
                    print(real)
            elif print_policy == "gf":
                print(real)
        
        phitcount = 0
        for pred in pr:
            if pred.hit:
                phitcount+=1
                if print_policy == "pt":
                    print(pred)
            elif print_policy == "pf":
                print(pred)           



    print("ALL:",file=sys.stderr)
    print("GT: {}/{}={}".format(rhitcount,len(gt),rhitcount/len(gt)),file=sys.stderr)
    print("PR: {}/{}={}".format(phitcount,len(pr),phitcount/len(pr)),file=sys.stderr)
