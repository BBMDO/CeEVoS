import sys
AA3 = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L',
       'MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y'}
pdb = sys.argv[1]
chain = sys.argv[2] if len(sys.argv)>2 else "A"
seq=[]
seen=set()
with open(pdb, errors="ignore") as f:
    for line in f:
        if not line.startswith("ATOM"): 
            continue
        if line[21].strip()!=chain:
            continue
        if line[12:16].strip()!="CA":
            continue
        resn=line[17:20].strip()
        resi=line[22:26].strip()
        if not resi.strip().isdigit():
            continue
        key=(chain,int(resi))
        if key in seen:
            continue
        seen.add(key)
        seq.append(AA3.get(resn,'X'))
print(">"+pdb.split("/")[-1]+"_"+chain)
print("".join(seq))
