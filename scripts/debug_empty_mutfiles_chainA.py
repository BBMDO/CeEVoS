import os, re, glob

PROJECT  = os.getcwd()
REPAIRED = os.path.join(PROJECT, "07_functional_inference/foldx/repaired_renum")
IN_DIR   = os.path.join(PROJECT, "07_functional_inference/foldx/mutations")          # originais
OUT_DIR  = os.path.join(PROJECT, "07_functional_inference/foldx/mutations_chainA")   # gerados

AA3_to_AA1 = {
 "ALA":"A","CYS":"C","ASP":"D","GLU":"E","PHE":"F","GLY":"G","HIS":"H","ILE":"I","LYS":"K","LEU":"L",
 "MET":"M","ASN":"N","PRO":"P","GLN":"Q","ARG":"R","SER":"S","THR":"T","VAL":"V","TRP":"W","TYR":"Y"
}

def residue_at(pdb_path, chain, pos):
    with open(pdb_path, errors="ignore") as f:
        for line in f:
            if not line.startswith("ATOM"): 
                continue
            atom = line[12:16].strip()
            ch   = line[21].strip()
            resn = line[17:20].strip()
            resi = line[22:26].strip()
            if atom == "CA" and ch == chain and resi.strip().isdigit() and int(resi) == pos:
                return AA3_to_AA1.get(resn)
    return None

pos_re = re.compile(r"[A-Z](\d+)[A-Z]")

def parse_positions(path):
    positions=[]
    with open(path) as f:
        for line in f:
            s=line.strip()
            if not s: 
                continue
            m = pos_re.search(s.replace(";",""))
            if m:
                positions.append(int(m.group(1)))
    return sorted(set(positions))

def parse_from_to(stem):
    # stem: OGxxxx__from_to_to__pp0.95
    parts = stem.split("__")
    og = parts[0]
    rest = "__".join(parts[1:])
    fromtok = rest.split("_to_")[0]
    totok   = rest.split("_to_")[1].split("__")[0]
    return og, fromtok, totok

def is_nonempty(path):
    return os.path.isfile(path) and os.path.getsize(path) > 0

print("file\tpositions\tfound_in_from\tfound_in_to\tequal_sites\tdiff_sites\tstatus")

for src in sorted(glob.glob(os.path.join(IN_DIR, "*.txt"))):
    base = os.path.basename(src)
    stem = base[:-4]
    og, fromtok, totok = parse_from_to(stem)

    from_pdb = os.path.join(REPAIRED, f"{og}__{fromtok}_Repair.pdb")
    to_pdb   = os.path.join(REPAIRED, f"{og}__{totok}_Repair.pdb")

    positions = parse_positions(src)
    if not positions:
        print(f"{base}\t0\t0\t0\t0\t0\tEMPTY_INPUT")
        continue

    found_from=0
    found_to=0
    equal_sites=0
    diff_sites=0

    for pos in positions:
        a1 = residue_at(from_pdb, "A", pos)
        a2 = residue_at(to_pdb,   "A", pos)
        if a1: found_from += 1
        if a2: found_to += 1
        if a1 and a2:
            if a1 == a2:
                equal_sites += 1
            else:
                diff_sites += 1

    out = os.path.join(OUT_DIR, base)
    if is_nonempty(out):
        status="NONEMPTY_OK"
    else:
        if found_from==0 or found_to==0:
            status="POS_NOT_FOUND_IN_PDB"
        elif diff_sites==0:
            status="ALL_EQUAL_NO_MUT"
        else:
            status="UNKNOWN_CHECK"

    print(f"{base}\t{len(positions)}\t{found_from}\t{found_to}\t{equal_sites}\t{diff_sites}\t{status}")
