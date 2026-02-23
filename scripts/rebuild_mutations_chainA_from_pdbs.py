import os, re, glob

PROJECT  = os.getcwd()
REPAIRED = os.path.join(PROJECT, "07_functional_inference/foldx/repaired_renum")
IN_DIR   = os.path.join(PROJECT, "07_functional_inference/foldx/mutations")          # mutfiles originais
OUT_DIR  = os.path.join(PROJECT, "07_functional_inference/foldx/mutations_chainA")   # mutfiles finais FoldX

os.makedirs(OUT_DIR, exist_ok=True)

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

# extrai posições de linhas tipo "G14K;" ou "A105V;" (só importa o número)
pos_re = re.compile(r"[A-Z](\d+)[A-Z]")

count_in = 0
count_out = 0
empty = 0
missing_pdb = 0

for src in sorted(glob.glob(os.path.join(IN_DIR, "*.txt"))):
    base = os.path.basename(src)          # OGxxxx__from_to_to__pp0.95.txt
    stem = base[:-4]
    parts = stem.split("__")
    if len(parts) < 2:
        continue

    og = parts[0]
    rest = "__".join(parts[1:])           # cereus_19790_to_fffd_23291__pp0.95
    if "_to_" not in rest:
        continue

    fromtok = rest.split("_to_")[0]
    totok   = rest.split("_to_")[1].split("__")[0]

    from_pdb = os.path.join(REPAIRED, f"{og}__{fromtok}_Repair.pdb")
    to_pdb   = os.path.join(REPAIRED, f"{og}__{totok}_Repair.pdb")

    # coletar posições do arquivo src
    positions = []
    with open(src) as f:
        for line in f:
            line=line.strip()
            if not line:
                continue
            m = pos_re.search(line.replace(";",""))
            if m:
                positions.append(int(m.group(1)))
    positions = sorted(set(positions))

    count_in += 1
    muts = []

    if not os.path.isfile(from_pdb) or not os.path.isfile(to_pdb):
        missing_pdb += 1
    else:
        for pos in positions:
            aa_from = residue_at(from_pdb, "A", pos)
            aa_to   = residue_at(to_pdb,   "A", pos)
            if not aa_from or not aa_to:
                continue
            if aa_from == aa_to:
                continue
            muts.append(f"{aa_from}A{pos}{aa_to};")

    out = os.path.join(OUT_DIR, base)
    with open(out, "w") as w:
        if muts:
            w.write("\n".join(muts) + "\n")
            count_out += 1
        else:
            empty += 1

print(f"[OK] inputs={count_in} outputs_nonempty={count_out} empty={empty} missing_pdb_pairs={missing_pdb}")
print(f"[OUTDIR] {OUT_DIR}")
