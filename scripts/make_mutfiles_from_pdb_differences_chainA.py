import os, glob
from collections import defaultdict

PROJECT  = os.getcwd()
REPAIRED = os.path.join(PROJECT, "07_functional_inference/foldx/repaired_renum")
OUTDIR   = os.path.join(PROJECT, "07_functional_inference/foldx/mutations_diff_chainA")
os.makedirs(OUTDIR, exist_ok=True)

AA3_to_AA1 = {
 "ALA":"A","CYS":"C","ASP":"D","GLU":"E","PHE":"F","GLY":"G","HIS":"H","ILE":"I","LYS":"K","LEU":"L",
 "MET":"M","ASN":"N","PRO":"P","GLN":"Q","ARG":"R","SER":"S","THR":"T","VAL":"V","TRP":"W","TYR":"Y"
}

def load_chain_seq_and_plddt(pdb_path, chain="A"):
    # retorna dict pos->(aa1, plddt_float)
    d = {}
    seen=set()
    with open(pdb_path, errors="ignore") as f:
        for line in f:
            if not line.startswith("ATOM"): 
                continue
            if line[21].strip() != chain:
                continue
            if line[12:16].strip() != "CA":
                continue
            resn = line[17:20].strip()
            resi = line[22:26].strip()
            if not resi.strip().isdigit():
                continue
            pos = int(resi)
            if pos in seen:
                continue
            seen.add(pos)
            aa1 = AA3_to_AA1.get(resn, "X")
            # AlphaFold usa B-factor como pLDDT frequentemente
            try:
                plddt = float(line[60:66].strip())
            except:
                plddt = None
            d[pos] = (aa1, plddt)
    return d

# ajuste aqui se quiser filtrar por confiança do AlphaFold
PLDDT_MIN = 0.0   # ex: 70.0 para filtrar

pdbs = sorted(glob.glob(os.path.join(REPAIRED, "OG*__*_Repair.pdb")))
by_og = defaultdict(list)
for p in pdbs:
    og = os.path.basename(p).split("__")[0]
    by_og[og].append(p)

wrote=0
skipped=0

for og, lst in sorted(by_og.items()):
    # acha o cereus como referência (primeiro que casar)
    refs = [p for p in lst if "__cereus_" in os.path.basename(p)]
    if not refs:
        continue
    ref = sorted(refs)[0]
    ref_tok = os.path.basename(ref).split("__")[1].replace("_Repair.pdb","")
    ref_data = load_chain_seq_and_plddt(ref, "A")

    for tgt in sorted(lst):
        if tgt == ref:
            continue
        tgt_tok = os.path.basename(tgt).split("__")[1].replace("_Repair.pdb","")
        tgt_data = load_chain_seq_and_plddt(tgt, "A")

        muts=[]
        common = sorted(set(ref_data.keys()) & set(tgt_data.keys()))
        for pos in common:
            a1, p1 = ref_data[pos]
            a2, p2 = tgt_data[pos]
            if a1 == "X" or a2 == "X":
                continue
            # filtro pLDDT (usa o menor dos dois)
            if p1 is not None and p2 is not None:
                if min(p1, p2) < PLDDT_MIN:
                    continue
            if a1 != a2:
                muts.append(f"{a1}A{pos}{a2};")

        out = os.path.join(OUTDIR, f"{og}__{ref_tok}_to_{tgt_tok}__diff.txt")
        if muts:
            with open(out,"w") as f:
                f.write("\n".join(muts)+"\n")
            wrote += 1
        else:
            # ainda escreve vazio? melhor não.
            skipped += 1

print(f"[OK] wrote_nonempty={wrote} skipped_empty={skipped} -> {OUTDIR}")
