from pathlib import Path
from Bio import AlignIO
from Bio.PDB import PDBParser, PPBuilder
import csv

pdb_dir = Path("07_functional_inference/foldx/pdb_flat")
aln_dir = Path("03_alignments/raw_nr")
out = Path("07_functional_inference/foldx/token_to_alignment_id.tsv")
out.parent.mkdir(parents=True, exist_ok=True)

parser = PDBParser(QUIET=True)
ppb = PPBuilder()

def seq_from_pdb(pdb_path: Path) -> str:
    structure = parser.get_structure("X", str(pdb_path))
    # tenta pegar o maior polipeptídeo (mais comum)
    peptides = []
    for model in structure:
        for chain in model:
            for pp in ppb.build_peptides(chain):
                peptides.append(str(pp.get_sequence()))
    if not peptides:
        return ""
    return max(peptides, key=len)

def identity(a: str, b: str) -> float:
    # identidade simples no prefixo do menor comprimento
    n = min(len(a), len(b))
    if n == 0:
        return 0.0
    matches = sum(1 for i in range(n) if a[i] == b[i])
    return matches / n

rows = []
for pdb in sorted(pdb_dir.glob("OG*__*.pdb")):
    name = pdb.stem
    og, token = name.split("__", 1)
    aln_path = aln_dir / f"{og}.fa.nr.fa.aln.fa"
    if not aln_path.exists():
        continue

    pdb_seq = seq_from_pdb(pdb)
    if not pdb_seq:
        rows.append([og, token, "", 0.0, "NO_PDB_SEQ"])
        continue

    aln = AlignIO.read(str(aln_path), "fasta")
    best = ("", 0.0)
    for rec in aln:
        aln_seq = str(rec.seq).replace("-", "")
        score = identity(pdb_seq, aln_seq)
        if score > best[1]:
            best = (rec.id, score)

    status = "OK" if best[1] >= 0.90 else "LOW_MATCH"
    rows.append([og, token, best[0], best[1], status])

with open(out, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["orthogroup","token","alignment_id","identity","status"])
    w.writerows(rows)

print(f"[OK] wrote {len(rows)} rows -> {out}")
