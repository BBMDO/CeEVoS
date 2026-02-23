import csv
from collections import defaultdict
from pathlib import Path

token_map = Path("07_functional_inference/foldx/token_to_alignment_id.tsv")
out_tsv   = Path("07_functional_inference/foldx/comparisons.tsv")

by_og = defaultdict(list)
with open(token_map) as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        if row["status"] not in ("OK","LOW_MATCH"):
            continue
        by_og[row["orthogroup"]].append(row)

lines = ["#OG\taln_fasta_path\tfrom_alignment_id\tto_alignment_id\tpp_min\tfrom_token\tto_token"]
for og, items in sorted(by_og.items()):
    # escolhe cereus_* como referência, se existir
    cereus = [x for x in items if x["token"].startswith("cereus_")]
    if not cereus:
        continue
    ref = cereus[0]
    aln = f"03_alignments/raw_nr/{og}.fa.nr.fa.aln.fa"
    for x in items:
        if x["token"] == ref["token"]:
            continue
        lines.append(f"{og}\t{aln}\t{ref['alignment_id']}\t{x['alignment_id']}\t0.95\t{ref['token']}\t{x['token']}")

out_tsv.parent.mkdir(parents=True, exist_ok=True)
out_tsv.write_text("\n".join(lines) + "\n")
print(f"[OK] wrote {out_tsv} with {len(lines)-1} comparisons")
