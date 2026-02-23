import csv, shutil
from pathlib import Path

map_tsv = Path("07_functional_inference/foldx/pdb_sample_map.tsv")
outdir = Path("07_functional_inference/foldx/pdb_flat")
outdir.mkdir(parents=True, exist_ok=True)

with open(map_tsv) as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        og = row["orthogroup"]
        pdb = Path(row["pdb_path"])
        token = row["token"]
        dest = outdir / f"{og}__{token}.pdb"
        shutil.copy2(pdb, dest)

print(f"[OK] wrote flat PDBs -> {outdir}")
