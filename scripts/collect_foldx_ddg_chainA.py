import glob, os, re, csv

ROOT = "07_functional_inference/foldx/buildmodel_nonempty"
OUT  = "07_functional_inference/foldx/results/foldx_ddg_summary.csv"

mut_re = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]A[0-9]+[ACDEFGHIKLMNPQRSTVWY]$")

def last_float(tokens):
    for t in reversed(tokens):
        try:
            return float(t)
        except:
            pass
    return None

rows = []
for dif in glob.glob(os.path.join(ROOT, "**", "Dif_*.fxout"), recursive=True):
    parts = dif.split(os.sep)
    mutation_set = parts[4] if len(parts) > 4 else "NA"
    pdb_run = parts[5] if len(parts) > 5 else "NA"
    orthogroup = mutation_set.split("__")[0] if mutation_set.startswith("OG") else "NA"

    with open(dif, errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if "FoldX" in s or s.startswith(("*", "#")):
                continue
            toks = re.split(r"[,\t; ]+", s)
            toks = [t for t in toks if t]
            if not toks:
                continue
            mut = toks[0]
            if not mut_re.match(mut):
                continue
            ddg = last_float(toks)
            if ddg is None:
                continue
            rows.append([orthogroup, mutation_set, pdb_run, os.path.basename(dif), mut, ddg])

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["orthogroup","mutation_set","pdb_run","dif_file","mutation","ddG_kcal_mol"])
    w.writerows(rows)

print(f"[OK] {len(rows)} rows -> {OUT}")
