import glob, os, csv

ROOT = "07_functional_inference/foldx/buildmodel_nonempty"
OUT  = "07_functional_inference/foldx/results/foldx_total_energy_from_Dif_fxout.csv"

rows = []

for dif in glob.glob(os.path.join(ROOT, "**", "Dif_*.fxout"), recursive=True):
    parts = dif.split(os.sep)
    mutation_set = parts[4] if len(parts) > 4 else "NA"
    pdb_run = parts[5] if len(parts) > 5 else "NA"
    orthogroup = mutation_set.split("__")[0] if mutation_set.startswith("OG") else "NA"

    in_table = False
    with open(dif, errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue

            # Detecta o header da tabela
            if s.startswith("Pdb"):
                in_table = True
                continue

            if not in_table:
                continue

            # Linhas da tabela: Pdb <total_energy> <...>
            cols = s.split()
            if len(cols) < 2:
                continue

            pdb_mut = cols[0]
            try:
                total_energy = float(cols[1])
            except:
                continue

            rows.append([orthogroup, mutation_set, pdb_run, os.path.basename(dif), pdb_mut, total_energy])

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["orthogroup","mutation_set","pdb_run","dif_file","pdb_mutant","total_energy_kcal_mol"])
    w.writerows(rows)

print(f"[OK] {len(rows)} rows -> {OUT}")
