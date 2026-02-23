import glob, os, csv, re

ROOT = "07_functional_inference/foldx/buildmodel_nonempty"
OUT  = "07_functional_inference/foldx/results/foldx_ddg_from_Dif_fxout.csv"

def parse_dif_fxout(path):
    """
    Parser robusto:
    - ignora linhas vazias/comentários
    - tenta achar um header com 'mutant'/'mutation' e coluna numérica final
    - se não achar header, usa heurística: primeira coluna = mut, último campo numérico = ddG
    """
    rows = []
    lines = open(path, errors="ignore").read().splitlines()
    for ln in lines:
        s = ln.strip()
        if not s:
            continue
        # pula linhas de banner/asteriscos
        if s.startswith(("*", "#")):
            continue
        # normaliza separador
        cols = re.split(r"[,\t; ]+", s)
        cols = [c for c in cols if c != ""]
        if len(cols) < 2:
            continue

        # tenta pegar último float como ddG
        ddg = None
        for i in range(len(cols)-1, -1, -1):
            try:
                ddg = float(cols[i])
                ddg_idx = i
                break
            except:
                pass
        if ddg is None:
            continue

        # mutante: geralmente 1ª coluna (ex GA14K) ou algo no começo
        mutant = cols[0]
        # descarta headers se mutant parece texto e não mutação
        if mutant.lower() in ("mutant", "mutation", "pdb", "wt", "wildtype"):
            continue

        rows.append((mutant, ddg))
    return rows

os.makedirs(os.path.dirname(OUT), exist_ok=True)

all_rows = []
for dif in glob.glob(os.path.join(ROOT, "**", "Dif_*.fxout"), recursive=True):
    parts = dif.split(os.sep)
    # .../buildmodel_nonempty/<mutation_set>/<pdbbase>/Dif_<pdbname>.fxout
    mutation_set = parts[4] if len(parts) > 4 else "NA"
    pdbbase = parts[5] if len(parts) > 5 else "NA"
    orthogroup = mutation_set.split("__")[0] if mutation_set.startswith("OG") else "NA"

    for mutant, ddg in parse_dif_fxout(dif):
        all_rows.append((orthogroup, mutation_set, pdbbase, os.path.basename(dif), mutant, ddg))

with open(OUT, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["orthogroup","mutation_set","pdb_run","dif_file","mutant","ddG_kcal_mol"])
    w.writerows(all_rows)

print(f"[OK] {len(all_rows)} rows -> {OUT}")
