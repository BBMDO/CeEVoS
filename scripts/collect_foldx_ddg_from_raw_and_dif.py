import glob, os, csv, re

ROOT = "07_functional_inference/foldx/buildmodel_nonempty"
OUT  = "07_functional_inference/foldx/results/foldx_ddg_final.csv"

def parse_table_fxout(path):
    """
    Parse genérico de fxout com tabela iniciando em linha que começa com 'Pdb'
    Retorna dict {pdb_name: total_energy}
    """
    d = {}
    in_table = False
    with open(path, errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("Pdb"):
                in_table = True
                continue
            if not in_table:
                continue
            cols = s.split()
            if len(cols) < 2:
                continue
            pdb = cols[0]
            try:
                te = float(cols[1])
            except:
                continue
            d[pdb] = te
    return d

rows = []

for rundir in glob.glob(os.path.join(ROOT, "*", "*"), recursive=False):
    # rundir esperado: buildmodel_nonempty/<mutation_set>/<pdb_run>/
    if not os.path.isdir(rundir):
        continue

    mutation_set = os.path.basename(os.path.dirname(rundir))
    pdb_run = os.path.basename(rundir)
    orthogroup = mutation_set.split("__")[0] if mutation_set.startswith("OG") else "NA"

    difs = glob.glob(os.path.join(rundir, "Dif_*.fxout"))
    raws = glob.glob(os.path.join(rundir, "Raw_*.fxout"))
    mutfile = os.path.join(rundir, "individual_list.txt")

    if not difs or not raws:
        continue

    dif_path = difs[0]
    raw_path = raws[0]

    dif_map = parse_table_fxout(dif_path)
    raw_map = parse_table_fxout(raw_path)

    # tenta achar WT no Raw:
    # prioridade: algo começando com WT_
    wt_candidates = [k for k in raw_map.keys() if k.startswith("WT_")]
    if wt_candidates:
        wt_pdb = wt_candidates[0]
    else:
        # fallback: o próprio nome do pdb_run com .pdb (ex: OG..._Repair.pdb)
        wt_guess = pdb_run + ".pdb"
        wt_pdb = wt_guess if wt_guess in raw_map else (list(raw_map.keys())[0] if raw_map else None)

    if wt_pdb is None or wt_pdb not in raw_map:
        continue

    wt_energy = raw_map[wt_pdb]

    # pega mutações do individual_list (se existir)
    muts = ""
    if os.path.isfile(mutfile):
        muts = open(mutfile).read().strip().replace("\n","|")

    for mut_pdb, mut_energy in dif_map.items():
        ddg = mut_energy - wt_energy
        rows.append([orthogroup, mutation_set, pdb_run, os.path.basename(raw_path), wt_pdb, wt_energy,
                     os.path.basename(dif_path), mut_pdb, mut_energy, ddg, muts])

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["orthogroup","mutation_set","pdb_run",
                "raw_file","wt_pdb","wt_total_energy",
                "dif_file","mut_pdb","mut_total_energy",
                "ddG_mut_minus_WT",
                "mutations_in_run"])
    w.writerows(rows)

print(f"[OK] {len(rows)} rows -> {OUT}")
