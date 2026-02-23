import os, glob, csv, re, math
from collections import defaultdict

ROOT = "07_functional_inference/foldx/buildmodel_diff_chainA"
OUT_REP = "07_functional_inference/foldx/results/foldx_ddg_replicates.csv"
OUT_SUM = "07_functional_inference/foldx/results/foldx_ddg_summary_mean_sd.csv"

# exemplo de mutação no seu arquivo: "GA14L;"
def read_mutations(rundir):
    p = os.path.join(rundir, "individual_list.txt")
    if not os.path.isfile(p):
        return ""
    txt = open(p, errors="ignore").read().strip()
    return txt.replace("\n", " ").strip()

def parse_raw_fxout(path):
    """
    Retorna dict pdb_name -> total_energy
    (lê tabela a partir da linha que começa com 'Pdb')
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

def mean_sd(vals):
    n = len(vals)
    if n == 0:
        return (None, None, 0)
    m = sum(vals) / n
    if n == 1:
        return (m, 0.0, 1)
    var = sum((x-m)**2 for x in vals) / (n-1)
    return (m, math.sqrt(var), n)

rows_rep = []
group_ddg = defaultdict(list)

raw_files = glob.glob(os.path.join(ROOT, "**", "Raw_*.fxout"), recursive=True)

for raw in sorted(raw_files):
    rundir = os.path.dirname(raw)

    # estrutura esperada: buildmodel_diff_chainA/<mutation_set>/<pdb_run>/
    parts = rundir.split(os.sep)
    # pega os 2 últimos componentes relevantes
    pdb_run = os.path.basename(rundir)
    mutation_set = os.path.basename(os.path.dirname(rundir))
    orthogroup = mutation_set.split("__")[0] if mutation_set.startswith("OG") else "NA"

    muts = read_mutations(rundir)

    d = parse_raw_fxout(raw)

    # para cada mutante, buscar WT correspondente
    # mutante: OGxxxx..._Repair_1.pdb
    # wt: WT_OGxxxx..._Repair_1.pdb
    for pdb_mut, e_mut in d.items():
        if pdb_mut.startswith("WT_"):
            continue
        wt_name = "WT_" + pdb_mut
        if wt_name not in d:
            # se faltar, ignora
            continue
        e_wt = d[wt_name]
        ddg = e_mut - e_wt

        # id da replica: sufixo _<n>.pdb
        m = re.search(r"_(\d+)\.pdb$", pdb_mut)
        rep = m.group(1) if m else "NA"

        rows_rep.append([
            orthogroup, mutation_set, pdb_run,
            os.path.basename(raw),
            pdb_mut, wt_name, rep,
            e_mut, e_wt, ddg,
            muts
        ])
        group_ddg[(orthogroup, mutation_set, pdb_run, muts)].append(ddg)

os.makedirs(os.path.dirname(OUT_REP), exist_ok=True)

with open(OUT_REP, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow([
        "orthogroup","mutation_set","pdb_run","raw_file",
        "mut_pdb","wt_pdb","replicate",
        "mut_total_energy","wt_total_energy","ddG_mut_minus_WT",
        "mutations_in_run"
    ])
    w.writerows(rows_rep)

# summary
rows_sum = []
for key, vals in group_ddg.items():
    orthogroup, mutation_set, pdb_run, muts = key
    m, sd, n = mean_sd(vals)
    rows_sum.append([orthogroup, mutation_set, pdb_run, n, m, sd, muts])

rows_sum.sort(key=lambda x: (x[0], x[1], x[2]))

with open(OUT_SUM, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["orthogroup","mutation_set","pdb_run","N_replicates","mean_ddG","sd_ddG","mutations_in_run"])
    w.writerows(rows_sum)

print(f"[OK] replicates: {len(rows_rep)} rows -> {OUT_REP}")
print(f"[OK] summary:    {len(rows_sum)} rows -> {OUT_SUM}")
