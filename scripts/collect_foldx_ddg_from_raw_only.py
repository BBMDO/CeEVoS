import glob, os, csv

ROOT = "07_functional_inference/foldx/buildmodel_nonempty"
OUT  = "07_functional_inference/foldx/results/foldx_ddg_from_RAW.csv"

def parse_raw(path):
    in_table = False
    mut = None
    wt = None
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
            if pdb.startswith("WT_"):
                wt = (pdb, te)
            else:
                mut = (pdb, te)
    return mut, wt

rows = []
for rundir in glob.glob(os.path.join(ROOT, "*", "*"), recursive=False):
    if not os.path.isdir(rundir):
        continue
    mutation_set = os.path.basename(os.path.dirname(rundir))
    pdb_run = os.path.basename(rundir)
    orthogroup = mutation_set.split("__")[0] if mutation_set.startswith("OG") else "NA"

    raws = glob.glob(os.path.join(rundir, "Raw_*.fxout"))
    if not raws:
        continue
    raw_path = raws[0]

    mutfile = os.path.join(rundir, "individual_list.txt")
    muts = ""
    if os.path.isfile(mutfile):
        muts = open(mutfile).read().strip().replace("\n","|")

    mut, wt = parse_raw(raw_path)
    if mut is None or wt is None:
        continue

    mut_pdb, mut_e = mut
    wt_pdb, wt_e = wt
    ddg = mut_e - wt_e

    rows.append([orthogroup, mutation_set, pdb_run,
                 os.path.basename(raw_path),
                 wt_pdb, wt_e,
                 mut_pdb, mut_e,
                 ddg, muts])

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["orthogroup","mutation_set","pdb_run",
                "raw_file",
                "wt_pdb","wt_total_energy",
                "mut_pdb","mut_total_energy",
                "ddG_mut_minus_WT",
                "mutations_in_run"])
    w.writerows(rows)

print(f"[OK] {len(rows)} rows -> {OUT}")
