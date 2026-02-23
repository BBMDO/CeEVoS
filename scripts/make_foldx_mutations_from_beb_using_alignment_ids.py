import csv, os
from Bio import AlignIO

beb_csv  = "01_psg_sites/BEB_sites_prevpaper.csv"
comp_tsv = "07_functional_inference/foldx/comparisons.tsv"
outdir   = "07_functional_inference/foldx/mutations"
os.makedirs(outdir, exist_ok=True)

# carrega BEB por OG
beb_by_og = {}
with open(beb_csv) as f:
    r = csv.DictReader(f)
    for row in r:
        og = row["orthogroup"]
        beb_by_og.setdefault(og, []).append((int(row["site"]), float(row["pp"])))

def aln_col_from_ungapped(seq, pos1):
    c = 0
    for i, aa in enumerate(str(seq)):
        if aa != "-":
            c += 1
        if c == pos1:
            return i
    return None

with open(comp_tsv) as f:
    for line in f:
        line=line.strip()
        if not line or line.startswith("#"):
            continue
        og, aln_path, from_id, to_id, ppmin, from_tok, to_tok = line.split("\t")
        ppmin = float(ppmin)

        if og not in beb_by_og:
            continue

        aln = AlignIO.read(aln_path, "fasta")

        A = None
        B = None
        for rec in aln:
            if rec.id == from_id:
                A = rec
            if rec.id == to_id:
                B = rec
        if A is None or B is None:
            print(f"[WARN] {og}: nao achei from/to no alinhamento ({from_id}, {to_id})")
            continue

        muts=[]
        for site, pp in beb_by_og[og]:
            if pp < ppmin:
                continue
            col = aln_col_from_ungapped(A.seq, site)
            if col is None:
                continue
            aaA = str(A.seq[col])
            aaB = str(B.seq[col])
            if aaA in "-X" or aaB in "-X":
                continue
            if aaA != aaB:
                muts.append(f"{aaA}{site}{aaB};")

        outfile = os.path.join(outdir, f"{og}__{from_tok}_to_{to_tok}__pp{ppmin}.txt")
        with open(outfile, "w") as out:
            out.write("\n".join(sorted(set(muts))) + ("\n" if muts else ""))

        print(f"[OK] {og} {from_tok}->{to_tok}: {len(set(muts))} mutacoes -> {outfile}")
