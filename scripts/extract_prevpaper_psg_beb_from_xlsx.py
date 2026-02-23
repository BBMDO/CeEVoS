import re
import pandas as pd
from pathlib import Path

xlsx = Path("00_inputs/selection_prevpaper/425_2024_4442_MOESM1_ESM.xlsx")
out_psg = Path("01_psg_sites/PSG_master_prevpaper.csv")
out_beb = Path("01_psg_sites/BEB_sites_prevpaper.csv")

tables = ["TableS5","TableS6","TableS7","TableS8","TableS9","TableS10","TableS11"]

def parse_branch_and_title(sheet_df):
    # linha 0 col 0 contém o título tipo: "Table S8. ... (Branch 4)."
    title = str(sheet_df.iloc[0,0])
    m = re.search(r"\(Branch\s+(\d+)\)", title)
    branch = int(m.group(1)) if m else None
    return branch, title

def read_sheet_with_header(xlsx, sheet):
    # No seu arquivo, a linha 1 (index=1) é o cabeçalho real.
    df = pd.read_excel(xlsx, sheet_name=sheet, header=None)
    branch, title = parse_branch_and_title(df)
    header = df.iloc[1].tolist()
    body = df.iloc[2:].copy()
    body.columns = header
    body["branch"] = branch
    body["table"] = sheet
    body["table_title"] = title
    return body

def find_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    raise KeyError(f"Não achei coluna entre: {candidates}. Colunas: {list(df.columns)}")

# regex para capturar entradas tipo: "117 R 0.980*" (podem haver várias na mesma célula)
beb_re = re.compile(r"(\d+)\s+([A-Z])\s+([01]\.\d+)\*?")

psg_rows = []
beb_rows = []

for t in tables:
    df = read_sheet_with_header(xlsx, t)

    col_gene = find_col(df, ["Gene","gene","OG","Orthogroup"])
    col_lrt  = find_col(df, ["ΔLRT","ΔLRT ", "DeltaLRT", "Δ LRT", "ΔLRT\t"])
    col_p    = find_col(df, ["p-value","p value","pvalue","P-value"])
    # coluna com BEB sites tem um nome longo no XLSX
    col_beb  = None
    for c in df.columns:
        if isinstance(c, str) and "BEB" in c:
            col_beb = c
            break

    col_ann  = None
    for cand in ["Annotation","annotation","Function","functional annotation"]:
        if cand in df.columns:
            col_ann = cand
            break

    col_ipr  = "InterPro" if "InterPro" in df.columns else None
    col_go   = "GO" if "GO" in df.columns else None

    # PSG master (um por gene)
    for _, r in df.iterrows():
        og = r.get(col_gene, None)
        if pd.isna(og):
            continue
        psg_rows.append({
            "branch": int(r["branch"]) if pd.notna(r["branch"]) else None,
            "orthogroup": str(og).strip(),
            "LRT": float(r[col_lrt]) if pd.notna(r[col_lrt]) else None,
            "p_value": float(r[col_p]) if pd.notna(r[col_p]) else None,
            "annotation": (str(r[col_ann]).strip() if col_ann and pd.notna(r[col_ann]) else ""),
            "interpro": (str(r[col_ipr]).strip() if col_ipr and pd.notna(r[col_ipr]) else ""),
            "go": (str(r[col_go]).strip() if col_go and pd.notna(r[col_go]) else ""),
            "source_table": t
        })

        # BEB sites (0..N por gene)
        if col_beb and pd.notna(r.get(col_beb, None)):
            s = str(r[col_beb])
            for m in beb_re.finditer(s):
                site = int(m.group(1))
                aa   = m.group(2)
                pp   = float(m.group(3))
                beb_rows.append({
                    "branch": int(r["branch"]) if pd.notna(r["branch"]) else None,
                    "orthogroup": str(og).strip(),
                    "site": site,
                    "aa_reported": aa,
                    "pp": pp,
                    "source_table": t
                })

psg = pd.DataFrame(psg_rows).drop_duplicates(subset=["branch","orthogroup","source_table"])
# FDR por branch (mais correto do que global misturado)
psg["FDR"] = psg.groupby("branch")["p_value"].transform(lambda x: pd.Series(pd.Series(x).rank(method="min")).values)  # placeholder
# vamos fazer BH de verdade:
def bh(p):
    p = pd.Series(p, dtype=float)
    n = p.notna().sum()
    if n == 0:
        return p
    order = p.sort_values().index
    ranks = pd.Series(range(1, n+1), index=order)
    q = p.loc[order] * n / ranks
    q = q[::-1].cummin()[::-1]
    out = pd.Series(index=p.index, dtype=float)
    out.loc[order] = q.values
    return out
psg["FDR"] = psg.groupby("branch")["p_value"].transform(bh)

psg.to_csv(out_psg, index=False)
beb = pd.DataFrame(beb_rows).sort_values(["branch","orthogroup","site"])
beb.to_csv(out_beb, index=False)

print(f"[OK] PSG master: {out_psg}  (rows={len(psg)})")
print(f"[OK] BEB sites:  {out_beb}  (rows={len(beb)})")
