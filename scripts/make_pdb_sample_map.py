import os, re, glob, csv, sys
from pathlib import Path
from Bio import SeqIO

# Ajuste estes caminhos se necessário
pdb_root = Path("04_structures/alphafold_best_pdb")
fasta_root = Path("04_structures/fasta")

# Heurística: o FASTA costuma ter header com sample (S115 etc.) ou algum identificador.
# Vamos procurar padrões tipo S\d+ dentro do header. Se não achar, cai para o ID inteiro.
sample_re = re.compile(r"(S\d{2,4})")

rows = []
for og_dir in sorted(pdb_root.glob("OG*/")):
    og = og_dir.name
    # Carrega todos FASTA do OG para construir um dicionário: seq_id_base -> sample
    # Ex.: arquivo: GeneX__C_insularis_S115.fasta ou header contém S115
    fasta_files = list((fasta_root / og).glob("*.fasta")) + list((fasta_root / og).glob("*.faa")) + list((fasta_root / og).glob("*.fa"))
    seq_to_sample = {}

    for ff in fasta_files:
        for rec in SeqIO.parse(str(ff), "fasta"):
            hdr = rec.id + " " + rec.description
            m = sample_re.search(hdr)
            sample = m.group(1) if m else rec.id
            # tenta reduzir rec.id para um "token" que pode bater com o foldname (cereus_19790 etc.)
            # Ex.: se rec.id contiver "cereus_19790" guardamos isso.
            for tok in ["cereus_", "fffd_", "fffe_", "ffff_", "fffg_", "fffh_", "fffi_", "fff3_"]:
                if tok in hdr:
                    # pega trecho tipo cereus_19790
                    mm = re.search(rf"({tok}\d+)", hdr)
                    if mm:
                        seq_to_sample[mm.group(1)] = sample

    # Agora percorre PDBs e tenta extrair token do filename
    for pdb in og_dir.glob("*.pdb"):
        m = re.search(r"fold_([a-z0-9]+_\d+)\.best_model_0\.pdb$", pdb.name)
        if not m:
            continue
        token = m.group(1)  # ex: cereus_19790
        sample = seq_to_sample.get(token, None)
        rows.append([og, pdb.as_posix(), token, sample if sample else "UNKNOWN"])

out = Path("07_functional_inference/foldx/pdb_sample_map.tsv")
out.parent.mkdir(parents=True, exist_ok=True)

with open(out, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["orthogroup","pdb_path","token","sample"])
    w.writerows(rows)

unknown = sum(1 for r in rows if r[3] == "UNKNOWN")
print(f"[OK] wrote {len(rows)} rows -> {out} | UNKNOWN={unknown}")
