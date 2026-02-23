import sys
from pathlib import Path

inp = Path(sys.argv[1])
out = Path(sys.argv[2])

new_resseq = 0
prev_key = None

with inp.open() as f, out.open("w") as w:
    for line in f:
        if line.startswith(("ATOM", "HETATM")):
            # força cadeia A (coluna 22, index 21)
            line = line[:21] + "A" + line[22:]

            resname = line[17:20]
            chain = line[21]
            resseq = line[22:26]
            icode  = line[26]

            key = (resname, chain, resseq, icode)
            if key != prev_key:
                new_resseq += 1
                prev_key = key

            # escreve novo resseq (colunas 23-26)
            line = line[:22] + f"{new_resseq:4d}" + line[26:]
            w.write(line)
        elif line.startswith("TER"):
            continue
        elif line.startswith("END"):
            continue
    w.write("END\n")
