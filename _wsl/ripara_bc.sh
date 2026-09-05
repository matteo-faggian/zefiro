#!/bin/bash
# I campi sono scritti in BINARIO: l'intestazione e il blocco boundaryField sono
# ASCII ma in mezzo c'e' un blocco di byte grezzi. Si lavora quindi sui BYTE,
# non sul testo: aprire un file binario come utf-8 fallisce alla prima
# occorrenza di 0xb2, che e' esattamente quello che e' successo.
pkill -f reactingFoam; sleep 3
C=/root/zefiro/runs/cfd/run3/caso
/root/zef/bin/python - "$C" <<'PY'
import re, sys, pathlib
C = pathlib.Path(sys.argv[1])
pc = 636100.0
nuovo = ("""    uscita
    {
        type            waveTransmissive;
        field           p;
        psi             thermo:psi;
        gamma           1.35;
        fieldInf        %.1f;
        lInf            0.05;
        value           uniform %.1f;
    }""" % (pc, pc)).encode()
schema = re.compile(rb"    uscita\s*\n    \{[^}]*\}")
tocchi, istante = 0, None
for proc in sorted(C.glob("processor*")):
    tempi = sorted((d for d in proc.iterdir() if d.is_dir()
                    and re.fullmatch(r"[0-9.eE+-]+", d.name)),
                   key=lambda d: float(d.name))
    if not tempi:
        continue
    istante = tempi[-1].name
    f = tempi[-1] / "p"
    if not f.exists():
        continue
    b = f.read_bytes()
    b2, n = schema.subn(nuovo, b, count=1)
    if n:
        f.write_bytes(b2)
        tocchi += 1
print("istante:", istante, " file di p corretti:", tocchi)
PY
echo "--- verifica ---"
strings $(ls -d $C/processor0/[0-9]* | sort -g | tail -1)/p | grep -A7 'uscita' | head -9
