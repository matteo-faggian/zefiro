#!/bin/bash
set +e
source /mnt/i/AA_ENGINE/_wsl/ambiente.sh
cp /mnt/i/AA_ENGINE/scripts/cfd/*.py /root/zefiro/scripts/cfd/
sed -i 's/\r$//' /root/zefiro/scripts/cfd/*.py
export PYTHONPATH=/root/zefiro/src
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
echo "--- dov'e' la mesh gia' fatta ---"
find /root/zefiro/runs/cfd -maxdepth 2 -name '*.msh' -exec ls -la {} \;
cd /root/zefiro
rm -rf runs/cfd/run4; mkdir -p runs/cfd/run4
#: LA MESH E' QUELLA DI RUN3, per nome esplicito. `find | head -1` ha pescato
#: runs/cfd/iniettore.msh, che e' la mesh grossolana della prova di fumo (11 MB
#: contro 57): la corsa sarebbe partita su un'altra mesh e il confronto con la
#: 3 non avrebbe voluto dire niente. Un file scelto per ordine alfabetico non
#: e' un file scelto.
M=/root/zefiro/runs/cfd/run3/iniettore.msh
[ -s "$M" ] || { echo "manca la mesh di run3: $M"; exit 2; }
echo "uso $M  ($(stat -c%s "$M") byte)"
cp "$M" runs/cfd/run4/iniettore.msh || exit 2
/root/zef/bin/python scripts/cfd/caso_mescolamento.py --caso runs/cfd/run4/caso \
    --solo-caso --fine 7.0e-5 --grossa 3.5e-4 --tempo 2.5e-4 \
    --scrittura 1.0e-5 --campionamento 5.0e-7 --proc 24 --binario 2>&1 | tail -20
echo "--- 0/U ingressi ---"
sed -n '/boundaryField/,/uscita/p' runs/cfd/run4/caso/0/U
echo "--- fvSolution PIMPLE ---"
sed -n '/^PIMPLE/,/^}/p' runs/cfd/run4/caso/system/fvSolution | grep -v '^ *//'
echo "--- fvOptions attivi ---"
grep -v '^ *//' runs/cfd/run4/caso/constant/fvOptions | grep -vE '^\s*$'
