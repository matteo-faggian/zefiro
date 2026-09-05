#!/bin/bash
# La pressione agli ingressi SALE, si assesta, o e' gia' ferma? A meta' del
# riempimento del dominio un salto d'iniezione alto non dice niente: dice solo
# che il dominio si sta caricando. La differenza fra "transitorio" e "il modello
# 0-D e' ottimista" la fa la DERIVATA, non il valore.
set +e
source /mnt/i/AA_ENGINE/_wsl/ambiente.sh
C=/root/zefiro/runs/cfd/run4/caso
cd "$C" || exit 2
for T in $(ls -d processor0/[0-9]* | xargs -n1 basename | sort -g); do
  [ "$T" = "0" ] && continue
  mpirun -np 24 postProcess -parallel -time "$T" -fields "(p T U phi)" \
      -dict system/misura_contorni > log.c 2>&1
  pa=$(grep 'areaAverage(ingresso_aria) of p' log.c | tail -1 | awk '{print $NF}')
  pg=$(grep 'areaAverage(ingresso_gpl) of p' log.c | tail -1 | awk '{print $NF}')
  pu=$(grep 'areaAverage(uscita) of p' log.c | tail -1 | awk '{print $NF}')
  qa=$(grep 'sum(ingresso_aria) of phi' log.c | tail -1 | awk '{print $NF}')
  qu=$(grep 'sum(uscita) of phi' log.c | tail -1 | awk '{print $NF}')
  printf "t=%-12s  p_aria=%-12s p_gpl=%-12s p_usc=%-12s  mdot_in=%-14s mdot_out=%s\n" \
         "$T" "$pa" "$pg" "$pu" "$qa" "$qu"
done
