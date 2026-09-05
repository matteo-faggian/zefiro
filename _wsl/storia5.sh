#!/bin/bash
# Storia completa di pressioni e portate ai bordi del run 4.
# Unisce tutti i riavvii di surfaceFieldValue in una sola serie temporale e la
# stampa in forma leggibile, con il residuo di massa (in - out)/in.
source /root/of/etc/profile 2>/dev/null || true
/root/zef/bin/python - <<'PY'
import numpy as np, pathlib

P = pathlib.Path("/root/zefiro/runs/cfd/run4/caso/postProcessing")

def serie(nome, colonna):
    """Serie (t, valore) unita fra tutti i riavvii, senza doppioni."""
    righe = {}
    for f in P.joinpath(nome).rglob("surfaceFieldValue.dat"):
        for l in f.read_text().splitlines():
            if l.startswith("#"):
                continue
            c = l.split()
            if len(c) > colonna:
                righe[float(c[0])] = float(c[colonna])
    t = np.array(sorted(righe))
    return t, np.array([righe[x] for x in t])

t_pa, pa = serie("p_ingresso_aria", 1)
t_pg, pg = serie("p_ingresso_gpl", 1)
t_pu, pu = serie("p_uscita", 1)
t_ma, ma = serie("portata_ingresso_aria", 1)
t_mg, mg = serie("portata_ingresso_gpl", 1)
t_mu, mu = serie("portata_uscita", 1)

print(f"campioni: p_aria {len(t_pa)}, p_gpl {len(t_pg)}, p_usc {len(t_pu)}, "
      f"m_aria {len(t_ma)}, m_gpl {len(t_mg)}, m_usc {len(t_mu)}")
print(f"finestra temporale: {t_pu.min()*1e6:.1f} - {t_pu.max()*1e6:.1f} us")
print()
print(f"{'t [us]':>9} {'p_aria':>9} {'p_gpl':>9} {'p_usc':>9} "
      f"{'m_in':>11} {'m_out':>11} {'sbil %':>8}")
comune = sorted(set(t_pu) & set(t_mu) & set(t_pa) & set(t_ma))
ia = {x: i for i, x in enumerate(t_pa)}
ig = {x: i for i, x in enumerate(t_pg)}
iu = {x: i for i, x in enumerate(t_pu)}
ja = {x: i for i, x in enumerate(t_ma)}
jg = {x: i for i, x in enumerate(t_mg)}
ju = {x: i for i, x in enumerate(t_mu)}
passo = max(1, len(comune)//40)
for k, x in enumerate(comune):
    if k % passo and k != len(comune)-1:
        continue
    m_in = ja and (ma[ja[x]] + (mg[jg[x]] if x in jg else 0.0))
    m_out = mu[ju[x]]
    sb = 100.0*(abs(m_in) - abs(m_out))/abs(m_in) if m_in else float('nan')
    print(f"{x*1e6:9.1f} {pa[ia[x]]/1e5:9.4f} "
          f"{(pg[ig[x]]/1e5 if x in ig else float('nan')):9.4f} "
          f"{pu[iu[x]]/1e5:9.4f} {m_in:11.6f} {m_out:11.6f} {sb:8.2f}")
PY
