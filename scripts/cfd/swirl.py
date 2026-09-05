#!/usr/bin/env python3
"""Inclinare i fori del gas crea swirl utile? I conti che decidono.

LA DISTINZIONE CHE FA TUTTA LA DIFFERENZA. "Swirl" indica due cose che si
comportano in modo opposto rispetto al mescolamento:

  ROTAZIONE RIGIDA - tutto il fluido gira insieme. Il disegno delle quattro
  strisce di combustibile RUOTA, e basta: due punti che stanno sulla stessa
  circonferenza restano alla stessa distanza per sempre. Non mescola nulla.

  ROTAZIONE DIFFERENZIALE (taglio azimutale) - la velocita' tangenziale
  cambia col raggio. Due elementi a raggi diversi si separano in azimut, la
  striscia di combustibile viene STIRATA in una spirale sempre piu' sottile, e
  quando lo spessore e' abbastanza piccolo la diffusione finisce il lavoro.
  E' questo che mescola, ed e' un meccanismo di DEFORMAZIONE, non di rotazione.

Quindi la domanda giusta non e' "quanto swirl", ma "quanto TAGLIO azimutale", e
quanta deformazione totale accumula il fluido nel tempo che ha.

IL CONTO DELLA DEFORMAZIONE. Sotto taglio semplice una lamella di spessore
iniziale s0 si assottiglia come s = s0 / (1 + gamma), con gamma la deformazione
accumulata. La diffusione chiude il resto quando il tempo diffusivo residuo
sta dentro il tempo che rimane:  s^2 / nu_t < t_res. Da qui la deformazione
necessaria, che e' il numero da confrontare con quello che una geometria puo'
davvero produrre.

E CHI PUO' PRODURLA. La quantita' di moto tangenziale che si puo' mettere nel
flusso e' limitata da chi la porta. Il combustibile e' il 5 % della portata: se
anche lo si sparasse TUTTO di lato, la velocita' tangenziale della miscela
sarebbe la sua, diluita di quel 5 %. E' il conto che dice se inclinare i fori
del gas possa creare swirl, o solo spostare un po' la striscia.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))


def dati():
    from mesh_iniettore import quote_dal_progetto
    d = quote_dal_progetto()
    from zefiro.sintesi import Impianto, Requisiti, progetta
    p = progetta(Requisiti(
        spinta=50.0, durata=5.0,
        impianto=Impianto(combustibile={"C3H8": 1.0}, T_bombola_min=288.15,
                          T_ossidante_iniezione=293.0,
                          T_combustibile_iniezione=283.0,
                          cd_ossidante=0.75, cd_combustibile=0.75)),
        geometria=False)
    return d, dict(p.geometria.params.derived)


def main() -> int:
    d, g = dati()
    n = d["n_getti"]
    Ri, Re, H = d["R_getti"], d["R_anello"], d["H"]
    R_c, r_cb = d["R_c"], d["r_centerbody"]
    L_c = g["L_c"]
    m_a, m_g = d["mdot_aria"], d["mdot_gpl"]
    V_a, V_g = d["V_aria"], d["V_gpl"]
    rho_a = d["rho_aria"]

    print("=== geometria e portate ===")
    print(f"  anello  Ri {Ri*1e3:.3f}  Re {Re*1e3:.3f}  H {H*1e3:.3f} mm")
    print(f"  camera  R_c {R_c*1e3:.3f}  r_corpo {r_cb*1e3:.3f}  L_c {L_c*1e3:.2f} mm")
    print(f"  portate aria {m_a*1e3:.3f} g/s   GPL {m_g*1e3:.3f} g/s "
          f"(GPL = {m_g/(m_a+m_g)*100:.2f} % della massa)")
    print()

    #: --- 1. tempi disponibili ------------------------------------------- #
    A_cam = math.pi * (R_c ** 2 - r_cb ** 2)
    V_cam = (m_a + m_g) / (rho_a * A_cam)
    t_anello = -d["x_getti"] / V_a
    t_camera = L_c / V_cam
    print("=== 1. quanto tempo c'e' per mescolare ===")
    print(f"  nell'anello, dai getti alla fine  : {t_anello*1e6:8.1f} us "
          f"({-d['x_getti']*1e3:.2f} mm a {V_a:.0f} m/s)")
    print(f"  in CAMERA, fino a fine camera     : {t_camera*1e6:8.1f} us "
          f"({L_c*1e3:.1f} mm a {V_cam:.1f} m/s)")
    print(f"  -> la camera offre {t_camera/t_anello:.0f} volte il tempo "
          "dell'anello: e' li' che il mescolamento puo' avvenire.")
    print()

    #: --- 2. quanto taglio serve ------------------------------------------ #
    #: nu_t dalla CFD della corsa 5 (media pesata sulla densita' nel condotto)
    nu_t_anello = 6.0e-4
    nu_t_camera = 1.0e-4     # ordine di grandezza: velocita' 10 volte minori
    print("=== 2. quanta DEFORMAZIONE serve (non quanta rotazione) ===")
    for dove, s0, t_res, nu_t in (
            ("anello  (mezzo passo azimutale)", math.pi * Ri / n, t_anello, nu_t_anello),
            ("camera  (mezzo passo azimutale)", math.pi * R_c / n, t_camera, nu_t_camera)):
        s_fin = math.sqrt(nu_t * t_res)      # spessore che la diffusione chiude
        gamma = s0 / s_fin - 1.0
        print(f"  {dove}: s0 {s0*1e3:6.3f} mm -> serve s {s_fin*1e3:6.3f} mm")
        print(f"      deformazione necessaria gamma = {gamma:6.1f}")
        #: gamma = (dV_t/dr) * t ; il raggio su cui si distribuisce il taglio
        dr = H if "anello" in dove else (R_c - r_cb)
        dV = gamma * dr / t_res
        print(f"      -> differenza di velocita' tangenziale su {dr*1e3:.2f} mm: "
              f"{dV:6.1f} m/s")
    print()

    #: --- 3. chi puo' produrla -------------------------------------------- #
    print("=== 3. i fori del gas possono produrla? ===")
    print("  Se il GPL fosse sparato COMPLETAMENTE di lato, la quantita' di")
    print("  moto tangenziale sarebbe m_gpl * V_gpl, spalmata su tutta la massa:")
    Vt_mix = m_g * V_g / (m_a + m_g)
    print(f"      V_t miscela = {m_g*1e3:.3f} g/s x {V_g:.1f} m/s / "
          f"{(m_a+m_g)*1e3:.3f} g/s = {Vt_mix:.1f} m/s")
    print("  ed e' il LIMITE SUPERIORE assoluto: a 90 gradi il getto non")
    print("  penetrerebbe piu' nell'anello e non ci sarebbe piu' iniezione.")
    print()
    print("  Con l'aria invece (94.6 % della massa) la velocita' tangenziale")
    print("  la si sceglie: e' l'aria che porta la quantita' di moto.")
    print()

    #: --- 4. che cosa compra davvero inclinare i fori ---------------------- #
    print("=== 4. che cosa compra inclinare i fori del gas (angolo composto) ===")
    S = 2.0 * math.pi * Ri / n
    print(f"  mezzo passo azimutale da coprire: {S/2*1e3:.3f} mm")
    print(f"{'beta':>6} {'V_t getto':>11} {'tiro laterale':>15} {'% del mezzo passo':>19} "
          f"{'J_eff/J':>9} {'C_eff':>8}")
    for beta in (0, 10, 15, 20, 25, 30, 40):
        b = math.radians(beta)
        Vt = V_g * math.sin(b)
        #: il tiro laterale: la componente tangenziale decade mentre il getto
        #: viene inglobato. Si stima con un decadimento esponenziale sul tempo
        #: di piegatura del getto (~ d/V_aria * sqrt(J)), che integrato da'
        #: circa META' del tiro "congelato". E' una stima, e come tale va letta.
        tiro = 0.5 * Vt * t_anello
        #: la penetrazione radiale usa solo la componente NEL piano meridiano:
        #: J efficace scala come cos^2(beta), e C = (S/H) sqrt(J) come cos(beta)
        print(f"{beta:6.0f} {Vt:11.1f} {tiro*1e3:14.3f} mm {tiro/(S/2)*100:17.1f} % "
              f"{math.cos(b)**2:9.3f} {d['C']*math.cos(b):8.3f}")
    print()

    #: --- 5. il prezzo dello swirl nell'aria ------------------------------- #
    print("=== 5. il prezzo, se lo swirl lo si mette nell'aria ===")
    dp_inj = g.get("dp_iniezione", None)
    q_ax = 0.5 * rho_a * V_a ** 2
    print(f"  pressione dinamica assiale gia' spesa: {q_ax/1e5:.3f} bar "
          f"(aria a {V_a:.0f} m/s)")
    for Vt in (50, 75, 100, 150):
        q = 0.5 * rho_a * Vt ** 2
        print(f"  V_t = {Vt:3d} m/s  ->  testa dinamica aggiuntiva "
              f"{q/1e5:.3f} bar   (swirl number ~ {Vt/V_a:.2f})")
    print()
    print("  Da confrontare con il salto d'iniezione di progetto, che e' il")
    print("  budget di pressione disponibile per TUTTO l'iniettore.")
    print()

    #: --- 6. l'effetto contrario ------------------------------------------- #
    print("=== 6. attenzione: lo swirl SEPARA anche ===")
    print(f"  rho GPL {d['rho_gpl']:.2f} contro rho aria {rho_a:.2f} kg/m3: "
          f"il combustibile e' {d['rho_gpl']/rho_a:.2f} volte piu' denso.")
    for Vt in (75, 150):
        a_c = Vt ** 2 / Ri
        print(f"  a V_t = {Vt} m/s l'accelerazione centripeta all'anello vale "
              f"{a_c:.3e} m/s2 = {a_c/9.81:.0f} g")
    print("  In un campo cosi' la fase piu' densa migra verso l'esterno: lo")
    print("  swirl uniforma in azimut ma STRATIFICA in raggio. E' un effetto")
    print("  da misurare, non da assumere trascurabile.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
