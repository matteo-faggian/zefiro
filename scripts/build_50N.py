#!/usr/bin/env python3
"""Costruisce il motore da 50 N in geometria implicita ed esporta lo STL.

Tutti i numeri del raffreddamento vengono da `heat_map_50N.py` e dal
dimensionamento idraulico: sono riportati qui accanto a ciascun campo perche'
si possa risalire al motivo senza cercarlo altrove.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from zefiro.config import load_operating_point                             # noqa: E402
from zefiro.geometry.parameters import default_design_vector, derive       # noqa: E402
from zefiro.l0.cycle import evaluate_l0                                    # noqa: E402
from zefiro.schemas import FuelSpec, OperatingPoint                        # noqa: E402
from zefiro.sdf.clearances import MIN_WALL_SLM                          # noqa: E402
from zefiro.injector import (                                              # noqa: E402
    BSPP, SOTTOMISURA_STAMPA, perdita_raccordo, progetta,
)
from zefiro.sdf.engine import (                                        # noqa: E402
    PORTE_PER_COLLETTORE, CircuitoRaffreddamento, TestaIniezione, costruisci,
)
from zefiro.sdf.meshing import (                                          # noqa: E402
    is_watertight, isosurface, mesh_area, mesh_volume, write_stl,
)

#: Le pressioni di alimentazione NON sono scritte qui: si derivano dalla
#: temperatura di progetto della bombola (config/operating_point.yaml,
#: fuel.bottle.T_design_K = 288.15 K) tramite feed.supply_pressures.
#: A 15 C il propano satura a 7.315 bar, quindi P_C = 6.361 bar.
#: Prima della revisione 0.4.0 qui c'erano 8.0e5 scritti a mano, che sono la
#: tensione di vapore a 18.3 C: il motore era progettato per una bombola tiepida
#: senza che nulla lo dicesse, e la raffica durava 5.26 s invece di 7.00.
_OP_CFG = load_operating_point()
P_SUP = _OP_CFG.p_air_supply
P_C = P_SUP / 1.15          # entrambi i vincoli di Dp si chiudono qui


#: Proprieta' dei due gas all'iniezione. Massa molare esatta dalla composizione;
#: gamma da CoolProp alle condizioni di monte (aria 293 K, propano 283 K).
MW_ARIA, GAMMA_ARIA = 28.9649, 1.400
MW_GPL, GAMMA_GPL = 44.0956, 1.130

#: Quote della testa che il dimensionamento NON fissa. Ognuna ha una ragione.
#:  x_getti     i getti stanno 3.5 mm a monte dello sbocco, contro i 3.1 mm di
#:              lunghezza di mescolamento richiesti: 13 % di margine.
#:  L_condotto  4.0 mm, cioe' 0.5 mm di condotto anche a monte dei getti, per
#:              non metterli sullo spigolo d'ingresso dove il profilo di
#:              velocita' non e' ancora quello del condotto.
#:  L_contr     3.0 mm per passare da (6.0, 9.5) a (3.10, 4.64) mm: la
#:              contrazione e' circa 6:1 e farla con pareti continue invece che
#:              con uno spigolo porta il Cd da 0.6 a circa 0.95.
#:  L_plenum    14.5 mm, e QUESTO NUMERO NON LO DECIDE L'ARIA: lo decide il
#:              raccordo. Il G3/8 vuole una guarnizione a legare da 21.2 mm di
#:              diametro esterno, quindi una faccia piana da 23; la faccia sta
#:              su una bozza radiale, e la bozza deve appoggiare sulla testa,
#:              che quindi non puo' essere piu' corta di 23 mm in tutto.
#:              Idraulicamente ne bastavano 6.5. E' il raccordo a dimensionare
#:              la testa, non il contrario, ed e' una cosa che si scopre solo
#:              provando a montarci un tubo vero.
#:  R_plenum    9.5 mm, testa a 11.0. Il raccordo conico verso il mantello
#:              (16.12 mm) va supportato in stampa: sono supporti ESTERNI e
#:              accessibili. Una testa dello stesso diametro del mantello non
#:              avrebbe sbalzi e peserebbe il doppio, per una tensione di
#:              cerchio di 3 MPa.
#:  R_corpo_plenum  6.0 mm: deve contenere il preforo del G1/8 (8.8 mm) con
#:              1.6 mm di parete.
#:  r_bore_gpl  2.1 mm: lascia 1.0 mm di parete verso l'aria (il doppio del
#:              minimo SLM) e un plenum 8 volte l'area dei getti.
#:  x_rampa     5.0 mm dentro la camera per riportare il corpo centrale da 3.10
#:              a 3.93 mm: semiangolo 9 gradi, nessun distacco.
X_GETTI = -3.5e-3
L_CONDOTTO = 4.0e-3
L_CONTRAZIONE = 3.0e-3
R_PLENUM = 9.5e-3
L_PLENUM = 14.5e-3
R_CORPO_PLENUM = 6.0e-3
T_MONTE = 1.5e-3
T_TESTA = 1.5e-3
R_BORE_GPL = 2.1e-3
X_RAMPA = 5.0e-3
L_CONO = 4.0e-3

#: SCELTA DELLE FILETTATURE, e non e' una scelta di catalogo: e' un bilancio
#: di pressione. Il salto d'iniezione dell'aria vale 0.954 bar ed e' TUTTO
#: quello che c'e', perche' p_c = p_sat/1.15 non lascia margini nascosti.
#: Quanto di quel salto se lo mangia il raccordo (injector.perdita_raccordo):
#:    G1/8  passaggio  5.0 mm   196 m/s   1.68 bar   impossibile
#:    G1/4  passaggio  7.5 mm    87 m/s   0.332 bar   35 % del salto
#:    G3/8  passaggio 10.0 mm    49 m/s   0.105 bar   11 %
#:    G1/2  passaggio 13.0 mm    29 m/s   0.037 bar    4 %
#: Si sceglie G3/8: il G1/4 e' fuori discussione, il G1/2 costerebbe una bozza
#: da 28 mm che allungherebbe la testa di altri 5 mm per 8 g di guadagno su un
#: bilancio in cui il tubo pesera' comunque piu' del raccordo (un tubo da 1/2"
#: lungo 2 m vale gia' 0.12 bar: e' il PRIMO pezzo da misurare, TODO A/F).
#: Sul GPL invece 1.93 g/s a 13.7 kg/m3 sono 0.14 l/s: un G1/8 costa 0.0035
#: bar, cioe' lo 0.4 %. Non serve niente di piu' grande.
FILETTO_ARIA = "G3/8"
FILETTO_GPL = "G1/8"
#: La bozza dell'aria e' centrata in modo da appoggiare tutta sulla testa.
X_PORTA_ARIA = -11.5e-3
R_BOSS_ARIA = 19.0e-3


def dimensiona_iniettore(op, params, l0):
    """L'iniettore a getti trasversali, dal dimensionamento alla geometria."""
    inj = progetta(
        mdot_aria=l0.mdot_air, mdot_gpl=l0.mdot_fuel_core,
        p_c=params.free["p_c"], p_aria_monte=op.p_air_supply,
        p_gpl_monte=op.p_fuel_supply, T_aria=float(op.T_air_in),
        T_gpl=float(op.T_fuel_in), MW_aria=MW_ARIA, MW_gpl=MW_GPL,
        gamma_aria=GAMMA_ARIA, gamma_gpl=GAMMA_GPL,
        cd=float(op.cd_injector_ox), R_iniezione=R_INIEZIONE)
    fa, fg = BSPP[FILETTO_ARIA], BSPP[FILETTO_GPL]
    testa = TestaIniezione(
        R_getti=inj.R_iniezione, R_anello=inj.R_medio, n_getti=inj.n_getti,
        d_getto=inj.d_getto, x_getti=X_GETTI, L_condotto=L_CONDOTTO,
        L_contrazione=L_CONTRAZIONE, R_plenum=R_PLENUM, L_plenum=L_PLENUM,
        R_corpo_plenum=R_CORPO_PLENUM,
        t_monte=T_MONTE, t_testa=T_TESTA, r_bore_gpl=R_BORE_GPL,
        x_rampa=X_RAMPA, L_cono=L_CONO,
        #: I fori si STAMPANO SOTTOMISURA e si portano a quota con la punta:
        #: un foro SLM esce piu' piccolo del nominale e con la parete rugosa,
        #: e maschiare dentro un foro rugoso significa un filetto storto.
        d_porta_aria=fa["punta"] - SOTTOMISURA_STAMPA,
        R_boss_aria=R_BOSS_ARIA,
        d_boss_aria=fa["d_guarnizione"] + 1.8e-3,
        x_porta_aria_fissa=X_PORTA_ARIA,
        d_filetto_gpl=fg["punta"] - SOTTOMISURA_STAMPA,
        L_filetto_gpl=fg["avvitamento"] + 2.0e-3)
    return inj, testa


#: Raggio di iniezione scelto in scripts/injector_50N.py: e' il valore che
#: rende C = 2.51 contro l'ottimo 2.5 di Holdeman CON il massimo numero di
#: getti compatibile con il 10 % di margine sul diametro minimo stampabile.
R_INIEZIONE = 3.10e-3


#: Setto fra il ramo di andata e quello di ritorno del circuito a U.
#:
#: Portato da 0.6 a 1.0 mm nella revisione 0.4.0, per due ragioni indipendenti.
#:
#: 1. INGEGNERISTICA. E' l'unica parete del motore la cui rottura e' SILENZIOSA.
#:    Se cede, l'acqua non esce: passa direttamente dall'andata al ritorno
#:    cortocircuitando il labbro, cioe' smette di raffreddare il punto piu'
#:    caldo del motore mentre portata, pressione e temperatura all'uscita
#:    restano quasi normali. Non c'e' strumento sul banco che lo veda. Il
#:    carico meccanico e' irrisorio (0.4 bar, meta' della caduta del circuito),
#:    quindi il margine non serve alla resistenza: serve alla POROSITA'. Un
#:    muro SLM in 316L sotto i 0.8 mm si stampa in una o due passate di laser e
#:    la sua tenuta al gas non e' garantita. Le pareti che separano acqua da
#:    gas (0.8 mm) hanno un guasto rumoroso e possono stare piu' sottili;
#:    questa no.
#:
#: 2. DI VERIFICABILITA'. A 0.6 mm il setto e' quattro voxel su una griglia da
#:    0.15 mm, e il genere topologico del pezzo NON convergeva: cambiava con la
#:    risoluzione, cioe' il modello non sapeva dire quanti fori avesse. Un
#:    dettaglio che il modello non riesce a risolvere e' un dettaglio che la
#:    macchina fatica a stampare: la non convergenza era essa stessa il
#:    risultato. A 1.0 mm sono sette voxel e la misura si stabilizza.
#:
#: Costo: 0.4 mm di raggio esterno e circa 3 g di acciaio.
SETTO_RITORNO = 1.0e-3

#: I due circuiti NON sono piu' definiti con quote assolute.
#:
#: Lo erano fino alla revisione 0.4.0 (camera 1.5-26.5 mm, gola 30.5-39.2 mm) e
#: quelle quote erano legate a un labbro che stava a x = 39.7 mm. Correggere la
#: pressione di bombola ha spostato il labbro a 42.09 mm: il circuito di gola
#: sarebbe finito 2.9 mm PRIMA della gola, cioe' avrebbe raffreddato tutto
#: tranne il punto piu' caldo del motore, e nessun controllo se ne sarebbe
#: accorto perche' la geometria resta perfettamente valida e stampabile.
#:
#: E' la stessa classe di difetto del raccordo da 1.5 mm che strozzava la gola:
#: un numero giusto per una versione precedente della geometria. La cura e'
#: strutturale, non numerica: le quote si derivano da (L_c, L_conv).
def circuiti(d: dict) -> list[CircuitoRaffreddamento]:
    """I due circuiti di raffreddamento, posizionati sulla geometria corrente.

    Ripartizione del carico (scripts/heat_map_50N.py, 3.72 kW totali):
      camera        2.32 kW su grande area, q = 0.90 MW/m2
      convergente   1.05 kW su 8 mm di corsa, q fino a 4.28 MW/m2
      plug          0.35 kW, e non e' raffreddabile: e' un corpo isolato

    In PARALLELO e non in serie: e' cio' che porta la perdita di carico
    complessiva da 2.5 a 0.77 bar, perche' la portata si divide e la caduta va
    come il quadrato della velocita'.
    """
    L_c, L_conv = d["L_c"], d["L_conv"]
    x_lip = L_c + L_conv

    #: Il circuito di gola deve arrivare al labbro ma non oltre: dopo il labbro
    #: il mantello non esiste piu' e i canali finirebbero nel vuoto.
    GOLA_MARGINE_LABBRO = 0.5e-3
    #: e comincia poco prima della fine del cilindro, cosi' copre tutto il
    #: convergente piu' un tratto di camera dove q sta gia' salendo.
    GOLA_ANTICIPO = 1.4e-3
    #: Fra i due circuiti resta pieno: due reticoli elicoidali che si
    #: incrociano producono pareti piu' sottili del passo della griglia, e la
    #: superficie esce aperta. Nella prima versione si accavallavano e uscivano
    #: 56 spigoli con un solo triangolo.
    SEPARAZIONE = 4.0e-3
    #: la piastra di testa resta piena: e' dove vanno gli iniettori.
    CAMERA_INIZIO = 1.5e-3

    x_gola_0 = L_c - GOLA_ANTICIPO
    x_gola_1 = x_lip - GOLA_MARGINE_LABBRO
    x_cam_1 = x_gola_0 - SEPARAZIONE
    if x_cam_1 - CAMERA_INIZIO < 10.0e-3:
        raise ValueError(
            f"camera cilindrica troppo corta ({(x_cam_1-CAMERA_INIZIO)*1e3:.1f} mm "
            "di corsa utile per i canali): la ripartizione dei due circuiti va "
            "ripensata, non allungata a forza."
        )

    #: Circuito CAMERA. 20 canali perche' sotto quel numero il setto fra canali
    #: supera i 120 K di salto (l'aletta conduce lateralmente: dT = q L^2/2kt).
    #: 1.0 mm e 3.5 m/s danno h = 12.6 kW/m2K: la parete lato acqua sta a 93 C,
    #: 50 K sotto l'ebollizione, e la perdita di carico e' 0.1 bar.
    #: passo 125 mm: autosostentamento vuole passo >= 2 pi R tan(alpha), e a
    #: R = 14.3 mm con 52 gradi di margine servono 115 mm.
    camera = CircuitoRaffreddamento(
        nome="camera", n_canali=20, lato=1.0e-3,
        parete_calda=1.2e-3, parete_fredda=1.2e-3,
        x_inizio=CAMERA_INIZIO, x_fine=x_cam_1, passo_elica=0.125, velocita=3.5,
    )
    #: Circuito GOLA + LABBRO. Canali piu' piccoli e veloci per h = 42.8
    #: kW/m2K, parete calda 0.8 mm per abbassare il salto di conduzione.
    #: A U: l'acqua scende al labbro, gira, e torna su uno strato piu' esterno,
    #: cosi' entrambi gli attacchi restano sulla parte cilindrica dove un foro
    #: radiale ha senso. Verso il labbro la parete e' conica e nessuna
    #: profondita' di foro raggiunge il collettore senza bucare il gas.
    gola = CircuitoRaffreddamento(
        nome="gola", n_canali=18, lato=0.6e-3,
        #: parete fredda 1.0 e non 1.4: il ramo a U occupa
        #: 0.8 + 0.6 + 1.0 + 0.6 + parete_fredda, e ogni decimo in piu' qui
        #: ispessisce il mantello su TUTTA la lunghezza del motore, perche' il
        #: mantello ha spessore unico. 1.0 mm e' il doppio del minimo SLM su
        #: una parete che separa acqua da aria: se cede si vede subito.
        parete_calda=0.8e-3, parete_fredda=1.0e-3,
        x_inizio=x_gola_0, x_fine=x_gola_1, passo_elica=0.125, velocita=12.0,
        ritorno=True, setto_ritorno=SETTO_RITORNO, sfalsamento_ritorno=3.0e-3,
    )
    return [camera, gola]


def punto_operativo():
    fuel = FuelSpec(composition={"C3H8": 1.0}, phase_at_injection="gas",
                    thermo_source="gri30.yaml")
    def op(m):
        return OperatingPoint(p_amb=101325.0, p_air_supply=P_SUP, p_fuel_supply=P_SUP,
                              fuel=fuel, T_air_in=293.0, T_fuel_in=283.0,
                              mdot_air_max=m, cd_injector_ox=0.75, cd_injector_fuel=0.75)
    m = brentq(lambda m: evaluate_l0(op(m), p_c=P_C, phi_core=0.9, f_film=0.0,
                                     thermal_severity=False).thrust - 50.0,
               1e-3, 0.2, xtol=1e-10)
    o = op(m)
    x = default_design_vector(p_c=P_C, phi_core=0.9, f_film=0.0, N_inj=12.0,
                              Dc_over_Dt=2.6, Lc_over_Dc=1.4, plug_trunc=0.8,
                              t_wall=0.0012, d_ox_ratio=0.05, fuel_vel_ratio=1.1,
                              conv_half_angle=0.65)
    return (o, *derive(x, o))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--passo", type=float, default=1.2e-4, help="passo griglia [m]")
    #: 0.5 mm e non 1.5. Un raccordo e' locale solo se i due corpi distano piu'
    #: del suo raggio: plug e labbro distano 1.78 mm, e con 1.5 mm il raccordo
    #: li univa ATTRAVERSO la gola, strozzandola del 25 %. Nessun altro
    #: controllo se ne accorgeva. Vedi engine.area_di_gola.
    ap.add_argument("--raccordo", type=float, default=5.0e-4)
    ap.add_argument("--out", type=Path, default=Path("runs/motore50"))
    ap.add_argument("--backend", default="numpy")
    a = ap.parse_args()

    from zefiro.sdf.core import backend
    xp = backend(a.backend)

    op, params, l0 = punto_operativo()
    d = params.derived
    print(f"MOTORE 50 N   p_c {P_C/1e5:.2f} bar   spinta {l0.thrust:.2f} N   "
          f"Isp {l0.Isp_s:.1f} s   eps {l0.epsilon:.3f}")
    print(f"  camera D {2*d['R_c']*1e3:.2f} mm  L {d['L_c']*1e3:.1f} mm   "
          f"labbro R {d['R_lip']*1e3:.2f} mm   gola {l0.A_t*1e6:.1f} mm2")

    t0 = time.perf_counter()
    inj, testa = dimensiona_iniettore(op, params, l0)
    print(f"  iniettore a getti trasversali: {inj.n_getti} getti da "
          f"{inj.d_getto*1e3:.3f} mm da R {inj.R_iniezione*1e3:.2f} mm, "
          f"condotto d'aria H {inj.altezza_anello*1e3:.3f} mm, "
          f"J {inj.J:.3f}, C {inj.C:.3f} (ottimo {inj.C_obiettivo})")
    for av in inj.avvertenze:
        print(f"    ATTENZIONE: {av}")
    print(f"  attacchi: aria {FILETTO_ARIA} (preforo stampato "
          f"{testa.d_porta_aria*1e3:.1f}, punta {BSPP[FILETTO_ARIA]['punta']*1e3:.2f}, "
          f"bozza {testa.d_boss_aria*1e3:.0f} mm a R {testa.R_boss_aria*1e3:.0f}, "
          f"avvitamento {(testa.R_boss_aria-testa.R_plenum)*1e3:.1f} mm)")
    print(f"            GPL {FILETTO_GPL} sulla faccia di monte (preforo "
          f"{testa.d_filetto_gpl*1e3:.1f}, punta {BSPP[FILETTO_GPL]['punta']*1e3:.2f}, "
          f"profondita' {testa.L_filetto_gpl*1e3:.0f} mm)")
    for nome, taglia, mdot, rho in (
            ("aria", FILETTO_ARIA, l0.mdot_air, 8.70),
            ("GPL", FILETTO_GPL, l0.mdot_fuel_core, 13.71)):
        v, dp = perdita_raccordo(mdot, rho, BSPP[taglia]["punta"] - 5.0e-3)
        print(f"    perdita nel raccordo {nome}: {v:.0f} m/s, {dp/1e5:.3f} bar "
              f"({dp/(0.15*params.free['p_c'])*100:.0f} % del salto d'iniezione)")

    m = costruisci(d, params.plug_contour_x, params.plug_contour_r,
                   circuiti(d), a.passo, raccordo=a.raccordo, xp=xp, testa=testa)
    print(f"\ngriglia {m.grid.shape} = {m.grid.n_voxels/1e6:.1f} Mvoxel, "
          f"{m.grid.memoria_mb():.0f} MB per campo, passo {a.passo*1e3:.3f} mm "
          f"({time.perf_counter()-t0:.1f} s)")
    for n in m.note:
        print(f"  ATTENZIONE: {n}")

    v_solido = m.solido.volume()
    v_canali = m.canali.volume()
    print(f"\nVOLUMI  solido {v_solido*1e6:8.3f} cm3   massa 316L "
          f"{v_solido*7990*1e3:6.1f} g")
    print(f"        canali {v_canali*1e6:8.3f} cm3   (cavita' del refrigerante)")

    print("\nesportazione...", flush=True)
    v, f = isosurface(m.solido)
    chiuso = is_watertight(f)
    v_mesh = mesh_volume(v, f)
    a.out.mkdir(parents=True, exist_ok=True)
    sha = write_stl(v, f, a.out / "motore50N.stl")
    print(f"  {len(f)} triangoli   chiuso: {chiuso}")
    print(f"  volume dalla mesh {v_mesh*1e6:.3f} cm3 contro {v_solido*1e6:.3f} "
          f"dal campo  (scarto {abs(v_mesh-v_solido)/v_solido:.2%})")
    print(f"  superficie {mesh_area(v, f)*1e4:.1f} cm2")
    print(f"  scritto {a.out/'motore50N.stl'}  sha {sha[:12]}")

    vc, fc = isosurface(m.canali)
    write_stl(vc, fc, a.out / "canali50N.stl")
    print(f"  scritto {a.out/'canali50N.stl'} ({len(fc)} triangoli): la sola rete "
          "di raffreddamento, per guardarla da sola")

    print("\nSTAMPABILITA'  (costruzione lungo l'asse, in piedi sulla piastra)")
    from zefiro.sdf.engine import area_di_gola
    from zefiro.sdf.printability import (
        analizza_sbalzi, polvere_evacuabile, sbalzi_interni,
    )
    rap = analizza_sbalzi(v, f, (1.0, 0.0, 0.0))
    rint = sbalzi_interni(v, f, m.cavita_gas, m.canali, (1.0, 0.0, 0.0))
    print(f"  da supportare {rap.area_da_supportare*1e4:.2f} cm2 su "
          f"{rap.area_totale*1e4:.1f} ({rap.frazione_da_supportare:.2%}); "
          f"appoggiati sulla piastra {rap.area_sulla_piastra*1e4:.2f} cm2")
    print("  area rivolta in basso per fascia di angolo [cm2]: "
          + str({k: round(x*1e4, 2) for k, x in rap.istogramma.items()}))
    if rint is not None:
        print(f"  DENTRO i canali, dove nessun supporto e' rimovibile: "
              f"{rint.area_da_supportare*1e4:.2f} cm2")
    from zefiro.sdf.clearances import rapporto_spessori
    from zefiro.sdf.meshing import connected_components
    comp = connected_components(v, f)
    # ATTENZIONE A COME SI LEGGE QUESTO NUMERO. Da quando i circuiti sono
    # davvero stagni, ognuno di essi e' una CAVITA' CHIUSA, e marching cubes
    # ne restituisce la superficie come componente a se'. Un pezzo sano ha
    # quindi 1 superficie esterna + una per ogni cavita' sigillata, non 1.
    # Il confronto giusto e' con il numero di domini stagni che la prova di
    # tenuta trova piu' avanti.
    print(f"  {len(comp)} superfici chiuse (1 esterna + le cavita' sigillate); "
          f"superficie chiusa: {chiuso}")
    ok_polvere, n_comp, sacche = polvere_evacuabile(m.canali, m.solido, m.grid)
    print(f"  polvere evacuabile: {ok_polvere} ({sacche} sacche chiuse su "
          f"{n_comp} componenti di vuoto)")
    x_lip = d["L_c"] + d["L_conv"]
    A, At = area_di_gola(m, x_lip, params.plug_contour_r[0], d["R_lip"])
    print(f"  AREA DI GOLA misurata sul solido {A*1e6:.2f} mm2 contro "
          f"{At*1e6:.2f} teorici ({A/At-1:+.1%})")

    print("\nPROVA DI TENUTA (flood-fill sul vuoto, attacchi tappati)")
    from zefiro.sdf.engine import sonde_motore, tappi
    from zefiro.sdf.tenuta import Sonda, verifica_tenuta
    tappato = m.solido.union(tappi(m, d))
    sonde = sonde_motore(m, d)
    # Cio' che DEVE comunicare: il condotto d'aria sbocca in camera, i getti
    # attraversano il condotto, e la camera e' aperta all'esterno dal labbro.
    attese = {("aria", "gas"), ("aria", "gpl"), ("gas", "gpl"),
              ("esterno", "gas"), ("aria", "esterno"), ("esterno", "gpl")}
    # I due strati di un circuito a U SONO lo stesso circuito: si uniscono
    # nell'inversione. Se NON comunicano il circuito e' interrotto, ed e' un
    # guasto grave quanto una perdita: l'acqua entrerebbe e non uscirebbe.
    for c in circuiti(d):
        if c.ritorno:
            attese.add(tuple(sorted((f"acqua_{c.nome}_andata",
                                     f"acqua_{c.nome}_ritorno"))))
    rap = verifica_tenuta(tappato, m.grid, sonde, attese=attese)
    print(rap.riassunto())
    print(f"  stagno: {rap.stagno}")

    print("\nPROVA DI TENUTA DEL SOLO CIRCUITO GPL (getti tappati)")
    print("  Con i cinque getti chiusi il GPL non deve toccare piu' NIENTE.")
    print("  E' la prova che verifica la parete da 1.0 mm fra il plenum del")
    print("  GPL e il condotto d'aria, che nessun'altra prova vede.")
    from zefiro.sdf.engine import _getti_gpl
    # I tappi dei getti vanno un voxel PIU' GRANDI del getto, non piu'
    # piccoli: con offset negativo restava un anello di vuoto spesso un voxel
    # tutto attorno, il getto non era tappato e la prova non provava niente.
    solo_gpl = tappato.union(_getti_gpl(m.grid, m.testa, xp).offset(m.grid.spacing))
    attese2 = {c for c in attese if "gpl" not in c}
    rap2 = verifica_tenuta(solo_gpl, m.grid, sonde, attese=attese2)
    print(rap2.riassunto())
    print(f"  stagno: {rap2.stagno}")

    print("\nSPESSORI DI PARETE fra vuoti che NON devono comunicare")
    import itertools
    # Le coppie che DEVONO toccarsi. Dichiararle e' obbligatorio: se non lo
    # fossero, il motore funzionante risulterebbe difettoso e si smetterebbe di
    # guardare i risultati - che e' il modo piu' comune in cui un controllo
    # automatico muore.
    #   aria-gas   il condotto d'aria sbocca in camera: e' l'iniezione.
    #   gpl-aria   i cinque getti attraversano il condotto: e' il mescolamento.
    #   gpl-gas    a valle dello sbocco sono lo stesso fluido.
    voluti = {("fori_iniezione", "gas"), ("aria", "gas"), ("aria", "gpl"),
              ("gas", "gpl")}
    for c in m.circuiti:
        for k in (0, 1):
            strato = f"{c.nome}_ritorno" if (c.ritorno and k == 1) else f"{c.nome}_andata"
            for j in range(PORTE_PER_COLLETTORE):
                voluti.add(tuple(sorted((f"{c.nome}_attacco_{k}_{j}", strato))))
    nomi = sorted(m.parti)
    coppie = [(a, b, MIN_WALL_SLM) for a, b in itertools.combinations(nomi, 2)
              if tuple(sorted((a, b))) not in voluti]
    # Il controllo costa una trasformata di distanza per parte: su una griglia
    # da venti milioni di voxel sono minuti. E' una verifica GEOMETRICA, e a
    # 0.15 mm risolve gia' tutto quello che c'e' da risolvere: sopra quella
    # soglia si dice di rifarla piu' grossa invece di far aspettare.
    if m.grid.n_voxels > 30_000_000:
        print(f"    saltato: griglia da {m.grid.n_voxels/1e6:.0f} Mvoxel. "
              "Rilancia con un passo piu' grosso per il controllo degli spessori.")
        return 0
    esiti = [x for x in rapporto_spessori(m, coppie) if x.minima != float("inf")]
    guai = [x for x in esiti if x.esito != "ok"]
    print(f"  {len(coppie)} coppie controllate, minimo richiesto "
          f"{MIN_WALL_SLM*1e3:.2f} mm (tolleranza = passo {m.grid.spacing*1e3:.2f} mm)")
    if guai:
        for x in sorted(guai, key=lambda z: z.minima)[:8]:
            print(f"    {x.a:<24}{x.b:<24}{x.minima*1e3:7.3f} mm   {x.esito}")
    else:
        print("    nessuna coppia sotto il minimo")
    print(f"  parete piu' sottile fra due vuoti: "
          f"{min(x.minima for x in esiti)*1e3:.3f} mm")

    print("\nRAFFREDDAMENTO (dal dimensionamento, vedi heat_map_50N.py)")
    tot = 0.0
    for c in m.circuiti:
        Q = c.n_canali * c.lato**2 * c.velocita
        tot += Q
        print(f"  {c.nome:<8}{c.n_canali:3d} canali {c.lato*1e3:.1f} mm a "
              f"{c.velocita:4.1f} m/s -> {Q*6e4:5.2f} l/min   "
              f"parete calda {c.parete_calda*1e3:.1f} mm")
    print(f"  {'totale':<8}{tot*6e4:33.2f} l/min")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
