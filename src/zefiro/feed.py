"""Impianto di alimentazione: compressore, serbatoio d'aria, bombola GPL.

Questo modulo esiste perche' con Zefiro l'impianto NON e' un dettaglio a monte:
e' il vincolo che dimensiona il motore. La portata d'aria disponibile fissa
l'area di gola, e la bombola GPL fissa il tetto di pressione di camera.

Provenienza dei dati:
  * proprieta' di propano e n-butano: CoolProp, equazioni di stato di
    riferimento (propano: Lemmon, McLinden & Wagner, JCED 2009; n-butano:
    Buecker & Wagner, JPCRD 2006). Nessun valore tabellato a mano.
  * aria: gas ideale con cp da CoolProp alle condizioni di aspirazione.

Cio' che NON e' calcolabile e resta un dato da fornire e' segnato come
argomento obbligatorio senza default.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from CoolProp.CoolProp import PropsSI

R_AIR = 287.052874       # J/(kg K), costante specifica dell'aria secca
HP_TO_W = 745.699872     # 1 cavallo vapore britannico (HP)
CV_TO_W = 735.49875      # 1 cavallo vapore metrico (CV) - in Italia spesso e' questo


# --------------------------------------------------------------------------- #
# Compressore
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class CompressorBounds:
    """Limiti TERMODINAMICI di portata, non prestazioni di una macchina reale."""
    mdot_isothermal: float   # kg/s, limite superiore assoluto (intercooling perfetto)
    mdot_adiabatic: float    # kg/s, monostadio adiabatico ideale
    w_isothermal: float      # J/kg
    w_adiabatic: float       # J/kg
    shaft_power: float       # W


def compressor_bounds(
    shaft_power: float, p_in: float, p_out: float, T_in: float
) -> CompressorBounds:
    """Portata massima compatibile con la potenza all'albero.

    Due limiti, entrambi esatti e non discutibili:

      * **isotermo**: w = R T ln(p2/p1). E' il lavoro MINIMO possibile per
        comprimere un gas ideale, raggiungibile solo con raffreddamento perfetto
        e infiniti stadi. Da' quindi la portata MASSIMA immaginabile.
      * **adiabatico monostadio**: w = cp T [(p2/p1)^((k-1)/k) - 1]. E' il
        lavoro di un compressore ideale senza alcun raffreddamento.

    Una macchina reale sta in mezzo, con un rendimento isotermo che per un
    compressore a pistoni di questa taglia e' tipicamente ben sotto l'unita'.
    Questa funzione NON stima quel rendimento: da' i due estremi e lascia che
    sia la targhetta (o una misura) a dire dove sta la macchina.
    """
    if p_out <= p_in:
        raise ValueError("p_out deve essere maggiore di p_in.")
    cp = PropsSI("C", "T", T_in, "P", p_in, "Air")
    cv = PropsSI("CVMASS", "T", T_in, "P", p_in, "Air")
    k = cp / cv
    pr = p_out / p_in

    w_iso = R_AIR * T_in * math.log(pr)
    w_ad = cp * T_in * (pr ** ((k - 1.0) / k) - 1.0)
    return CompressorBounds(
        mdot_isothermal=shaft_power / w_iso,
        mdot_adiabatic=shaft_power / w_ad,
        w_isothermal=w_iso,
        w_adiabatic=w_ad,
        shaft_power=shaft_power,
    )


def mdot_from_fad(fad_litres_per_min: float, T_ref: float = 293.15,
                  p_ref: float = 101325.0) -> float:
    """Converte la portata dichiarata in targa ("aria resa", l/min) in kg/s.

    La FAD e' un volume di aria ASPIRATA, quindi va valutata alle condizioni di
    aspirazione, non a quelle di mandata. E' l'errore piu' comune su questo
    numero e vale un fattore pari al rapporto di compressione.
    """
    rho = p_ref / (R_AIR * T_ref)
    return fad_litres_per_min * 1.0e-3 / 60.0 * rho


# --------------------------------------------------------------------------- #
# Serbatoio d'aria: funzionamento a raffica
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class TankBlowdown:
    mass_initial: float      # kg
    mass_usable_isothermal: float
    mass_usable_adiabatic: float
    T_final_adiabatic: float  # K
    volume: float
    p_initial: float
    p_final: float


def tank_blowdown(volume: float, p_initial: float, p_final: float,
                  T_initial: float) -> TankBlowdown:
    """Massa d'aria estraibile da un serbatoio fra due pressioni.

    Due limiti, di nuovo entrambi esatti:

      * **isotermo**: lo svuotamento e' lento e le pareti del serbatoio
        reintegrano il calore. m_usabile = (p1 - p2) V / (R T).
      * **adiabatico**: lo svuotamento e' rapido, il gas che resta si espande
        e si raffredda: T2 = T1 (p2/p1)^((k-1)/k). Resta piu' massa dentro,
        quindi se ne estrae MENO.

    Una raffica di pochi secondi da un serbatoio d'acciaio sta vicino
    all'adiabatico: usare il limite isotermo sovrastima la durata.
    """
    if p_final >= p_initial:
        raise ValueError("p_final deve essere minore di p_initial.")
    cp = PropsSI("C", "T", T_initial, "P", p_initial, "Air")
    cv = PropsSI("CVMASS", "T", T_initial, "P", p_initial, "Air")
    k = cp / cv

    m1 = p_initial * volume / (R_AIR * T_initial)
    m2_iso = p_final * volume / (R_AIR * T_initial)
    T2 = T_initial * (p_final / p_initial) ** ((k - 1.0) / k)
    m2_ad = p_final * volume / (R_AIR * T2)
    return TankBlowdown(
        mass_initial=m1,
        mass_usable_isothermal=m1 - m2_iso,
        mass_usable_adiabatic=m1 - m2_ad,
        T_final_adiabatic=T2,
        volume=volume,
        p_initial=p_initial,
        p_final=p_final,
    )


def burst_duration(blowdown: TankBlowdown, mdot: float,
                   mdot_recharge: float = 0.0) -> tuple[float, float]:
    """(durata_adiabatica, durata_isoterma) in secondi.

    `mdot_recharge` e' la portata che il compressore reintegra durante la
    raffica: allunga la durata ma non la rende infinita finche' e' inferiore a
    `mdot`.
    """
    net = mdot - mdot_recharge
    if net <= 0.0:
        return (math.inf, math.inf)
    return (blowdown.mass_usable_adiabatic / net, blowdown.mass_usable_isothermal / net)


# --------------------------------------------------------------------------- #
# Bombola GPL
# --------------------------------------------------------------------------- #
COOLPROP_NAME = {"C3H8": "Propane", "C4H10": "n-Butane", "IC4H10": "IsoButane"}


def saturation_pressure(composition: dict[str, float], T: float) -> float:
    """Pressione della bombola per legge di Raoult [Pa].

        p = SOMMA_i  x_i * p_sat,i(T)

    Ipotesi: miscela liquida IDEALE. Propano e butano sono idrocarburi
    omologhi con struttura simile, quindi Raoult e' un'approssimazione
    ragionevole; lo scostamento reale e' dell'ordine di qualche punto
    percentuale e va tenuto presente se si usa questa relazione per dedurre la
    composizione (vedi `composition_from_pressure`).

    Nota fisica che rende utile questa funzione: la pressione di bombola NON
    dipende da quanto liquido e' rimasto. Finche' c'e' liquido, dipende solo
    dalla temperatura e dalla composizione. Una bombola quasi vuota ha la
    stessa pressione di una piena.
    """
    return sum(
        x * PropsSI("P", "T", T, "Q", 0, COOLPROP_NAME[s])
        for s, x in composition.items()
    )


def composition_from_pressure(p: float, T: float,
                              heavy: str = "C4H10") -> float:
    """Frazione molare di propano dedotta da (pressione, temperatura) misurate.

    Inversione di Raoult per una miscela binaria propano / butano:

        x_C3H8 = (p - p_sat,butano(T)) / (p_sat,propano(T) - p_sat,butano(T))

    E' il modo piu' economico che hai per sapere che cosa c'e' davvero nella
    bombola: **un manometro e un termometro**, senza analisi di laboratorio.
    Misura entrambi nello stesso momento, con la bombola in equilibrio termico
    (non subito dopo un prelievo).

    Ritorna un valore fuori da [0, 1] se la misura non e' compatibile con una
    binaria propano/butano: e' un'informazione, non un errore da nascondere.
    """
    p_light = PropsSI("P", "T", T, "Q", 0, "Propane")
    p_heavy = PropsSI("P", "T", T, "Q", 0, COOLPROP_NAME[heavy])
    return (p - p_heavy) / (p_light - p_heavy)


@dataclass(frozen=True)
class AutorefrigerationResult:
    heat_rate: float          # W sottratti dal prelievo di vapore
    thermal_capacity: float   # J/K della bombola (liquido + involucro)
    dT_dt: float              # K/s
    dT_over_burn: float       # K nella durata indicata
    p_start: float            # Pa
    p_end: float              # Pa


def autorefrigeration(
    mdot_fuel: float,
    T: float,
    composition: dict[str, float],
    liquid_mass: float,
    shell_mass: float,
    shell_specific_heat: float,
    burn_time: float,
) -> AutorefrigerationResult:
    """Raffreddamento della bombola dovuto al prelievo di vapore.

    Il vapore prelevato deve essere prodotto vaporizzando liquido, e il calore
    latente lo si prende dalla bombola stessa. Bilancio (trascurando lo scambio
    con l'ambiente, che a queste potenze e' piccolo e va nel verso favorevole):

        (m_liq c_liq + m_involucro c_involucro) dT/dt = - mdot h_fg(T)

    Trascurare lo scambio con l'ambiente rende il risultato CONSERVATIVO: la
    bombola reale si raffredda un po' meno di cosi'.

    `liquid_mass` e `shell_mass` non sono inventabili ma sono facilissimi da
    ottenere: la tara e' stampigliata sul collare della bombola, e il peso
    totale lo dai con una bilancia. m_liquido = peso_totale - tara.
    """
    h_fg = sum(
        x * (PropsSI("H", "T", T, "Q", 1, COOLPROP_NAME[s])
             - PropsSI("H", "T", T, "Q", 0, COOLPROP_NAME[s]))
        for s, x in composition.items()
    )
    c_liq = sum(
        x * PropsSI("C", "T", T, "Q", 0, COOLPROP_NAME[s])
        for s, x in composition.items()
    )
    q = mdot_fuel * h_fg
    capacity = liquid_mass * c_liq + shell_mass * shell_specific_heat
    dT_dt = q / capacity
    dT = dT_dt * burn_time
    return AutorefrigerationResult(
        heat_rate=q,
        thermal_capacity=capacity,
        dT_dt=dT_dt,
        dT_over_burn=dT,
        p_start=saturation_pressure(composition, T),
        p_end=saturation_pressure(composition, T - dT),
    )


# --------------------------------------------------------------------------- #
# Pressioni di alimentazione: DERIVATE, non scelte
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class SupplyPressures:
    """Le due pressioni a monte dei riduttori, alla fine della raffica."""
    p_fuel_supply: float     # Pa assoluti
    p_air_supply: float      # Pa assoluti, fondo scarico del serbatoio
    p_c_max: float           # Pa, tetto di pressione di camera
    T_bottle: float          # K, temperatura di progetto della bombola
    binding: str             # "bombola" oppure "serbatoio"


def supply_pressures(
    composition: dict[str, float],
    T_bottle: float,
    p_tank_max: float,
    dp_fraction: float,
) -> SupplyPressures:
    """Deriva le pressioni di alimentazione da (composizione, temperatura).

    Il punto che questa funzione esiste per rendere impossibile da sbagliare:
    **la pressione del GPL non e' un parametro di progetto, e' una proprieta'
    termodinamica**. Vale p_sat(T) e nient'altro. Scriverla a mano in un file di
    configurazione significa poter scrivere un numero che la fisica non produce
    a nessuna temperatura raggiungibile, e non accorgersene.

    Da qui discendono, in cascata:

      * il tetto di camera  p_c,max = p_sat / (1 + f_Dp) : sopra, il GPL non ha
        abbastanza salto per entrare in modo stabile;
      * il fondo scarico del serbatoio d'aria, che si pone UGUALE a p_sat. Piu'
        in alto si spreca aria (e quindi durata) senza guadagnare nulla, perche'
        il tetto lo fissa comunque il GPL; piu' in basso il riduttore dell'aria
        perde il controllo prima di quello del GPL, e la miscela va ricca in
        modo incontrollato proprio a fine raffica.

    Se la bombola e' piu' calda del serbatoio (p_sat > p_tank_max) il vincolo si
    scambia: comanda il serbatoio. Il campo `binding` dice quale dei due.

    `dp_fraction` e' la frazione di p_c richiesta come salto d'iniezione
    (0.15 in Zefiro, vedi opt.objectives.MIN_INJECTOR_DP_FRACTION): sotto quella
    soglia l'iniettore non disaccoppia piu' l'alimentazione dalla camera e il
    sistema puo' entrare in oscillazione di tipo chugging.
    """
    p_sat = saturation_pressure(composition, T_bottle)
    if p_sat >= p_tank_max:
        p_sup = p_tank_max
        binding = "serbatoio"
    else:
        p_sup = p_sat
        binding = "bombola"
    return SupplyPressures(
        p_fuel_supply=p_sat,
        p_air_supply=p_sup,
        p_c_max=p_sup / (1.0 + dp_fraction),
        T_bottle=T_bottle,
        binding=binding,
    )
