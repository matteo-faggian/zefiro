"""Aerospike a espansione esterna, coppia gas-gas, raffreddamento ad acqua.

E' l'architettura di Zefiro, ma qui NON e' un progetto: e' una regola per ogni
quota. La differenza si vede dal fatto che spinta, propellenti e temperatura di
esercizio sono liberi, e trenta numeri che prima erano scritti a mano ora
escono da un conto o da una decisione dichiarata.

COME LEGGERE QUESTO MODULO. Ogni quota passa da `reg.aggiungi(Scelta(...))` con
la sua origine. Cercare `Origine.DECISA` da' l'elenco esatto di cio' che in
questo motore e' ancora giudizio umano: sono quelle le righe da attaccare per
migliorare il modello, e sono quelle che vanno rimesse in discussione quando i
requisiti cambiano molto.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

from scipy.optimize import brentq

from zefiro import cooling
from zefiro.carico import mappa_di_calore, potenza_per_zona, q_massimo_per_zona
from zefiro.feed import supply_pressures, tank_blowdown
from zefiro.geometry.parameters import DESIGN_BOUNDS, default_design_vector, derive
from zefiro.injector import (
    BSPP, SOTTOMISURA_STAMPA, Iniettore, perdita_raccordo, progetta as progetta_iniettore,
)
from zefiro.l0.cycle import evaluate_l0
from zefiro.schemas import FuelSpec, MissingDatum, OperatingPoint
from zefiro.sdf import engine
from zefiro import injector as injector_module
from zefiro.sdf.engine import (
    CircuitoRaffreddamento, CollettoreAcqua, FessuraFilm, InterfacciaMontaggio,
    Portagomma, TestaIniezione)
from zefiro.sdf.printability import passo_elica_minimo
from zefiro.sintesi.esito import Inviluppo
from zefiro.sintesi.provenienza import Origine, Registro, Scelta
from zefiro.sintesi.requisiti import Requisiti

NOME = "aerospike-gas-gas"
DESCRIZIONE = ("aerospike a espansione esterna, iniezione a getti trasversali, "
               "raffreddamento ad acqua in due circuiti paralleli")

#: --- costanti FISICHE e di NORMA, non decisioni ------------------------------
MW_ARIA = 28.9649          # kg/kmol, coerente con units.AIR_MOLE_FRACTIONS
GAMMA_ARIA = 1.400
#: Conducibilita' dell'AISI 316L a temperatura di parete [W/(m K)]. E' un dato
#: di materiale, non una scelta; resta un TODO confermarlo sulla lega del
#: fornitore (TODO J/6): l'intervallo di letteratura e' 14-19 a 500-800 K.
K_316L = 15.0

#: --- soglie di PROGETTO, decise, con la loro ragione -------------------------
#: Rapporto fra area di passaggio in camera e area di gola. Con 6 il numero di
#: Mach in camera sta sotto 0.1: la caduta di pressione totale fra testa e gola
#: e' allora sotto l'1 %, e il riscaldamento della parete di camera resta
#: governato dall'area e non dalla velocita'. Sotto 4 la camera comincia a
#: comportarsi come un condotto e Bartz in camera sale in fretta.
RAPPORTO_CONTRAZIONE = 6.0
#: Salto massimo di temperatura ammesso ATTRAVERSO la parete calda [K].
#: Non e' una soglia di resistenza: e' la tensione termica. Con
#: sigma ~ E alpha dT / (1-nu), 250 K in 316L danno gia' oltre 1 GPa, dieci
#: volte lo snervamento. La parete plasticizza a ogni accensione (ratcheting):
#: e' cio' che uccide le camere rigenerative, ed e' la ragione per cui la
#: parete calda si fa la piu' sottile che il processo consente e non la piu'
#: spessa che la pressione consentirebbe.
DELTA_T_PARETE_MAX = 250.0
#: Margine all'ebollizione lato acqua [K].
MARGINE_EBOLLIZIONE = 40.0
#: Salto di temperatura dell'acqua fra ingresso e uscita [K]. Deciso: piu'
#: piccolo vuole piu' portata (e piu' perdita di carico), piu' grande avvicina
#: l'uscita alla saturazione.
DELTA_T_ACQUA = 30.0
#: Perdita di carico ammessa sul circuito dell'acqua [Pa]. La rete del banco
#: da' circa 4 bar: si spende meno di un quarto.
DP_ACQUA_MAX = 1.0e5

#: Quanto della caduta totale ammessa puo' consumare UN SOLO tratto di canali.
#: Serve a scegliere fra canali larghi (poca caduta, ma collettori enormi per
#: distribuirvi l'acqua in modo uniforme) e canali stretti (piu' caduta, ma
#: collettori piccoli). Non e' una preferenza: e' il modo di spendere il budget
#: dove compra qualcosa.
FRAZIONE_DP_RAMO = 0.35

#: Quanto puo' valere la caduta lungo un anello di collettore rispetto a quella
#: di un canale del tratto che alimenta. E' il vincolo che dimensiona gli anelli:
#: con 0.10, fra il canale piu' favorito e il piu' sfavorito la portata
#: differisce del 5 % e h del 4 %. Non e' una velocita' scelta a mano: e' la
#: conseguenza che si vuole evitare, scritta come numero.
FRAZIONE_MALDISTRIBUZIONE = 0.10

#: Altezza radiale minima di un vano anulare [m]. NON e' un vincolo di processo:
#: l'SLM fa vani piu' bassi. E' un vincolo di VERIFICABILITA'. Un vano da 1.3 mm
#: con pareti da 1 mm, campionato su una griglia da 0.35 mm, sono tre voxel e
#: mezzo di vuoto fra tre di materiale: la prova di tenuta lo trovava aperto e
#: non si poteva distinguere un difetto vero dalla discretizzazione. E costa
#: quasi niente, perche' a sezione fissata la massa del collare va come L per h,
#: cioe' come l'area: alzare il vano e accorciarlo non pesa.
H_MIN_ANELLO = 2.5e-3
#: Margine sul diametro minimo stampabile quando quel dato NON e' misurato.
MARGINE_FORO = 1.10


@dataclass
class Dimensionamento:
    """Tutto cio' che serve per costruire la geometria, piu' come ci si e'
    arrivati."""
    registro: Registro
    requisiti: Requisiti
    op: OperatingPoint
    x: Any
    params: Any
    l0: Any
    iniettore: Iniettore
    testa: TestaIniezione
    circuiti: list[CircuitoRaffreddamento]
    inviluppo: Inviluppo
    collettore: CollettoreAcqua | None = None
    montaggio: InterfacciaMontaggio | None = None
    fessura: FessuraFilm | None = None
    carico: dict[str, float] = field(default_factory=dict)
    plug: Any = None            # zefiro.film.TransitorioPlug
    film: Any = None            # zefiro.film.Film, None se non serve
    note: list[str] = field(default_factory=list)


# --------------------------------------------------------------------------- #
# 0. applicabilita'
# --------------------------------------------------------------------------- #
def applicabile(req: Requisiti) -> list[str]:
    """Perche' questa architettura NON va bene per questi requisiti.

    Lista vuota = va bene. Il sintetizzatore la usa per scegliere, e per dire
    all'utente perche' ha rifiutato invece di limitarsi a rifiutare.
    """
    m = []
    if req.impianto.ossidante != "aria":
        m.append(f"ossidante {req.impianto.ossidante!r}: questa famiglia e' gas-gas ad aria")
    if not 5.0 <= req.spinta <= 2000.0:
        m.append(
            f"spinta {req.spinta:.0f} N fuori dall'intervallo in cui l'architettura e' "
            "stata verificata (5-2000 N). Sotto, la gola diventa piu' piccola del "
            "minimo stampabile; sopra, il raffreddamento ad acqua di banco non basta "
            "e serve un rigenerativo vero")
    return m


# --------------------------------------------------------------------------- #
# 1. punto operativo: le pressioni sono termodinamica, non progetto
# --------------------------------------------------------------------------- #
def _punto_operativo(req: Requisiti, reg: Registro,
                     f_film: float = 0.0) -> tuple[OperatingPoint, float, float]:
    imp = req.impianto
    sup = supply_pressures(dict(imp.combustibile), imp.T_bombola_min,
                           imp.p_serbatoio_max, req.frazione_dp_iniezione)
    reg.aggiungi(Scelta(
        "p_combustibile", sup.p_fuel_supply, "Pa", Origine.DERIVATA,
        f"tensione di vapore della miscela a {imp.T_bombola_min-273.15:.1f} C. "
        "Non e' un parametro: una bombola eroga p_sat(T) e nient'altro",
        fonte="feed.supply_pressures, CoolProp"))
    reg.aggiungi(Scelta(
        "p_ossidante", sup.p_air_supply, "Pa", Origine.DERIVATA,
        f"fondo scarico del serbatoio, posto uguale a p_sat perche' il vincolo lo "
        f"detta {sup.binding}: piu' in alto si spreca aria senza alzare p_c, piu' in "
        "basso il riduttore dell'aria perde il controllo prima di quello del GPL",
        fonte="feed.supply_pressures"))
    p_c = sup.p_c_max
    reg.aggiungi(Scelta(
        "p_camera", p_c, "Pa", Origine.DERIVATA,
        f"p_alimentazione / (1 + {req.frazione_dp_iniezione:.2f}): sopra questa "
        "pressione l'iniettore non ha piu' il salto che disaccoppia camera e "
        "alimentazione, e il sistema puo' entrare in chugging"))

    fuel = FuelSpec(composition=dict(imp.combustibile), phase_at_injection="gas",
                    thermo_source=imp.meccanismo)

    def op_con(m: float) -> OperatingPoint:
        return OperatingPoint(
            p_amb=req.p_ambiente, p_air_supply=sup.p_air_supply,
            p_fuel_supply=sup.p_fuel_supply, fuel=fuel,
            T_air_in=imp.T_ossidante_iniezione, T_fuel_in=imp.T_combustibile_iniezione,
            mdot_air_max=m, cd_injector_ox=imp.cd_ossidante,
            cd_injector_fuel=imp.cd_combustibile)

    def scarto(m: float) -> float:
        return evaluate_l0(op_con(m), p_c=p_c, phi_core=req.phi, f_film=f_film,
                           thermal_severity=False).thrust - req.spinta

    lo, hi = 1.0e-5, 5.0
    if scarto(lo) > 0.0 or scarto(hi) < 0.0:
        raise MissingDatum(
            f"nessuna portata d'aria fra {lo*1e3:.2f} e {hi*1e3:.0f} g/s da' "
            f"{req.spinta:.1f} N a p_c = {p_c/1e5:.2f} bar.")
    mdot = brentq(scarto, lo, hi, xtol=1.0e-12)
    reg.aggiungi(Scelta(
        "mdot_ossidante", mdot, "kg/s", Origine.DERIVATA,
        f"la portata che da' esattamente {req.spinta:.1f} N a p_c = {p_c/1e5:.3f} bar "
        f"con phi = {req.phi} e {f_film:.1%} di combustibile al film",
        fonte="l0.cycle.evaluate_l0"))
    return op_con(mdot), p_c, mdot


def _durata(req: Requisiti, reg: Registro, mdot: float) -> float:
    imp = req.impianto
    bd = tank_blowdown(imp.volume_serbatoio, imp.p_serbatoio_max,
                       reg.valore("p_ossidante"), imp.T_serbatoio)
    #: LIMITE ADIABATICO, non isotermo. Una raffica di pochi secondi da un
    #: serbatoio d'acciaio sta vicino all'adiabatico: il gas che resta dentro si
    #: espande e si raffredda, quindi ne resta di piu' e se ne estrae di meno.
    #: Usare l'isotermo sovrastima la durata del 36 %.
    massa = bd.mass_usable_adiabatic
    reg.aggiungi(Scelta(
        "massa_ossidante_utile", massa, "kg", Origine.DERIVATA,
        f"svuotamento adiabatico del serbatoio da {imp.p_serbatoio_max/1e5:.1f} bar "
        f"al fondo scarico. La ricarica del compressore NON e' contata: e' un dato "
        "non ancora misurato, e quando arrivera' sara' margine",
        fonte="feed.tank_blowdown"))
    durata = massa / mdot
    reg.aggiungi(Scelta(
        "durata_raffica", durata, "s", Origine.DERIVATA,
        "massa utile diviso portata. Cresce se la bombola e' piu' FREDDA, perche' "
        "il fondo scarico scende con p_sat e resta piu' aria estraibile"))
    return durata


# --------------------------------------------------------------------------- #
# 2. vettore di progetto: quattro quote derivate, quattro decise, quattro inerti
# --------------------------------------------------------------------------- #
def _vettore_di_progetto(req: Requisiti, reg: Registro, op: OperatingPoint,
                         p_c: float, mdot: float, f_film: float = 0.0):
    proc = req.processo

    #: PARETE STRUTTURALE. Derivata dal massimo fra due condizioni, e la
    #: pressione perde di quattro ordini di grandezza: a 6.4 bar su raggio 12 mm
    #: la tensione di cerchio e' 6 MPa contro i 110 dello snervamento a 600 C.
    #: Comanda il processo.
    t_press = p_c * 0.05 / (110.0e6 / 1.5)      # stima grossolana e generosa
    t_wall = max(2.0 * proc.parete_minima, t_press)
    t_wall = min(max(t_wall, DESIGN_BOUNDS["t_wall"][0]), DESIGN_BOUNDS["t_wall"][1])
    reg.aggiungi(Scelta(
        "t_parete_strutturale", t_wall, "m", Origine.DERIVATA,
        f"il doppio della parete minima di processo ({proc.parete_minima*1e3:.2f} mm). "
        f"La pressione ne chiederebbe {t_press*1e6:.1f} um: su questo motore la "
        "struttura non e' un problema, la fabbricazione si'"))

    def costruisci(dc_over_dt: float, lc_over_dc: float):
        x = default_design_vector(
            p_c=p_c, phi_core=req.phi, f_film=f_film,
            Dc_over_Dt=dc_over_dt, Lc_over_Dc=lc_over_dc,
            N_inj=12.0, theta_swirl=0.0, d_ox_ratio=0.05, fuel_vel_ratio=1.1,
            conv_half_angle=CONV_HALF_ANGLE, plug_trunc=PLUG_TRUNC, t_wall=t_wall)
        return derive(x, op, mdot_air=mdot)

    #: RAPPORTO DI CONTRAZIONE. Dc_over_Dt non si sceglie: si sceglie il
    #: rapporto fra area di passaggio ANULARE in camera e area di gola, e
    #: Dc_over_Dt e' quello che lo realizza. La differenza conta perche' la
    #: camera e' anulare: pi R_c^2 non e' l'area di passaggio, lo e'
    #: pi (R_c^2 - r_corpo^2), e il corpo centrale non e' piccolo.
    def scarto_contrazione(dc: float) -> float:
        params, l0 = costruisci(dc, 1.4)
        d = params.derived
        area = math.pi * (d["R_c"] ** 2 - d["r_centerbody"] ** 2)
        return area / l0.A_t - RAPPORTO_CONTRAZIONE

    lo, hi = DESIGN_BOUNDS["Dc_over_Dt"]
    if scarto_contrazione(lo) * scarto_contrazione(hi) > 0.0:
        raise MissingDatum(
            f"nessun rapporto Dc/Dt fra {lo} e {hi} da' una contrazione di "
            f"{RAPPORTO_CONTRAZIONE}: la camera anulare non ci sta.")
    dc_over_dt = brentq(scarto_contrazione, lo, hi, xtol=1.0e-10)
    reg.aggiungi(Scelta(
        "Dc_su_Dt", dc_over_dt, "-", Origine.DERIVATA,
        f"il valore che porta l'area ANULARE di camera a {RAPPORTO_CONTRAZIONE} volte "
        "quella di gola, cioe' Mach 0.1 in camera. Con l'area del cerchio invece che "
        "della corona si sbaglierebbe di un fattore 1.1"))

    #: LUNGHEZZA DI CAMERA. Si sceglie L*, non L_c: L* = V_c/A_t e' la grandezza
    #: che governa il tempo di permanenza, e L_c e' quella che lo realizza data
    #: la contrazione. Il valore minimo viene da un confronto di TEMPI, che e'
    #: la ragione fisica: la permanenza deve stare molto sopra la somma di
    #: mescolamento e chimica.
    lc = LC_SU_DC
    params, l0 = costruisci(dc_over_dt, lc)
    reg.aggiungi(Scelta(
        "Lc_su_Dc", lc, "-", Origine.DECISA,
        f"L* che ne risulta = {params.derived['L_star']*1e2:.0f} cm. Il criterio "
        "vero e' il rapporto fra tempo di permanenza e tempi di mescolamento e "
        "chimica, che viene verificato a valle: se scende sotto 3 il progetto e' "
        "rifiutato",
        scartate=("piu' corta: meno massa ma margine di permanenza sotto 3",
                  "piu' lunga: piu' superficie da raffreddare e piu' massa")))
    reg.aggiungi(Scelta("L_star", params.derived["L_star"], "m", Origine.DERIVATA,
                        "volume di camera diviso area di gola, conseguenza di "
                        "contrazione e lunghezza"))
    for nome, val, motivo in (
        ("phi_camera", req.phi, "requisito: leggermente magro, quasi al massimo di "
         "T_ad senza lasciare incombusto che esce dall'ugello"),
        ("f_film", f_film, "frazione di combustibile deviata al velo che protegge "
         "il corpo centrale. NON e' film in testa - quello proteggeva una camera "
         "che non ne ha bisogno - e' una fessura alla BASE DEL PLUG, l'unica parte "
         "del motore che nessun circuito puo' raggiungere"),
        ("semiangolo_convergente", CONV_HALF_ANGLE, "37 gradi: sopra i 45 il flusso "
         "si stacca, sotto i 25 il convergente si allunga e porta massa e superficie "
         "da raffreddare"),
        ("troncamento_plug", PLUG_TRUNC, "80 % del plug ideale: la coda contribuisce "
         "pochissima spinta e molta massa e superficie. E' il compromesso classico "
         "dell'aerospike"),
    ):
        reg.aggiungi(Scelta(nome, val, "-" if nome != "semiangolo_convergente" else "rad",
                            Origine.DECISA, motivo))
    return params, l0


CONV_HALF_ANGLE = 0.65      # rad
PLUG_TRUNC = 0.80
LC_SU_DC = 1.4


# --------------------------------------------------------------------------- #
# 3. raffreddamento: dai flussi termici ai canali, non viceversa
# --------------------------------------------------------------------------- #
#: Velocita' massima ammessa nell'acqua [m/s]. Sopra, l'erosione e il rumore
#: diventano un problema e la perdita di carico esplode (va come v^2).
V_ACQUA_MAX = 18.0
#: Lati di canale che si prendono in considerazione, dal piu' piccolo. Il
#: piu' piccolo che soddisfa il bilancio di carico vince: un canale piccolo
#: da' h alto a parita' di velocita' (h ~ D^-0.2 a v fissata) e lascia piu'
#: materiale fra un canale e l'altro.
LATI_CANDIDATI = (0.6e-3, 0.8e-3, 1.0e-3, 1.2e-3, 1.5e-3, 2.0e-3)
#: Salto di temperatura massimo ammesso lungo il SETTO fra due canali [K].
#: Il setto e' un'aletta: conduce lateralmente e il suo apice sta piu' caldo
#: della base. dT = q L^2 / (2 k t) con L semipasso e t spessore.
DELTA_T_SETTO_MAX = 120.0


@dataclass
class RamoRaffreddamento:
    nome: str
    potenza: float
    q_max: float
    mdot: float
    lato: float
    velocita: float
    n_canali: int
    parete_calda: float
    parete_fredda: float
    dp: float
    h: float
    T_parete_gas: float
    margine_ebollizione: float


def _parete_calda(q_max: float, proc, reg: Registro, nome: str) -> float:
    """Lo spessore della parete calda lo detta la TENSIONE TERMICA, non la
    pressione: si fa la piu' sottile che il processo consente, non la piu'
    spessa che la struttura permetterebbe."""
    #: SI PARTE DAL PIU' SOTTILE, non dal piu' spesso ammesso. Alla prima
    #: stesura questa funzione restituiva il massimo spessore compatibile col
    #: limite termico, e per la camera dava 4 mm: un errore di verso. Una parete
    #: calda spessa non serve a niente e fa danno tre volte - piu' salto di
    #: temperatura, piu' tensione termica, piu' massa - e la pressione non la
    #: chiede (6 MPa di cerchio contro 110 di snervamento a caldo).
    t_termico = DELTA_T_PARETE_MAX * K_316L / max(q_max, 1.0)
    t_processo = 2.0 * proc.parete_minima
    if t_termico < proc.parete_minima:
        raise MissingDatum(
            f"il ramo '{nome}' vede {q_max/1e6:.2f} MW/m2: per restare entro "
            f"{DELTA_T_PARETE_MAX:.0f} K di salto la parete dovrebbe essere "
            f"{t_termico*1e3:.2f} mm, sotto il minimo che il processo stampa a "
            f"tenuta ({proc.parete_minima*1e3:.2f} mm). Non c'e' parete possibile.")
    t = min(t_processo, t_termico)
    if t < t_processo:
        motivo = (f"il limite e' TERMICO: q = {q_max/1e6:.2f} MW/m2 con k = {K_316L} "
                  f"W/mK da' {DELTA_T_PARETE_MAX:.0f} K di salto a {t*1e3:.2f} mm. "
                  "Piu' spessa plasticizza a ogni accensione")
    else:
        motivo = (f"il doppio della parete minima di processo. Il limite termico ne "
                  f"consentirebbe {t_termico*1e3:.2f} mm, ma piu' spessa non serve: "
                  "la pressione chiede 6 MPa contro 110 di snervamento a caldo")
    reg.aggiungi(Scelta(f"parete_calda_{nome}", t, "m", Origine.DERIVATA, motivo))
    return t


def _dimensiona_ramo(nome: str, potenza: float, q_max: float,
                     lunghezza: float, raggio: float, req: Requisiti,
                     reg: Registro, T_acqua: float, p_acqua: float,
                     mdot_imposta: float | None = None) -> RamoRaffreddamento:
    """Dal flusso termico ai canali. L'ORDINE DEI PASSI E' IL PUNTO.

    Alla prima stesura avevo derivato la portata dal salto di temperatura
    dell'acqua e da li' il numero di canali. E' sbagliato, e il modello se n'e'
    accorto rifiutando il progetto: **in gola la portata non la decide il salto
    entalpico, la decide il coefficiente di scambio**. Servono 51 kW/m2K, la
    velocita' che li produce e' quella, l'area di passaggio segue, e il salto di
    temperatura che ne risulta e' di pochi gradi. Derivare la portata dai 30 K
    voluti dava un decimo dell'acqua necessaria.

    L'ordine giusto:
      1. parete calda dalla TENSIONE TERMICA (la piu' sottile che il processo fa);
      2. h necessario perche' la parete lato acqua resti sotto saturazione;
      3. numero di canali dal SETTO, che e' un'aletta: il suo apice sta piu'
         caldo della base di q L^2 / (2 k t), e con pochi canali larghi quel
         salto vale centinaia di gradi;
      4. velocita' da h, area da velocita', portata come RISULTATO;
      5. perdita di carico e salto di temperatura come verifiche, non come dati.

    CON `mdot_imposta` L'ORDINE SI ROVESCIA, ed e' voluto. I due rami sono in
    SERIE: la portata non e' piu' un risultato di questo ramo, e' un dato che
    arriva dall'altro. Restano liberi il lato e il numero di canali, e la
    velocita' diventa il risultato: v = mdot / (rho n lato^2). Si sceglie allora
    la combinazione che raggiunge h_min con la caduta piu' bassa, perche' in
    serie le cadute si SOMMANO e non si dividono.
    """
    proc = req.processo
    parete_calda = _parete_calda(q_max, proc, reg, nome)
    parete_fredda = max(2.0 * proc.parete_minima, proc.parete_minima)
    reg.aggiungi(Scelta(
        f"parete_fredda_{nome}", parete_fredda, "m", Origine.DERIVATA,
        "il doppio della parete minima: separa acqua da atmosfera, la pressione "
        "e' irrisoria e comanda la tenuta del processo"))

    h_min = cooling.required_h_for_no_boiling(
        q_max, T_acqua, p_acqua, safety_margin_K=MARGINE_EBOLLIZIONE)
    rho = cooling.PropsSI("D", "T", T_acqua, "P", p_acqua, "Water")
    cp = cooling.PropsSI("C", "T", T_acqua, "P", p_acqua, "Water")
    circonferenza = 2.0 * math.pi * raggio

    def numero_di_canali(lato: float):
        """Il MINIMO numero di canali che tiene l'apice del setto entro il salto
        ammesso. Il setto fra due canali e' un'aletta: meno canali significa un
        setto piu' largo, e il salto va come il QUADRATO della larghezza."""
        n_max = int(circonferenza / (lato + proc.parete_minima))
        for n in range(4, max(n_max, 5)):
            setto = circonferenza / n - lato
            if setto < 2.0 * proc.parete_minima:
                return None
            dT = q_max * (0.5 * setto) ** 2 / (2.0 * K_316L * parete_calda)
            if dT <= DELTA_T_SETTO_MAX:
                return n, setto, dT
        return None

    scelto = None
    scartati: list[str] = []
    candidati: list[tuple] = []
    for lato in LATI_CANDIDATI:
        if lato < proc.margine_foro(MARGINE_FORO):
            continue
        n_ok = numero_di_canali(lato)
        if n_ok is None:
            scartati.append(f"{lato*1e3:.1f} mm: nessun numero di canali tiene il "
                            f"setto sotto {DELTA_T_SETTO_MAX:.0f} K con "
                            f"{parete_calda*1e3:.2f} mm di parete calda")
            continue
        n, setto, dT_setto = n_ok

        if mdot_imposta is None:
            def scarto(v: float) -> float:
                return cooling.dittus_boelter_h(lato, v, T_acqua, p_acqua) - h_min

            if scarto(V_ACQUA_MAX) < 0.0:
                scartati.append(f"{lato*1e3:.1f} mm: nemmeno a {V_ACQUA_MAX:.0f} m/s "
                                f"si arriva a {h_min/1e3:.0f} kW/m2K")
                continue
            v = brentq(scarto, 0.05, V_ACQUA_MAX, xtol=1.0e-9)
            mdot = rho * v * n * lato ** 2
        else:
            #: In serie la portata e' un dato: la velocita' e' cio' che ne esce.
            #: Con n scelto dal setto, l'unico modo di alzare h e' STRINGERE il
            #: canale, e l'unico modo di abbassare la caduta e' allargarlo. Il
            #: compromesso lo decide il ciclo, non una regola scritta.
            mdot = mdot_imposta
            v = mdot / (rho * n * lato ** 2)
            if v > V_ACQUA_MAX:
                scartati.append(f"{lato*1e3:.1f} mm x {n}: con la portata della "
                                f"serie ({mdot:.4f} kg/s) la velocita' sarebbe "
                                f"{v:.1f} m/s, sopra {V_ACQUA_MAX:.0f}")
                continue
            if cooling.dittus_boelter_h(lato, v, T_acqua, p_acqua) < h_min:
                scartati.append(f"{lato*1e3:.1f} mm x {n}: a {v:.1f} m/s da' "
                                f"{cooling.dittus_boelter_h(lato, v, T_acqua, p_acqua)/1e3:.1f} "
                                f"kW/m2K contro i {h_min/1e3:.1f} richiesti")
                continue

        h = cooling.dittus_boelter_h(lato, v, T_acqua, p_acqua)
        dp = cooling.channel_pressure_drop(lato, lunghezza, v, T_acqua, p_acqua)
        dT_acqua = potenza / (mdot * cp)
        if dp > DP_ACQUA_MAX:
            scartati.append(f"{lato*1e3:.1f} mm: caduta {dp/1e5:.2f} bar sopra il "
                            f"limite di {DP_ACQUA_MAX/1e5:.1f}")
            continue
        candidati.append((dp, lato, v, n, h, setto, dT_setto, mdot, dT_acqua))
        if mdot_imposta is None:
            break

    if candidati:
        #: A PORTATA IMPOSTA SI SCEGLIE IL CANALE PIU' PICCOLO, non quello che
        #: cade di meno. E' l'opposto di quello che avevo scritto per primo, e
        #: il motivo e' che la caduta non e' l'unica cosa che si paga.
        #:
        #: Un canale largo cade poco, e sembra gratis. Ma la caduta di un canale
        #: e' anche il RIFERIMENTO rispetto a cui si giudica la caduta lungo
        #: l'anello del collettore che lo alimenta: se il canale cade 0.02 bar,
        #: l'anello deve cadere meno di 0.002, e per riuscirci deve avere una
        #: sezione enorme - cioe' un collare esterno che pesa piu' di cento
        #: grammi. Con un canale piu' stretto la stessa uniformita' di
        #: distribuzione si ottiene con un anello molto piu' piccolo.
        #:
        #: Si prende quindi il piu' piccolo che sta nel budget di caduta, non il
        #: piu' comodo. Il budget e' quello che c'e': in serie le cadute si
        #: sommano.
        ammessi = [c for c in candidati if c[0] <= FRAZIONE_DP_RAMO * DP_ACQUA_MAX]
        scelta_c = min(ammessi or candidati, key=lambda c: (c[1], c[0]))
        dp, lato, v, n, h, setto, dT_setto, mdot, dT_acqua = scelta_c
        scelto = (lato, v, n, h, dp, setto, dT_setto, mdot, dT_acqua)

    if scelto is None:
        raise MissingDatum(
            f"nessun canale raffredda il ramo '{nome}' (q = {q_max/1e6:.2f} MW/m2, "
            f"serve h = {h_min/1e3:.0f} kW/m2K). Scartati: " + "; ".join(scartati))
    lato, v, n, h, dp, setto, dT_setto, mdot, dT_acqua = scelto

    reg.aggiungi(Scelta(
        f"lato_canale_{nome}", lato, "m", Origine.DERIVATA,
        f"il piu' piccolo lato stampabile che raggiunge h = {h/1e3:.1f} kW/m2K "
        f"(ne servono {h_min/1e3:.1f}) restando sotto {DP_ACQUA_MAX/1e5:.1f} bar. "
        + ("Scartati: " + "; ".join(scartati) if scartati else "")))
    reg.aggiungi(Scelta(
        f"n_canali_{nome}", float(n), "-", Origine.DERIVATA,
        f"il minimo che tiene l'apice del setto entro {DELTA_T_SETTO_MAX:.0f} K: "
        f"setto {setto*1e3:.2f} mm, salto {dT_setto:.0f} K. Con meta' dei canali il "
        f"salto quadruplicherebbe"))
    if mdot_imposta is None:
        reg.aggiungi(Scelta(
            f"velocita_acqua_{nome}", v, "m/s", Origine.DERIVATA,
            f"la velocita' che porta Dittus-Boelter a {h_min/1e3:.0f} kW/m2K su un "
            f"canale da {lato*1e3:.1f} mm"))
        reg.aggiungi(Scelta(
            f"portata_acqua_{nome}", mdot, "kg/s", Origine.DERIVATA,
            f"RISULTATO di velocita' e sezione, non dato: il salto di temperatura che "
            f"ne esce e' {dT_acqua:.1f} K, cioe' l'acqua esce quasi fredda. In gola "
            "comanda lo scambio, non l'entalpia"))
    else:
        reg.aggiungi(Scelta(
            f"velocita_acqua_{nome}", v, "m/s", Origine.DERIVATA,
            f"RISULTATO, non scelta: in serie la portata e' quella del ramo che la "
            f"detta ({mdot:.4f} kg/s), e {n} canali da {lato*1e3:.1f} mm la fanno "
            f"passare a questa velocita'. Ne esce h = {h/1e3:.1f} kW/m2K contro i "
            f"{h_min/1e3:.1f} richiesti"))

    w = cooling.steady_wall_temperatures(q_max, parete_calda, K_316L, h,
                                         T_acqua, p_acqua)
    return RamoRaffreddamento(
        nome=nome, potenza=potenza, q_max=q_max, mdot=mdot, lato=lato, velocita=v,
        n_canali=n, parete_calda=parete_calda, parete_fredda=parete_fredda,
        dp=dp, h=h, T_parete_gas=w.T_gas_side, margine_ebollizione=w.margin_to_boiling)


def _raffreddamento(req: Requisiti, reg: Registro, op, params, l0, p_c: float):
    """Due tratti, camera e gola, IN SERIE su una sola acqua.

    PERCHE' IN SERIE, dopo che la prima versione li metteva in parallelo. In
    parallelo servivano due ingressi e due uscite (l'impianto ne ha uno solo) e
    soprattutto una strozzatura di bilanciamento: due rami in parallelo hanno
    per definizione la stessa caduta, quindi con cadute di progetto diverse
    l'acqua se ne va quasi tutta nel ramo che oppone meno resistenza, cioe' via
    dalla gola. Quel foro calibrato e' un pezzo che DEVE essere giusto, e se
    sbagliato non si vede: portata e temperatura all'uscita restano normali
    mentre il labbro non e' piu' raffreddato.

    In serie non c'e' niente da bilanciare: la stessa acqua passa da tutto, per
    costruzione. Il prezzo sarebbe la caduta, e il conto che aveva bocciato la
    serie alla prima stesura era fatto mandando la portata piena nei canali
    della camera DIMENSIONATI PER LA PORTATA RIDOTTA del parallelo. Ridimensionati
    per la portata vera, la camera cade pochi millesimi di bar e in piu' esce
    dal regime di transizione (Re ~ 1150) in cui Dittus-Boelter e' meno
    affidabile.

    ORDINE: gola per prima. Non per ragioni termiche - l'acqua si scalda in
    tutto di pochi gradi, quindi l'ordine e' termicamente indifferente - ma
    perche' e' l'ordine in cui nessun condotto deve scavalcare il collettore di
    un altro tratto."""
    d = params.derived
    righe = mappa_di_calore(op, params, l0, p_c)
    pot = potenza_per_zona(righe, params)
    qmax = q_massimo_per_zona(righe)
    Q_tot = sum(pot.values())
    reg.aggiungi(Scelta("potenza_termica", Q_tot, "W", Origine.DERIVATA,
                        "integrale di q sulla superficie bagnata: "
                        + ", ".join(f"{k} {v:.0f} W" for k, v in sorted(pot.items())),
                        fonte="carico.mappa_di_calore (Bartz)"))

    T_acqua, p_acqua = 293.15, 4.0e5

    #: Il ramo di gola prende il convergente E il plug: sono la stessa zona
    #: termica, quella dove q supera i MW/m2 su pochi millimetri di corsa.
    Q_camera = pot.get("camera", 0.0)
    Q_gola = Q_tot - Q_camera
    q_camera = qmax.get("camera", 1.0)
    q_gola = max(qmax.get("convergente", 0.0), qmax.get("plug", 0.0))

    L_c, L_conv = d["L_c"], d["L_conv"]
    x_lip = L_c + L_conv
    #: Quote assiali dei due rami, derivate dalla geometria e non scritte:
    #: erano scritte, ed erano giuste per un labbro che stava 2.4 mm piu' in
    #: qua. Il circuito di gola sarebbe finito PRIMA della gola.
    MARGINE_LABBRO = 2.0e-3       # >= meta' larghezza del collettore di inversione
    ANTICIPO_GOLA = 1.4e-3
    SEPARAZIONE = 4.0e-3
    CAMERA_INIZIO = 1.5e-3
    x_gola_0 = L_c - ANTICIPO_GOLA
    x_gola_1 = x_lip - 0.5e-3
    x_cam_1 = x_gola_0 - SEPARAZIONE

    #: LA GOLA DETTA LA PORTATA. E' il tratto in cui la portata non e' libera:
    #: serve h = 51 kW/m2K perche' la parete lato acqua non arrivi a saturazione,
    #: e la velocita' che lo produce fissa l'area e quindi la portata. La camera
    #: la subisce.
    ramo_gola = _dimensiona_ramo("gola", Q_gola, q_gola, abs(x_gola_1 - x_gola_0),
                                 0.5 * (d["R_c"] + d["R_lip"]), req, reg,
                                 T_acqua, p_acqua)
    ramo_camera = _dimensiona_ramo("camera", Q_camera, q_camera,
                                   abs(x_cam_1 - CAMERA_INIZIO), d["R_c"], req, reg,
                                   T_acqua, p_acqua, mdot_imposta=ramo_gola.mdot)
    rami = [ramo_camera, ramo_gola]
    reg.aggiungi(Scelta(
        "portata_acqua_totale", ramo_gola.mdot, "kg/s", Origine.DERIVATA,
        "la portata del ramo di gola, che in serie e' anche quella della camera e "
        "quella della canna. E' un requisito verso l'impianto, non un parametro del "
        "motore: se la rete non la garantisce, non c'e' geometria che rimedi"))
    cp_w = cooling.PropsSI("C", "T", T_acqua, "P", p_acqua, "Water")
    dT_serie = Q_tot / (ramo_gola.mdot * cp_w)
    reg.aggiungi(Scelta(
        "salto_temperatura_acqua", dT_serie, "K", Origine.DERIVATA,
        f"tutta la potenza termica ({Q_tot:.0f} W) su tutta la portata: e' il "
        "riscaldamento dell'acqua dall'ingresso all'uscita. E' anche la misura di "
        "quanto conti l'ordine dei due tratti, cioe' pochissimo"))

    circuiti = []
    for r, (x0, x1) in zip(rami, ((CAMERA_INIZIO, x_cam_1), (x_gola_0, x_gola_1))):
        ritorno = r.nome == "gola"
        setto = 2.0 * req.processo.parete_minima if ritorno else 0.6e-3
        if ritorno:
            reg.aggiungi(Scelta(
                "setto_ritorno_gola", setto, "m", Origine.DERIVATA,
                "il doppio della parete minima. E' l'UNICA parete del motore la cui "
                "rottura e' silenziosa: l'acqua non esce, cortocircuita il labbro e "
                "smette di raffreddare il punto piu' caldo mentre portata e "
                "temperatura all'uscita restano normali"))
        #: Passo dell'elica: l'autosostentamento vuole passo >= 2 pi R tan(alpha).
        #: Deriva dal raggio del canale e dall'angolo di progetto, non si sceglie.
        prof = r.parete_calda + r.lato + (setto + r.lato if ritorno else 0.0)
        raggio_canale = max(_raggio_parete_locale(d, x0),
                            _raggio_parete_locale(d, x1)) + prof
        passo = 1.15 * passo_elica_minimo(raggio_canale, req.processo.angolo_progetto)
        reg.aggiungi(Scelta(
            f"passo_elica_{r.nome}", passo, "m", Origine.DERIVATA,
            f"2 pi R tan({req.processo.angolo_progetto:.0f} gradi) a R = "
            f"{raggio_canale*1e3:.1f} mm, piu' il 15 % di margine: sotto, il tetto "
            "del canale e' uno sbalzo, e dentro un canale nessun supporto e' rimovibile"))
        circuiti.append(CircuitoRaffreddamento(
            nome=r.nome, n_canali=r.n_canali, lato=r.lato,
            parete_calda=r.parete_calda, parete_fredda=r.parete_fredda,
            x_inizio=x0, x_fine=x1, passo_elica=passo, velocita=r.velocita,
            ritorno=ritorno, setto_ritorno=setto, sfalsamento_ritorno=3.0e-3))
    return circuiti, rami, pot, (T_acqua, p_acqua)


#: Velocita' massima nei condotti del collettore [m/s]. Non e' il limite dei
#: canali (18 m/s: li' la velocita' SERVE, compra scambio termico). Qui la
#: velocita' non compra niente e costa come il quadrato, quindi il limite si
#: DERIVA dal budget di caduta invece di sceglierlo: si ammette che le perdite
#: concentrate del collettore valgano al piu' `FRAZIONE_DP_COLLETTORE` del
#: budget, e da li' esce la velocita' e quindi il diametro.
FRAZIONE_DP_COLLETTORE = 0.35

#: Perdite concentrate, in pressioni dinamiche riferite a ciascun tratto.
K_PORTAGOMMA = 1.0     # imbocco dalla canna nel foro
K_FORI = 1.5           # ingresso in un fascio di fori + sbocco
K_ANELLO = 1.0         # allargamento nell'anello e restringimento all'uscita

#: Diametri di foro ammessi [m]. Non e' un catalogo: sono i passi con cui ha
#: senso quotare un foro stampato.
DIAMETRI_FORO = tuple(x * 1.0e-3 for x in (2.0, 2.5, 3.0, 3.2, 4.0, 5.0, 6.0,
                                           7.0, 8.0, 9.0, 10.0, 12.0))

#: Diametro massimo di un foro ORIZZONTALE stampato senza supporti [m]. Sopra,
#: la calotta superiore e' uno sbalzo che cede. E' la stessa soglia che ha
#: prodotto gli attacchi da 3.2 mm invece di uno da 5.6.
D_FORO_ORIZZONTALE_MAX = 3.2e-3


#: Fattore fra la spinta e il carico di progetto dell'interfaccia. NON e' un
#: coefficiente di sicurezza sulla spinta: la spinta e' 50 N e quattro bulloni
#: qualsiasi la reggono mille volte. E' il riconoscimento che il carico vero non
#: e' la spinta, sono le mani, i tubi che tirano, il banco che vibra e il colpo
#: d'ariete all'apertura delle valvole - e nessuno di questi e' calcolato.
#: Dimensionare sulla sola spinta darebbe M1.6: si dichiara il fattore invece di
#: fingere che 50 N siano il caso peggiore.
FATTORE_CARICO_MONTAGGIO = 20.0

#: Filettature metriche a passo grosso: (nominale, nocciolo da maschiare).
#: I due valori sono normativi (ISO 261/ISO 965), non scelti.
METRICHE = ((3.0e-3, 2.5e-3), (4.0e-3, 3.3e-3), (5.0e-3, 4.2e-3),
            (6.0e-3, 5.0e-3), (8.0e-3, 6.8e-3))

#: Tensione ammissibile su un bullone di classe 8.8 [Pa], con il coefficiente di
#: sicurezza gia' dentro: 640 MPa di snervamento diviso 4.
SIGMA_AMM_BULLONE = 160.0e6


#: Inclinazione della fessura di film sull'asse del motore. Piu' e' piccola,
#: piu' il film esce tangente alla parete che deve proteggere - che e' quello
#: che serve. Il limite e' la stampa dalla parte opposta: la fessura e' un
#: condotto, e un condotto inclinato di `angolo` sull'asse di costruzione e'
#: autosostentato solo se `angolo` sta entro (90 - angolo di progetto). Si sta
#: quindi il piu' basso possibile, cioe' proprio a quel limite.
def _fessura_film(req: Requisiti, reg: Registro, d, params, testa, l0,
                  film, p_c: float) -> FessuraFilm:
    """La fessura del film e il condotto che la alimenta.

    IL DIAMETRO DEL FORO DI DOSATURA NON SI SCEGLIE: e' il piu' piccolo che il
    processo sa fare. E' la frazione di film a seguirlo, non il contrario. La
    ragione sta in `alimentabilita_film`: la fessura vede lo stesso salto di
    pressione dei getti principali, quindi se fosse lei a dosare si porterebbe
    via una frazione confrontabile della portata. Chi dosa e' un foro, e sotto
    il minimo di processo i fori non esistono.
    """
    proc = req.processo
    #: `plug_contour_r[0]` e' il raggio del plug ALLA SUA BASE, che e' dove la
    #: fessura deve stare: e' la stessa quota che `dimensiona_film` usa per
    #: l'area di efflusso, e le due devono essere la stessa cosa o il modello
    #: del film descrive una fessura diversa da quella costruita.
    r_slot = params.plug_contour_r[0]
    x_lip = d["L_c"] + d["L_conv"]
    x_uscita = x_lip + params.plug_contour_x[0]
    angolo = math.radians(90.0 - proc.angolo_progetto)
    d_dos = proc.margine_foro(MARGINE_FORO)
    #: il condotto: il piu' grande che lascia parete dentro il corpo centrale
    #: nel suo punto piu' STRETTO, che e' in testa (R_getti), non alla base del
    #: plug. Un condotto dimensionato sul punto largo buca il punto stretto.
    r_condotto = max(0.5 * d_dos,
                     min(0.5 * r_slot, testa.R_getti - 2.0 * proc.parete_minima
                         - 0.5 * d_dos))
    spessore_labbro = 2.0 * proc.parete_minima
    #: IL FORO DI DOSATURA COMINCIA DENTRO IL PLENUM DEL GPL, non a x = 0. Il
    #: plenum finisce poco a valle dei getti principali: un foro che comincia
    #: alla faccia d'iniezione non lo tocca, e tutto il condotto del film resta
    #: un vuoto cieco - 217 mm3 di polvere che non esce e nessun film.
    x_dos = testa.x_getti
    x_cond0 = x_dos + 3.0 * d_dos

    #: verifica che la fessura ci stia: deve scendere dal labbro al condotto
    #: restando dentro il tratto cilindrico del corpo centrale
    h_r = film.altezza_fessura / math.cos(angolo)
    L = (r_slot - h_r - r_condotto) / math.sin(angolo)
    x0 = x_uscita - L * math.cos(angolo)
    avvertenze: list[str] = []
    if x0 <= d["L_c"]:
        avvertenze.append(
            f"la fessura comincia a x = {x0*1e3:.1f} mm, cioe' dentro la camera "
            f"(che finisce a {d['L_c']*1e3:.1f}): li' il corpo centrale non e' piu' "
            "cilindrico e il condotto esce dalla parete")
    if r_condotto <= 0.5 * d_dos:
        avvertenze.append(
            "il condotto del film e' stretto quanto il foro di dosatura: non e' "
            "piu' il foro a decidere la portata ma l'attrito nel condotto, che "
            "dipende dalla rugosita' di stampa e non e' un dato che abbiamo")

    for nome, val, unita, motivo in (
        ("d_dosatura_film", d_dos, "m",
         "il piu' piccolo foro che il processo sa fare. NON e' una scelta: e' il "
         "vincolo, e la frazione di film e' la sua conseguenza. Se la macchina "
         "facesse fori piu' piccoli, il film sarebbe piu' magro"),
        ("r_condotto_film", r_condotto, "m",
         f"il piu' grande che lascia due pareti minime dentro il corpo centrale nel "
         f"suo punto piu' stretto, che e' in testa ({testa.R_getti*1e3:.2f} mm di "
         "raggio) e non alla base del plug"),
        ("angolo_fessura_film", math.degrees(angolo), "gradi",
         f"90 meno l'angolo di progetto ({proc.angolo_progetto:.0f}): la fessura e' "
         "un condotto, e sopra questa inclinazione ha un soffitto. Piu' bassa "
         "sarebbe meglio per il film, che vuole uscire tangente"),
        ("spessore_labbro_film", spessore_labbro, "m",
         f"il labbro raggiunge questo spessore gia' "
         f"{spessore_labbro/math.tan(angolo)*1e3:.2f} mm prima del bordo, perche' la "
         "fessura e' inclinata. Gli ultimi decimi usciranno arrotondati dalla "
         "stampa, e per un labbro di fessura va bene"),
    ):
        reg.aggiungi(Scelta(nome, val, unita, Origine.DERIVATA, motivo))

    return FessuraFilm(
        x_uscita=x_uscita, r_uscita=r_slot, altezza=film.altezza_fessura,
        angolo=angolo, r_condotto=r_condotto, x_condotto_0=x_cond0,
        d_dosatura=d_dos, x_dosatura=x_dos, spessore_labbro=spessore_labbro,
        avvertenze=tuple(avvertenze))


def _montaggio(req: Requisiti, reg: Registro, d, testa, spinta: float,
               collettore) -> InterfacciaMontaggio:
    """L'interfaccia di montaggio, sulla faccia di monte.

    Il numero di bracci non e' quattro per abitudine: gli assi a 0 e 180 gradi
    sono gia' occupati (il bocchettone dell'aria e le canne dell'acqua) e ne
    servono almeno tre per definire un piano. Quattro a 45, 135, 225 e 315 gradi
    danno un vincolo simmetrico e nessuna interferenza.

    IL RAGGIO DEL CERCHIO DEI BULLONI NON LO DECIDE LA RESISTENZA. Lo decide la
    CHIAVE: i fori sono passanti e il dado sta dietro la flangia, quindi dietro
    ci deve essere spazio libero per una chiave da `ingombro_chiave`. Il corpo
    della testa e' un cilindro da 10.5 mm di raggio, e il dado deve stargli
    fuori.
    """
    proc = req.processo
    n = 4
    carico = FATTORE_CARICO_MONTAGGIO * spinta / n
    area_serve = carico / SIGMA_AMM_BULLONE
    d_min_resistenza = math.sqrt(4.0 * area_serve / math.pi)
    scelta = next(((dn, dc) for dn, dc in METRICHE
                   if math.pi / 4.0 * dc ** 2 >= area_serve), None)
    if scelta is None:
        raise MissingDatum("nessuna metrica dell'elenco regge il carico di montaggio")
    #: fra le metriche che REGGONO si prende comunque la M5, e non e' una
    #: ragione strutturale: sotto M5 la chiave e' minuscola, la rondella non c'e'
    #: e un bullone da banco piu' piccolo non si trova in cassetta. E' una
    #: scelta di officina, e va dichiarata come tale invece di farla passare per
    #: un calcolo.
    d_nom = max(scelta[0], 5.0e-3)
    d_foro = d_nom + 0.5e-3                     # foro di passaggio, ISO 273 media
    #: apertura di chiave di un dado esagonale metrico (ISO 4032): 8 mm per M5,
    #: 7 per M4, 10 per M6. Il cerchio circoscritto e' 2/sqrt(3) volte tanto.
    chiave = {3.0e-3: 5.5e-3, 4.0e-3: 7.0e-3, 5.0e-3: 8.0e-3,
              6.0e-3: 10.0e-3, 8.0e-3: 13.0e-3}[d_nom]
    ingombro = chiave * 2.0 / math.sqrt(3.0)

    #: RAGGIO. Tre vincoli, e vince il piu' esterno:
    #:   * fuori dalla sede filettata del GPL, che sta al centro della stessa
    #:     faccia;
    #:   * il DADO deve stare fuori dal corpo della testa, altrimenti la chiave
    #:     non entra;
    #:   * dentro il bocchettone dell'aria, cosi' una chiave arriva anche a
    #:     quello.
    R_gpl = (0.5 * testa.d_filetto_gpl + 2.0 * proc.parete_minima
             + 0.5 * d_foro + 2.0 * proc.parete_minima)
    R_dado = testa.R_testa + 0.5 * ingombro + 2.0 * proc.parete_minima
    R_b = max(R_gpl, R_dado)
    R_max = testa.R_boss_aria + 0.5 * ingombro
    d_pad = d_foro + 4.0 * proc.parete_minima
    larghezza = d_pad
    spessore = max(4.0 * proc.parete_minima, 3.0e-3)

    avvertenze: list[str] = []
    if R_b > R_max:
        avvertenze.append(
            f"il cerchio dei bulloni sta a {R_b*1e3:.1f} mm contro un massimo "
            f"consigliato di {R_max*1e3:.1f}: le piazzole sporgono oltre il "
            "bocchettone dell'aria e possono ostacolare la chiave del raccordo")
    avvertenze.append(
        "la piastra di banco deve avere un foro centrale: il raccordo del GPL e' "
        "assiale al centro della stessa faccia di montaggio")
    avvertenze.append(
        f"i fori sono PASSANTI da {d_foro*1e3:.1f} mm: la filettatura sta nel dado, "
        f"dietro la flangia, e non nel pezzo stampato. Serve {ingombro*1e3:.0f} mm "
        "di spazio libero dietro ogni piazzola")

    for nome, val, unita, motivo in (
        ("n_bulloni_montaggio", float(n), "-",
         "quattro: gli assi a 0 e 180 gradi sono occupati dal bocchettone "
         "dell'aria e dalle canne dell'acqua, e ne servono almeno tre per "
         "definire un piano"),
        ("d_bulloni_montaggio", d_nom, "m",
         f"M{d_nom*1e3:.0f}. Il carico di progetto e' {carico:.0f} N per bullone "
         f"({FATTORE_CARICO_MONTAGGIO:.0f} volte la spinta divisa per {n}), che "
         f"chiederebbe {d_min_resistenza*1e3:.2f} mm di nocciolo: la resistenza NON "
         "e' il criterio, lo e' la maneggiabilita' in officina"),
        ("R_bulloni_montaggio", R_b, "m",
         f"il maggiore fra {R_gpl*1e3:.1f} mm (fuori dalla sede del GPL) e "
         f"{R_dado*1e3:.1f} (il dado deve stare fuori dal corpo della testa, che ha "
         f"raggio {testa.R_testa*1e3:.1f}): lo decide la chiave, non il carico"),
        ("spessore_flangia_montaggio", spessore, "m",
         "quattro pareti minime. Sta tutto nel piano della piastra di stampa, "
         "quindi non ha sbalzi: una flangia anulare continua peserebbe 35 g, "
         "quattro bracci ne pesano una decina e reggono lo stesso carico assiale"),
    ):
        reg.aggiungi(Scelta(nome, val, unita, Origine.DERIVATA, motivo))
    reg.aggiungi(Scelta(
        "apertura_chiave_montaggio", chiave, "m", Origine.NORMATIVA,
        f"ISO 4032: dado esagonale M{d_nom*1e3:.0f}"))

    return InterfacciaMontaggio(
        n_bracci=n, angolo_offset=0.25 * math.pi, R_bulloni=R_b, d_foro=d_foro,
        d_piazzola=d_pad, larghezza_braccio=larghezza, spessore=spessore,
        r_attacco=0.5 * testa.R_testa, d_filetto=d_nom, ingombro_chiave=ingombro,
        carico_per_bullone=carico, avvertenze=tuple(avvertenze))


def _collettore_acqua(req: Requisiti, reg: Registro, d, circuiti, rami, testa,
                      T_acqua: float, p_acqua: float) -> CollettoreAcqua:
    """Il collettore a un ingresso e una uscita, dimensionato e non disegnato.

    Ogni quota esce da un bilancio: i diametri dalla velocita' che il budget di
    caduta consente, il numero di fori dalla velocita' dentro il collettore
    anulare che alimentano, l'altezza dei vani anulari dalla portata che ci
    passa, le pendenze dall'autosostentamento.
    """
    rho = cooling.PropsSI("D", "T", T_acqua, "P", p_acqua, "Water")
    mdot = max(r.mdot for r in rami)
    Q_v = mdot / rho
    proc = req.processo
    parete = 2.0 * proc.parete_minima
    pendenza = math.tan(math.radians(proc.angolo_progetto))

    cam = next(c for c in circuiti if c.nome == "camera")
    gola = next(c for c in circuiti if c.nome == "gola")
    t_max = max(engine.spessore_mantello(c) for c in circuiti)
    x_lip = d["L_c"] + d["L_conv"]

    K_tot = 2.0 * K_PORTAGOMMA + 2.0 * K_FORI + 2.0 * K_ANELLO
    v_max = math.sqrt(2.0 * FRAZIONE_DP_COLLETTORE * DP_ACQUA_MAX / (rho * K_tot))
    reg.aggiungi(Scelta(
        "velocita_condotti_acqua", v_max, "m/s", Origine.DERIVATA,
        f"non scelta: e' la velocita' per cui le perdite concentrate del collettore "
        f"(K totale {K_tot:.1f}) valgono il {FRAZIONE_DP_COLLETTORE*100:.0f} % del "
        f"budget di {DP_ACQUA_MAX/1e5:.1f} bar. Nei canali il limite e' "
        f"{V_ACQUA_MAX:.0f} m/s perche' li' la velocita' compra scambio termico; qui "
        "non compra niente e costa come il quadrato"))
    area_serve = Q_v / v_max

    def fori_per(area, d_max, R_corona, nome, n_max=24):
        """Il minimo numero di fori che da' l'area voluta, col diametro piu'
        grande ammesso. Fra due fori adiacenti deve restare una parete: un
        fascio di fori che si compenetrano non e' un fascio di fori, e' una
        fessura anulare, e la parete che doveva reggere non c'e' piu'."""
        for dd in reversed([x for x in DIAMETRI_FORO if x <= d_max + 1e-12]):
            n = max(3, math.ceil(area / (math.pi / 4.0 * dd ** 2)))
            if n > n_max:
                continue
            if 2.0 * math.pi * R_corona / n - dd < proc.parete_minima:
                continue
            return n, dd
        raise MissingDatum(
            f"nessuna combinazione di fori sotto {d_max*1e3:.1f} mm e {n_max} pezzi "
            f"su una corona di raggio {R_corona*1e3:.1f} mm da' i {area*1e6:.1f} mm2 "
            f"che servono ({nome})")

    def giro(L, h, k, R):
        """(area, velocita', caduta) lungo un anello di sezione (L, h).

        La portata cala linearmente con l'angolo mentre l'acqua si distribuisce,
        quindi l'integrale di v^2 lungo il giro vale un TERZO di quello a
        portata piena.
        """
        area = L * h - 0.5 * k * pendenza * h * h
        if area <= 1.0e-9:
            return None
        v = 0.5 * Q_v / area
        Dh = 4.0 * area / (2.0 * (L + h))
        Re = rho * v * Dh / 1.0e-3
        f = (0.316 / Re ** 0.25) if Re > 2300 else (64.0 / max(Re, 1.0))
        return area, v, f * (math.pi * R / Dh) * 0.5 * rho * v ** 2 / 3.0

    def anello(x_v0, x_max, pavimento_inclinato, R_med, dp_ammessa, x_c0_min=None):
        """Sezione e posizione di un vano anulare, dal vincolo di DISTRIBUZIONE.

        NON si dimensiona su una velocita' scelta. Il compito dell'anello e' che
        tutti i canali che alimenta ricevano la stessa acqua, e la misura di
        quanto ci riesce e' il rapporto fra la caduta lungo il giro e la caduta
        di un canale: se l'anello cade quanto un canale, quelli lontani
        dall'attacco non ricevono niente. `dp_ammessa` e' quella soglia.

        LA SEZIONE NON E' UN RETTANGOLO. Il tetto e' rivolto verso la piastra e
        deve arretrare di `pendenza` per unita' di altezza: sezione a trapezio,
        L h - p h^2 / 2. Se ANCHE il pavimento e' inclinato - perche' sotto c'e'
        la rampa con cui il collare nasce dal corpo e non c'e' spazio a monte
        per finirla prima - diventa un rombo, L h - p h^2: meta' area a parita'
        di lunghezza.

        FRA LE COPPIE (L, h) CHE VANNO BENE SI PRENDE LA PIU' BASSA, non la piu'
        corta. La massa del collare va come h per L per il raggio, ma h alza
        anche il raggio esterno di tutto il pezzo: un anello lungo e basso pesa
        meno di uno corto e alto a parita' di sezione.
        """
        k = 2.0 if pavimento_inclinato else 1.0
        L_max = x_max - x_v0
        if L_max <= 1.0e-3:
            return ((x_v0, x_v0 + 1.0e-3, 1.5e-3, x_v0 - parete),
                    1.0e-9, 1.0e-3, (0.0, float("inf")))
        migliore = None
        n_L = 40
        for i in range(1, n_L + 1):
            L = L_max * i / n_L
            h = H_MIN_ANELLO
            while h <= L / (k * pendenza) + 1.0e-12:
                r = giro(L, h, k, R_med)
                if r is not None:
                    area, v, dp_g = r
                    if v <= V_ACQUA_MAX and dp_g <= dp_ammessa:
                        costo = h * (L + pendenza * h)
                        if migliore is None or costo < migliore[0]:
                            migliore = (costo, L, h, area, v, dp_g)
                        break
                h += 1.0e-4
        if migliore is None:
            #: nessuna sezione rispetta il vincolo: si prende la piu' grande che
            #: ci sta e lo si DICHIARA, invece di far finta di niente
            L = L_max
            h = L / (k * pendenza)
            area, v, dp_g = giro(L, h, k, R_med)
            migliore = (h * L, L, h, area, v, dp_g)
        _, L, h, area, v, dp_g = migliore
        x_c0 = x_v0 - parete - (0.0 if pavimento_inclinato else pendenza * h)
        if x_c0_min is not None and x_c0 < x_c0_min:
            x_c0 = x_c0_min
        return (x_v0, x_v0 + L, h, x_c0), area, L, (v, dp_g)

    def sede_portagomma(an, R_par):
        """Dove ci sta l'attacco, e quanto grosso.

        DUE VINCOLI DI SEGNO OPPOSTO, ed e' la loro combinazione a decidere il
        diametro - non la portata.

        In alto: il foro e' a goccia, e l'apice sta `pendenza * raggio` piu' su
        del centro. L'apice deve restare dentro il collare, quindi
        x <= x_v1 - parete - pendenza * d/2.

        In basso: sotto il centro il foro attraversa il tratto di collare che
        sta fra la parete esterna del mantello e il raggio a cui comincia il
        gambo, e li' il collare non e' ancora alto - nasce dal corpo con una
        rampa. Quindi x >= x_c0 + pendenza * (h/2) + parete + d/2.

        I due si incontrano solo se il vano e' abbastanza LUNGO. La prima volta
        avevo imposto solo il primo, e l'apice usciva dal collare; la seconda
        solo il secondo, e usciva il fondo. Scriverli insieme da' anche il
        diametro massimo che ci sta, che e' l'informazione che serviva.
        """
        x_v0, x_v1, h, x_c0 = an
        basso = x_c0 + pendenza * 0.5 * h + parete
        alto = x_v1 - parete
        d_fit = 2.0 * (alto - basso) / (1.0 + pendenza)
        return basso, alto, max(d_fit, 0.0)

    # --- ingresso: anello esterno all'inizio della camera --------------------
    #: IL COLLARE COMINCIA A x = 0. Prima non c'e' corpo a cui attaccarsi: la
    #: testa non e' un disco del diametro del mantello, e' un cilindro da
    #: 10.5 mm di raggio che si allarga a cono solo nell'ultimo tratto. Un
    #: collare piazzato piu' a monte starebbe in aria - ed e' esattamente quello
    #: che faceva il plenum della prima stesura.
    dp_cam0 = next(r.dp for r in rami if r.nome == "camera")
    dp_gola0 = next(r.dp for r in rami if r.nome == "gola")
    an_in, area_in, L_in, (v_an_in, dp_giro_in) = anello(
        parete, cam.x_fine - 4.0e-3, True,
        d["R_c"] + t_max + 1.5e-3, FRAZIONE_MALDISTRIBUZIONE * dp_cam0,
        x_c0_min=0.0)
    prof_cam = cam.parete_calda + 0.5 * cam.lato
    R_par_in = d["R_c"]
    n_in, d_in = fori_per(area_serve, D_FORO_ORIZZONTALE_MAX,
                          R_par_in + t_max, "fori d'ingresso")
    area_col_cam = engine.COLLETTORE_LARGHEZZA * (cam.lato + 2.0 * engine.COLLETTORE_MARGINE)
    n_in = max(n_in, math.ceil(Q_v / (2.0 * area_col_cam * v_max)))
    v_fori_in = Q_v / (n_in * math.pi / 4.0 * d_in ** 2)

    # --- raccordo camera -> gola: vano anulare INTERNO -----------------------
    #: I due collettori stanno a profondita' che si SOVRAPPONGONO: basta un vano
    #: anulare interno al mantello che li unisca su tutta la circonferenza.
    #: Interno e non esterno, e la differenza non e' estetica: un vano dentro
    #: materiale pieno ha il pavimento appoggiato, e solo il tetto e' uno sbalzo.
    #: E costa zero - nessun collare, nessuna massa aggiunta.
    prof_gola = gola.parete_calda + 0.5 * gola.lato
    semi_cam = 0.5 * cam.lato + engine.COLLETTORE_MARGINE
    semi_gola = 0.5 * gola.lato + engine.COLLETTORE_MARGINE
    prof_r0 = min(prof_cam - semi_cam, prof_gola - semi_gola)
    prof_r1 = max(prof_cam + semi_cam, prof_gola + semi_gola)
    x_r0 = cam.x_fine - 0.5 * engine.COLLETTORE_LARGHEZZA
    x_r1 = gola.x_inizio + 0.5 * engine.COLLETTORE_LARGHEZZA
    R_racc = _raggio_parete_locale(d, 0.5 * (x_r0 + x_r1)) + 0.5 * (prof_r0 + prof_r1)
    v_raccordo = Q_v / (2.0 * math.pi * R_racc * (prof_r1 - prof_r0))

    # --- uscita: anello esterno fra il ritorno di gola e il labbro -----------
    x_gu = gola.x_inizio + gola.sfalsamento_ritorno
    prof_rit = gola.parete_calda + gola.lato + gola.setto_ritorno + 0.5 * gola.lato
    R_par_out = _raggio_parete_locale(d, x_gu)
    n_out, d_out = fori_per(area_serve, D_FORO_ORIZZONTALE_MAX,
                            R_par_out + t_max, "fori d'uscita")
    area_col_gola = engine.COLLETTORE_LARGHEZZA * (gola.lato + 2.0 * engine.COLLETTORE_MARGINE)
    n_out = max(n_out, math.ceil(Q_v / (2.0 * area_col_gola * v_max)))
    v_fori_out = Q_v / (n_out * math.pi / 4.0 * d_out ** 2)
    #: il collare puo' cominciare molto prima del vano: a monte c'e' mantello
    #: in abbondanza, quindi la rampa fa in tempo a finire e il pavimento resta
    #: piatto. Il vano comincia mezza larghezza di collettore prima dei fori.
    an_out, area_out, L_out, (v_an_out, dp_giro_out) = anello(
        x_gu - 0.5 * engine.COLLETTORE_LARGHEZZA, x_lip - parete, False,
        R_par_out + t_max + 1.5e-3, FRAZIONE_MALDISTRIBUZIONE * dp_gola0)

    # --- portagomma ----------------------------------------------------------
    #: DIAMETRO DEL FORO DEL PORTAGOMMA. Lo detta il piu' stretto fra tre
    #: vincoli: la velocita' ammessa, il posto che c'e' nell'anello d'ingresso e
    #: quello che c'e' nell'anello d'uscita.
    d_velocita = engine.CANNA_DIAMETRO_INTERNO - 2.0 * parete
    for dd in DIAMETRI_FORO:
        if Q_v / (math.pi / 4.0 * dd ** 2) <= v_max:
            d_velocita = dd
            break
    basso_in, alto_in, d_fit_in = sede_portagomma(an_in, R_par_in)
    basso_out, alto_out, d_fit_out = sede_portagomma(an_out, R_par_out)
    d_barbo = min(d_velocita, engine.CANNA_DIAMETRO_INTERNO - 2.0 * parete,
                  d_fit_in, d_fit_out)
    #: si arrotonda per difetto a un decimo di millimetro
    d_barbo = max(math.floor(d_barbo * 1.0e4) / 1.0e4, 1.0e-3)
    v_barbo = Q_v / (math.pi / 4.0 * d_barbo ** 2)

    avvertenze_pg: list[str] = []

    def portagomma(an, basso, alto):
        x_b = 0.5 * ((basso + 0.5 * d_barbo)
                     + (alto - pendenza * 0.5 * d_barbo))
        h = an[2]
        return Portagomma(
            x_base=x_b, verso=0.0,
            R_posizione=_raggio_parete_locale(d, x_b) + t_max + 0.5 * h,
            angolo=math.pi, d_canna=engine.CANNA_DIAMETRO_INTERNO, d_foro=d_barbo,
            lunghezza=engine.PORTAGOMMA_LUNGHEZZA, n_creste=engine.PORTAGOMMA_CRESTE,
            passo_creste=engine.PORTAGOMMA_PASSO,
            interferenza=engine.PORTAGOMMA_INTERFERENZA)

    pg_in = portagomma(an_in, basso_in, alto_in)
    pg_out = portagomma(an_out, basso_out, alto_out)
    for nome_c, pg, an, d_fit in (("ingresso", pg_in, an_in, d_fit_in),
                                  ("uscita", pg_out, an_out, d_fit_out)):
        if pg.x_base + 0.5 * d_barbo < an[0] or pg.x_base - 0.5 * d_barbo > an[1]:
            avvertenze_pg.append(
                f"il portagomma d'{nome_c} e' centrato a x = {pg.x_base*1e3:.1f} mm e "
                f"l'anello va da {an[0]*1e3:.1f} a {an[1]*1e3:.1f}: il foro non "
                "incrocia l'anello e il circuito non e' collegato")
        v_b = Q_v / (math.pi / 4.0 * d_barbo ** 2)
        if d_fit < d_velocita - 1e-9 and v_b > V_ACQUA_MAX:
            avvertenze_pg.append(
                f"nell'anello d'{nome_c} ci sta un foro da {d_fit*1e3:.1f} mm contro i "
                f"{d_velocita*1e3:.1f} che vorrebbe la portata, e l'acqua ci passa a "
                f"{v_b:.1f} m/s: sopra il limite dei canali")

    def conc(K, v):
        return K * 0.5 * rho * v ** 2

    dp = (2.0 * conc(K_PORTAGOMMA, v_barbo)
          + conc(K_FORI, v_fori_in) + conc(K_FORI, v_fori_out)
          + conc(K_ANELLO, v_an_in) + conc(K_ANELLO, v_an_out)
          + conc(K_ANELLO, v_raccordo) + dp_giro_in + dp_giro_out)

    avvertenze: list[str] = list(avvertenze_pg)
    dp_cam = next(r.dp for r in rami if r.nome == "camera")
    dp_gola = next(r.dp for r in rami if r.nome == "gola")
    for nome_c, dp_giro, dp_can, n_can in (("ingresso", dp_giro_in, dp_cam, cam.n_canali),
                                           ("uscita", dp_giro_out, dp_gola, gola.n_canali)):
        #: MALDISTRIBUZIONE. I canali sboccano nell'anello in punti diversi, e
        #: l'anello ha una caduta: il canale vicino all'attacco vede una
        #: contropressione minore e prende piu' acqua.
        if dp_giro >= dp_can:
            avvertenze.append(
                f"la caduta lungo l'anello d'{nome_c} ({dp_giro/1e5:.3f} bar) e' "
                f"maggiore di quella di un canale ({dp_can/1e5:.3f}): i canali lontani "
                "dall'attacco non ricevono acqua")
        else:
            spread = 1.0 - math.sqrt(1.0 - dp_giro / dp_can)
            if spread > 0.10:
                avvertenze.append(
                    f"all'anello d'{nome_c} la portata fra il canale piu' favorito e "
                    f"il piu' sfavorito differisce del {spread*100:.0f} %: h cala del "
                    f"{(1-(1-spread)**0.8)*100:.0f} % sul canale peggiore fra i "
                    f"{n_can}, e il margine all'ebollizione va verificato su quello")
    for nome_c, v in (("ingresso", v_an_in), ("uscita", v_an_out)):
        if v > V_ACQUA_MAX:
            avvertenze.append(
                f"nell'anello d'{nome_c} l'acqua va a {v:.1f} m/s, sopra il limite dei "
                f"canali ({V_ACQUA_MAX:.0f}): non c'e' corsa assiale per fare la "
                "sezione che serve")
    if an_in[1] > cam.x_fine - 4.0e-3:
        avvertenze.append("l'anello d'ingresso arriva sotto il collettore d'uscita "
                          "della camera")

    for nome, val, unita, motivo in (
        ("n_fori_ingresso_acqua", float(n_in), "-",
         f"{n_in} fori radiali da {d_in*1e3:.1f} mm a {v_fori_in:.1f} m/s. Radiali per "
         f"forza, e per questo non superano {D_FORO_ORIZZONTALE_MAX*1e3:.1f} mm: un "
         "foro orizzontale piu' grande ha una calotta che in SLM cede"),
        ("n_fori_uscita_acqua", float(n_out), "-",
         f"{n_out} fori radiali da {d_out*1e3:.1f} mm a {v_fori_out:.1f} m/s"),
        ("altezza_anello_ingresso", an_in[2], "m",
         f"sezione a ROMBO ({area_in*1e6:.1f} mm2 su {L_in*1e3:.1f} mm di corsa): il "
         f"collare comincia dove comincia il mantello, quindi anche il pavimento del "
         f"vano segue la rampa. Acqua a {v_an_in:.1f} m/s, caduta lungo l'anello "
         f"{dp_giro_in/1e5:.3f} bar"),
        ("altezza_anello_uscita", an_out[2], "m",
         f"sezione a TRAPEZIO ({area_out*1e6:.1f} mm2 su {L_out*1e3:.1f} mm): qui il "
         f"collare comincia 6 mm prima del vano, la rampa e' finita e il pavimento e' "
         f"piatto. Acqua a {v_an_out:.1f} m/s, caduta {dp_giro_out/1e5:.3f} bar"),
        ("altezza_anello_raccordo", prof_r1 - prof_r0, "m",
         f"il vano che unisce i due collettori su tutta la circonferenza, DENTRO il "
         f"mantello: l'acqua ci passa a {v_raccordo:.2f} m/s, quindi non c'e' nessuna "
         "distribuzione da garantire e nessun foro da dimensionare"),
        ("d_portagomma", d_barbo, "m",
         f"foro di passaggio: {v_barbo:.1f} m/s. Il gambo esterno lo detta la canna "
         f"({engine.CANNA_DIAMETRO_INTERNO*1e3:.1f} mm), il foro la portata"),
        ("caduta_collettore_acqua", dp, "Pa",
         f"due portagomma {2*conc(K_PORTAGOMMA, v_barbo)/1e5:.3f}, fori "
         f"{(conc(K_FORI, v_fori_in)+conc(K_FORI, v_fori_out))/1e5:.3f}, anelli "
         f"{(conc(K_ANELLO, v_an_in)+conc(K_ANELLO, v_an_out)+conc(K_ANELLO, v_raccordo))/1e5:.3f}, "
         f"giri {(dp_giro_in+dp_giro_out)/1e5:.3f} bar"),
    ):
        reg.aggiungi(Scelta(nome, val, unita, Origine.DERIVATA, motivo))
    reg.aggiungi(Scelta(
        "d_canna_acqua", engine.CANNA_DIAMETRO_INTERNO, "m", Origine.DECISA,
        "1/2 pollice, il formato da giardinaggio piu' diffuso. NON E' STATO "
        "MISURATO: se la canna e' da 5/8 o 3/4 cambia il portagomma e nient'altro"))

    return CollettoreAcqua(
        ingresso=pg_in, anello_ingresso=an_in,
        x_collettore_ingresso=cam.x_inizio, n_fori_ingresso=n_in, d_foro_ingresso=d_in,
        raccordo=(x_r0, x_r1, prof_r0, prof_r1),
        uscita=pg_out, anello_uscita=an_out,
        x_collettore_uscita=x_gu, n_fori_uscita=n_out, d_foro_uscita=d_out,
        pendenza_tetto=pendenza, parete=parete, dp_totale=dp,
        avvertenze=tuple(avvertenze))


def _raggio_parete_locale(d, x: float) -> float:
    R_c, R_lip, L_c, L_conv = d["R_c"], d["R_lip"], d["L_c"], d["L_conv"]
    if x <= L_c:
        return R_c
    return R_c + min((x - L_c) / max(L_conv, 1e-12), 1.0) * (R_lip - R_c)


# --------------------------------------------------------------------------- #
# 4. iniettore a getti trasversali
# --------------------------------------------------------------------------- #
#: Frazione massima del salto d'iniezione che si accetta di spendere in un
#: raccordo. Non e' generosita': il salto d'iniezione e' TUTTO quello che c'e',
#: perche' p_c e' fissata a p_alimentazione/1.15 senza margini nascosti, e il
#: tubo a monte pesera' piu' del raccordo.
FRAZIONE_DP_RACCORDO = 0.15
#: Passaggio libero di un raccordo, stimato come diametro della punta meno 5 mm
#: (il codolo di un portagomma o di un innesto rapido non e' il nocciolo della
#: filettatura). E' una stima CONSERVATIVA e va confermata sul raccordo vero.
STRETTOIA_RACCORDO = 5.0e-3


def _iniettore(req: Requisiti, reg: Registro, op, params, l0, p_c: float) -> Iniettore:
    d = params.derived
    imp = req.impianto
    MW_gpl = _massa_molare(op)
    gamma_gpl = _gamma_gpl(op)

    def prova(R: float, n=None):
        return progetta_iniettore(
            mdot_aria=l0.mdot_air, mdot_gpl=l0.mdot_fuel_core, p_c=p_c,
            p_aria_monte=op.p_air_supply, p_gpl_monte=op.p_fuel_supply,
            T_aria=float(op.T_air_in), T_gpl=float(op.T_fuel_in),
            MW_aria=MW_ARIA, MW_gpl=MW_gpl, gamma_aria=GAMMA_ARIA,
            gamma_gpl=gamma_gpl, cd=float(op.cd_injector_ox),
            R_iniezione=R, n_getti=n)

    #: REGOLA DI SCELTA DEL RAGGIO, dichiarata prima di guardare i numeri.
    #: A portate e salti fissati vale d_getto = k H con k indipendente dal
    #: raggio, e H = f(A_aria, R): abbassare il raggio ALZA il diametro dei
    #: fori. Il corpo centrale si puo' strizzare, ed e' il solo modo di
    #: guadagnare diametro senza toccare nessuna pressione.
    #: Fra i punti che la correlazione di Holdeman giudica equivalenti non si
    #: sceglie il minimo di |C/C*-1|: la correlazione ha una dispersione molto
    #: piu' larga dell'1 %. Si sceglie il MASSIMO numero di getti che conserva
    #: il margine sul diametro minimo stampabile (piu' getti = un foro otturato
    #: fa meno danno), e fra quelli il raggio con C piu' vicino all'ottimo.
    d_min = req.processo.margine_foro(MARGINE_FORO)
    candidati = []
    r_cb = d["r_centerbody"]
    n_passi = 200
    for k in range(n_passi + 1):
        R = r_cb * (0.55 + 0.55 * k / n_passi)
        try:
            i = prova(R)
        except (ValueError, ZeroDivisionError):
            continue
        if i.d_getto >= d_min and i.n_getti >= 4 and abs(i.scarto_C) <= 0.10:
            candidati.append((R, i))
    if not candidati:
        #: Il diametro dei getti va come la radice della portata di
        #: combustibile, che va come la spinta: da qui si ricava la spinta
        #: MINIMA a cui questa architettura sa fare fori stampabili, ed e'
        #: un'informazione molto piu' utile di "non ci riesco".
        migliore = max(
            (prova(r_cb * (0.55 + 0.55 * k / n_passi)) for k in range(n_passi + 1)),
            key=lambda i: i.d_getto, default=None)
        rapporto = (d_min / migliore.d_getto) ** 2 if migliore else float("inf")
        raise MissingDatum(
            f"nessun raggio d'iniezione da' getti sopra {d_min*1e3:.2f} mm con almeno "
            "4 fori e C entro il 10 % dall'ottimo di Holdeman: il massimo ottenibile e' "
            f"{migliore.d_getto*1e3:.3f} mm. Con questa portata di combustibile i fori "
            "sarebbero piu' piccoli di cio' che la macchina fa. Servirebbe una spinta "
            f"di almeno {req.spinta*rapporto:.0f} N, oppure fori realizzati per "
            "foratura o elettroerosione dopo la stampa.")
    n_max = max(i.n_getti for _, i in candidati)
    R, inj = min((c for c in candidati if c[1].n_getti == n_max),
                 key=lambda c: abs(c[1].scarto_C))

    reg.aggiungi(Scelta(
        "R_iniezione", R, "m", Origine.DERIVATA,
        f"il corpo centrale si strizza da {r_cb*1e3:.2f} a {R*1e3:.2f} mm: e' il "
        f"raggio che da' il massimo numero di getti ({inj.n_getti}) sopra il minimo "
        f"stampabile con margine, e fra quelli C = {inj.C:.3f} contro l'ottimo "
        f"{inj.C_obiettivo} di Holdeman"))
    reg.aggiungi(Scelta("d_getto", inj.d_getto, "m", Origine.DERIVATA,
                        f"portata di GPL divisa in {inj.n_getti} getti al salto "
                        "disponibile"))
    reg.aggiungi(Scelta("n_getti", float(inj.n_getti), "-", Origine.DERIVATA,
                        "passo di Holdeman S = C H / sqrt(J) sulla circonferenza "
                        "d'iniezione"))
    reg.aggiungi(Scelta("H_condotto_aria", inj.altezza_anello, "m", Origine.DERIVATA,
                        "altezza della corona d'aria: A = pi (Re^2 - Ri^2), esatta e "
                        "non 2 pi R H, che a questo H/R sbaglia del 25 %"))
    reg.aggiungi(Scelta("J_iniezione", inj.J, "-", Origine.DERIVATA,
                        "rapporto dei flussi di quantita' di moto. NON e' un "
                        "parametro: salti uguali in frazione di p_c e densita' vicine "
                        "lo portano vicino a 1 da soli", fonte="Holdeman 1993"))
    reg.aggiungi(Scelta("C_Holdeman", inj.C, "-", Origine.NORMATIVA,
                        f"(S/H) sqrt(J); l'ottimo per getti da un lato e' "
                        f"{inj.C_obiettivo}", fonte="Prog. Energy Combust. Sci. 19 (1993) 31-70"))
    return inj


def _massa_molare(op) -> float:
    from zefiro.l0.mixture import MixtureModel
    g = MixtureModel.from_fuel(op.fuel).gas
    g.X = {k: v for k, v in op.fuel.composition.items()}
    return float(g.mean_molecular_weight)


def _gamma_gpl(op) -> float:
    from zefiro.l0.mixture import MixtureModel
    g = MixtureModel.from_fuel(op.fuel).gas
    g.TPX = float(op.T_fuel_in), 1.0e5, {k: v for k, v in op.fuel.composition.items()}
    return float(g.cp_mass / g.cv_mass)


def _filettatura(nome: str, mdot: float, rho: float, dp_disponibile: float,
                 reg: Registro) -> str:
    """La piu' PICCOLA filettatura che non si mangia il salto d'iniezione.

    E' un bilancio di pressione, non una scelta di catalogo. Un raccordo non e'
    un tubo: e' una strozzatura seguita da un allargamento brusco, e
    l'allargamento brusco non recupera niente.
    """
    righe = []
    for taglia in ("G1/8", "G1/4", "G3/8", "G1/2", "G3/4", "G1"):
        v, dp = perdita_raccordo(mdot, rho, BSPP[taglia]["punta"] - STRETTOIA_RACCORDO)
        righe.append(f"{taglia} {v:.0f} m/s {dp/1e5:.3f} bar")
        if dp <= FRAZIONE_DP_RACCORDO * dp_disponibile:
            reg.aggiungi(Scelta(
                f"filettatura_{nome}", BSPP[taglia]["d_esterno"], "m", Origine.NORMATIVA,
                f"{taglia}: la piu' piccola che spende meno del "
                f"{FRAZIONE_DP_RACCORDO:.0%} del salto d'iniezione "
                f"({dp/1e5:.3f} bar su {dp_disponibile/1e5:.3f}). "
                + "Scartate: " + "; ".join(righe[:-1] or ["nessuna"]),
                fonte="ISO 228-1"))
            return taglia
    raise MissingDatum(
        f"nessuna filettatura fino a G1 porta {mdot*1e3:.1f} g/s di {nome} senza "
        f"spendere piu' del {FRAZIONE_DP_RACCORDO:.0%} del salto d'iniezione: "
        + ", ".join(righe))


# --------------------------------------------------------------------------- #
# 5. la testa: la dimensiona il RACCORDO, non l'idraulica
# --------------------------------------------------------------------------- #
def _testa(req: Requisiti, reg: Registro, op, params, l0, inj: Iniettore,
           p_c: float) -> TestaIniezione:
    proc = req.processo
    d = params.derived
    dp_iniezione = req.frazione_dp_iniezione * p_c

    rho_aria = op.p_air_supply / (8314.462618 / MW_ARIA * float(op.T_air_in))
    rho_gpl = op.p_fuel_supply / (8314.462618 / _massa_molare(op) * float(op.T_fuel_in))
    t_aria = _filettatura("aria", l0.mdot_air, rho_aria, dp_iniezione, reg)
    t_gpl = _filettatura("GPL", l0.mdot_fuel_core, rho_gpl, dp_iniezione, reg)
    fa, fg = BSPP[t_aria], BSPP[t_gpl]

    d_filetto_gpl = fg["punta"] - proc.sottomisura_foro
    d_porta_aria = fa["punta"] - proc.sottomisura_foro
    reg.aggiungi(Scelta(
        "preforo_gpl", d_filetto_gpl, "m", Origine.DERIVATA,
        f"punta da {fg['punta']*1e3:.2f} mm meno la sottomisura di stampa: un foro "
        "SLM esce piu' piccolo del nominale e rugoso, e maschiare in un foro rugoso "
        "da' un filetto storto"))
    reg.aggiungi(Scelta("preforo_aria", d_porta_aria, "m", Origine.DERIVATA,
                        f"punta da {fa['punta']*1e3:.2f} mm meno la sottomisura"))

    #: Il corpo centrale nel plenum deve contenere la filettatura del GPL.
    R_corpo = 0.5 * fg["punta"] + 2.0 * proc.parete_minima
    reg.aggiungi(Scelta(
        "R_corpo_plenum", R_corpo, "m", Origine.DERIVATA,
        f"meta' della punta del {t_gpl} piu' due pareti minime. E' piu' grosso del "
        f"raggio d'iniezione ({inj.R_iniezione*1e3:.2f} mm): il corpo centrale si "
        "assottiglia attraverso la contrazione, dove il raggio piccolo serve"))

    #: Plenum del GPL dentro il corpo centrale: almeno quattro volte l'area dei
    #: getti (altrimenti i getti lontani dall'ingresso ricevono meno) e almeno
    #: due pareti minime dall'aria.
    r_bore_min = math.sqrt(4.0 * inj.n_getti * math.pi / 4.0 * inj.d_getto ** 2 / math.pi)
    r_bore_max = inj.R_iniezione - 2.0 * proc.parete_minima
    if r_bore_min > r_bore_max:
        raise MissingDatum(
            f"il plenum del GPL vorrebbe {r_bore_min*1e3:.2f} mm di raggio ma il "
            f"corpo centrale strizzato ne lascia {r_bore_max*1e3:.2f}.")
    r_bore = r_bore_max
    reg.aggiungi(Scelta(
        "r_plenum_gpl", r_bore, "m", Origine.DERIVATA,
        f"il massimo che lascia due pareti minime verso l'aria. Ne servivano almeno "
        f"{r_bore_min*1e3:.2f} mm perche' il plenum sia 4 volte l'area dei getti"))

    #: LA TESTA LA DIMENSIONA IL RACCORDO. La guarnizione a legare del raccordo
    #: dell'aria ha un diametro esterno, la faccia piana deve contenerlo, e la
    #: faccia sta su una bozza che deve appoggiare TUTTA sulla testa.
    d_boss = fa["d_guarnizione"] + 1.8e-3
    L_condotto = max(2.0 * inj.lunghezza_mescolamento, 4.0e-3)
    L_contrazione = 1.0 * (7.0e-3)     # ricalcolata sotto, serve un primo valore
    t_monte = 3.0 * proc.parete_minima
    t_testa = 3.0 * proc.parete_minima

    #: L_plenum e' il minimo che (a) fa appoggiare la bozza e (b) tiene il
    #: plenum piu' capiente del condotto.
    #: Il criterio e' "passaggio del plenum >= 2 volte quello del condotto", e
    #: 2 * L * (Rp - Rc) = 2 * A lo soddisfa ESATTAMENTE, cioe' lo sfiora: al
    #: primo motore diverso da Zefiro la verifica falliva per un millimetro
    #: quadrato di arrotondamento. Un criterio si soddisfa con margine, non al
    #: pelo: qui il 10 %.
    MARGINE_PLENUM = 1.10
    L_plenum = max(d_boss - t_monte - L_condotto, 4.0e-3)
    R_plenum = R_corpo + MARGINE_PLENUM * inj.area_aria / L_plenum
    L_contrazione = R_plenum - inj.R_medio      # cono a 45 gradi
    L_plenum = max(d_boss - t_monte - L_condotto - L_contrazione, 4.0e-3)
    R_plenum = R_corpo + MARGINE_PLENUM * inj.area_aria / L_plenum
    L_contrazione = max(R_plenum - inj.R_medio, 2.0e-3)
    reg.aggiungi(Scelta(
        "L_plenum_aria", L_plenum, "m", Origine.DERIVATA,
        f"il minimo che fa appoggiare tutta la bozza da {d_boss*1e3:.0f} mm del "
        f"raccordo {t_aria} sulla testa. Idraulicamente ne basterebbero meno: e' il "
        "raccordo a dimensionare la testa, non il contrario"))
    reg.aggiungi(Scelta(
        "R_plenum_aria", R_plenum, "m", Origine.DERIVATA,
        "raggio che rende la sezione di passaggio del plenum doppia di quella del "
        "condotto: sotto, sarebbe il plenum a distribuire male e i getti vedrebbero "
        "un flusso trasversale diverso da un settore all'altro"))
    reg.aggiungi(Scelta("L_contrazione", L_contrazione, "m", Origine.DERIVATA,
                        "cono a 45 gradi dal plenum al condotto: una parete continua "
                        "invece di uno spigolo porta il Cd da 0.6 a circa 0.95"))
    reg.aggiungi(Scelta(
        "L_condotto_aria", L_condotto, "m", Origine.DERIVATA,
        f"il doppio della lunghezza di mescolamento ({inj.lunghezza_mescolamento*1e3:.2f} "
        "mm): meta' a valle dei getti per il mescolamento, meta' a monte perche' i "
        "getti non stiano sullo spigolo d'ingresso"))

    x_getti = -0.5 * L_condotto
    reg.aggiungi(Scelta(
        "x_getti", x_getti, "m", Origine.DERIVATA,
        f"a meta' del condotto: restano {-x_getti*1e3:.2f} mm confinati a valle "
        f"contro i {inj.lunghezza_mescolamento*1e3:.2f} richiesti dal mescolamento"))

    R_testa = R_plenum + t_testa
    R_out = d["R_c"] + _spessore_mantello_max(req, reg)
    #: Raccordo esterno fra testa e mantello. Costruendo lungo +x e' una
    #: superficie rivolta in basso: perche' sia autosostentata serve
    #: dx >= dr / tan(90 - alpha), cioe' dx >= dr * tan(alpha) con alpha misurato
    #: dal piano. Farlo autosostentato costa lunghezza ma non costa supporti.
    dr = max(R_out - R_testa, 0.0)
    L_cono = dr * math.tan(math.radians(req.processo.angolo_progetto))
    reg.aggiungi(Scelta(
        "L_cono_testa", L_cono, "m", Origine.DERIVATA,
        f"passa da {R_testa*1e3:.1f} a {R_out*1e3:.1f} mm restando autosostentato a "
        f"{req.processo.angolo_progetto:.0f} gradi: cosi' non servono supporti "
        "esterni su quella faccia"))

    #: Rampa del corpo centrale dal raggio strizzato a quello di camera, dentro
    #: la camera, con semiangolo blando: e' un allargamento, e un allargamento
    #: troppo brusco stacca.
    SEMIANGOLO_RAMPA = math.radians(12.0)
    x_rampa = max((d["r_centerbody"] - inj.R_iniezione) / math.tan(SEMIANGOLO_RAMPA),
                  1.0e-3)
    reg.aggiungi(Scelta(
        "x_rampa_corpo_centrale", x_rampa, "m", Origine.DERIVATA,
        f"riporta il corpo centrale da {inj.R_iniezione*1e3:.2f} a "
        f"{d['r_centerbody']*1e3:.2f} mm con 12 gradi di semiangolo, dentro la camera "
        "dove di area ce n'e' da vendere"))

    #: Avvitamento del raccordo dell'aria: la bozza deve stare abbastanza fuori
    #: perche' il naso del raccordo, che entra oltre la filettatura, non arrivi
    #: a toccare il corpo centrale.
    gioco = 1.0e-3
    R_boss = max(R_plenum + fa["avvitamento"] - (R_plenum - R_corpo - gioco),
                 R_plenum + 6.0e-3)
    reg.aggiungi(Scelta(
        "R_bozza_aria", R_boss, "m", Origine.DERIVATA,
        f"il naso del raccordo {t_aria} entra oltre la filettatura: la bozza sta "
        f"abbastanza fuori perche' non tocchi il corpo centrale, con {gioco*1e3:.1f} "
        "mm di gioco"))

    x_monte = -(L_condotto + L_contrazione + L_plenum)
    x_porta = 0.5 * (x_monte - t_monte + 0.0)
    x_porta = min(max(x_porta, x_monte + 0.5e-3), -L_condotto - L_contrazione - 0.5e-3)

    return TestaIniezione(
        R_getti=inj.R_iniezione, R_anello=inj.R_medio, n_getti=inj.n_getti,
        d_getto=inj.d_getto, x_getti=x_getti, L_condotto=L_condotto,
        L_contrazione=L_contrazione, R_plenum=R_plenum, L_plenum=L_plenum,
        R_corpo_plenum=R_corpo, t_monte=t_monte, t_testa=t_testa,
        r_bore_gpl=r_bore, x_rampa=x_rampa, L_cono=L_cono,
        d_porta_aria=d_porta_aria, R_boss_aria=R_boss, d_boss_aria=d_boss,
        x_porta_aria_fissa=x_porta, d_filetto_gpl=d_filetto_gpl,
        L_filetto_gpl=fg["avvitamento"] + 2.0e-3)


def _spessore_mantello_max(req: Requisiti, reg: Registro) -> float:
    """Serve prima che i circuiti esistano: si usa il valore gia' registrato."""
    return reg.valore("spessore_mantello")


# --------------------------------------------------------------------------- #
# 6. orchestrazione
# --------------------------------------------------------------------------- #
def dimensiona(req: Requisiti) -> Dimensionamento:
    """Dai requisiti a tutte le quote, con la provenienza di ognuna.

    L'ordine non e' arbitrario: ogni passo consuma solo cio' che i precedenti
    hanno gia' dimostrato. Non ci sono iterazioni fra termochimica e geometria
    perche' su un motore a pressione costante il ciclo apparente si spezza: la
    termochimica non dipende dalla geometria.
    """
    reg = Registro()
    note: list[str] = []

    #: DUE PASSATE, e servono entrambe. Il film devia combustibile dal nucleo,
    #: quindi sposta il rapporto di equivalenza, quindi la c*, quindi la portata
    #: che da' la spinta richiesta. Con una sola passata la spinta usciva 50.30
    #: invece di 50.00: piccolo, ma e' un requisito, non un'aspirazione.
    op0, p_c, mdot0 = _punto_operativo(req, Registro())
    f_film = _frazione_di_film(req, op0, p_c, mdot0)
    #: PUNTO FISSO, e serve davvero. La frazione di film cambia il rapporto di
    #: equivalenza del nucleo, quindi T_ad, quindi il flusso termico sul plug,
    #: quindi la frazione di film che serve. La bisezione risolve il problema
    #: INTERNO; qui si verifica quello vero, sul dimensionamento completo, e si
    #: sale di mezzo punto finche' non torna. Senza, il progetto usciva a 1253 K
    #: contro un limite di 1250: tre gradi, ma dalla parte sbagliata.
    for _ in range(8):
        op, p_c, mdot = _punto_operativo(req, Registro(), f_film)
        params_p, l0_p = _vettore_di_progetto(req, Registro(), op, p_c, mdot, f_film)
        if stato_plug(req, op, params_p, l0_p, p_c, f_film)[0].T_finale <= T_LIMITE_PLUG:
            break
        f_film = round(f_film + 0.005, 4)
        if f_film > DESIGN_BOUNDS["f_film"][1]:
            raise MissingDatum(
                "la frazione di film necessaria supera il massimo ammesso dal "
                "vettore di progetto: il plug non e' salvabile a questa durata.")
    op, p_c, mdot = _punto_operativo(req, reg, f_film)
    durata = _durata(req, reg, mdot)
    reg.aggiungi(Scelta(
        "frazione_film_plug", f_film, "-", Origine.DERIVATA,
        "la MINIMA frazione di combustibile che tiene il corpo centrale sotto "
        f"{T_LIMITE_PLUG:.0f} K per tutta la raffica ANCHE SE il film rende "
        f"il {FATTORE_EFFICACIA_FILM:.0%} di quanto promette la correlazione. Senza film il plug fonde a "
        "2.8 s: non e' raffreddabile in nessun altro modo, e le alternative "
        "(gas dall'interno, acqua nel corpo centrale, troncare, ingrandire) sono "
        "state provate e scartate con dei numeri - vedi zefiro.film"))
    params, l0 = _vettore_di_progetto(req, reg, op, p_c, mdot, f_film)

    #: Lo spessore del mantello serve alla testa prima che i circuiti esistano.
    #: Si registra un valore provvisorio e si verifica a valle che i circuiti ci
    #: stiano dentro: se non ci stanno, il progetto viene rifiutato invece di
    #: essere tagliato in silenzio - che e' esattamente il difetto che apriva il
    #: circuito della gola sull'atmosfera.
    from zefiro.sdf.engine import spessore_mantello
    circuiti, rami, potenze, condizioni_acqua = _raffreddamento(req, reg, op, params, l0, p_c)
    t_max = max(spessore_mantello(c) for c in circuiti)
    reg.aggiungi(Scelta(
        "spessore_mantello", t_max, "m", Origine.DERIVATA,
        "il massimo fra i rami di parete_calda + canale (+ setto + canale di "
        "ritorno) + parete_fredda. Contare il ritorno NON e' un dettaglio: senza, "
        "lo strato di ritorno veniva tagliato sulla pelle e il circuito della gola "
        "era aperto sull'atmosfera"))

    inj = _iniettore(req, reg, op, params, l0, p_c)
    testa = _testa(req, reg, op, params, l0, inj, p_c)
    #: il collettore dell'acqua si dimensiona DOPO la testa: il suo plenum sta
    #: dentro la piastra, e deve sapere dove finisce quello dell'aria.
    collettore = _collettore_acqua(req, reg, params.derived, circuiti, rami, testa,
                                   *condizioni_acqua)
    #: L'INTERFACCIA DI MONTAGGIO. Fino a questa revisione non c'era: il motore
    #: non si attaccava a niente, ed era un pezzo bellissimo da tenere in mano.
    montaggio = _montaggio(req, reg, params.derived, testa, l0.thrust, collettore)
    for a in montaggio.avvertenze:
        note.append(f"montaggio: {a}")

    #: Tempi: il confronto che dice se il motore ha senso.
    rho_c = p_c * l0.MW_c / (8314.462618 * l0.T_ad)
    d = params.derived
    A_ann = math.pi * (d["R_c"] ** 2 - d["r_centerbody"] ** 2)
    u_c = (l0.mdot_air + l0.mdot_fuel_core) / (rho_c * A_ann)
    t_perm = d["L_c"] / u_c
    t_mix = inj.lunghezza_mescolamento / inj.V_aria
    reg.aggiungi(Scelta("tempo_permanenza", t_perm, "s", Origine.DERIVATA,
                        f"L_c / velocita' di camera ({u_c:.1f} m/s)"))
    reg.aggiungi(Scelta("tempo_mescolamento", t_mix, "s", Origine.DERIVATA,
                        "lunghezza di mescolamento di Holdeman diviso la velocita' "
                        "del flusso trasversale"))

    fuori = [
        f"bombola sotto {req.impianto.T_bombola_min-273.15:.0f} C: p_sat scende sotto "
        f"{reg.valore('p_combustibile')/1e5:.2f} bar e il GPL non ha piu' il salto per "
        "entrare in camera; la miscela smagrisce fino allo spegnimento",
        f"bombola piu' calda: il riduttore parzializza e p_c resta {p_c/1e5:.2f} bar, "
        "ma il fondo scarico del serbatoio non cambia e la durata resta quella",
        f"serbatoio sotto {reg.valore('p_ossidante')/1e5:.2f} bar: il riduttore "
        "dell'aria perde il controllo PRIMA di quello del GPL e la miscela va ricca "
        "in modo incontrollato, con il motore gia' caldo",
    ]
    inviluppo = Inviluppo(
        T_bombola_minima=req.impianto.T_bombola_min, durata_raffica=durata,
        spinta=l0.thrust, isp=l0.Isp_s, p_camera=p_c, epsilon=l0.epsilon,
        fuori_inviluppo=tuple(fuori))

    for a in inj.avvertenze:
        note.append(f"iniettore: {a}")
    note.extend(testa.verifica(inj.lunghezza_mescolamento))
    for r in rami:
        note.append(
            f"ramo {r.nome}: {r.n_canali} canali da {r.lato*1e3:.1f} mm a "
            f"{r.velocita:.1f} m/s, h {r.h/1e3:.1f} kW/m2K, parete gas "
            f"{r.T_parete_gas:.0f} K, margine all'ebollizione "
            f"{r.margine_ebollizione:.0f} K, caduta {r.dp/1e5:.3f} bar")
    #: BILANCIO DI CARICO DELLA SERIE. In parallelo le cadute dei due rami si
    #: devono PAREGGIARE, e qui c'era una strozzatura calibrata a farlo. In
    #: serie si SOMMANO, e non c'e' niente da calibrare: il conto o torna o non
    #: torna, e se non torna e' il progetto a essere sbagliato, non un foro.
    dp_serie = sum(r.dp for r in rami) + collettore.dp_totale
    reg.aggiungi(Scelta(
        "caduta_acqua_totale", dp_serie, "Pa", Origine.DERIVATA,
        "somma delle cadute dei due tratti in serie piu' quelle del collettore: "
        + ", ".join(f"{r.nome} {r.dp/1e5:.3f} bar" for r in rami)
        + f", collettore {collettore.dp_totale/1e5:.3f} bar. In serie si sommano, "
        "e questo e' cio' che la canna deve vincere"))
    if dp_serie > DP_ACQUA_MAX:
        note.append(
            f"la serie chiede {dp_serie/1e5:.2f} bar contro un limite di "
            f"{DP_ACQUA_MAX/1e5:.1f}: la rete del banco deve garantirli, oppure "
            "vanno allargati i canali del tratto che cade di piu'")
    note.append(
        f"raffreddamento IN SERIE: canna -> anello d'ingresso -> camera -> raccordo "
        f"-> gola -> inversione al labbro -> ritorno -> anello d'uscita -> canna. "
        f"Un solo ingresso e una sola uscita, per canna da "
        f"{collettore.ingresso.d_canna*1e3:.1f} mm, e nessuna strozzatura di "
        f"bilanciamento da azzeccare")
    for a in collettore.avvertenze:
        note.append(f"collettore acqua: {a}")

    #: IL FILM SI PUO' ALIMENTARE? Domanda diversa da "il film funziona", e
    #: finora non se la poneva nessuno.
    area_film, area_min_film, note_film = alimentabilita_film(
        req, op, l0, p_c, _DP_INIEZIONE_FRAZIONE * p_c, params.derived)
    note.extend(note_film)
    _, film_prova = stato_plug(req, op, params, l0, p_c, f_film)
    if film_prova is not None:
        for a in film_prova.avvertenze:
            note.append(f"film: {a}")
    if l0.mdot_fuel_film > 0.0:
        reg.aggiungi(Scelta(
            "area_metering_film", area_film, "m2", Origine.DERIVATA,
            f"sezione di efflusso che fa passare esattamente {l0.mdot_fuel_film*1e3:.2f} "
            f"g/s di film sotto il salto d'iniezione: corrisponde a un foro da "
            f"{math.sqrt(4.0*area_film/math.pi)*1e3:.2f} mm, contro un minimo di "
            f"processo di {math.sqrt(4.0*area_min_film/math.pi)*1e3:.2f} mm che pero' "
            "e' a sua volta un valore ASSUNTO e non misurato (TODO)"))

    plug, film = stato_plug(req, op, params, l0, p_c, f_film)
    fessura = None
    if film is not None:
        fessura = _fessura_film(req, reg, params.derived, params, testa, l0,
                                film, p_c)
        for a in fessura.avvertenze:
            note.append(f"fessura di film: {a}")
    reg.aggiungi(Scelta(
        "T_plug_fine_raffica", plug.T_finale, "K", Origine.DERIVATA,
        f"transitorio del corpo centrale trattato come termicamente sottile "
        f"(penetrazione 7.3 mm contro 3.9 di raggio), con il flusso che cala "
        f"mentre la parete si scalda: {plug.potenza_iniziale:.0f} W all'inizio, "
        f"{plug.potenza_finale:.0f} W alla fine"))
    return Dimensionamento(
        registro=reg, requisiti=req, op=op, x=None, params=params, l0=l0,
        iniettore=inj, testa=testa, circuiti=circuiti, inviluppo=inviluppo,
        collettore=collettore, montaggio=montaggio, fessura=fessura,
        carico=potenze, plug=plug, film=film, note=note)


def costruisci_geometria(dim: Dimensionamento, passo: float, raccordo: float = 5.0e-4,
                         xp=None):
    from zefiro.sdf.core import backend
    from zefiro.sdf.engine import costruisci

    p = dim.params
    return costruisci(p.derived, p.plug_contour_x, p.plug_contour_r, dim.circuiti,
                      passo, raccordo=raccordo, xp=xp or backend("numpy"),
                      testa=dim.testa, collettore=dim.collettore,
                      montaggio=dim.montaggio, film=dim.fessura)


# --------------------------------------------------------------------------- #
# 7. il film che salva il corpo centrale
# --------------------------------------------------------------------------- #
#: Temperatura massima ammessa per il plug [K]. Non e' la fusione (1673): e' il
#: punto oltre il quale il 316L non tiene piu' la forma sotto il carico di
#: pressione del gas, e una punta di aerospike che si deforma cambia l'ugello.
T_LIMITE_PLUG = 1250.0
#: Altezza della fessura anulare del film [m]. Il risultato NON dipende da
#: questo numero: nella correlazione xi va come x/(M s) e M va come 1/s, quindi
#: il prodotto M*s e' fissato dalla portata. Conta solo la frazione di
#: combustibile. 0.5 mm e' scelto perche' e' la piu' stretta che lascia uscire
#: la polvere.
ALTEZZA_FESSURA_FILM = 5.0e-4
#: Quanto si crede alla correlazione del film. 0.5 = il progetto deve reggere
#: anche se il film rende la meta' di quanto Stollery & El-Ehwany promettono.
#: La correlazione viene da lastra piana con gas simili; qui la parete e'
#: conica, il gas sta 1800 K sopra il refrigerante e il numero di Mach non e'
#: quello degli esperimenti. Dimensionare a 1.0 vorrebbe dire giocarsi il pezzo
#: sulla fedelta' di una formula del 1965.
FATTORE_EFFICACIA_FILM = 0.5
#: Temperatura del GPL alla fessura: ha percorso tutto il corpo centrale, che e'
#: caldo. 400 K e' una stima prudente in eccesso (piu' caldo = film meno
#: efficace); va confermata da una CFD coniugata.
T_FILM = 400.0


def _stazioni_plug(op, params, l0, p_c):
    """(x dalla fessura, area, h_gas, T_aw) per ogni tratto del plug."""
    d = params.derived
    x_lip = d["L_c"] + d["L_conv"]
    righe = [r for r in mappa_di_calore(op, params, l0, p_c) if r[2] == "plug"]
    staz = []
    for a, b in zip(righe[:-1], righe[1:]):
        ds = math.hypot(b[0] - a[0], b[1] - a[1])
        staz.append((0.5 * (a[0] + b[0]) - x_lip,
                     math.pi * (a[1] + b[1]) * ds,
                     0.5 * (a[6] + b[6]), 0.5 * (a[5] + b[5])))
    return staz


def _massa_plug(params) -> float:
    cx, cr = params.plug_contour_x, params.plug_contour_r
    v = 0.0
    for (x0, r0), (x1, r1) in zip(zip(cx[:-1], cr[:-1]), zip(cx[1:], cr[1:])):
        v += math.pi / 3.0 * abs(x1 - x0) * (r0 * r0 + r0 * r1 + r1 * r1)
    return RHO_316L_ARCH * v


RHO_316L_ARCH, CP_316L_ARCH, T_FUS_ARCH = 7990.0, 550.0, 1673.0


#: Coefficiente di efflusso di un foro stampato, lo stesso usato per gli
#: iniettori. Serve qui per chiedersi se il film si puo' ALIMENTARE, non solo se
#: funziona una volta alimentato.
CD_FORO_FILM = 0.75

#: Frazione della pressione di camera spesa nel salto d'iniezione. E' la stessa
#: usata per gli iniettori: il film esce nella stessa camera e vede lo stesso
#: salto, ed e' proprio questo il motivo per cui la sua sezione non e' libera.
_DP_INIEZIONE_FRAZIONE = 0.15


def alimentabilita_film(req: Requisiti, op, l0, p_c: float,
                        dp_iniezione: float, d) -> tuple[float, float, list[str]]:
    """Il film si puo' alimentare? Restituisce (area richiesta, area minima
    stampabile, avvertenze).

    LA DOMANDA CHE MANCAVA. Il modello del film dimensiona la FESSURA da cui
    l'aria protettiva esce, e verifica che protegga. Non si chiedeva da dove
    arrivasse. E la risposta non e' scontata, perche' il film esce nella stessa
    camera in cui esce il combustibile principale, e quindi vede LO STESSO salto
    di pressione: se ha una sezione confrontabile con quella dei getti, si porta
    via una frazione confrontabile della portata - non il 9 %.

    Il conto: alla portata di film voluta corrisponde una sezione di efflusso

        A = mdot / (Cd sqrt(2 rho dp))

    e se quella sezione e' piu' piccola del piu' piccolo foro che il processo sa
    fare, il film NON e' alimentabile da un foro. Va allora strozzato per
    ATTRITO (un condotto lungo e stretto) oppure alimentato da una linea propria
    a pressione piu' bassa - due architetture diverse, non un dettaglio.
    """
    rho = p_c * (1.0 + dp_iniezione / p_c) * _massa_molare(op) / (8314.462618 * 293.0)
    area = l0.mdot_fuel_film / (CD_FORO_FILM * math.sqrt(2.0 * rho * dp_iniezione))
    d_min = req.processo.margine_foro(MARGINE_FORO)
    area_min = math.pi / 4.0 * d_min ** 2
    note: list[str] = []
    if l0.mdot_fuel_film <= 0.0:
        return 0.0, area_min, note
    if area < area_min:
        note.append(
            f"il film chiede {area*1e6:.3f} mm2 di sezione di efflusso per passare "
            f"{l0.mdot_fuel_film*1e3:.2f} g/s sotto {dp_iniezione/1e5:.2f} bar, contro "
            f"i {area_min*1e6:.3f} mm2 del piu' piccolo foro stampabile "
            f"({d_min*1e3:.2f} mm): un solo foro minimo ne farebbe passare "
            f"{area_min/area:.0f} volte tanto. Il film NON e' alimentabile da un foro "
            "sul plenum del combustibile: o lo si strozza per attrito con un condotto "
            "lungo, o arriva da una linea propria a pressione piu' bassa. Finche' non "
            "e' deciso, la frazione di film e' un desiderio e non una portata")
    return area, area_min, note


def stato_plug(req: Requisiti, op, params, l0, p_c: float, f_film: float):
    """Temperatura del plug a fine raffica, con o senza film."""
    from zefiro.film import dimensiona_film, transitorio_plug

    d = params.derived
    staz = _stazioni_plug(op, params, l0, p_c)
    massa = _massa_plug(params)
    film = None
    if f_film > 0.0:
        rho_g = p_c * l0.MW_c / (8314.462618 * l0.T_ad)
        A_ann = math.pi * (d["R_c"] ** 2 - d["r_centerbody"] ** 2)
        u_g = (l0.mdot_air + l0.mdot_fuel_core + l0.mdot_fuel_film) / (rho_g * A_ann)
        rho_c = p_c / (8314.462618 / _massa_molare(op) * T_FILM)
        film = dimensiona_film(l0.mdot_fuel_film, ALTEZZA_FESSURA_FILM,
                               d["r_centerbody"], rho_c, 8.0e-6, rho_g, u_g,
                               T_FILM, f_film)
    return transitorio_plug(staz, massa, req.durata, film, CP_316L_ARCH,
                            T_FUS_ARCH,
                            fattore_efficacia=FATTORE_EFFICACIA_FILM), film


def _frazione_di_film(req: Requisiti, op, p_c: float, mdot: float) -> float:
    """La MINIMA frazione di combustibile che tiene il plug sotto la soglia.

    Bisezione, non tentativi: la frazione entra nel ciclo L0 (sposta phi del
    nucleo e quindi T_ad) e nella geometria, quindi a ogni tentativo si rifa'
    tutto il dimensionamento. Sono un centinaio di millisecondi a colpo.
    """
    reg_finto = Registro()

    def temperatura(f: float) -> float:
        params, l0 = _vettore_di_progetto(req, Registro(), op, p_c, mdot, f)
        return stato_plug(req, op, params, l0, p_c, f)[0].T_finale

    if temperatura(0.0) <= T_LIMITE_PLUG:
        return 0.0
    lo, hi = 0.0, min(DESIGN_BOUNDS["f_film"][1], 0.35)
    if temperatura(hi) > T_LIMITE_PLUG:
        raise MissingDatum(
            f"nemmeno deviando il {hi:.0%} del combustibile al film il corpo "
            f"centrale resta sotto {T_LIMITE_PLUG:.0f} K per {req.durata:.1f} s. "
            "Con questa architettura la raffica non puo' durare cosi' a lungo.")
    for _ in range(14):
        mid = 0.5 * (lo + hi)
        if temperatura(mid) > T_LIMITE_PLUG:
            lo = mid
        else:
            hi = mid
    #: si prende il capo ALTO della bisezione: e' quello che soddisfa il
    #: vincolo, e su una correlazione con questa dispersione arrotondare per
    #: difetto sarebbe ottimismo.
    return math.ceil(hi * 200.0) / 200.0        # arrotondato allo 0.5 %
