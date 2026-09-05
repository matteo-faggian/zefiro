#!/usr/bin/env python3
"""Disuniformita' della miscela sui piani campionati dalla CFD.

PERCHE' NON SI RILEGGE IL VOLUME. La prima versione tagliava una fetta di
spessore finito dal campo di volume e ne faceva la statistica. Due difetti:
la fetta mescola stazioni diverse (a valle dei getti il gradiente assiale e'
forte, quindi la fetta APPIATTISCE proprio la quantita' che si misura), e i
campi di volume si possono scrivere solo poche volte perche' pesano. I piani
campionati da OpenFOAM durante il calcolo sono tagli ESATTI e costano poco:
si scrivono fitti nel tempo e dicono anche QUANDO il regime e' arrivato.

LA MISURA.

    U(x) = sqrt( < (Y - <Y>)^2 > ) / <Y>

con medie pesate sulla PORTATA rho*u_x*dA e non sull'area: conta come il
combustibile e' distribuito nel flusso, non nello spazio. Un ristagno pieno di
GPL fermo non alimenta niente e non deve contare come miscela.

VERIFICA INDIPENDENTE. Su ogni piano a valle dei getti la media pesata sulla
portata deve valere mdot_gpl/(mdot_gpl+mdot_aria). Se non torna, il conto e'
sbagliato o il transitorio non e' finito: e' un controllo che non costa nulla
e che smaschera l'errore piu' probabile, cioe' pesare male.
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))


def leggi_vtp(path: Path):
    """(punti, poligoni, campi) da un .vtp scritto dal functionObject."""
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy

    r = vtk.vtkXMLPolyDataReader()
    r.SetFileName(str(path))
    r.Update()
    d = r.GetOutput()
    punti = vtk_to_numpy(d.GetPoints().GetData()).astype(np.float64)
    celle = d.GetPolys().GetData()
    conn = vtk_to_numpy(celle).astype(np.int64) if celle.GetNumberOfTuples() else np.empty(0, np.int64)
    poligoni, i = [], 0
    while i < conn.size:
        n = conn[i]
        poligoni.append(conn[i + 1:i + 1 + n])
        i += 1 + n
    pd = d.GetPointData()
    campi = {}
    for k in range(pd.GetNumberOfArrays()):
        a = pd.GetArray(k)
        if a is None:
            continue
        campi[a.GetName()] = vtk_to_numpy(a).astype(np.float64)
    return punti, poligoni, campi


def per_faccia(punti, poligoni, campi, nomi):
    """Area vettoriale e valori medi di faccia.

    I poligoni del taglio non sono triangoli: si triangolano a ventaglio dal
    primo vertice. L'area vettoriale di un poligono piano e' la somma delle
    aree vettoriali dei triangoli del ventaglio, esattamente.
    """
    A = np.zeros((len(poligoni), 3))
    val = {n: np.zeros((len(poligoni),) + campi[n].shape[1:]) for n in nomi}
    for i, p in enumerate(poligoni):
        P = punti[p]
        a = np.zeros(3)
        for j in range(1, len(p) - 1):
            a += 0.5 * np.cross(P[j] - P[0], P[j + 1] - P[0])
        A[i] = a
        for n in nomi:
            val[n][i] = campi[n][p].mean(axis=0)
    return A, val


def disuniformita(path: Path):
    """(U, Y_medio, portata) sul piano campionato in `path`."""
    punti, poligoni, campi = leggi_vtp(path)
    if not poligoni:
        return None
    nomi = [n for n in ("C3H8", "U", "rho") if n in campi]
    if "C3H8" not in nomi or "U" not in nomi:
        return None
    if "T" in campi:
        nomi.append("T")
    A, val = per_faccia(punti, poligoni, campi, nomi)
    #: il piano ha normale x: il flusso di massa attraverso la faccia e'
    #: rho * (u . A). Si prende il valore assoluto perche' un ricircolo che
    #: torna indietro trasporta comunque combustibile.
    rho = val["rho"] if "rho" in val else np.ones(len(poligoni))
    w = np.abs(rho * np.einsum("ij,ij->i", val["U"], A))
    if w.sum() <= 0:
        return None
    Y = val["C3H8"]
    Ym = float(np.average(Y, weights=w))
    if Ym <= 1e-12:
        return None
    var = float(np.average((Y - Ym) ** 2, weights=w))
    #: QUANTA DELLA PORTATA ATTRAVERSA CELLE FREDDE. Il limite numerico
    #: dichiarato in `fvOptions` taglia la temperatura a 250 K, e in una zona
    #: attorno al getto ci arriva davvero: sono celle in cui la densita' - e
    #: quindi il peso con cui entrano in questa media - non e' quella vera.
    #: Se la loro quota di portata e' trascurabile sul piano di misura, il
    #: risultato non ne risente; se non lo e', va detto prima del numero.
    freddo = 0.0
    if "T" in val:
        freddo = float(w[val["T"] < 260.0].sum() / w.sum())
    return math.sqrt(var) / Ym, Ym, float(w.sum()), freddo


def tempo_di_lavaggio(d, x_target: float, n: int = 4000) -> float:
    """Tempo per rinnovare TUTTO il fluido a monte della stazione `x` [s].

    E' il volume a monte diviso la portata volumetrica, che con densita'
    costante si scrive come integrale lungo l'asse

        t(x) = integrale da x_in a x di  A(s)/Q ds  =  integrale di ds/V(s)

    e su questa geometria e' esatto, non stimato: A(s) e' nota in forma chiusa.

    NON E' IL TEMPO DI ARRIVO DEL COMBUSTIBILE, e la differenza conta. Il
    fronte di GPL viaggia sul nucleo veloce e arriva molto prima: nelle
    immagini a 30 us era gia' a x = +5 mm, dove questo conto da' 500 us. Il
    fronte e' la punta; questo e' il tempo perche' anche l'ultimo angolo di
    ricircolo sia stato sostituito. Le STATISTICHE (media e scarto della
    frazione di massa su un piano) dipendono da tutto il fluido a monte, non
    dalla punta: e' questo il numero che dice quando una statistica smette di
    ricordarsi delle condizioni iniziali. Ed e' un LIMITE INFERIORE, perche' un
    ricircolo si svuota piu' lentamente della media.

    PERCHE' NE SERVE UNO PER PIANO. Il dominio ha due tempi caratteristici che
    differiscono di venti volte:

        anello d'iniezione   36.1 mm2   120 m/s      55 us
        camera              410.0 mm2   10.6 m/s   1132 us

    La camera ha 11.4 volte l'area dell'anello, quindi la stessa portata ci
    scorre undici volte piu' lenta. Un verdetto unico "il transitorio e'
    finito" per tutti i piani sarebbe sbagliato in partenza: a 250 us i piani
    dentro l'anello sono lavati quattro volte e quelli in camera non hanno
    ancora visto passare un quarto del loro volume.
    """
    mdot = d["mdot_gpl"] + d["mdot_aria"]
    rho = d["rho_aria"]
    x0 = d["x_ingresso"]
    if x_target <= x0:
        return 0.0
    xs = np.linspace(x0, x_target, n)
    r_i = np.where(xs < 0.0, d["R_getti"],
                   d["R_getti"] + np.clip(xs / d["x_rampa"], 0.0, 1.0)
                   * (d["r_centerbody"] - d["R_getti"]))
    R_e = np.where(xs < 0.0, d["R_anello"], d["R_c"])
    A = math.pi * (R_e ** 2 - r_i ** 2)
    V = mdot / (rho * A)
    integra = getattr(np, "trapezoid", None) or np.trapz
    return float(integra(1.0 / V, xs))


#: Quante volte il tempo di lavaggio bisogna aver simulato prima di credere a
#: una statistica su quel piano. Non e' un numero da manuale: e' la regola
#: pratica dei tre-quattro lavaggi, e vale come SOGLIA, non come garanzia -
#: la stazionarieta' vera la decide il confronto fra due finestre temporali,
#: che e' quello che fa `--storia`. Questo serve solo a non perdere tempo a
#: guardare numeri che di sicuro non sono pronti.
LAVAGGI_MINIMI = 3.0


def main() -> int:
    from mesh_iniettore import quote_dal_progetto

    ap = argparse.ArgumentParser()
    ap.add_argument("--caso", type=Path, required=True)
    ap.add_argument("--tempo", default=None, help="istante (default: l'ultimo)")
    ap.add_argument("--storia", action="store_true",
                    help="stampa l'andamento nel tempo invece del solo istante")
    ap.add_argument("--csv", type=Path, default=None)
    ap.add_argument("--campionamento", default="superfici",
                    help="cartella di postProcessing da leggere")
    ap.add_argument("--tutti", action="store_true",
                    help="stampa ogni istante invece delle sole statistiche")
    a = ap.parse_args()

    base = a.caso / "postProcessing" / a.campionamento
    if not base.exists():
        raise SystemExit(f"{base} non esiste")
    tempi = sorted((p for p in base.iterdir() if p.is_dir()),
                   key=lambda p: float(p.name))
    if not tempi:
        raise SystemExit("nessun istante campionato")

    d = quote_dal_progetto()
    xg, Lm = d["x_getti"], d["L_mescolamento"]
    piani = d["piani_x"]
    Y_atteso = d["mdot_gpl"] / (d["mdot_gpl"] + d["mdot_aria"])

    print(f"getti a x = {xg*1e3:.2f} mm;  la correlazione di Holdeman promette")
    print(f"miscela uniforme entro {Lm*1e3:.2f} mm, cioe' a x = {(xg+Lm)*1e3:.2f} mm")
    print(f"frazione di massa di GPL attesa a valle: {Y_atteso:.5f}\n")

    #: il tempo che ogni piano richiede prima che il suo numero voglia dire
    #: qualcosa. Si stampa ACCANTO al numero, non in una nota a pie' di pagina:
    #: una tabella che mette sulla stessa riga un valore maturo e uno che non
    #: lo e' senza dirlo invita a leggerli allo stesso modo.
    lavaggio = [tempo_di_lavaggio(d, x) for x in piani]

    def tabella(t: Path):
        t_ora = float(t.name)
        print(f"istante {t.name} s  ({t_ora*1e6:.1f} us)")
        print(f"  {'piano':>7} {'x [mm]':>9} {'(x-xg)/Lmix':>12} "
              f"{'Y medio':>10} {'scarto %':>9} {'U':>8} {'freddo':>8} "
              f"{'t lavaggio':>11} {'lavaggi':>8}  stato")
        ultimo = None
        for i, x in enumerate(piani):
            f = t / f"sez_{i:02d}.vtp"
            if not f.exists():
                continue
            r = disuniformita(f)
            if r is None:
                continue
            U, Ym, _, freddo = r
            sc = 100.0 * (Ym - Y_atteso) / Y_atteso if x > xg else float("nan")
            ta = lavaggio[i]
            lav = t_ora / ta if ta > 0 else float("inf")
            if lav >= LAVAGGI_MINIMI:
                stato = "maturo"
            elif lav >= 1.0:
                stato = "acerbo"
            else:
                stato = "NON ANCORA RAGGIUNTO"
            nota = ""
            if abs(x - (xg + Lm)) < 1e-9:
                nota = "  <-- fine lunghezza di mescolamento"
            if abs(x) < 1e-9:
                nota = "  <-- sbocco in camera"
            print(f"  sez_{i:02d} {x*1e3:9.2f} {(x-xg)/Lm:12.2f} {Ym:10.5f} "
                  f"{sc:9.1f} {U:8.3f} {freddo:7.1%} {ta*1e6:9.1f} "
                  f"{lav:8.1f}  {stato}{nota}")
            if abs(x - (xg + Lm)) < 1e-9 and lav >= LAVAGGI_MINIMI:
                ultimo = U
        print(f"\n  \"lavaggi\" = tempo simulato / tempo di lavaggio del volume a "
              f"monte di quel piano.\n  Sotto {LAVAGGI_MINIMI:.0f} il numero "
              "non e' un risultato: e' una fotografia del transitorio.\n"
              "  Il FRONTE di combustibile arriva molto prima - ma una media e "
              "uno scarto\n  dipendono da tutto il fluido a monte, non dalla "
              "punta.")
        return ultimo

    if a.storia:
        idx = next((i for i, x in enumerate(piani)
                    if abs(x - (xg + Lm)) < 1e-9), None)
        righe = []
        for t in tempi:
            f = t / f"sez_{idx:02d}.vtp"
            if idx is None or not f.exists():
                continue
            r = disuniformita(f)
            if r is None:
                continue
            righe.append((float(t.name), r[0], r[1], r[3]))
        if not righe:
            print("nessun istante utile sul piano di fine mescolamento\n")
        else:
            if a.tutti:
                print(f"  {'t [us]':>9} {'U':>10} {'Y medio':>10} {'freddo':>8}")
                for t, u, y, f in righe:
                    print(f"  {t*1e6:9.2f} {u:10.3f} {y:10.5f} {f:7.1%}")
            #: UN GETTO TRASVERSALE NON E' STAZIONARIO. Si stacca, oscilla, e la
            #: disuniformita' letta a un istante qualunque puo' valere il doppio
            #: di quella letta un microsecondo dopo. Leggere l'ultimo fotogramma
            #: e chiamarlo risultato e' un modo di ottenere qualunque risposta si
            #: voglia. Si media allora sull'ultimo terzo della corsa, e si
            #: confronta con il terzo precedente: se le due medie non
            #: coincidono, il transitorio non e' finito e il numero non vale.
            n = len(righe)
            terzo = max(n // 3, 1)
            def statistiche(blocco):
                u = np.array([r[1] for r in blocco])
                y = np.array([r[2] for r in blocco])
                f = np.array([r[3] for r in blocco])
                return u.mean(), u.std(), y.mean(), f.mean()
            u2, s2, y2, f2 = statistiche(righe[-terzo:])
            u1, s1, y1, f1 = statistiche(righe[-2 * terzo:-terzo] or righe[:terzo])
            deriva = abs(u2 - u1) / max(u2, 1e-12)
            print(f"  finestra    {'U medio':>10} {'scarto tipo':>12} {'Y medio':>10} "
                  f"{'portata fredda':>15}")
            print(f"  penultimo terzo {u1:10.3f} {s1:12.3f} {y1:10.5f} {f1:14.1%}")
            print(f"  ultimo terzo    {u2:10.3f} {s2:12.3f} {y2:10.5f} {f2:14.1%}")
            print(f"  deriva fra le due finestre: {deriva*100:.1f} %"
                  + ("   -> STAZIONARIO" if deriva < 0.10 else
                     "   -> TRANSITORIO NON FINITO: il numero non vale"))
            print(f"  oscillazione entro l'ultimo terzo: {s2/max(u2,1e-12)*100:.0f} %"
                  " (il getto e' instabile, non e' rumore numerico)")
        if a.csv:
            a.csv.write_text("t,U,Ymedio,portata_fredda\n" + "\n".join(
                f"{t},{u},{y},{f}" for t, u, y, f in righe), encoding="utf-8")
        print()

    t = tempi[-1] if a.tempo is None else next(
        p for p in tempi if p.name == a.tempo)
    U = tabella(t)
    if U is None:
        #: DUE MOTIVI DIVERSI, e confonderli sarebbe la cosa peggiore: uno e'
        #: "manca il dato", l'altro e' "il dato c'e' ma non e' ancora un
        #: risultato". Il secondo, detto come il primo, invita a rilanciare un
        #: calcolo che sta gia' andando bene.
        i = next((k for k, x in enumerate(piani)
                  if abs(x - (xg + Lm)) < 1e-9), None)
        esiste = i is not None and (t / f"sez_{i:02d}.vtp").exists()
        if not esiste:
            print("\nNESSUN VERDETTO: il piano di fine mescolamento non e'")
            print("stato campionato.")
        else:
            ta = lavaggio[i]
            print(f"\nNESSUN VERDETTO ANCORA. Il piano di fine mescolamento e'")
            print(f"campionato, ma a {float(t.name)*1e6:.1f} us sono passati solo "
                  f"{float(t.name)/ta:.1f} tempi di lavaggio")
            print(f"({ta*1e6:.1f} us l'uno) contro i {LAVAGGI_MINIMI:.0f} "
                  "richiesti. Il numero in tabella e' vero,")
            print(f"ma e' il transitorio, non il regime: servono almeno "
                  f"{LAVAGGI_MINIMI*ta*1e6:.0f} us di calcolo.")
        return 0
    print(f"\nVERDETTO: disuniformita' = {U:.3f} alla fine della lunghezza di")
    print("mescolamento promessa dalla correlazione.")
    if U < 0.10:
        print("  La correlazione REGGE: il combustibile e' distribuito.")
    elif U < 0.25:
        print("  Mescolamento parziale: la correlazione e' ottimista, ma il")
        print("  combustibile e' comunque distribuito prima dello sbocco.")
    else:
        print("  La correlazione NON regge su questa geometria: i getti non")
        print("  penetrano come previsto e l'iniettore va ridisegnato.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
