#!/usr/bin/env python3
"""Mesh del settore d'iniezione per la CFD di mescolamento.

CHE COSA SI VUOLE SAPERE. Tutto l'iniettore poggia su una correlazione:
C = (S/H) sqrt(J) = 2.5 dice che i getti penetrano il condotto e che la miscela
e' uniforme entro un'altezza di condotto a valle. Se non e' vero, il motore e'
un coassiale a taglio travestito. Questa mesh serve a chiederlo a Navier-Stokes
invece che a una formula del 1993.

PERCHE' UN SETTORE E NON TUTTO IL MOTORE. I getti sono `n` identici e
regolarmente spaziati: il piano a meta' strada fra due getti e' un piano di
SIMMETRIA esatto, non un'approssimazione. Si modella quindi un settore di
360/n gradi con il getto sulla bisettrice e simmetria sui due fianchi. Con 4
getti sono 90 gradi, cioe' un quarto del lavoro.

La simmetria e' preferita alla periodicita' ciclica perche' qui e' ESATTA (il
getto sta al centro del settore) e perche' una ciclica mal accoppiata e' uno
dei modi piu' comuni di ottenere un risultato sbagliato che sembra giusto.

DOMINIO. Da monte del condotto d'aria fino a un tratto di camera a valle dello
sbocco: serve vedere sia il mescolamento dentro il condotto (dove la
correlazione vale) sia che cosa ne resta dopo l'allargamento, dove non vale
piu' e dove il motore brucia davvero.
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import gmsh

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))


#: Di quanto il cilindro del getto deve sconfinare oltre le due pareti che
#: attraversa. Serve a evitare facce coincidenti nell'operazione booleana, ed
#: e' il solo motivo per cui esiste.
SCONFINAMENTO_GETTO = 3.0e-4


def tratto_del_getto(r_bore: float, R_getti: float, R_anello: float):
    """(raggio di partenza, lunghezza) del cilindro che scava il foro del getto.

    IL DIFETTO CHE QUESTA FUNZIONE CHIUDE. La prima stesura faceva

        lung = (R_anello - r_bore) + 2 mm

    partendo da r_bore: il cilindro arrivava a r = 6.300 mm, cioe' **2.000 mm
    oltre la parete esterna del condotto**, che sta a 4.300. Il "+2 mm" e' il
    trucco standard per non lasciare facce coincidenti in un'operazione
    booleana, ed e' giusto per un TAGLIO - li' l'eccedenza finisce dentro il
    solido e sparisce. Qui pero' il cilindro viene FUSO con il dominio fluido,
    e l'eccedenza non sparisce: **diventa fluido**. Il risultato era una tasca
    cieca di diametro d_getto e profonda 2 mm scavata nel metallo davanti a
    ogni getto, che nel pezzo vero non esiste.

    Misurata sulla corsa 4 a regime, quella tasca conteneva Y_GPL = 0.0000
    (massimo compreso, su 927 punti): il flusso trasversale piega il getto
    prima che arrivi a quel raggio, e la tasca resta un ricircolo d'aria a
    25 m/s. Non falsava dunque il mescolamento - ma era geometria che il motore
    non ha, e con getti piu' penetranti o inclinati potrebbe non restare
    innocua. Si toglie perche' e' sbagliata, non perche' fa danno.

    La lunghezza giusta e' la sola che serve: **collegare il plenum all'anello**.
    Tutto cio' che sta oltre R_getti e' gia' fluido, quindi la fusione non
    aggiunge niente; basta sconfinare quel tanto che evita le facce coincidenti
    sulla parete attraversata.

    LO SCONFINAMENTO E' DA UNA PARTE SOLA, e la prima correzione sbagliava
    anche questo. Metterlo anche all'estremita' interna (r_bore - 0.3 mm)
    sembrava simmetrico e innocuo, e faceva due danni:

    (1) IL FORO DIVENTAVA PIU' LUNGO DEL VERO. Il condotto reale va dalla
        parete del plenum alla parete dell'anello, cioe' R_getti - r_bore =
        1.000 mm, che su d = 0.731 da' L/d = 1.37 - il valore da cui dipende
        se il getto riattacca dopo la vena contracta, che e' gia' una voce
        aperta del progetto. Con lo sconfinamento interno L/d saliva a 2.19:
        un altro iniettore.

    (2) SPARIVA LA PATCH D'INGRESSO DEL GPL. La base interna del cilindro E'
        la faccia d'ingresso del combustibile, e `_marca_bordi` la riconosce
        perche' sta a raggio r_bore. Spostandola a r_bore - 0.3 mm non veniva
        piu' riconosciuta, finiva nel gruppo "pareti", e la mesh usciva con
        CINQUE bordi invece di sei: il combustibile non entrava da nessuna
        parte. Trovato guardando `constant/polyMesh/boundary` prima di
        lanciare, non dopo.

    All'interno non c'e' niente con cui la faccia possa coincidere - il raggio
    minore di R_getti e' stato tolto dal dominio, li' non c'e' solido da
    intersecare - quindi lo sconfinamento interno non serviva nemmeno allo
    scopo per cui esiste.
    """
    if not (r_bore < R_getti < R_anello):
        raise ValueError(
            f"raggi incoerenti: r_bore {r_bore*1e3:.3f} < R_getti "
            f"{R_getti*1e3:.3f} < R_anello {R_anello*1e3:.3f} mm")
    r0 = r_bore
    r1 = R_getti + SCONFINAMENTO_GETTO
    if r1 >= R_anello:
        raise ValueError(
            f"il cilindro del getto arriverebbe a {r1*1e3:.3f} mm, cioe' oltre "
            f"la parete esterna a {R_anello*1e3:.3f}: tornerebbe a scavare una "
            "tasca cieca nel metallo.")
    return r0, r1 - r0


def costruisci(d, out: Path, fine: float, grossa: float, L_camera: float = 12.0e-3):
    """`d` e' il dizionario di quote prodotto da `quote_dal_progetto`."""
    n = d["n_getti"]
    settore = 2.0 * math.pi / n
    Ri, Re = d["R_getti"], d["R_anello"]
    x_in, x_getti = d["x_ingresso"], d["x_getti"]
    R_c, r_cb, x_rampa = d["R_c"], d["r_centerbody"], d["x_rampa"]
    r_bore, d_getto = d["r_bore"], d["d_getto"]

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("iniettore")
    occ = gmsh.model.occ

    # --- condotto d'aria: corona fra Ri e Re, da x_in a 0 -------------------
    L = -x_in
    esterno = occ.addCylinder(x_in, 0, 0, L, 0, 0, Re)
    interno = occ.addCylinder(x_in - 1e-3, 0, 0, L + 2e-3, 0, 0, Ri)
    #: removeTool=True. Con False il cilindro interno RESTA nel modello, e le
    #: sue due basi diventano due entita' superficiali con baricentro sull'asse:
    #: la classificazione per angolo le prendeva per fianchi del settore e
    #: OpenFOAM si fermava dicendo che il piano di simmetria non e' piano.
    #: Aveva ragione lui.
    condotto, _ = occ.cut([(3, esterno)], [(3, interno)], removeObject=True,
                          removeTool=True)

    # --- camera: da 0 a L_camera, corona fra il corpo centrale e R_c --------
    cam = occ.addCylinder(0, 0, 0, L_camera, 0, 0, R_c)
    #: il corpo centrale si allarga da Ri a r_cb entro x_rampa e poi resta
    #: cilindrico: e' la rampa che riporta la strizione al raggio di camera.
    rampa = occ.addCone(0, 0, 0, x_rampa, 0, 0, Ri, r_cb)
    resto = occ.addCylinder(x_rampa, 0, 0, L_camera - x_rampa, 0, 0, r_cb)
    corpo, _ = occ.fuse([(3, rampa)], [(3, resto)])
    camera, _ = occ.cut([(3, cam)], corpo, removeObject=True, removeTool=True)

    # --- il getto: cilindro radiale sulla BISETTRICE del settore ------------
    #: sulla bisettrice, cosi' i due fianchi del settore sono piani di
    #: simmetria esatti. Parte dentro il corpo centrale (r_bore) e arriva
    #: oltre la parete esterna del condotto.
    ang = 0.5 * settore
    dy, dz = math.cos(ang), math.sin(ang)
    r0, lung = tratto_del_getto(r_bore, Ri, Re)
    #: ATTENZIONE ALL'ORDINE DEGLI ARGOMENTI. addCylinder(x,y,z, dx,dy,dz, r):
    #: l'asse e' (dx,dy,dz). Il getto e' RADIALE, quindi dx = 0 e la direzione
    #: sta tutta nel piano y-z. Alla prima stesura avevo scritto
    #: (lung*dy, lung*dz, 0), cioe' un getto diagonale nel piano x-y: la mesh
    #: si costruiva lo stesso, il caso sarebbe girato lo stesso, e il risultato
    #: sarebbe stato quello di un iniettore che non esiste. L'ha trovato il
    #: confronto fra l'area della faccia d'ingresso e quella analitica.
    getto = occ.addCylinder(x_getti, r0 * dy, r0 * dz,
                            0.0, lung * dy, lung * dz, 0.5 * d_getto)

    fluido, _ = occ.fuse(condotto + camera, [(3, getto)])

    # --- taglio del settore: due semispazi ---------------------------------
    #: un cuneo costruito come prisma, non come box ruotato: il box lascerebbe
    #: spigoli non allineati e OCC produrrebbe facce sottilissime.
    BIG = 4.0 * R_c
    x0, x1 = x_in - 2e-3, L_camera + 2e-3
    p = [occ.addPoint(x0, 0, 0),
         occ.addPoint(x0, BIG, 0),
         occ.addPoint(x0, BIG * math.cos(settore), BIG * math.sin(settore))]
    l = [occ.addLine(p[0], p[1]), occ.addLine(p[1], p[2]), occ.addLine(p[2], p[0])]
    cl = occ.addCurveLoop(l)
    sf = occ.addPlaneSurface([cl])
    cuneo = occ.extrude([(2, sf)], x1 - x0, 0, 0)
    vol_cuneo = [t for t in cuneo if t[0] == 3]
    dominio, _ = occ.intersect(fluido, vol_cuneo)
    occ.synchronize()

    _marca_bordi(d, dominio, settore, x_in, L_camera)
    _dimensioni(d, fine, grossa)

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.optimize("Netgen")
    out.parent.mkdir(parents=True, exist_ok=True)
    gmsh.write(str(out))
    n_el = len(gmsh.model.mesh.getElementsByType(4)[0])
    gmsh.finalize()
    return n_el


def _marca_bordi(d, dominio, settore, x_in, L_camera):
    """Assegna i gruppi fisici guardando dove sta il baricentro di ogni faccia.

    Si marca per POSIZIONE e non per numero d'ordine: l'ordine delle facce che
    OCC produce dopo un'intersezione non e' stabile fra una versione e l'altra,
    e una condizione al contorno finita sulla faccia sbagliata e' un errore che
    non si vede finche' non si guardano i risultati.
    """
    #: Si classificano SOLO le facce che bordano davvero il dominio finale.
    #: Iterare su tutte le entita' del modello raccoglie anche i residui delle
    #: operazioni booleane: esistono ancora nel kernel ma non sono superfici del
    #: fluido, e finiscono in un gruppo a caso. E' successo con le due basi del
    #: cilindro usato per scavare la corona.
    facce_dominio = {abs(t) for _, t in gmsh.model.getBoundary(dominio,
                                                               oriented=False)}
    #: DUE gruppi di simmetria, non uno. OpenFOAM accetta `symmetryPlane` solo
    #: su una patch PIANA: mettendo i due fianchi del settore nello stesso
    #: gruppo, ne calcola la normale media - che non e' la normale di nessuno
    #: dei due - e si ferma. E' un errore che si vede subito perche' il
    #: solutore non parte; l'alternativa (`symmetry`, che accetta patch non
    #: piane) partirebbe e riflettereb-be sul piano sbagliato.
    gruppi = {"ingresso_aria": [], "ingresso_gpl": [], "uscita": [],
              "simmetria_1": [], "simmetria_2": [], "pareti": []}
    for dim, tag in gmsh.model.getEntities(2):
        if tag not in facce_dominio:
            continue
        cm = gmsh.model.occ.getCenterOfMass(dim, tag)
        x, y, z = cm
        n = gmsh.model.getNormal(tag, [0.5, 0.5])
        r = math.hypot(y, z)
        ang = math.atan2(z, y)
        if abs(x - x_in) < 1e-5:
            gruppi["ingresso_aria"].append(tag)
        elif abs(x - L_camera) < 1e-5:
            gruppi["uscita"].append(tag)
        elif abs(r - d["r_bore"]) < 5e-5 and abs(x - d["x_getti"]) < 2e-3:
            gruppi["ingresso_gpl"].append(tag)
        elif abs(ang) < 1e-4:
            gruppi["simmetria_1"].append(tag)
        elif abs(ang - settore) < 1e-4:
            gruppi["simmetria_2"].append(tag)
        else:
            gruppi["pareti"].append(tag)
    #: UN GRUPPO VUOTO E' UN ERRORE, NON UN CASO DA SALTARE. Il `if tag:` da
    #: solo faceva questo: una faccia non riconosciuta finiva in "pareti", il
    #: gruppo che le avrebbe dovute contenere restava vuoto, e la mesh usciva
    #: con cinque bordi invece di sei senza un rigo di avviso. E' successo
    #: davvero: spostando la base del getto di 0.3 mm l'ingresso del GPL non
    #: e' piu' stato riconosciuto, ed e' diventato parete. Il combustibile non
    #: entrava piu' nel motore, e la mesh risultava "OK" a checkMesh - perche'
    #: geometricamente lo era. Se ne sarebbe accorto solo il solutore,
    #: fermandosi su una condizione al contorno riferita a una patch assente.
    vuoti = [n for n, t in gruppi.items() if not t]
    if vuoti:
        raise RuntimeError(
            f"gruppi di bordo rimasti vuoti: {', '.join(vuoti)}. Le facce ci "
            "sono ma la classificazione per posizione non le riconosce: sono "
            f"finite in 'pareti' ({len(gruppi['pareti'])} facce). Controllare "
            "che i raggi e le quote usati per classificare siano ancora quelli "
            "con cui la geometria viene costruita.")
    for nome, tag in gruppi.items():
        g = gmsh.model.addPhysicalGroup(2, tag)
        gmsh.model.setPhysicalName(2, g, nome)
    vol = [t for _, t in gmsh.model.getEntities(3)]
    g = gmsh.model.addPhysicalGroup(3, vol)
    gmsh.model.setPhysicalName(3, g, "fluido")
    return gruppi


def _dimensioni(d, fine: float, grossa: float):
    """Raffinamento: fine attorno al getto e nel condotto, grossa in camera."""
    campo_getto = gmsh.model.mesh.field.add("Ball")
    gmsh.model.mesh.field.setNumber(campo_getto, "Radius", 6.0 * d["d_getto"])
    gmsh.model.mesh.field.setNumber(campo_getto, "Thickness", 4.0 * d["d_getto"])
    ang = math.pi / d["n_getti"]
    rm = 0.5 * (d["R_getti"] + d["R_anello"])
    gmsh.model.mesh.field.setNumber(campo_getto, "XCenter", d["x_getti"])
    gmsh.model.mesh.field.setNumber(campo_getto, "YCenter", rm * math.cos(ang))
    gmsh.model.mesh.field.setNumber(campo_getto, "ZCenter", rm * math.sin(ang))
    gmsh.model.mesh.field.setNumber(campo_getto, "VIn", fine)
    gmsh.model.mesh.field.setNumber(campo_getto, "VOut", grossa)

    campo_condotto = gmsh.model.mesh.field.add("Box")
    gmsh.model.mesh.field.setNumber(campo_condotto, "XMin", d["x_ingresso"])
    gmsh.model.mesh.field.setNumber(campo_condotto, "XMax", 2.0e-3)
    gmsh.model.mesh.field.setNumber(campo_condotto, "YMin", -1.0)
    gmsh.model.mesh.field.setNumber(campo_condotto, "YMax", 1.0)
    gmsh.model.mesh.field.setNumber(campo_condotto, "ZMin", -1.0)
    gmsh.model.mesh.field.setNumber(campo_condotto, "ZMax", 1.0)
    gmsh.model.mesh.field.setNumber(campo_condotto, "VIn", 2.0 * fine)
    gmsh.model.mesh.field.setNumber(campo_condotto, "VOut", grossa)
    gmsh.model.mesh.field.setNumber(campo_condotto, "Thickness", 2.0e-3)

    minimo = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(minimo, "FieldsList",
                                     [campo_getto, campo_condotto])
    gmsh.model.mesh.field.setAsBackgroundMesh(minimo)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)      # HXT, veloce e robusto


def _senza_ripetizioni(xs, tolleranza: float = 1.0e-9):
    """Le stazioni di misura, in ordine e senza doppioni entro `tolleranza`."""
    fuori = []
    for x in xs:
        if not any(abs(x - y) <= tolleranza for y in fuori):
            fuori.append(x)
    return fuori


def quote_dal_progetto(spinta=50.0, durata=5.0, T_bombola=288.15):
    """Le quote NON si scrivono qui: si chiedono al sintetizzatore.

    E' il punto di tutta l'infrastruttura: la CFD deve verificare il motore che
    il modello produce, non una sua copia scritta a mano che diverge al primo
    cambiamento.
    """
    from zefiro.sintesi import Impianto, Requisiti, progetta

    p = progetta(Requisiti(
        spinta=spinta, durata=durata,
        impianto=Impianto(combustibile={"C3H8": 1.0}, T_bombola_min=T_bombola,
                          T_ossidante_iniezione=293.0, T_combustibile_iniezione=283.0,
                          cd_ossidante=0.75, cd_combustibile=0.75)), geometria=False)
    dim, t, inj = p.geometria, p.geometria.testa, p.geometria.iniettore
    dd = dim.params.derived
    return {
        "n_getti": t.n_getti, "R_getti": t.R_getti, "R_anello": t.R_anello,
        "d_getto": t.d_getto, "x_getti": t.x_getti, "r_bore": t.r_bore_gpl,
        "x_ingresso": t.x_condotto, "R_c": dd["R_c"], "r_centerbody": dd["r_centerbody"],
        "x_rampa": t.x_rampa,
        # condizioni al contorno
        "V_aria": inj.V_aria, "rho_aria": inj.rho_aria, "T_aria": 293.0,
        "V_gpl": inj.V_getto, "rho_gpl": inj.rho_getto, "T_gpl": 283.0,
        "p_c": dim.registro.valore("p_camera"),
        "mdot_aria": dim.l0.mdot_air, "mdot_gpl": dim.l0.mdot_fuel_core,
        "H": inj.altezza_anello, "C": inj.C, "J": inj.J,
        "L_mescolamento": inj.lunghezza_mescolamento,
        #: Piani su cui misurare la disuniformita': dai getti fino a valle
        #: dello sbocco, passando per la lunghezza di mescolamento promessa.
        #: Se la correlazione ha ragione, la disuniformita' deve essere gia'
        #: piccola al piano che sta a x_getti + L_mescolamento.
        #: I PIANI SI DEDUPLICANO. `x_getti` vale esattamente
        #: -lunghezza_mescolamento, quindi il piano a f = 1.0 CADE su x = 0 ed
        #: era lo stesso di quello messo dopo in coda: due `sez_XX` con le
        #: stesse coordinate, campionati due volte, stampati come se fossero
        #: due misure indipendenti. Non falsava nessun numero - erano identici -
        #: ma invitava a leggere una conferma dove c'era una copia.
        "piani_x": _senza_ripetizioni(
            [t.x_getti + f * inj.lunghezza_mescolamento
             for f in (0.25, 0.5, 1.0, 1.5)] + [0.0, 3.0e-3, 8.0e-3]),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fine", type=float, default=6.0e-5, help="cella al getto [m]")
    ap.add_argument("--grossa", type=float, default=4.0e-4)
    ap.add_argument("--out", type=Path, default=Path("runs/cfd/iniettore.msh"))
    a = ap.parse_args()
    d = quote_dal_progetto()
    print(f"getti {d['n_getti']} x {d['d_getto']*1e3:.3f} mm  "
          f"condotto {d['R_getti']*1e3:.2f}-{d['R_anello']*1e3:.2f} mm "
          f"(H {d['H']*1e3:.3f})  C {d['C']:.3f}  J {d['J']:.3f}")
    print(f"aria {d['V_aria']:.1f} m/s   GPL {d['V_gpl']:.1f} m/s   "
          f"lunghezza di mescolamento attesa {d['L_mescolamento']*1e3:.2f} mm")
    n = costruisci(d, a.out, a.fine, a.grossa)
    print(f"mesh scritta in {a.out}: {n} tetraedri")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
