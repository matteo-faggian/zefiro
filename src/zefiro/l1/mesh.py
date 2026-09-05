"""Mesh del volume FLUIDO come settore periodico, con Gmsh (kernel OCC).

Tre cose sono delicate e vengono verificate invece che sperate:

1. **La periodicita' deve essere CONFORME.** La condizione `cyclic` di
   OpenFOAM accoppia le facce a coppie: se le due facce del settore non hanno
   la stessa mesh, l'accoppiamento e' interpolato, e su un flusso reagente
   supersonico l'interpolazione a cavallo del piano periodico introduce una
   sorgente numerica di massa. Si usa `gmsh.model.mesh.setPeriodic`, e poi si
   VERIFICA nodo per nodo (`verify_periodicity`), perche' l'API accetta
   volentieri una trasformazione sbagliata senza protestare.

2. **I nomi delle frontiere.** Vengono dal poligono meridiano
   (`geometry.profile.fluid_polygon`), dove nascono in codice puro e testabile,
   e vengono propagati alle superfici di rivoluzione. Non si indovinano dai
   baricentri: una condizione al contorno sulla faccia sbagliata non fa
   fallire niente e produce solo un risultato falso.

3. **Gli iniettori sono IMPRINTATI, non estrusi.** I fori si tagliano sulla
   faccia d'iniezione come superfici separate, senza modellare il condotto a
   monte. Non e' una scorciatoia: a monte del foro la condizione e' comunque
   una portata imposta, quindi modellare il condotto aggiungerebbe celle senza
   aggiungere fisica. Cio' che conta - il getto discreto e la sua penetrazione
   nella corrente d'aria - resta risolto.
"""
from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass
from pathlib import Path

from zefiro.geometry.profile import fluid_polygon
from zefiro.schemas import (
    CANONICAL_BOUNDARIES,
    GeometryParams,
    MeshArtifact,
    ZefiroError,
)

#: Frontiere che DEVONO esserci. Le altre di `CANONICAL_BOUNDARIES` possono
#: mancare legittimamente: `inlet_fuel_film` non esiste se f_film = 0, e
#: `wall_throat`/`wall_cowl` non esistono su un plug a espansione esterna con
#: labbro di spessore nullo (la gola e' delimitata da uno spigolo, non da una
#: superficie). Vedi la docstring di `fluid_polygon`.
REQUIRED_BOUNDARIES: frozenset[str] = frozenset({
    "inlet_air", "inlet_fuel_core", "wall_chamber", "wall_convergent",
    "wall_plug", "outlet_far", "periodic_a", "periodic_b",
})

#: Tolleranza relativa (sul raggio del labbro) per considerare corrispondenti
#: due nodi dopo la rotazione periodica.
#:
#: Il valore viene da cosa serve a valle, non da un'idea di precisione: la
#: `matchTolerance` di default delle patch cyclic di OpenFOAM e' 1e-4
#: relativa. Qui si tiene dieci volte piu' stretto. Piu' stretto ancora non
#: avrebbe senso: OCC valuta le spline con una tolleranza sua, e sul contorno
#: del plug lascia scarti di qualche decina di nanometri (3.6e-6 relativi)
#: che non sono un difetto di periodicita' ma il rumore del kernel geometrico.
PERIODIC_TOL_REL = 1.0e-5

#: Soglie di qualita', quelle su cui decide `checkMesh` di OpenFOAM.
#:
#: Sono trattate diversamente perche' contano diversamente. La skewness oltre 4
#: e' un ERRORE: i flussi convettivi diventano inconsistenti e non c'e' niente
#: nel solutore che lo compensi. La non-ortogonalita' oltre 70 gradi e' un
#: AVVISO: OpenFOAM la corregge con i `nNonOrthogonalCorrectors`, al prezzo di
#: qualche iterazione in piu'. Fermare la generazione per la seconda sarebbe
#: rifiutare mesh perfettamente utilizzabili.
MAX_SKEWNESS = 4.0
MAX_NON_ORTHOGONALITY = 70.0

#: Setto minimo fra due fori adiacenti, in frazione del diametro maggiore.
#: Non e' un vincolo di mesh ma di FABBRICAZIONE: sotto questa soglia il ponte
#: di materiale fra due fori in SLM non regge, ed e' lo stesso ragionamento
#: che governa `film_land` in geometry/parameters.
MIN_HOLE_GAP_FRACTION = 0.25

#: Segmenti minimi lungo la CIRCONFERENZA di un foro d'iniezione.
#:
#: Un cerchio meshato con n segmenti diventa un poligono inscritto, la cui area
#: e' (n / 2 pi) sin(2 pi / n) volte quella del cerchio: con 7 segmenti manca
#: il 13 %. E l'area del foro non e' un dettaglio estetico - e' cio' che, a
#: portata imposta, fissa la VELOCITA' di iniezione, quindi il rapporto delle
#: quantita' di moto, quindi la miscelazione: esattamente la grandezza per cui
#: questa mesh esiste. Con 16 segmenti l'errore d'area scende all'1.3 %.
#: Misurato, non stimato: il test sull'area delle patch lo verifica.
MIN_SEGMENTS_PER_HOLE = 16


class MeshError(ZefiroError):
    """La mesh non soddisfa una precondizione dichiarata."""


@dataclass(frozen=True)
class MeshOptions:
    #: Tetto alla dimensione di cella sui fori. La dimensione EFFETTIVA di
    #: ciascun foro e' il minimo fra questo e pi*d/MIN_SEGMENTS_PER_HOLE,
    #: cioe' ogni foro e' risolto in proporzione al proprio diametro.
    injector_size: float = 1.5e-4      # m
    wall_size: float = 4.0e-4          # m, sulle pareti calde
    farfield_size: float = 8.0e-3      # m, lontano da tutto
    #: Entro questa distanza da una parete la cella resta della dimensione di
    #: parete; oltre, cresce linearmente fino a `farfield_size` a
    #: `growth_distance`. Due numeri e non uno perche' servono due cose
    #: diverse: uno strato di celle fini che risolva il gradiente, e una
    #: crescita abbastanza dolce da non creare celle storte.
    refine_distance: float = 1.5e-3    # m, alone fine attorno alle pareti
    #: Alone fine attorno a un iniettore, in DIAMETRI del foro: e' il campo
    #: vicino del getto, e la sua scala e' il foro, non un millimetro assoluto.
    injector_halo_diameters: float = 2.0
    growth_distance: float = 6.0e-2    # m
    order: int = 1
    algorithm3d: int = 1               # 1 = Delaunay, robusto su domini sottili
    #: Passate dell'ottimizzatore Netgen. ZERO di default, ed e' una scelta
    #: MISURATA: una passata alza la skewness da 1.43 a 2.15 senza migliorare
    #: la non-ortogonalita'. Le due strade istintive per le schegge residue -
    #: raffinare vicino all'asse e ottimizzare - peggiorano entrambe la mesh su
    #: questa geometria. La manopola resta perche' su una geometria diversa
    #: puo' servire, ma il default e' quello che i numeri dicono.
    optimize_passes: int = 0
    write_msh: bool = True


def generate_mesh(
    params: GeometryParams,
    out_dir: Path,
    run_id: str,
    opts: MeshOptions = MeshOptions(),
) -> MeshArtifact:
    """Costruisce e mesha il settore periodico 2 pi / N del dominio fluido."""
    import gmsh

    d = params.derived
    n_inj = int(round(d["N_inj"]))
    sector = 2.0 * math.pi / n_inj
    pts2d, nomi, curve = fluid_polygon(d, params.plug_contour_x,
                                       params.plug_contour_r)

    out_dir.mkdir(parents=True, exist_ok=True)
    msh_path = out_dir / f"{run_id}_fluid.msh"

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add(run_id)
        occ = gmsh.model.occ

        # --- profilo meridiano nel piano z = 0 ---------------------------- #
        vtx = [occ.addPoint(x, r, 0.0) for x, r in pts2d]
        n = len(vtx)

        # I tratti dichiarati lisci da `fluid_polygon` diventano UNA spline.
        # Il contorno del plug ha ~135 punti su una decina di millimetri: come
        # spezzata, ogni punto sarebbe un nodo obbligato e forzerebbe celle da
        # 0.07 mm accanto a celle da 0.8 mm. Misurato prima della correzione:
        # non-ortogonalita' 89 gradi e skewness 4.8, cioe' mesh da buttare.
        linee: list[int] = []
        nomi_curva: list[str] = []
        vertici_curva: list[list[int]] = []      # quali vertici copre ogni curva
        i = 0
        while i < n:
            corrente = [(a, b) for a, b in curve if a == i]
            if corrente:
                a, b = corrente[0]
                linee.append(occ.addSpline(vtx[a:b + 1]))
                nomi_curva.append(nomi[a])
                vertici_curva.append(list(range(a, b + 1)))
                i = b
            else:
                linee.append(occ.addLine(vtx[i], vtx[(i + 1) % n]))
                nomi_curva.append(nomi[i])
                vertici_curva.append([i, (i + 1) % n])
                i += 1
        anello = occ.addCurveLoop(linee)
        faccia = occ.addPlaneSurface([anello])

        # --- rivoluzione del settore attorno all'asse x -------------------- #
        # Si ruota di -sector/2 e si estrude di +sector, cosi' il settore e'
        # simmetrico rispetto al piano z = 0: le due facce periodiche stanno a
        # +-sector/2 e la rotazione che le lega e' esattamente `sector`.
        occ.rotate([(2, faccia)], 0, 0, 0, 1, 0, 0, -0.5 * sector)
        rivoluzione = occ.revolve([(2, faccia)], 0, 0, 0, 1, 0, 0, sector)
        occ.synchronize()

        volumi = [t for dim, t in rivoluzione if dim == 3]
        if len(volumi) != 1:
            raise MeshError(f"la rivoluzione ha prodotto {len(volumi)} volumi, ne serve 1")
        vol = volumi[0]

        # Corrispondenza lato -> superficie, per TOPOLOGIA e non per ordine.
        #
        # La tentazione e' di fidarsi dell'ordine con cui `revolve` restituisce
        # le superfici laterali. Non regge: un lato che GIACE sull'asse (il
        # tratto `axis`, r = 0 a entrambi gli estremi) ruotando non genera
        # nessuna superficie, e l'ordine scala di uno senza avvisare. Il
        # risultato sarebbe ogni condizione al contorno spostata di una
        # superficie: silenzioso e disastroso.
        #
        # Invece: ogni superficie laterale ha fra le proprie curve di bordo
        # ESATTAMENTE il lato meridiano che l'ha generata, il cui tag non e'
        # cambiato. Si interroga la topologia e non si indovina niente.
        tutte_2d = [t for dim, t in rivoluzione if dim == 2]
        faccia_a = faccia          # a -sector/2, quella di partenza
        faccia_b = tutte_2d[0]     # a +sector/2, quella di arrivo
        laterali = tutte_2d[1:]

        da_linea = dict(zip(linee, nomi_curva))
        per_nome: dict[str, list[int]] = {}
        assegnate = 0
        for sup in laterali:
            bordo = {abs(t) for _, t in gmsh.model.getBoundary([(2, sup)],
                                                               oriented=False)}
            candidati = sorted(bordo & set(da_linea))
            if len(candidati) != 1:
                raise MeshError(
                    f"la superficie {sup} tocca {len(candidati)} lati meridiani "
                    f"({candidati}): la corrispondenza non e' univoca"
                )
            per_nome.setdefault(da_linea[candidati[0]], []).append(sup)
            assegnate += 1

        # Solo i lati sull'asse possono non aver generato superficie.
        # Quali curve giacciono sull'asse si decide dal POLIGONO, dove r vale
        # esattamente 0.0, e non dal bounding box di OCC, che porta con se' una
        # tolleranza di 1e-7 m: con un raggio di base di 0.76 mm quella
        # tolleranza non e' innocua, e comunque una soglia arbitraria qui non
        # serve a niente quando il dato esatto ce l'abbiamo gia'.
        senza_superficie = [nome for k, nome in enumerate(nomi_curva)
                            if all(pts2d[j][1] == 0.0 for j in vertici_curva[k])]
        attese = len(linee) - len(senza_superficie)
        if assegnate != attese:
            raise MeshError(
                f"{assegnate} superfici laterali per {attese} curve non degeneri "
                f"(su {len(linee)} totali, {len(senza_superficie)} sull'asse)"
            )
        if set(senza_superficie) - {"axis"}:
            raise MeshError(
                f"lati degeneri non sull'asse: {sorted(set(senza_superficie))}"
            )

        per_nome["periodic_a"] = [faccia_a]
        per_nome["periodic_b"] = [faccia_b]

        # --- fori d'iniezione imprintati sulla faccia --------------------- #
        per_nome = _imprint_injectors(gmsh, per_nome, d, sector)
        occ.synchronize()

        # --- physical groups ---------------------------------------------- #
        for nome, tags in sorted(per_nome.items()):
            if nome not in CANONICAL_BOUNDARIES:
                raise MeshError(
                    f"frontiera {nome!r} non e' in CANONICAL_BOUNDARIES "
                    f"{CANONICAL_BOUNDARIES}: un refuso qui diventa una "
                    "condizione al contorno sbagliata in silenzio"
                )
            gmsh.model.addPhysicalGroup(2, tags, name=nome)
        gmsh.model.addPhysicalGroup(3, [vol], name="fluid")

        mancanti = REQUIRED_BOUNDARIES - set(per_nome)
        if mancanti:
            raise MeshError(f"frontiere obbligatorie mancanti: {sorted(mancanti)}")

        # --- periodicita' conforme ----------------------------------------- #
        # matrice 4x4, riga per riga: rotazione di `sector` attorno a x
        c, s = math.cos(sector), math.sin(sector)
        affine = [1, 0, 0, 0,
                  0, c, -s, 0,
                  0, s, c, 0,
                  0, 0, 0, 1]
        gmsh.model.mesh.setPeriodic(2, [faccia_b], [faccia_a], affine)

        # Ogni foro ha la SUA dimensione di cella, proporzionale al SUO
        # diametro. Usarne una sola, presa dal foro piu' piccolo, era un
        # errore costoso: su un punto reale il foro del film e' 0.095 mm e
        # imponeva celle da 19 um anche attorno al foro d'aria da 1.8 mm, che
        # di quella finezza non ha nessun bisogno. La mesh non finiva piu'.
        diametri = {"inlet_air": d["d_ox"], "inlet_fuel_core": d["d_fuel"]}
        if d.get("d_film", 0.0) > 0.0:
            diametri["inlet_fuel_film"] = d["d_film"]
        _size_fields(gmsh, per_nome, opts, diametri)
        gmsh.option.setNumber("Mesh.Algorithm3D", opts.algorithm3d)
        gmsh.option.setNumber("Mesh.ElementOrder", opts.order)
        gmsh.model.mesh.generate(3)
        for _ in range(opts.optimize_passes):
            gmsh.model.mesh.optimize("Netgen")

        n_celle = len(gmsh.model.mesh.getElementsByType(4)[0])
        if n_celle == 0:
            raise MeshError("nessun tetraedro generato")

        periodica, scarto = verify_periodicity(gmsh, faccia_a, faccia_b, sector,
                                               d["R_lip"])
        non_orto, skew, frazione = mesh_quality(gmsh)

        avvisi: list[str] = []
        if not periodica:
            raise MeshError(
                f"le facce periodiche non corrispondono: scarto massimo "
                f"{scarto:.2e} volte il raggio del labbro, contro una tolleranza "
                f"di {PERIODIC_TOL_REL:.0e}. Con patch cyclic non conformi "
                "OpenFOAM interpola attraverso il piano periodico e introduce "
                "una sorgente numerica di massa."
            )
        if skew > MAX_SKEWNESS:
            raise MeshError(
                f"skewness massima {skew:.2f} oltre il limite {MAX_SKEWNESS}: "
                "i flussi convettivi diventano inconsistenti e il solutore non "
                "ha modo di compensare. Non e' una mesh utilizzabile."
            )
        if non_orto > MAX_NON_ORTHOGONALITY:
            avvisi.append(
                f"non-ortogonalita' massima {non_orto:.1f} gradi oltre "
                f"{MAX_NON_ORTHOGONALITY:.0f}, su {frazione:.3%} delle facce "
                "interne. Utilizzabile, ma il caso OpenFOAM deve usare almeno "
                "un nNonOrthogonalCorrector."
            )

        if opts.write_msh:
            gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)   # letta da OpenFOAM
            gmsh.write(str(msh_path))
    finally:
        gmsh.finalize()

    sha = hashlib.sha256(msh_path.read_bytes()).hexdigest() if msh_path.exists() else ""
    return MeshArtifact(
        run_id=run_id,
        msh_path=msh_path,
        n_cells=int(n_celle),
        boundary_names=tuple(sorted(per_nome)),
        is_periodic=bool(periodica),
        sector_angle=sector,
        max_non_orthogonality=float(non_orto),
        max_skewness=float(skew),
        sha256=sha,
        frac_non_orthogonal=float(frazione),
        periodic_mismatch=float(scarto),
        warnings=tuple(avvisi),
    )


def _imprint_injectors(gmsh, per_nome, d, sector):
    """Taglia i fori sulla faccia d'iniezione come superfici separate.

    Un foro per settore per famiglia, che e' esattamente la periodicita' del
    motore: N fori su 2 pi, quindi uno su 2 pi / N. I fori stanno al centro del
    settore (z = 0 dopo la rotazione), che e' l'unica posizione compatibile con
    la periodicita': un foro fuori centro verrebbe tagliato in due dai piani
    periodici.
    """
    occ = gmsh.model.occ
    facce = per_nome.pop("wall_faceplate")
    if len(facce) != 1:
        raise MeshError(f"{len(facce)} facce d'iniezione, ne serve 1")
    faccia = facce[0]

    R_inj = d["R_inj"]
    larghezza = R_inj * sector
    fori: list[tuple[str, float, float]] = [
        ("inlet_air", d["d_ox"], R_inj),
        ("inlet_fuel_core", d["d_fuel"], R_inj),
    ]
    if d.get("d_film", 0.0) > 0.0:
        fori.append(("inlet_fuel_film", d["d_film"], d["R_film"]))

    # --- disposizione angolare dei fori ------------------------------------ #
    # Il primo tentativo distribuiva i centri a passo uniforme senza guardare i
    # DIAMETRI. Su un punto reale il foro dell'aria e' 4.2 mm e il passo era
    # 1.5 mm: i fori si sovrapponevano, il `fragment` li fondeva, e le patch
    # uscivano con l'area sbagliata. La portata imposta avrebbe dato la
    # velocita' sbagliata, e nessun residuo se ne sarebbe accorto.
    # Il test sull'AREA delle patch e' cio' che lo ha trovato.
    somma = sum(dia for _, dia, _ in fori)
    gioco = (larghezza - somma) / (len(fori) + 1)
    if gioco <= MIN_HOLE_GAP_FRACTION * max(dia for _, dia, _ in fori):
        raise MeshError(
            f"i fori non ci stanno nel settore: {somma*1e3:.2f} mm di diametri "
            f"su {larghezza*1e3:.2f} mm di arco a R_inj lasciano {gioco*1e3:.2f} mm "
            f"di setto, sotto il minimo richiesto. Alza N_inj o riduci d_ox_ratio."
        )

    posti = []
    percorso = -0.5 * larghezza
    for nome, dia, raggio in fori:
        percorso += gioco + 0.5 * dia
        ang = percorso / raggio
        posti.append((nome, dia, raggio * math.cos(ang), raggio * math.sin(ang), ang))
        percorso += 0.5 * dia

    # verifica esplicita di non sovrapposizione, in 3D e non per costruzione
    for i in range(len(posti)):
        for j in range(i + 1, len(posti)):
            (ni, di, yi, zi, _), (nj, dj, yj, zj, _) = posti[i], posti[j]
            distanza = math.hypot(yi - yj, zi - zj)
            if distanza < 0.5 * (di + dj):
                raise MeshError(
                    f"i fori {ni} e {nj} si sovrappongono: interasse "
                    f"{distanza*1e3:.3f} mm contro {(di+dj)/2*1e3:.3f} mm richiesti"
                )

    # e che stiano dentro l'anello di camera, altrimenti il disco sporge dalla
    # faccia e il fragment lascia un pezzo che non e' un foro
    for nome, dia, raggio in fori:
        if raggio - 0.5 * dia < d["r_centerbody"] or raggio + 0.5 * dia > d["R_c"]:
            raise MeshError(
                f"il foro {nome} (d = {dia*1e3:.2f} mm a r = {raggio*1e3:.2f} mm) "
                f"sporge dall'anello di camera "
                f"[{d['r_centerbody']*1e3:.2f}, {d['R_c']*1e3:.2f}] mm"
            )

    dischi = []
    for nome, dia, y, z, _ in posti:
        cerchio = occ.addCircle(0.0, y, z, 0.5 * dia, zAxis=[1.0, 0.0, 0.0])
        loop = occ.addCurveLoop([cerchio])
        dischi.append((nome, occ.addPlaneSurface([loop])))

    frammenti, mappa = occ.fragment(
        [(2, faccia)], [(2, t) for _, t in dischi]
    )
    occ.synchronize()

    # `mappa` e' allineata agli input: mappa[0] sono i pezzi della faccia,
    # mappa[1:] i pezzi di ciascun disco. Un disco che sta dentro la faccia
    # produce un pezzo solo, ed e' quello il foro.
    per_nome["wall_faceplate"] = []
    fori_tags: dict[str, list[int]] = {}
    for (nome, _), pezzi in zip(dischi, mappa[1:]):
        fori_tags[nome] = [t for dim, t in pezzi if dim == 2]
    tutti_fori = {t for v in fori_tags.values() for t in v}
    resto = [t for dim, t in mappa[0] if dim == 2 and t not in tutti_fori]

    for nome, tags in fori_tags.items():
        if not tags:
            raise MeshError(f"il foro {nome} non ha lasciato traccia sulla faccia")
        per_nome[nome] = tags
    per_nome["wall_faceplate"] = resto
    return per_nome


def _size_fields(gmsh, per_nome, opts: MeshOptions,
                 diametri: dict[str, float]) -> None:
    """Dimensione delle celle: fine dove la fisica e' fine, grossa altrove.

    Le zone fini non sono scelte per gusto. Il foro del GPL e' il piu' piccolo
    (frazioni di millimetro) ed e' proprio la scala che decide la
    miscelazione, cioe' la grandezza per cui questa mesh esiste. Le pareti
    calde servono al flusso termico. Il campo lontano non serve a niente:
    li' le celle possono essere grosse.
    """
    f = gmsh.model.mesh.field
    campi = []

    def crescita(tags, dimensione, alone):
        """Cella `dimensione` fino a distanza `alone` dalla superficie, poi
        crescita lineare fino a `farfield_size` a `growth_distance`."""
        if not tags:
            return
        d_id = f.add("Distance")
        f.setNumbers(d_id, "SurfacesList", list(tags))
        t_id = f.add("Threshold")
        f.setNumber(t_id, "InField", d_id)
        f.setNumber(t_id, "SizeMin", dimensione)
        f.setNumber(t_id, "SizeMax", opts.farfield_size)
        f.setNumber(t_id, "DistMin", alone)
        f.setNumber(t_id, "DistMax", opts.growth_distance)
        campi.append(t_id)

    # NOTA su un secondo tentativo FALLITO, tenuta perche' e' istruttiva.
    # La prima versione metteva un campo di dimensione anche su `outlet_far`,
    # con l'idea di "infittire verso l'uscita". Fa l'opposto di quel che serve:
    # la frontiera del blocco di campo lontano E' outlet_far, quindi ogni punto
    # del campo lontano le e' vicino, e tutto il blocco veniva meshato fine.
    # Risultato: 130 000 celle di cui la gran parte sprecate nel pennacchio,
    # dove non c'e' niente da risolvere. Il campo lontano non va infittito:
    # va lasciato crescere partendo dalle pareti.
    # Un campo per foro, dimensionato sul DIAMETRO di quel foro. Anche l'alone
    # fine si misura in diametri e non in millimetri assoluti: la scala del
    # getto e' il foro che lo genera.
    #
    # Un solo campo per tutti gli iniettori, preso dal foro piu' piccolo, era
    # un errore costoso: su un punto reale il foro del film e' 0.095 mm e
    # imponeva celle da 19 um anche attorno al foro d'aria da 1.8 mm, che di
    # quella finezza non ha nessun bisogno. La mesh non finiva piu'.
    for nome, dia in diametri.items():
        tags = per_nome.get(nome, [])
        if tags:
            crescita(tags,
                     min(opts.injector_size, math.pi * dia / MIN_SEGMENTS_PER_HOLE),
                     opts.injector_halo_diameters * dia)
    crescita([t for n in ("wall_chamber", "wall_convergent", "wall_plug",
                          "wall_faceplate")
              for t in per_nome.get(n, [])], opts.wall_size, opts.refine_distance)

    minimo = f.add("Min")
    f.setNumbers(minimo, "FieldsList", campi)
    f.setAsBackgroundMesh(minimo)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)


# --- verifiche indipendenti -------------------------------------------------- #

def verify_periodicity(gmsh, faccia_a: int, faccia_b: int, sector: float,
                       scala: float) -> tuple[bool, float]:
    """Controlla NODO PER NODO che le due facce del settore corrispondano.

    Perche' non ci si fida di `setPeriodic`: l'API accetta una matrice affine
    sbagliata senza protestare, e produce una mesh che sembra a posto. Il
    difetto si manifesta solo dentro OpenFOAM, come uno squilibrio di massa a
    cavallo del piano periodico, e a quel punto e' difficilissimo da attribuire.

    Restituisce `(conforme, scarto_massimo_relativo)`. Lo scarto e'
    adimensionalizzato sul raggio del labbro: un numero in metri non direbbe
    se e' tanto o poco.
    """
    import numpy as np

    def nodi(sup: int) -> np.ndarray:
        _, coord, _ = gmsh.model.mesh.getNodes(2, sup, includeBoundary=True)
        return np.array(coord, dtype=float).reshape(-1, 3)

    A, B = nodi(faccia_a), nodi(faccia_b)
    if len(A) == 0 or len(B) == 0:
        return False, float("inf")
    if len(A) != len(B):
        return False, float("inf")

    c, s = math.cos(sector), math.sin(sector)
    R = np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])
    Ar = A @ R.T

    # per ogni nodo ruotato di A, il piu' vicino di B
    d2 = ((Ar[:, None, :] - B[None, :, :]) ** 2).sum(axis=2)
    scarto = float(np.sqrt(d2.min(axis=1)).max()) / scala
    return bool(scarto < PERIODIC_TOL_REL), scarto


def mesh_quality(gmsh) -> tuple[float, float, float]:
    """Non-ortogonalita' massima [gradi], skewness massima, e la FRAZIONE di
    facce interne oltre i 70 gradi. Calcolate qui, non da checkMesh.

    Sono le due metriche su cui `checkMesh` di OpenFOAM decide, ed e' il motivo
    per cui vengono calcolate adesso e non dopo: scoprire che la mesh non va
    bene DOPO aver scritto il caso, decomposto e lanciato il solutore costa
    un'ora; scoprirlo qui costa niente.

    Definizioni (le stesse di OpenFOAM):
      * non-ortogonalita': angolo fra la congiungente i baricentri di due celle
        adiacenti e la normale alla faccia che le separa;
      * skewness: distanza fra il centro della faccia e il punto in cui la
        congiungente la attraversa, rapportata alla congiungente stessa.
    """
    import numpy as np

    tipi, tag_el, nodi_el = gmsh.model.mesh.getElements(3)
    if 4 not in tipi:
        return float("inf"), float("inf"), 1.0
    i4 = list(tipi).index(4)
    conn = np.array(nodi_el[i4], dtype=np.int64).reshape(-1, 4)

    tag_nodi, coord, _ = gmsh.model.mesh.getNodes()
    xyz = np.array(coord, dtype=float).reshape(-1, 3)
    indice = np.zeros(int(np.max(tag_nodi)) + 1, dtype=np.int64)
    indice[np.array(tag_nodi, dtype=np.int64)] = np.arange(len(tag_nodi))
    celle = xyz[indice[conn]]                       # (n_celle, 4, 3)
    centri = celle.mean(axis=1)

    # le 4 facce di ogni tetraedro, come terne ORDINATE di nodi: due celle
    # adiacenti condividono la stessa terna ordinata
    FACCE = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))
    n_celle = len(conn)
    chiavi = np.empty((n_celle * 4, 3), dtype=np.int64)
    proprietaria = np.empty(n_celle * 4, dtype=np.int64)
    for k, f in enumerate(FACCE):
        blocco = slice(k * n_celle, (k + 1) * n_celle)
        chiavi[blocco] = np.sort(conn[:, f], axis=1)
        proprietaria[blocco] = np.arange(n_celle)

    ordine = np.lexsort((chiavi[:, 2], chiavi[:, 1], chiavi[:, 0]))
    ch, pr = chiavi[ordine], proprietaria[ordine]
    uguali = np.all(ch[:-1] == ch[1:], axis=1)
    (interne,) = np.nonzero(uguali)
    if len(interne) == 0:
        return float("inf"), float("inf"), 1.0

    a, b = pr[interne], pr[interne + 1]
    nodi_faccia = ch[interne]                       # (n_interne, 3)
    P = xyz[indice[nodi_faccia]]                    # (n_interne, 3, 3)
    centro_faccia = P.mean(axis=1)
    normale = np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0])
    normale /= np.maximum(np.linalg.norm(normale, axis=1, keepdims=True), 1e-300)

    dvec = centri[b] - centri[a]
    dnorm = np.maximum(np.linalg.norm(dvec, axis=1), 1e-300)
    coseno = np.abs((dvec * normale).sum(axis=1)) / dnorm
    angoli = np.degrees(np.arccos(np.clip(coseno, -1.0, 1.0)))
    non_orto = float(angoli.max())
    frazione = float((angoli > MAX_NON_ORTHOGONALITY).mean())

    # punto in cui la congiungente attraversa il piano della faccia
    denom = (dvec * normale).sum(axis=1)
    sicuro = np.abs(denom) > 1e-300
    t = np.zeros(len(a))
    t[sicuro] = (((centro_faccia - centri[a]) * normale).sum(axis=1)[sicuro]
                 / denom[sicuro])
    incrocio = centri[a] + t[:, None] * dvec
    skew = float((np.linalg.norm(incrocio - centro_faccia, axis=1) / dnorm).max())
    return non_orto, skew, frazione
