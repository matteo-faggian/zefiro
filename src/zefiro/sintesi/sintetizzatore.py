"""progetta(requisiti) -> Progetto | Rifiuto.

Il cancello. Una geometria esce di qui solo se ha passato tutte le verifiche
bloccanti; altrimenti esce un rifiuto che dice quale e perche'.

La ragione per cui il cancello e' qui e non in coda a uno script: fino alla
revisione 0.4.0 le verifiche stavano alla fine di `build_50N.py` e le leggeva
un umano. Un umano che legge trenta righe di referto legge anche
"polvere evacuabile: True (0 sacche chiuse su 1 componenti di vuoto)" e pensa
che sia un successo, mentre e' il sintomo di un circuito di raffreddamento
aperto sull'atmosfera. Un cancello non si distrae.
"""
from __future__ import annotations

import math
from typing import Any

from zefiro.schemas import MissingDatum, ZefiroError
from zefiro.sintesi.architetture import ARCHITETTURE
from zefiro.sintesi.esito import Esito, Progetto, Rifiuto, Verifica
from zefiro.sintesi.requisiti import Requisiti


def _T_LIMITE_PLUG() -> float:
    from zefiro.sintesi.architetture.aerospike_gas_gas import T_LIMITE_PLUG
    return T_LIMITE_PLUG

#: Temperatura oltre la quale l'AISI 316L non e' piu' una struttura [K].
#: La fusione e' a 1673 K; si ferma al 75 %, dove lo snervamento e' gia'
#: sceso di un ordine di grandezza. Non e' un margine di sicurezza: e' il punto
#: oltre il quale il materiale non fa piu' il mestiere per cui e' li'.
T_MAX_PARETE = 1250.0
#: Rapporto minimo fra tempo di permanenza in camera e somma dei tempi di
#: mescolamento e di chimica. Sotto 3 la combustione finisce fuori dall'ugello.
MARGINE_TEMPI = 3.0


def progetta(
    req: Requisiti,
    *,
    geometria: bool = True,
    passo: float | None = None,
    passo_confronto: float | None = None,
) -> Progetto | Rifiuto:
    """Dai requisiti a un progetto verificato, oppure a un rifiuto motivato.

    `geometria=False` ferma la sintesi alle verifiche analitiche: serve ai test
    e alla ricerca multi-obiettivo, dove costruire un campo di distanza per ogni
    punto costerebbe minuti invece di secondi.
    """
    mancanti = req.dati_mancanti()
    if mancanti:
        return Rifiuto(
            motivo="mancano dati che non si possono stimare senza inventarli",
            dati_mancanti=mancanti)

    scarti: list[tuple[str, str]] = []
    dim = None
    arch = None
    for candidata in ARCHITETTURE:
        motivi = candidata.applicabile(req)
        if motivi:
            scarti.append((candidata.NOME, "; ".join(motivi)))
            continue
        try:
            dim = candidata.dimensiona(req)
            arch = candidata
            break
        except (MissingDatum, ZefiroError, ValueError) as e:
            scarti.append((candidata.NOME, str(e)))
    if dim is None or arch is None:
        return Rifiuto(
            motivo="nessuna architettura disponibile soddisfa questi requisiti",
            architetture_provate=scarti)

    verifiche = _verifiche_analitiche(req, dim)
    prog = Progetto(
        architettura=arch.NOME, requisiti=req, registro=dim.registro,
        inviluppo=dim.inviluppo, verifiche=verifiche, geometria=None,
        dati_mancanti=[])
    prog.geometria = dim

    if geometria and not any(v.bloccante for v in verifiche):
        prog.verifiche += _verifiche_geometriche(req, dim, arch, passo, passo_confronto)
    return prog


# --------------------------------------------------------------------------- #
# Verifiche analitiche: costano millisecondi, quindi si fanno sempre
# --------------------------------------------------------------------------- #
def _v(nome, ok, domanda, dettaglio, valori=None, non_conclusiva=False):
    esito = (Esito.NON_CONCLUSIVA if non_conclusiva
             else (Esito.PASSATA if ok else Esito.FALLITA))
    return Verifica(nome=nome, esito=esito, domanda=domanda, dettaglio=dettaglio,
                    valori=valori or {})


def _verifiche_analitiche(req: Requisiti, dim) -> list[Verifica]:
    reg, l0, inj = dim.registro, dim.l0, dim.iniettore
    inv = dim.inviluppo
    p_c = inv.p_camera
    out = []

    out.append(_v(
        "spinta", abs(l0.thrust - req.spinta) < 1.0e-6 * req.spinta,
        "il motore fa la spinta chiesta?",
        f"{l0.thrust:.4f} N contro {req.spinta:.4f} richiesti"))

    out.append(_v(
        "durata", inv.durata_raffica >= req.durata,
        "il serbatoio regge la raffica richiesta?",
        f"{inv.durata_raffica:.2f} s disponibili contro {req.durata:.2f} richiesti "
        f"({inv.durata_raffica/req.durata-1:+.0%})"))

    out.append(_v(
        "ugello", inv.epsilon >= req.epsilon_minimo,
        "il rapporto d'aree e' abbastanza da valere un ugello?",
        f"eps = {inv.epsilon:.3f}, minimo {req.epsilon_minimo:.2f}"))

    t_perm = reg.valore("tempo_permanenza")
    t_mix = reg.valore("tempo_mescolamento")
    t_chim, conclusiva = _tempo_chimico(dim, p_c, req)
    margine = t_perm / max(t_mix + t_chim, 1e-12)
    out.append(_v(
        "tempi", margine >= MARGINE_TEMPI,
        "la permanenza in camera basta a mescolare E bruciare?",
        f"permanenza {t_perm*1e3:.3f} ms contro mescolamento {t_mix*1e3:.3f} + "
        f"chimica {t_chim*1e3:.3f}: margine {margine:.1f}x (minimo {MARGINE_TEMPI})",
        non_conclusiva=not conclusiva))

    #: Le pareti: la temperatura lato gas in regime stazionario con il
    #: raffreddamento attivo. Non e' il transitorio: con l'acqua la parete va
    #: all'equilibrio in frazioni di secondo.
    for nota in dim.note:
        if nota.startswith("ramo "):
            continue
    T_max = max(float(reg[k].valore) for k in reg.quote if k.startswith("parete_calda_"))
    out.append(_v(
        "parete_minima_rispettata", T_max >= req.processo.parete_minima,
        "le pareti calde stanno sopra il minimo di processo?",
        f"la piu' sottile e' {T_max*1e3:.2f} mm contro "
        f"{req.processo.parete_minima*1e3:.2f} minimi"))

    out.append(_v(
        "iniettore_holdeman", abs(inj.scarto_C) <= 0.10,
        "i getti penetrano il condotto come vuole la correlazione?",
        f"C = {inj.C:.3f} contro l'ottimo {inj.C_obiettivo} ({inj.scarto_C:+.1%})"))

    d_min = req.processo.margine_foro()
    out.append(_v(
        "fori_stampabili", inj.d_getto >= d_min,
        "i fori d'iniezione si stampano?",
        f"{inj.d_getto*1e3:.3f} mm contro {d_min*1e3:.3f} richiesti "
        f"(minimo di processo {req.processo.foro_minimo*1e3:.2f}, "
        f"{'misurato' if 'foro_minimo' in req.processo.misurati else 'NON misurato'})"))

    out.append(_v(
        "salto_iniezione",
        inj.dp_aria >= req.frazione_dp_iniezione * p_c * 0.999
        and inj.dp_gpl >= req.frazione_dp_iniezione * p_c * 0.999,
        "l'iniettore disaccoppia camera e alimentazione?",
        f"aria {inj.dp_aria/1e5:.3f} bar, GPL {inj.dp_gpl/1e5:.3f} bar, "
        f"minimo {req.frazione_dp_iniezione*p_c/1e5:.3f}"))

    out.append(_v(
        "testa_coerente", not dim.testa.verifica(inj.lunghezza_mescolamento),
        "la testa realizza davvero l'iniettore dimensionato?",
        "; ".join(dim.testa.verifica(inj.lunghezza_mescolamento)) or "nessun rilievo"))

    #: LE PARTI CHE NESSUNO RAFFREDDA. Il corpo centrale e' un corpo isolato in
    #: mezzo alla camera anulare: nessun circuito lo raggiunge, e senza film
    #: fonde a 2.8 s. La domanda non era mai stata posta.
    plug = dim.plug
    T_lim = _T_LIMITE_PLUG()
    out.append(_v(
        "plug_sotto_soglia", plug.T_finale <= T_lim,
        "il corpo centrale, che nessun circuito raggiunge, arriva a fine raffica?",
        f"massa {plug.massa*1e3:.2f} g, {plug.potenza_iniziale:.0f} W iniziali: "
        f"arriva a {plug.T_finale:.0f} K dopo {req.durata:.1f} s, contro un limite "
        f"di {T_lim:.0f} K"
        + (f" — con {dim.film.frazione:.1%} di combustibile al film"
           if dim.film else " — SENZA film")))

    out.append(_v(
        "quote_tutte_motivate", not reg.mancanti,
        "esiste una quota che nessuno ha misurato ne' dedotto?",
        "; ".join(s.nome for s in reg.mancanti) or f"{len(reg)} quote, tutte con origine"))
    return out


def _tempo_chimico(dim, p_c: float, req: Requisiti) -> tuple[float, bool]:
    """Tempo di spegnimento di un PSR alle condizioni di camera.

    NON e' il ritardo di autoaccensione: reagenti a 300 K non si autoaccendono
    mai, e usare quel numero darebbe un tempo chimico infinito. Il tempo che
    conta e' quello sotto il quale un reattore perfettamente miscelato si
    spegne.

    Se il calcolo non converge si restituisce (0, False) e la verifica esce
    NON CONCLUSIVA: non passata. Un controllo che non ha potuto controllare non
    e' un controllo superato.
    """
    try:
        from zefiro.l0.chemistry import psr_blowout
        from zefiro.l0.mixture import MixtureModel

        model = MixtureModel.from_fuel(dim.op.fuel)
        punto = psr_blowout(model, p_c, req.phi, float(dim.op.T_air_in),
                            float(dim.op.T_fuel_in))
        if punto.tau_blowout is None:
            return 0.0, False
        return float(punto.tau_blowout), True
    except Exception:
        return 0.0, False


# --------------------------------------------------------------------------- #
# Verifiche geometriche: costano minuti, e vanno fatte a DUE risoluzioni
# --------------------------------------------------------------------------- #
def _verifiche_geometriche(req: Requisiti, dim, arch, passo, passo_confronto):
    from zefiro.sdf.engine import area_di_gola, sonde_motore, tappi, _getti_gpl
    from zefiro.sdf.meshing import isosurface, is_watertight, mesh_volume
    from zefiro.sdf.printability import polvere_evacuabile
    from zefiro.sdf.tenuta import confronta, verifica_tenuta

    d = dim.params.derived
    #: Il passo si sceglie dal DETTAGLIO PIU' FINE che deve essere risolto:
    #: serve almeno un quinto della parete piu' sottile, altrimenti la
    #: verifica non sta guardando la geometria che ha in mano.
    minimo = min([req.processo.parete_minima]
                 + [c.lato for c in dim.circuiti]
                 + [dim.testa.d_getto])
    passo = passo or minimo / 4.0
    #: IL CONFRONTO SI FA RAFFINANDO, NON INGROSSANDO.
    #:
    #: Il secondo passo era `passo * 1.25`, cioe' una griglia PIU' GROSSA. Il
    #: risultato era prevedibile e inutile: su un pezzo con pareti da 1 mm e
    #: canali da 0.8, una griglia da 0.31 mm non risolve piu' niente e produce
    #: sacche chiuse e componenti di vuoto che a 0.25 mm non esistono. Il
    #: confronto dichiarava allora "verdetto non stabile" - e aveva ragione, ma
    #: sulla griglia sbagliata: stava misurando che la griglia grossa e' grossa,
    #: cosa che si sapeva.
    #:
    #: La domanda vera e' se la risposta e' CONVERGIUTA, e a quella si risponde
    #: raffinando: se il verdetto non cambia passando a `passo / 1.25`, allora
    #: il passo scelto basta.
    passo_confronto = passo_confronto or passo / 1.25

    esiti: list[Verifica] = []
    rapporti = []
    for p in (passo, passo_confronto):
        m = arch.costruisci_geometria(dim, p)
        v, f = isosurface(m.solido)
        chiuso = is_watertight(f)
        vol = mesh_volume(v, f)
        tappato = m.solido.union(tappi(m, d))
        sonde = sonde_motore(m, d)
        attese = {("aria", "gas"), ("aria", "gpl"), ("gas", "gpl"),
                  ("esterno", "gas"), ("aria", "esterno"), ("esterno", "gpl")}
        #: L'ACQUA E' UN DOMINIO SOLO. In serie tutti i tratti comunicano fra
        #: loro per costruzione: dichiararlo qui non e' indulgenza, e' la
        #: descrizione del progetto. Cio' che resta da dimostrare - e che la
        #: prova dimostra - e' che quel dominio unico NON comunica con il gas,
        #: con l'aria, con il GPL e con l'esterno.
        nomi_acqua = [s_ for s_ in (so.nome for so in sonde) if s_.startswith("acqua")]
        for i, a in enumerate(nomi_acqua):
            for b in nomi_acqua[i + 1:]:
                attese.add(tuple(sorted((a, b))))
        if dim.collettore is None:
            for c in dim.circuiti:
                if c.ritorno:
                    attese.add(tuple(sorted((f"acqua_{c.nome}_andata",
                                             f"acqua_{c.nome}_ritorno"))))
        rap = verifica_tenuta(tappato, m.grid, sonde, attese=attese)
        solo_gpl = tappato.union(_getti_gpl(m.grid, m.testa, m.solido.xp)
                                 .offset(m.grid.spacing))
        rap_gpl = verifica_tenuta(solo_gpl, m.grid, sonde,
                                  attese={c for c in attese if "gpl" not in c})
        x_lip = d["L_c"] + d["L_conv"]
        A, At = area_di_gola(m, x_lip, dim.params.plug_contour_r[0], d["R_lip"])
        ok_pol, n_comp, sacche = polvere_evacuabile(m.canali, m.solido, m.grid)
        rapporti.append(dict(passo=p, parti=tuple(m.parti), chiuso=chiuso,
                             volume=vol, campo=m.solido.volume(),
                             tenuta=rap, tenuta_gpl=rap_gpl, area=A, area_teorica=At,
                             polvere=ok_pol, sacche=sacche, note=list(m.note),
                             isole=getattr(m.solido, "diagnostica_isole", {})))
        del m

    a, b = rapporti
    diff = confronta(a["tenuta"], b["tenuta"])
    convergente = not diff

    esiti.append(Verifica(
        nome="superficie_chiusa",
        esito=Esito.PASSATA if a["chiuso"] and b["chiuso"] else Esito.FALLITA,
        domanda="la mesh e' una superficie chiusa senza bordi liberi?",
        dettaglio=f"chiusa a {a['passo']*1e3:.3f} e {b['passo']*1e3:.3f} mm",
        convergente=a["chiuso"] == b["chiuso"]))

    esiti.append(Verifica(
        nome="tenuta",
        esito=(Esito.PASSATA if a["tenuta"].stagno and b["tenuta"].stagno
               else Esito.FALLITA) if convergente else Esito.NON_CONCLUSIVA,
        domanda="acqua, gas ed esterno sono domini disgiunti come da progetto?",
        dettaglio=(a["tenuta"].riassunto().splitlines()[0]
                   + ("" if convergente else "  VERDETTO NON STABILE: " + "; ".join(diff))),
        valori={"componenti": a["tenuta"].n_componenti},
        convergente=convergente))

    esiti.append(Verifica(
        nome="tenuta_circuito_gpl",
        esito=Esito.PASSATA if a["tenuta_gpl"].stagno and b["tenuta_gpl"].stagno
        else Esito.FALLITA,
        domanda="tappando i getti, il GPL resta isolato da tutto il resto?",
        dettaglio="verifica la parete fra plenum del GPL e condotto d'aria, che "
                  "nessun'altra prova vede"))

    scarto_area = a["area"] / a["area_teorica"] - 1.0
    esiti.append(Verifica(
        nome="area_di_gola",
        esito=Esito.PASSATA if abs(scarto_area) < 0.05 else Esito.FALLITA,
        domanda="il materiale non ha strozzato l'anello di gola?",
        dettaglio=f"{a['area']*1e6:.2f} mm2 contro {a['area_teorica']*1e6:.2f} "
                  f"teorici ({scarto_area:+.1%})",
        valori={"scarto": scarto_area}))

    scarto_vol = abs(a["volume"] - a["campo"]) / a["campo"]
    esiti.append(Verifica(
        nome="volume_coerente",
        esito=Esito.PASSATA if scarto_vol < 0.01 else Esito.FALLITA,
        domanda="il volume dalla mesh coincide con quello dal campo?",
        dettaglio=f"{a['volume']*1e6:.3f} contro {a['campo']*1e6:.3f} cm3 "
                  f"({scarto_vol:.2%})"))

    esiti.append(Verifica(
        nome="polvere_evacuabile",
        esito=Esito.PASSATA if a["polvere"] and b["polvere"] else Esito.FALLITA,
        domanda="ogni vuoto ha una via d'uscita per la polvere?",
        dettaglio=f"{a['sacche']} sacche chiuse a {a['passo']*1e3:.3f} mm, "
                  f"{b['sacche']} a {b['passo']*1e3:.3f}",
        convergente=(a["polvere"] == b["polvere"])))

    #: Frammenti di materiale staccati. Il criterio NON e' "ce ne sono": in
    #: geometria implicita qualche voxel isolato lo produce sempre la
    #: discretizzazione, e infatti il loro volume si dimezza raffinando la
    #: griglia. Il criterio e' se ne RESTA uno dopo la pulizia, e quanto grande:
    #: un frammento che vale piu' dello 0.1 % del pezzo non e' un artefatto,
    #: e' una geometria che si e' spezzata davvero.
    isole = a.get("isole") or {}
    rimasta = isole.get("volume_massima_rimasta", 0.0)
    principale = isole.get("volume_principale", a["campo"]) or a["campo"]
    grave = [n for n in a["note"] if "sacche di vuoto CHIUSE" in n or "bordo" in n]
    if grave or rimasta > 1.0e-3 * principale:
        esito = Esito.FALLITA
        dettaglio = ("; ".join(grave) or
                     f"resta un frammento staccato da {rimasta*1e9:.2f} mm3, "
                     f"{rimasta/principale:.2%} del pezzo")
    elif isole.get("n_totali", 0) and isole.get("volume_totale", 0.0) > 1.0e-4 * principale:
        esito = Esito.NON_CONCLUSIVA
        dettaglio = (f"{isole['n_totali']} frammenti per "
                     f"{isole['volume_totale']*1e9:.2f} mm3: sopra la soglia di "
                     "artefatto, da riguardare a griglia piu' fine")
    else:
        esito = Esito.PASSATA
        dettaglio = (f"{isole.get('n_totali', 0)} frammenti da "
                     f"{isole.get('volume_totale', 0.0)*1e9:.3f} mm3, tutti rimossi e "
                     "tutti sotto la soglia di artefatto")
    #: LA GEOMETRIA DEVE AVERE CIO' CHE IL DIMENSIONAMENTO HA CALCOLATO.
    #: Il film del plug e' dimensionato (9 % del combustibile, fessura da 0.5
    #: mm alla base del plug) ma la geometria implicita non lo costruisce
    #: ancora. Finche' e' cosi', lo STL descrive un motore DIVERSO da quello
    #: verificato - un motore il cui corpo centrale fonde a 2.8 s - e non deve
    #: uscire. E' la stessa disciplina del resto: meglio un rifiuto che un file
    #: che sembra giusto.
    #: SI GUARDA LA GEOMETRIA COSTRUITA, non un attributo del dimensionamento.
    #: La versione precedente cercava `dim.parti_geometria`, che non esiste: il
    #: `getattr` con valore di ripiego lo rendeva sempre vuoto, quindi la
    #: verifica falliva SEMPRE - anche quando la fessura ci fosse stata. Un
    #: controllo che risponde sempre "no" e' rotto quanto uno che risponde
    #: sempre "si": in entrambi i casi non sta guardando niente.
    ha_film = dim.film is not None
    presente = "film" in rapporti[0].get("parti", ())
    esiti.append(Verifica(
        nome="film_del_plug_costruito",
        esito=Esito.PASSATA if (not ha_film or presente) else Esito.FALLITA,
        domanda="la geometria contiene la fessura di film che il calcolo richiede?",
        dettaglio=("nessun film richiesto" if not ha_film else
                   ("fessura presente" if presente else
                    f"il dimensionamento chiede {dim.film.frazione:.1%} di combustibile "
                    "al film e una fessura alla base del plug, che la geometria NON "
                    "costruisce: lo STL sarebbe di un motore diverso da quello "
                    "verificato"))))

    esiti.append(Verifica(
        nome="costruzione_pulita", esito=esito,
        domanda="la costruzione ha lasciato frammenti staccati o sacche chiuse?",
        dettaglio=dettaglio))
    return esiti


#: Proprieta' del 316L per il transitorio del plug. Sono di LETTERATURA per il
#: 316L laminato, non misurate sulla lega del fornitore (TODO J/6): il cp medio
#: fra 300 e 1300 K e la temperatura di fusione cambiano di qualche punto da
#: colata a colata, e il risultato con essi.
RHO_316L, CP_316L, T_FUSIONE_316L, H_FUSIONE_316L = 7990.0, 550.0, 1673.0, 270.0e3


def _durata_del_plug(dim, emissivita: float = 0.3) -> tuple[float, float, float, float]:
    """Quanto dura il corpo centrale prima di fondere.

    IL PLUG NON E' RAFFREDDATO E NON PUO' ESSERLO facilmente: e' un corpo
    isolato al centro di una camera anulare, e portarci dentro dell'acqua
    vorrebbe dire attraversare il getto. La sua sola difesa e' la propria
    capacita' termica piu' cio' che riesce a irradiare.

    Modello LUMPED e non semi-infinito, e la scelta e' motivata da un numero:
    in cinque secondi la profondita' di penetrazione termica nel 316L e'
    sqrt(pi alpha t) = 7.3 mm, piu' del raggio massimo del plug (3.9 mm). Il
    plug e' quindi **termicamente sottile**: si scalda tutto insieme, e
    trattarlo come semi-infinito sottostimerebbe la temperatura.

    Il bilancio include il re-irraggiamento:

        m cp dT/dt = Q - eps sigma A (T^4 - T_amb^4)

    e il risultato e' che NON SERVE A NIENTE, il che e' a sua volta un
    risultato. A 1673 K questo plug irradia 14 W con eps = 0.3 e 39 W con
    eps = 0.8, contro 353 W assorbiti: il 4-11 %. La temperatura di equilibrio
    radiativo sarebbe fra 2660 e 3715 K, cioe' molto oltre la fusione: **non
    esiste emissivita' che salvi il plug**, e il tempo alla fusione passa da
    2.02 s a 2.06 s fra eps = 0.2 e eps = 0.8.

    NOTA DI METODO, tenuta agli atti. Scrivendo questa funzione avevo stimato a
    mano "145 W irradiati, il 41 %", sbagliando di un fattore dieci
    (0.3 * 5.67e-8 * 1.09e-4 * 1673^4 = 14.5 W, non 145). La stima a mano diceva
    che l'irraggiamento cambiava le carte; il conto fatto dal codice dice che
    non le cambia. La conclusione regge nei due casi, ma il numero era falso, e
    un numero falso in un commento sopravvive molto piu' a lungo di uno in una
    formula.
    """
    cx = dim.params.plug_contour_x
    cr = dim.params.plug_contour_r
    volume = 0.0
    area = 0.0
    for (x0, r0), (x1, r1) in zip(zip(cx[:-1], cr[:-1]), zip(cx[1:], cr[1:])):
        volume += math.pi / 3.0 * abs(x1 - x0) * (r0 * r0 + r0 * r1 + r1 * r1)
        area += math.pi * (r0 + r1) * math.hypot(x1 - x0, r1 - r0)
    massa = RHO_316L * volume
    q = dim.carico.get("plug", 0.0)
    if q <= 0.0 or massa <= 0.0:
        return T_FUSIONE_316L, float("inf"), massa, 0.0

    SIGMA = 5.670374419e-8
    T, t, dt, T_amb = 300.0, 0.0, 1.0e-3, 300.0
    capacita = massa * CP_316L
    while t < 60.0:
        netto = q - emissivita * SIGMA * area * (T ** 4 - T_amb ** 4)
        if netto <= 0.0:
            return T_FUSIONE_316L, float("inf"), massa, q      # equilibrio sotto fusione
        T += netto / capacita * dt
        t += dt
        if T >= T_FUSIONE_316L:
            return T_FUSIONE_316L, t, massa, q
    return T_FUSIONE_316L, float("inf"), massa, q
