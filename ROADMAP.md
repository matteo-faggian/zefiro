# Zefiro — Roadmap

Ogni fase ha **criteri di chiusura verificabili**: o il criterio è un test che
gira, o è un numero confrontato con un riferimento indipendente. "Sembra
funzionare" non è un criterio.

Le fasi non sono strettamente sequenziali, ma la **fase 0 blocca tutto ciò che
produce numeri sul motore reale**: finché i dati mancanti non ci sono, il
software gira ma i suoi risultati descrivono un motore ipotetico.

---

## Fase 0 — Chiudere i dati mancanti  ⟵ **BLOCCANTE, e non dipende dal codice**

Elenco completo e motivato in `docs/architettura.md` §8.4. Stanno come `null` in
`config/`.

**La decisione che viene prima di tutte le altre: continuo o raffica.**
`scripts/plant_report.py` mostra che il compressore 3 HP dà 4–7 g/s continui,
mentre il serbatoio da 100 L permette 50–150 g/s per 2–8 secondi. Non è una
sfumatura: in continuo gli iniettori vengono da 0.2 mm e il motore **non è
stampabile in SLM**; a 50 g/s le quote entrano in un intervallo fabbricabile.

Poi, in ordine di impatto:

1. **Aria resa del compressore** (l/min di targa, non l'aria aspirata).
2. **Pressione e temperatura della bombola misurate insieme** → composizione.
3. **Portata massima del riduttore GPL** — a 100 g/s d'aria servono 27 kg/h,
   fuori portata per un riduttore da barbecue.
4. **Tara e peso della bombola** → massa di liquido.
5. **Temperature a valle dei riduttori.**
6. **AISI 316L da SLM**: σ_y(T), E(T), α(T), k(T), ρ, criterio di ammissibilità.
7. **Limiti di processo della stampante**: diametro minimo di foro, spessore
   minimo, angolo massimo di overhang.

**Criteri di chiusura**
- `config/operating_point.yaml` non contiene più alcun `null`.
- `config/materials/aisi316l.yaml` ha `source` compilato con una citazione
  verificabile per ogni proprietà.
- `check_manufacturability` non restituisce violazioni sul punto di progetto
  scelto.
- I test di `test_config.py` che oggi verificano *"i null sono ancora null"*
  diventano test sui range dei valori.

---

## Fase 1 — Fondamenta  ✅ **chiusa**

Contratti dati, scaffold, geometria end-to-end, L0 in Cantera.

**Criteri di chiusura — tutti soddisfatti**
- [x] `pip install -e .` riesce e `pytest` è verde (99 test).
- [x] Da un dizionario di parametri escono uno STEP e uno STL.
- [x] Tenuta verificata due volte con metodi indipendenti: `is_valid` di OCCT
      e manifold-ness della tassellazione controllata a mano.
- [x] Il settore periodico × N riproduce il volume del pezzo intero
      (rapporto 1.000000).
- [x] L0 riproduce un caso verificabile a mano: AFR = 15.6799 (identico al
      conto stechiometrico), T_ad = 2266.7 K contro ~2267 K di letteratura.
- [x] Il contorno dell'aerospike chiude **esattamente** sull'asse al Mach di
      progetto, e l'area di gola geometrica coincide con `πR²/ε`.
- [x] Il butano non viene sostituito di nascosto: solleva `MissingThermoData`.
- [x] Impianto di alimentazione (`feed.py`): limiti del compressore per via
      termodinamica, blowdown del serbatoio, composizione della bombola dedotta
      da (p, T) misurate, autorefrigerazione. 126 test totali.

---

## Fase 2 — Storia delle run e orchestrazione

Rendere ripetibile ciò che oggi è ripetibile solo a mano.

**Da fare**
- `store/db.py`: SQLite con una colonna scalare per ogni voce di `x`, di
  `L0Result` e di `Objectives`; campi spaziali in Parquet referenziati per path.
- Snakemake: DAG `design → L0 → geometria`, con `run_id` come wildcard.
- `scripts/export_winner.py`: copia gli artefatti del design vincente.

**Criteri di chiusura**
- Rilanciare la stessa run produce lo **stesso `run_id`** e file **byte-identici**
  (confronto SHA256), su una macchina appena installata da `environment.yml`.
- La query `SELECT run_id FROM runs WHERE plug_trunc BETWEEN 0.2 AND 0.4 AND
  Isp_s > 120` risponde su 10⁴ run in meno di 100 ms.
- Una run con working tree sporco finisce in `runs_dirty` e **non** compare
  nelle query di sintesi.

---

## Fase 3 — L1 fluidodinamica

**Da fare**
- `l1/mesh.py`: Gmsh sul **negativo** del solido (dominio fluido), settore
  periodico, strato limite risolto, physical groups = `CANONICAL_BOUNDARIES`.
- `l1/cfd.py`: `reactingFoam`, inizializzato dallo stato di equilibrio L0.
- Meccanismo ridotto per la CFD, ottenuto dal San Diego mech.

**Criteri di chiusura**
- **Grid convergence study** su almeno 3 livelli di mesh, con indice GCI di
  Roache calcolato: la spinta deve convergere entro il 2 %.
- I nomi delle boundary prodotte coincidono **esattamente** con
  `CANONICAL_BOUNDARIES` (test automatico).
- Bilancio di massa e di energia chiuso entro lo 0.1 % sul dominio.
- La CFD riproduce `c*` di L0 entro il 5 % su un caso a φ = 1 senza film
  cooling. Se lo scarto è maggiore, va **spiegato** (perdite di ristagno,
  combustione incompleta) prima di procedere, non tarato.
- Il meccanismo ridotto riproduce τ_ign e S_L del San Diego completo entro il
  10 % nel range di p e T di camera.

---

## Fase 4 — L1 strutturale e accoppiamento

**Da fare**
- `l1/mapping.py`: trasferimento one-way di p e q̇ dalla mesh fluida a quella
  strutturale.
- `l1/fem.py`: CalculiX, analisi termica transitoria seguita da meccanica, con
  proprietà del 316L SLM dipendenti dalla temperatura.

**Criteri di chiusura**
- **Conservazione dell'integrale nel mapping**: `|∫q̇dA_CFD − ∫q̇dA_FEM| /
  ∫q̇dA_CFD < 10⁻³`, verificata automaticamente. Se non conserva, il FEM sta
  risolvendo un problema diverso da quello posto dalla CFD.
- Verifica del FEM su un caso a soluzione analitica nota: cilindro in pressione
  (Lamé) e cilindro con gradiente termico radiale stazionario, entrambi entro
  l'1 %.
- Convergenza di mesh sul FEM: σ_von Mises di picco stabile entro il 5 %.
- Il caso di riferimento **con film cooling disattivato** deve dare una parete
  che supera la temperatura ammissibile. Se non la supera, il modello termico
  non sta caricando abbastanza e va indagato prima di fidarsi dei casi con film.

---

## Fase 5 — Obiettivi, vincoli, ottimizzazione multi-obiettivo

**Da fare**
- `opt/objectives.py`: completare `objectives_l1`, con margini strutturali e
  termici e vincoli di fabbricabilità SLM.
- `opt/driver.py`: NSGA-II su L0 (migliaia di valutazioni), con `seed`
  obbligatorio.

**Criteri di chiusura**
- Il fronte di Pareto è **stabile** rispetto al seme: due run con semi diversi
  danno fronti la cui distanza di Hausdorff normalizzata è sotto il 5 %.
- Nessun punto del fronte viola un vincolo di fabbricabilità (verifica a
  posteriori indipendente dall'ottimizzatore).
- Esiste almeno un punto del fronte che soddisfa **tutti** i vincoli hard. Se
  non esiste, il problema è sovravincolato e va detto — non è un fallimento
  dell'ottimizzatore.

---

## Fase 6 — Surrogato e ottimizzazione bayesiana

**Da fare**
- DOE Latin Hypercube su L1 sui candidati sopravvissuti a L0.
- Gaussian Process (BoTorch/GPyTorch) su (x → obiettivi L1).
- Acquisition multi-obiettivo (qNEHVI) che decide quali run L1 lanciare.

**Criteri di chiusura**
- Errore del surrogato in **leave-one-out** sotto il 10 % sugli obiettivi
  principali.
- Il surrogato non estrapola: ogni punto proposto ha incertezza predittiva
  sotto una soglia dichiarata, altrimenti si lancia la run L1.
- A parità di budget di run L1, il bayesiano trova un ipervolume superiore al
  DOE puro. Se non lo fa, il surrogato non serve e va detto.

---

## Fase 7 — L2 e validazione sperimentale

**Da fare**
- Poche run L2 di validazione ad alta fedeltà sul design vincente.
- Confronto con dati di banco: spinta, p_c, temperatura di parete.

**Criteri di chiusura**
- Spinta misurata entro il 10 % della previsione L2, oppure discrepanza
  **spiegata** con un meccanismo fisico identificato.
- Temperatura di parete misurata dentro la banda di incertezza predetta.
- Il modello viene aggiornato **una sola volta** con i dati sperimentali, e la
  versione pre-aggiornamento resta nel repo. Tarare iterativamente il modello
  sui dati finché non torna significa smettere di predire e iniziare a
  interpolare.

---

## Debiti tecnici già noti

Elencati qui perché tacerli li renderebbe invisibili, non inesistenti.

| # | Debito | Quando diventa bloccante |
|---|---|---|
| 1 | Il **supporto del plug** non è modellato: nella realtà servono razze o un pilone che attraversano il flusso. Alterano l'area di gola effettiva e generano scie. | Fase 3: la CFD sul modello senza razze sovrastima le prestazioni. |
| 2 | Il **labbro** è modellato come spigolo vivo con spessore radiale `t_wall`. Il metodo di Angelino assume espansione centrata, cioè labbro affilato rispetto alla scala della gola: qui non lo è. | Fase 3: attendersi uno scostamento reale del contorno. |
| 3 | Nessun modello di **pressione di base** per il plug troncato. `plug_trunc < 1` dà oggi un limite superiore di prestazione. | Fase 3: solo la CFD può dare la pressione di base. |
| 4 | Nessuno **strato limite** nel contorno: manca la correzione per lo spessore di spostamento δ*. | Fase 4: correggere il contorno con δ* dalla RANS. |
| 5 | Il **film cooling** entra a L0 solo come sottrazione di massa dal core, senza modello di efficienza. | Fase 3: senza un modello di efficienza, `f_film` non è ottimizzabile. |
| 6 | Il GPL è trattato come **gas ideale** all'iniezione. Vicino alla tensione di vapore non lo è. | Fase 0, punto 2: dipende dalla modalità di prelievo. |
| 7 | Il sistema di **accensione** non è modellato. | Fase 3, se serve simulare il transitorio di avvio. |
| 8 | Non c'è modello di **stabilità di combustione** (né bassa né alta frequenza). Il Δp iniettore è oggi solo un vincolo di soglia. | Fase 5: se il motore instabile passa i vincoli, il vincolo è sbagliato. |
