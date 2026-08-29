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

**Decisione chiusa: 5 s in regime stazionario, ṁ_aria = 71.8 g/s,
p_c = 5.22 bar, 104 N** (`docs/architettura.md` §8.4). Ne è seguito che il film
cooling in testa era nel posto sbagliato (§8.5), e `f_film` è a 0.

Resta da chiudere, in ordine di impatto:

1. **Aria resa del compressore** (l/min di targa, non l'aria aspirata). Oggi
   non è contata nel punto di progetto: quando la avrai, diventa margine.
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

## Fase 2 — Storia delle run e DOE  ✅ **chiusa**

`store/db.py` (SQLite + Parquet), `doe.py` (Latin Hypercube deterministico),
`opt/objectives.py` completato con vincoli reali, `scripts/sweep_l0.py`.

**Criteri di chiusura — tutti soddisfatti**
- [x] Due esportazioni della stessa geometria sono **byte-identiche** (SHA256
      di STEP e STL). Ha richiesto di normalizzare l'header STEP: OCCT ci
      scriveva l'ora di creazione, che rendeva il file diverso a ogni export.
      Ora l'header porta il `run_id` invece del timestamp, quindi il file si
      autoidentifica e resta riproducibile.
- [x] La query di riferimento
      `SELECT run_id FROM runs WHERE plug_trunc BETWEEN 0.2 AND 0.4 AND
      Isp_s > 120 AND feasible = 1` risponde su **10 000 run in 4.0 ms**
      (criterio: < 100 ms).
- [x] Le run con working tree sporco finiscono in `runs_dirty` per costruzione
      (routing dal suffisso del `run_id`, non da un flag) e **non** compaiono
      nelle query di sintesi.
- [x] LHS verificato come tale: ogni proiezione monodimensionale ha esattamente
      un punto per strato, e la discrepanza di Kolmogorov-Smirnov batte il
      campionamento casuale su 18 semi su 20.
- [x] Sentinella di coerenza: se `derive` alterasse una variabile libera,
      l'inserimento nel database fallisce. È l'unico punto del sistema che se
      ne accorgerebbe.
- [x] Sweep di 400 punti L0 eseguito senza scarti: 96 fattibili.

**Aggiunta fuori piano: interfaccia web** (`src/zefiro/web/`, `web/`). API
versionata con errori strutturati, coda per i lavori lunghi, spaccato quotato,
solido 3D con taglio, grafici termici e storico. Criteri soddisfatti: la GUI
produce lo **stesso `run_id`** della riga di comando (test dedicato), ogni
errore di dominio ha un codice stabile, e la pagina non lancia valutazioni
destinate a fallire quando mancano dati.

**Rinviato di proposito: Snakemake.** L'architettura lo prevede e resta la
scelta giusta, ma oggi non guadagnerebbe nulla: un DAG serve quando i job sono
costosi, falliscono a metà e vanno ripresi. Le valutazioni L0 costano
millisecondi e `sweep_l0.py` le fa già in modo deterministico e ripartibile
(`INSERT OR REPLACE` sul `run_id`). Snakemake entra in **fase 3**, quando i job
sono run OpenFOAM da ore. Scriverlo adesso significherebbe consegnare codice
non esercitato.

---

## Fase 3 — Mesh e CFD reattivo

**Stato: la MESH è chiusa** (2026-08-29). Il solutore resta aperto.

**Fatto**
- `geometry/profile.fluid_polygon`: dominio fluido come complemento del solido
  dentro un contenitore che arriva all'ambiente (l'aerospike espande
  all'esterno, quindi il contorno del getto non è una parete). I nomi delle
  frontiere nascono qui, in codice puro e testabile.
- `l1/mesh.generate_mesh`: settore 2π/N, periodicità conforme, fori
  d'iniezione **imprintati** sulla faccia come patch separate.
- `l1/mesh.verify_periodicity` e `mesh_quality`: verifiche indipendenti, nodo
  per nodo e faccia per faccia, senza fidarsi né di Gmsh né di `checkMesh`.
- `scripts/mesh_report.py`: la mesh si guarda, non solo si misura.

**Criteri di chiusura — verifica**

| criterio | esito |
|---|---|
| patch periodiche conformi | scarto **3.6e-6** raggi di labbro (34 nm); stesso numero di triangoli e stessa area a 8 cifre |
| skewness sotto il limite OpenFOAM (4) | **1.45** |
| non-ortogonalità | max **70.7°** su lo **0.0001 %** delle facce (una) |
| nomi delle frontiere corretti | verificati contro la geometria, non contro l'ordine di costruzione |
| volume della mesh = volume analitico | verificato in forma chiusa dal poligono meridiano |
| area delle patch d'ingresso = πd²/4 | entro il **3 %**, con risoluzione imposta dal mesher |

**Cosa ha trovato la mesh**

1. **Il ginocchio del fronte di Pareto non è fabbricabile.** 18 fori d'aria da
   2.62 mm su un arco di 4.11 mm lasciano 0.23 mm di materiale; a N = 24 il
   setto è negativo. `min_feature` non poteva vederlo perché guarda le *quote*,
   non le *distanze*. Aggiunto `injector_pitch`.
2. **Il foro del film è 0.095 mm**, sotto qualunque minimo SLM. Non c'è vincolo
   che lo fermi finché `min_feature_size` resta il TODO 7: è la prova che quel
   TODO non è burocrazia.
3. Due correzioni istintive alla qualità della mesh — raffinare vicino all'asse
   e ottimizzare con Netgen — **peggiorano** entrambe la mesh. La causa vera era
   il contorno del plug dato come 135 segmenti invece che come una spline.

**Ancora aperto**
- `l1/cfd.run_reacting_foam`: scrittura del caso OpenFOAM e lancio. OpenFOAM non
  è installato in questo ambiente, quindi va fatto e verificato in WSL.
- Da lì escono le due cose che bloccano il resto: la **metrica di miscelazione**
  (fase 5, ramo L1) e i **carichi non uniformi** per l'ottimizzazione topologica
  (fase 7).

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
- Verifica del FEM sui casi analitici già implementati in `zefiro/structural.py`
  (Lamé per il cilindro in pressione, parete impedita per la tensione termica),
  entrambi entro l'1 %.
- **Attenzione a che cosa si sta verificando**: a 5.22 bar la pressione produce
  6 MPa contro i ~460 del termico. Un FEM che azzecca Lamé e sbaglia il termico
  passerebbe la verifica e fallirebbe il progetto.
- Convergenza di mesh sul FEM: σ_von Mises di picco stabile entro il 5 %.
- Il caso di riferimento **con film cooling disattivato** deve dare una parete
  che supera la temperatura ammissibile. Se non la supera, il modello termico
  non sta caricando abbastanza e va indagato prima di fidarsi dei casi con film.

---

## Fase 4bis — Design generativo  ⟵ **motore costruito, in attesa dei carichi veri**

Quello costruito finora è design **parametrico**: 12 quote dentro una topologia
decisa a mano. Il design **generativo** è un'altra cosa — la forma è il
risultato, non l'ingresso. Sono tre problemi distinti, con matematiche diverse,
e conflaterli è il modo più rapido di non risolverne nessuno.

### (a) Ottimizzazione topologica della struttura ✅ **motore pronto**

`src/zefiro/topopt/`. FEM Q4 piano e assialsimmetrico, filtri di densità,
Heaviside, filtro di stampabilità SLM, SIMP, obiettivo di tensione con
variabile aggiunta, Lagrangiano aumentato.

Verificato dal basso: patch test a 3·10⁻¹⁶ · cilindro in pressione contro Lamé
allo 0.002 % · adjoint del filtro di stampabilità a 3·10⁻¹⁰ · sensibilità di
compliance e tensione contro differenze finite a 10⁻⁸ · trave MBB che riproduce
il traliccio noto.

**Due risultati emersi facendolo girare, entrambi non ovvi:**

1. **Minimizzare la compliance è sbagliato con carico termico.** Il carico è la
   dilatazione impedita, quindi è proporzionale al materiale: aggiungere
   materia aggiunge rigidezza *e* carico, e la compliance sale. Misurato: la
   sensibilità che ignora il termine `2 u' df/dx` sbaglia del 159 % e ha il
   **segno invertito in 19 elementi su 80**. Un ottimizzatore guidato da quella
   toglie materiale dove serve. Risolto passando a un obiettivo di tensione.
2. **Il criterio di ottimalità (OC) diverge sulla tensione.** OC è derivato per
   la compliance e presuppone sensibilità sempre negativa; sulla tensione cambia
   segno. Forzandolo, la tensione di picco è salita a 1.4·10⁷ MPa. Risolto con
   Lagrangiano aumentato + L-BFGS-B.

**Che cosa manca perché produca nervature vere.** Con carico uniforme lungo
l'asse la risposta ottima *è* spessore uniforme, e infatti esce quella — è
fisica corretta, non un difetto. Le nervature nascono dove il carico è
**disuniforme**: picchi locali di flusso termico, reazioni di flangia, la zona
di gola. Quei carichi oggi sono una mia idealizzazione; devono venire dalla CFD
(fase 3). **Il motore generativo è pronto e in attesa dei carichi veri.**

*Criterio di chiusura:* la topologia ottimizzata gira sui carichi mappati dalla
CFD, non su un gradiente lineare assunto, e il risultato passa i vincoli di
processo SLM reali (TODO K).

### (b) Forma libera del percorso fluido, accoppiata al mixing ❌ **bloccata dalla fase 3**

È la critica giusta: oggi il mixing entra solo attraverso il rapporto delle
quantità di moto `J`, che è un proxy 0D, e non è accoppiato con nulla di
aerodinamico. Accoppiarli davvero vuol dire un'unica funzione obiettivo che
vede insieme qualità di miscelazione e perdita di pressione totale, su un
contorno **libero** (FFD o spline) invece che sul mio contorno di Angelino.

Il metodo elegante sarebbe l'**aggiunto** del solutore reattivo. Non è
praticabile: l'aggiunto disponibile in OpenFOAM è incomprimibile e non reattivo,
e scriverne uno per un RANS reattivo comprimibile è un lavoro da tesi di
dottorato, non da fine settimana.

La strada praticabile, e che la tua macchina regge, è quella già prevista
dall'architettura: **ottimizzazione su surrogato**. 20–40 variabili di forma
FFD, DOE Latin Hypercube, Gaussian Process, e acquisizione multi-obiettivo che
decide quali run CFD lanciare. Ore o giorni di calcolo, ma reale.

*Prerequisito non aggirabile:* la CFD della fase 3. Ottimizzare la forma su una
correlazione ±30 % significa ottimizzare il rumore.

### (c) Topologia dell'iniezione ❌ **dopo la (b)**

Numero, disposizione e forma degli elementi. È in parte discreto e non si tratta
con SIMP; va con ottimizzazione mista intera sul surrogato, una volta che (b)
sa misurare il mixing.

---

## Fase 5 — Obiettivi, vincoli, ottimizzazione multi-obiettivo

**Stato: il ramo L0 è CHIUSO** (2026-08-29). Il ramo L1 resta aperto perché
dipende dalla fase 3.

**Fatto**
- `opt/objectives.py`: registro a 6 obiettivi e 10 vincoli. La scelta di *quali*
  mettere nel fronte è **misurata**, non assunta: `scripts/objective_screening.py`
  campiona 400 punti e calcola la correlazione di rango. Vedi
  `docs/architettura.md` §8ter.1.
- `l0/chemistry.py`: tempo chimico come **blowout di un PSR**, non come ritardo
  di autoaccensione (che a 300 K darebbe la risposta sbagliata). Tabulato su
  griglia, errore d'interpolazione misurato.
- `l0/mixture.solution`: cache della `ct.Solution`. 107 → 13 ms per valutazione.
- `geometry/profile.py`: volume del pezzo in **forma chiusa esatta** (rivoluzione
  del poligono meridiano). Sostituisce una stima di parete sottile che sbagliava
  del 7.7 %. Coincide con OCCT a ~1e-15 e costa microsecondi.
- `opt/driver.py`: NSGA-II, `seed` obbligatorio, interi riparati dentro
  l'algoritmo, `front_quality` per la robustezza al seme.
- `store/db.RunStore.migrate()`: i registri definiscono le colonne, quindi
  aggiungere un obiettivo è un cambio di schema. Migrazione per `ALTER TABLE`;
  un database più recente del codice viene rifiutato.
- `scripts/optimize_l0.py` e `scripts/pareto_report.py`.

**Criteri di chiusura — verifica**

| criterio | esito |
|---|---|
| fronte stabile rispetto al seme, sotto il 5 % | **0.64 %** di dispersione dell'ipervolume normalizzato su 3 semi |
| nessun punto viola un vincolo di fabbricabilità | verificato: `cv_min = cv_avg = 0` dalla 4ª generazione |
| esiste un punto che soddisfa tutti i vincoli hard | sì, l'intera popolazione finale |
| il driver ritrova un fronte noto | ZDT1, errore mediano 1.1e-4, con convergenza dimostrata al crescere del budget |

**Cosa ha trovato la prima run, che è il motivo per cui si fanno**

L'ottimizzatore ha portato `d_ox_ratio` al massimo del box per **azzerare il Δp
dell'iniettore d'aria** e potersi permettere p_c = 5.94 bar contro i 6.0 bar di
serbatoio. Non era un difetto numerico: era un vincolo **mancante**. Il Δp
minimo era imposto sul combustibile e non sull'aria, che è il 94 % della
portata. Aggiunto `ox_dp_stability`; il tetto di p_c torna a 5.22 bar, cioè
esattamente il valore ricavato a mano mesi prima per un'altra strada.

È il comportamento atteso di un ottimizzatore ed è il motivo per cui vale la
pena farlo girare: **trova le cose che non hai scritto.**

**Ancora aperto (dipende dalla fase 3)**
- `objectives_l1`: margini strutturali e termici dal FEM.
- Il vincolo di **mixing**. A L0 τ_chem è 30–100 µs, ordini di grandezza sotto
  qualunque tempo di miscelazione reale: Zefiro è limitato dalla miscelazione,
  non dalla chimica, e la miscelazione L0 non la vede.

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
| 5 | ~~Il raffreddamento della gola non ha soluzione~~ → **RISOLTO per il banco**: canale ad acqua da 0.5 mm a 12 m/s porta la gola a 446 °C (§8.6). Resta aperto per un eventuale motore volante, che non può portarsi la canna dell'acqua. | Non blocca più la fase 3. |
| 5b | Il **film cooling** entra a L0 solo come sottrazione di massa dal core, senza modello di efficienza, e la sua posizione (testa) è quella sbagliata. | Fase 3: serve un modello di efficienza e un'iniezione vicino alla gola. |
| 6 | Il GPL è trattato come **gas ideale** all'iniezione. Vicino alla tensione di vapore non lo è. | Fase 0, punto 2: dipende dalla modalità di prelievo. |
| 7 | Il sistema di **accensione** non è modellato. | Fase 3, se serve simulare il transitorio di avvio. |
| 8 | Non c'è modello di **stabilità di combustione** (né bassa né alta frequenza). Il Δp iniettore è oggi solo un vincolo di soglia. | Fase 5: se il motore instabile passa i vincoli, il vincolo è sbagliato. |
