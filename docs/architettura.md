# Zefiro — Architettura del software

Versione documento: `arch-0.1.0` · Schema dati: `zefiro-schema-0.1.0`

Questo documento definisce **struttura del repo**, **contratti dati fra i moduli**,
**convenzioni di naming** e **versionamento delle run**. È il contratto: il codice si
adegua a questo documento, non viceversa. Ogni modifica ai contratti incrementa
`SCHEMA_VERSION` e va annotata in fondo (§9).

---

## 1. Premessa fisica: che cosa è Zefiro

Dalle scelte chiuse in fase di intervista:

| Voce | Scelta | Conseguenza sull'architettura |
|---|---|---|
| Ciclo | Pressure-fed, **regolatore a valle** → `p_c` costante | Il punto operativo è **stazionario**: una valutazione per candidato, non una traiettoria. Il loop resta trattabile. |
| Ossidante | Aria da compressore 3 HP + serbatoio 100 L, `p_max = 10 bar` | `p_c` ha un tetto **hard** ben sotto i 10 bar (§1.1). Il serbatoio abilita il funzionamento a raffica (§8.1). |
| Combustibile | GPL da bombola da barbecue, indicata a 8 bar | Il Δp iniettore combustibile è la risorsa più scarsa dell'intero sistema. La pressione misurata identifica la composizione (§8.2). |
| Ugello | **Aerospike** (plug anulare, espansione esterna) | Contorno generato per via analitica (metodo di Angelino, §4.3), non tabellato. |
| Termica | **Film cooling** con GPL | `f_film` è una variabile di progetto, non un post-processing: entra già in L0. |
| Simmetria | **Settore periodico 3D**, `N_inj` elementi | La geometria deve saper emettere sia il solido completo sia il settore 1/N con facce periodiche marcate. |
| Cinetica | San Diego mech | Meccanismo **non vendorizzato**: scaricato e hashato (§6.3). |
| Fabbricazione | **SLM, AISI 316L** | Libertà topologica, ma vincoli su spessore minimo e angolo di overhang → sono **vincoli geometrici**, non estetica. |
| Obiettivi | "tutto" → **Pareto multi-obiettivo** | Il contratto `Objectives` porta un *vettore* di obiettivi + un vettore di vincoli, mai uno scalare pesato. |

### 1.1 Il vincolo che domina tutto il progetto

Non è una scelta di architettura, è aritmetica, e va scritta qui perché condiziona
ogni numero a valle.

La camera è alimentata da due sorgenti a pressione limitata. Perché il flusso
attraverso un iniettore sia **stabile** (non accoppiato all'acustica di camera)
serve un Δp dell'ordine del 15–25 % di `p_c`. Quindi:

```
p_c  <  min( p_aria_reg , p_gpl_reg ) / (1 + Δp/p_c)
p_c  <  8 bar / 1.20  ≈ 6.7 bar        (limite dal lato GPL, il più stringente)
```

e realisticamente, tenendo conto delle perdite di linea e del regolatore,
`p_c ∈ [3, 5] bar`. Il rapporto di espansione disponibile è quindi
`p_c/p_amb ≈ 3–5`, che con γ ≈ 1.25 dei prodotti dà `M_e ≈ 1.5–1.7` e
`ε = A_e/A_t ≈ 1.2–1.4`.

**Conseguenza da tenere presente:** a questo rapporto di pressione l'aerospike è
un ugello legittimo, ma il suo vantaggio caratteristico — la compensazione di
quota — è quasi nullo, perché il getto è appena supersonico e non c'è
sovraespansione da compensare. L'aerospike qui va giustificato come *oggetto di
studio e di validazione del metodo di progetto*, non come scelta di prestazione.
Se in futuro `p_c` salisse (bombole ad alta pressione al posto del compressore),
lo stesso modulo produce il contorno corretto senza modifiche: il codice è
parametrico su `p_c`.

Il modulo `zefiro.l0.nozzle` **emette un warning esplicito** quando
`p_c/p_amb < 8`, invece di lasciar credere che l'aerospike stia lavorando nel suo
regime di progetto.

### 1.2 La seconda cosa che va detta

Con aria come ossidante, il rapporto stechiometrico in massa è
`AFR ≈ 15.67` (propano puro, aria reale — valore calcolato, §7.1). Cioè il
**94 %** del flusso di massa è aria, e di questa il **76 %** è azoto inerte che
assorbe calore senza contribuire. Due conseguenze operative:

1. La portata d'aria è la risorsa dimensionante, e con questo impianto la
   scelta **continuo vs raffica** cambia il motore di un ordine di grandezza
   (§8.1).
2. La portata di GPL è piccolissima in assoluto. Il film cooling sottrae
   `f_film` a un flusso già minimo: prima di ottimizzare `f_film` bisogna
   verificare che il GPL disponibile basti a schermare la parete. Il modulo L0
   calcola il rapporto `ṁ_film / (carico termico)` e lo espone come vincolo, non
   lo assume risolto.

---

## 2. Struttura del repository

```
AA_ENGINE/
├── pyproject.toml              # package `zefiro`, dipendenze, config pytest/ruff
├── environment.yml             # env conda riproducibile (WSL2/Ubuntu)
├── README.md
├── ROADMAP.md
├── .gitignore
│
├── docs/
│   └── architettura.md         # questo file
│
├── config/                     # DATI DI INGRESSO, versionati in Git
│   ├── operating_point.yaml    # condizioni operative reali (contiene TODO espliciti)
│   ├── design_default.yaml     # vettore di progetto x di partenza + bounds
│   └── materials/
│       └── aisi316l.yaml       # proprietà materiale (contiene TODO espliciti)
│
├── src/zefiro/
│   ├── schemas.py              # ⟵ IL CONTRATTO. Tutte le dataclass scambiate.
│   ├── units.py                # costanti e unità: tutto SI, nessuna eccezione
│   ├── config.py               # caricamento YAML + validazione + errori su TODO
│   │
│   ├── geometry/               # x → solido → STEP/STL      (IMPLEMENTATO)
│   │   ├── parameters.py       # separazione liberi / derivati
│   │   ├── aerospike.py        # contorno plug, metodo di Angelino
│   │   ├── chamber.py          # profilo camera + convergente
│   │   ├── injector.py         # elementi di iniezione, settore periodico
│   │   └── build.py            # assemblaggio build123d + export + metriche
│   │
│   ├── l0/                     # termochimica 0D/1D          (IMPLEMENTATO)
│   │   ├── mixture.py          # composizione GPL, φ, AFR
│   │   ├── equilibrium.py      # T_ad, composizione, c*
│   │   ├── nozzle.py           # espansione 1D, C_F, I_sp, spinta
│   │   ├── kinetics.py         # τ_ign, S_L (richiede mech esterno)
│   │   └── cycle.py            # bilancio completo → L0Result
│   │
│   ├── l1/                     # RANS + FEM                  (STUB, fase 3–4)
│   │   ├── mesh.py             # gmsh
│   │   ├── cfd.py              # OpenFOAM
│   │   ├── fem.py              # CalculiX
│   │   └── mapping.py          # CFD → FEM, one-way
│   │
│   ├── opt/                    # obiettivi e ottimizzatore   (STUB, fase 5–6)
│   │   ├── objectives.py
│   │   ├── surrogate.py
│   │   └── driver.py
│   │
│   └── store/                  # storia dei design           (STUB, fase 2)
│       ├── runid.py            # hashing deterministico       (IMPLEMENTATO)
│       └── db.py               # SQLite + Parquet
│
├── scripts/
│   ├── fetch_mechanism.py      # scarica San Diego mech e ne registra lo SHA256
│   ├── build_geometry.py       # CLI: config → STEP + STL
│   └── run_l0.py               # CLI: config → L0Result (JSON)
│
├── tests/
└── runs/                       # ⟵ NON versionato. Output per run-id. Vedi §5.
```

### 2.1 Perché `src/` layout

Perché impedisce che i test importino il package dalla directory di lavoro invece
che dall'installato. Se `import zefiro` funziona nei test, allora funziona anche
quando Snakemake lancia un job in una working directory arbitraria. È l'unico
layout che rende l'errore "funziona da qui ma non da lì" impossibile per
costruzione.

### 2.2 Python su Windows o su WSL — la scelta e il perché

**Raccomandazione: tutto il layer Python dentro WSL2/Ubuntu. Una sola toolchain.**

Non è una preferenza, sono tre argomenti concreti:

1. **Il confine costa e non è trasparente.** OpenFOAM e CalculiX girano solo in
   WSL. Con Python su Windows, ogni valutazione L1 attraversa il confine
   `9p`/`drvfs` con mesh e campi da centinaia di MB. Non è solo lento: la
   traduzione dei path, i permessi POSIX inventati da drvfs e la conversione dei
   line ending sono tre sorgenti indipendenti di errore silenzioso in un sistema
   che deve essere bit-riproducibile.

2. **Riproducibilità binaria.** Cantera, OCCT e NumPy hanno wheel diverse su
   Windows e su Linux, compilate con compilatori e BLAS diversi. Le differenze
   sono nell'ultimo bit, ma un ottimizzatore che gira 10⁴ volte amplifica
   l'ultimo bit. Se il requisito è "stesso input → stesso output *verificabile*",
   avere un solo target di compilazione non è pignoleria: è la condizione perché
   l'affermazione sia falsificabile.

3. **Un solo lockfile.** `environment.yml` + `conda-lock` su Linux descrive
   l'intero ambiente numerico. Con due OS servono due lockfile che divergono.

**Il prezzo, e come si paga.** I visualizzatori CAD e Paraview stanno comodi su
Windows. Si risolve così: il **repo** vive su `I:\AA_ENGINE` (montato in WSL come
`/mnt/i/AA_ENGINE`) perché è codice, è piccolo, e ti serve vederlo da Windows.
La cartella **`runs/`** invece va nel filesystem ext4 di WSL
(es. `~/zefiro-runs`), con `RUNS_ROOT` da variabile d'ambiente, perché è I/O
pesante e su drvfs perdi un fattore ~5–10 su file piccoli e numerosi — ed è
esattamente quello che produce una mesh. Gli artefatti da guardare (STEP, STL,
`.foam`) si copiano su `I:` solo a fine run.

> Se in futuro volessi girare il layer Python su Windows per comodità di
> debugging in un IDE, va bene — ma allora **i risultati prodotti da Windows non
> entrano nel database delle run**. Il campo `env_fingerprint` (§5.2) serve
> proprio a rendere questa regola verificabile invece che sperata.

---

## 3. Il loop, e chi produce cosa

```
                    ┌─────────────────────────────────────────┐
                    │            DesignVector  x              │
                    └────────────────┬────────────────────────┘
                                     ▼
                    ┌─────────────────────────────────────────┐
      derivazione   │  geometry.parameters.derive(x, op)      │
      algebrica     │        → GeometryParams                 │
                    └────────────────┬────────────────────────┘
                            ┌────────┴────────┐
                            ▼                 ▼
         ┌──────────────────────────┐   ┌──────────────────────────┐
         │ geometry.build           │   │ l0.cycle.evaluate        │
         │   → GeometryArtifact     │   │   → L0Result             │
         │   (STEP, STL, metriche)  │   │   (ms, migliaia di volte) │
         └────────────┬─────────────┘   └────────────┬─────────────┘
                      │                              │
                      │      ┌───── screening L0 ────┘
                      │      │  (scarta il 95 % dei candidati)
                      ▼      ▼
         ┌─────────────────────────────────────────────────┐
         │ l1.mesh → MeshArtifact                          │
         │ l1.cfd  → CFDResult   (OpenFOAM, reattivo)      │   minuti/ore
         │ l1.mapping → LoadSet  (p, q̇ a parete)           │
         │ l1.fem  → FEMResult   (CalculiX, termo-mecc.)   │
         └────────────────────┬────────────────────────────┘
                              ▼
         ┌─────────────────────────────────────────────────┐
         │ opt.objectives.assemble → Objectives            │
         │   f = [f1..fm] vettore,  g = [g1..gk] vincoli   │
         └────────────────────┬────────────────────────────┘
                              ▼
         ┌─────────────────────────────────────────────────┐
         │ opt.surrogate (GP)  ⇄  opt.driver (pymoo/BoTorch)│
         └────────────────────┬────────────────────────────┘
                              └──────────► nuovo x  (ritorno in cima)
```

**La geometria è input dei solver.** Non esiste alcun percorso in cui un solver
produce geometria. L'unico export STEP/STL "finale" è una copia dell'artefatto
del design vincente, prodotta da `scripts/export_winner.py` a convergenza.

---

## 4. Contratti dati

Regole generali, valide per **tutti** i contratti:

- **Unità: SI stretto, sempre.** Pa, K, kg/s, m, N, s. Nessun bar, nessun mm,
  nessun grado. La conversione avviene **solo** al confine I/O (lettura YAML,
  scrittura report). `zefiro.units` contiene i fattori; il resto del codice non
  li vede.
- Ogni dataclass è **frozen** (immutabile) e ha `schema_version: str`.
- Ogni dataclass serializza a JSON con `to_dict()` / `from_dict()`; i float
  vanno in JSON con `repr()` (round-trip esatto), mai formattati.
- Nessun campo opzionale "che tanto di solito c'è". Se un valore può mancare, il
  tipo è `float | None` e chi legge deve gestirlo.
- **Nessun default fisico dentro il codice.** I default stanno nei YAML in
  `config/`, sotto controllo di versione, dove si vedono nel diff.

### 4.1 `DesignVector` — l'unica cosa che l'ottimizzatore muove

```python
@dataclass(frozen=True)
class DesignVector:
    schema_version: str          # "zefiro-schema-0.1.0"
    values: Mapping[str, float]  # nome → valore, SI
    bounds: Mapping[str, tuple[float, float]]
```

Invariante: `set(values) == set(bounds)` e ogni valore è dentro i suoi bounds.
Violarlo è un errore, non un warning: un ottimizzatore che esce dai bounds ha un
bug, e mascherarlo con un clamp rende il bug invisibile.

**Parametri liberi (12).** Scelti col criterio: un parametro è libero se e solo se
(a) è fisicamente realizzabile in modo indipendente dagli altri e (b) esiste un
motivo per cui il suo ottimo non è ovvio a priori.

| # | nome | simbolo | unità | intervallo proposto | perché è libero |
|---|---|---|---|---|---|
| 1 | `p_c` | p_c | Pa | 3e5 … 6e5 | trade-off: alza c\* e ε, ma erode il Δp iniettore |
| 2 | `phi_core` | φ_core | – | 0.7 … 1.15 | l'ottimo non è φ=1: T_ad max è a φ≈1.05, ma la parete preferisce magro |
| 3 | `f_film` | f_film | – | 0.0 … 0.35 | sottrae combustibile al core (perde spinta) per salvare la parete |
| 4 | `Dc_over_Dt` | D_c/D_t | – | 2.0 … 5.0 | rapporto di contrazione: alto = bassa velocità in camera = più tempo di residenza, ma più massa e più area da raffreddare |
| 5 | `Lc_over_Dc` | L_c/D_c | – | 0.8 … 4.0 | lunghezza camera in forma adimensionale: governa τ_res |
| 6 | `N_inj` | N | – (intero) | 6 … 24 | più elementi = mixing più fine, ma Δp e fabbricabilità peggiorano |
| 7 | `theta_swirl` | θ_s | rad | 0 … 1.05 (0–60°) | swirl: accorcia la fiamma, costa pressione totale |
| 8 | `d_ox_ratio` | d_ox/D_c | – | 0.02 … 0.12 | area di passaggio aria per elemento → velocità di iniezione |
| 9 | `fuel_vel_ratio` | u_f/u_ox | – | 0.3 … 3.0 | rapporto di quantità di moto: è **il** parametro che governa il mixing |
| 10 | `conv_half_angle` | β | rad | 0.35 … 0.79 (20–45°) | angolo del convergente: perdite vs lunghezza |
| 11 | `plug_trunc` | – | – | 0.2 … 1.0 | frazione di plug mantenuta (1.0 = plug completo). Tronca massa e lunghezza pagando prestazione |
| 12 | `t_wall` | t | m | 8e-4 … 4e-3 | spessore parete: vincolo strutturale vs inerzia termica vs SLM |

> `N_inj` è intero: l'ottimizzatore lo tratta come variabile discreta (pymoo lo
> supporta nativamente). Non arrotondare un float — produrrebbe due `x` diversi
> con lo stesso `run_id`, che rompe §5.

**Parametri derivati.** Calcolati da `derive(x, operating_point)`, mai forniti
dall'esterno, mai salvati come input. Se un derivato ti serve come input, allora
la parametrizzazione è sbagliata e va cambiata qui, non aggirata nel codice.

| nome | da cosa deriva |
|---|---|
| `mdot_air` | dalla gola: `ṁ = p_c·A_t/c*`, con `c*` da L0 → **iterazione interna** (§4.5) |
| `mdot_fuel` | `mdot_air/AFR_st · φ_core / (1 - f_film)` |
| `A_t`, `D_t` | dalla portata e da `c*` |
| `R_lip` | raggio di labbro dell'aerospike, da `A_t` e dalla geometria anulare |
| `epsilon` | `A_e/A_t` da `p_c/p_amb` isentropico con γ dei prodotti |
| `M_e`, `nu_e` | Mach di uscita e angolo di Prandtl-Meyer di progetto |
| `plug_contour` | array (x, r) dal metodo di Angelino |
| `D_c`, `L_c` | da `Dc_over_Dt`, `Lc_over_Dc`, `D_t` |
| `V_c`, `L_star` | volume di camera e lunghezza caratteristica `L* = V_c/A_t` |
| `tau_res` | `ρ_c V_c / ṁ_tot` |
| `A_ox_elem`, `A_f_elem` | aree di passaggio per elemento, da `d_ox_ratio` e `fuel_vel_ratio` |
| `dp_inj_ox`, `dp_inj_f` | Δp iniettore, da portata, area e `C_d` (`C_d` è un **dato di config**, §8) |
| `sector_angle` | `2π/N_inj` |

### 4.2 `OperatingPoint` — le condizioni al contorno

```python
@dataclass(frozen=True)
class OperatingPoint:
    schema_version: str
    p_amb: float          # Pa
    T_air_in: float       # K, a valle del regolatore
    T_fuel_in: float      # K
    p_air_supply: float   # Pa, monte regolatore  (1.0e6 = 10 bar)
    p_fuel_supply: float  # Pa                    (8.0e5 = 8 bar)
    mdot_air_max: float | None   # kg/s  ⟵ TODO: dato del compressore
    fuel: FuelSpec
    cd_injector_ox: float | None
    cd_injector_fuel: float | None
    burn_time: float | None      # s
```

```python
@dataclass(frozen=True)
class FuelSpec:
    schema_version: str
    composition: Mapping[str, float]  # frazioni MOLARI, es. {"C3H8": 1.0}
    phase_at_injection: str           # "gas" | "liquid"
    thermo_source: str                # nome del file meccanismo
```

`composition` è molare per una ragione precisa: le schede tecniche del GPL danno
percentuali volumetriche, che in fase gassosa **sono** frazioni molari. Convertire
a massa in lettura introdurrebbe un passaggio non tracciabile.

### 4.3 `GeometryParams` e `GeometryArtifact`

`GeometryParams` è il *record completo* (liberi + derivati) e viene serializzato
accanto al file CAD: è ciò che rende un STEP interpretabile fra sei mesi.

```python
@dataclass(frozen=True)
class GeometryArtifact:
    schema_version: str
    run_id: str
    step_path: Path
    stl_path: Path
    params: GeometryParams
    volume: float             # m³ del solido di parete
    wetted_area: float        # m² superficie bagnata dai gas
    is_valid_brep: bool       # esito OCCT BRepCheck
    is_watertight_mesh: bool  # ogni edge della mesh condiviso da esattamente 2 triangoli
    n_triangles: int
    mesh_tolerance: float     # m, tolleranza di tassellazione usata per l'STL
    sha256_step: str
    sha256_stl: str
```

Il water-tightness è verificato **due volte, con metodi indipendenti**: la
validità B-Rep secondo OCCT e la manifold-ness della tassellazione. La prima non
implica la seconda (un solido valido può tassellare male con tolleranza troppo
lasca) ed è esattamente quel caso a rompere il mesher a valle.

#### Contorno dell'aerospike — perché il metodo di Angelino e come si deriva

Non è una formula presa da una tabella: si ricava in cinque righe e il codice la
implementa esattamente così (`geometry/aerospike.py`).

Ipotesi: espansione **esterna** con ventaglio di Prandtl-Meyer **centrato sul
labbro** (raggio `R_e`, sull'asse `x = 0`). Ogni caratteristica che parte dal
labbro è retta, e lungo di essa il flusso è uniforme a Mach `M`. Il plug è la
linea di corrente che chiude il campo.

1. Deflessione del flusso a Mach `M`, misurata dall'asse:
   `θ(M) = ν(M_e) − ν(M)`, con ν angolo di Prandtl-Meyer.
   (Al labbro `M = 1`, `ν = 0`, quindi il flusso parte inclinato di `ν(M_e)`;
   all'uscita `M = M_e`, `θ = 0`, cioè assiale. Come deve essere.)
2. Angolo di Mach: `μ(M) = arcsin(1/M)`.
3. La caratteristica fa quindi con l'asse l'angolo `α = θ + μ`.
4. Il punto di plug su quella caratteristica è a distanza `ℓ` dal labbro:
   `x_p = ℓ cos α`, `r_p = R_e − ℓ sin α`.
5. **Conservazione della massa attraverso la caratteristica.** La superficie di
   rivoluzione generata dal segmento labbro→plug ha area laterale
   `S = π (R_e + r_p) ℓ`; la sua proiezione normale alla direzione del flusso
   vale `S sin μ`. Uguagliando a `(A/A*)(M) · A_t` con `A_t = π R_e²/ε`:

   ```
   π (R_e² − r_p²) sin μ / sin α  =  (A/A*)(M) · π R_e² / ε
   ```

   da cui, con `sin μ = 1/M`,

   ```
   r_p(M) / R_e = sqrt( 1 − M · (A/A*)(M) · sin α(M) / ε )
   ℓ            = (R_e − r_p) / sin α
   x_p          = ℓ cos α
   ```

Il contorno si ottiene marciando `M` da 1 a `M_e`. Nessun dato tabellato entra
nel calcolo: solo γ e il rapporto di pressione.

**Dove questa derivazione si rompe** (da tenere in mente quando la CFD non
tornerà):
- assume caratteristiche rette e flusso uniforme su di esse — vero solo se
  l'espansione è veramente centrata, cioè se il labbro è affilato rispetto alla
  scala della gola. Con `t_wall` di 1–4 mm e una gola piccola, **il labbro non è
  affilato**: aspettarsi uno scostamento reale;
- ignora completamente lo strato limite. Su un plug lungo e a bassa `Re`,
  lo spessore di spostamento non è trascurabile: il contorno andrà corretto in
  fase 4 con `δ*` preso dalla RANS;
- ignora la scia di base se `plug_trunc < 1`. La pressione di base è un
  parametro che **solo la CFD** può dare: `l0.nozzle` la tratta come incognita e
  restituisce una prestazione *limite superiore* per il plug troncato, marcata
  come tale nel risultato (`base_pressure_model = "none"`).

### 4.4 `L0Result`

```python
@dataclass(frozen=True)
class L0Result:
    schema_version: str
    # miscela
    phi_core: float; phi_global: float; AFR_stoich: float
    mdot_air: float; mdot_fuel_core: float; mdot_fuel_film: float
    # termochimica
    T_ad: float                       # K, equilibrio HP
    X_eq: Mapping[str, float]         # frazioni molari (specie > 1e-6)
    gamma_c: float; MW_c: float; cp_c: float
    # prestazione 1D
    c_star: float                     # m/s
    M_e: float; epsilon: float; C_F: float
    Isp_s: float                      # s, su massa TOTALE (aria + GPL)
    Isp_fuel_s: float                 # s, sul solo GPL (metrica air-breathing)
    thrust: float                     # N
    # tempi caratteristici
    tau_res: float | None             # s
    tau_chem: float | None            # s, None se il meccanismo non è disponibile
    damkohler: float | None
    # diagnostica onesta
    warnings: tuple[str, ...]
    assumptions: tuple[str, ...]      # ogni ipotesi attiva, in chiaro
```

`Isp_s` **e** `Isp_fuel_s` sono entrambe presenti perché con l'aria a bordo la
prima è la metrica di sistema e la seconda è la metrica di combustione: usarne
una sola sarebbe fuorviante in un verso o nell'altro.

`assumptions` non è decorativo: è il meccanismo con cui un risultato porta con sé
il proprio dominio di validità. Un `L0Result` con
`"frozen expansion"` dentro `assumptions` non è confrontabile con uno che ha
`"shifting equilibrium"`, e il codice di analisi lo può controllare.

### 4.5 Il punto fisso fra geometria e L0

C'è un accoppiamento circolare: `A_t` serve per la geometria, ma dipende da `c*`,
che dipende dalla termochimica, che dipende da `p_c` e φ — non da `A_t`.

Il ciclo si **spezza**, non si itera: `c*` dipende solo da `(p_c, φ, composizione,
T_in)`. Quindi l'ordine corretto è:

```
(p_c, φ, fuel, T_in) ──► l0.equilibrium ──► c*, γ, MW
                                             │
        ṁ_aria (da config)  ────────────────►├─► A_t = ṁ_tot · c* / p_c
                                             └─► ε, M_e ──► contorno plug
```

Nessuna iterazione, nessun solver non lineare, nessuna tolleranza da tarare.
Questo è il motivo per cui `derive()` è una funzione pura e deterministica, e per
cui il modulo geometria può essere testato senza CFD.

### 4.6 Contratti L1 (definiti ora, implementati in fase 3–4)

```python
@dataclass(frozen=True)
class MeshArtifact:
    schema_version: str; run_id: str
    msh_path: Path; n_cells: int
    boundary_names: tuple[str, ...]   # nomi CANONICI, §6.2
    is_periodic: bool; sector_angle: float
    max_non_orthogonality: float; max_skewness: float
    sha256: str

@dataclass(frozen=True)
class CFDResult:
    schema_version: str; run_id: str; mesh_id: str
    converged: bool; residual_final: Mapping[str, float]
    eta_combustion: float             # frazione di ΔH rilasciata a valle uscita camera
    p_wall: Path                      # campo di parete, Parquet (patch_id, x, y, z, p)
    q_wall: Path                      # idem, flusso termico
    thrust_cfd: float
    solver_version: str

@dataclass(frozen=True)
class FEMResult:
    schema_version: str; run_id: str; cfd_id: str
    T_wall_max: float; von_mises_max: float
    margin_yield: float               # σ_amm/σ_max − 1
    margin_temp: float                # T_amm/T_max − 1
    solver_version: str
```

L'accoppiamento CFD→FEM è **one-way** e passa da due file Parquet con schema
fisso `(patch_id: str, x: f8, y: f8, z: f8, value: f8)`. Il mapping da mesh
fluida a mesh strutturale è `l1.mapping`, e deve conservare l'integrale: la
verifica è `|∫q̇ dA|_CFD − |∫q̇ dA|_FEM| / |∫q̇ dA|_CFD < 1e-3`. Se non conserva,
il FEM sta risolvendo un problema diverso da quello che la CFD ha posto, e i
margini strutturali non significano nulla.

### 4.7 `Objectives` — vettoriale, mai scalarizzato

```python
@dataclass(frozen=True)
class Objectives:
    schema_version: str; run_id: str; fidelity: str   # "L0" | "L1" | "L2"
    f: Mapping[str, float]   # da MINIMIZZARE, sempre. Segno invertito a monte.
    g: Mapping[str, float]   # vincoli: FATTIBILE se g <= 0, sempre.
    feasible: bool
    source: Mapping[str, str]  # obiettivo → modulo che l'ha prodotto
```

Convenzione rigida: **minimizzare `f`, `g ≤ 0` fattibile**. Un obiettivo da
massimizzare entra come `-valore`. Chi produce l'obiettivo fa l'inversione, non
l'ottimizzatore: così cambiare ottimizzatore non cambia il significato dei numeri.

`source` esiste perché in un fronte di Pareto misto L0/L1 devi poter sapere quale
punto è stato valutato con quale fedeltà senza consultare il database.

---

## 5. Versionamento di una run

### 5.1 `run_id`

```
run_id = blake2b_16( canonical_json({
    "schema_version": SCHEMA_VERSION,
    "design":         x.values,          # chiavi ordinate, float via repr()
    "operating":      operating_point,   # idem
    "code_rev":       git_describe(),    # commit del repo, "-dirty" se sporco
}) )
```

Proprietà volute:
- **deterministico**: stesso input → stesso id, su qualunque macchina;
- **sensibile al codice**: cambiare un modulo cambia gli id, quindi non si
  confrontano mai risultati prodotti da versioni diverse credendoli uguali;
- **rifiuta il dirty**: se il working tree è sporco, `run_id` porta il suffisso
  `-dirty` e il record va in una tabella separata `runs_dirty`. Le run sporche
  servono per esplorare, non per concludere.

Formato directory: `runs/<YYYYMMDD>/<run_id>/` con dentro
`design.json`, `geometry.step`, `geometry.stl`, `l0.json`, `mesh/`, `cfd/`,
`fem/`, `objectives.json`, `env.json`, `log.txt`.

La data nel path è solo per comodità umana; **l'identità è il `run_id`**, e il
database ha un indice unico su di esso.

### 5.2 `env_fingerprint`

In `env.json`: versione Python, versioni di cantera / OCCT / numpy / gmsh /
OpenFOAM / CalculiX, `platform.uname()`, e SHA256 del meccanismo cinetico. È il
campo che rende falsificabile l'affermazione "riproducibile": due run con lo
stesso `run_id` e `env_fingerprint` diverso che danno numeri diversi non sono un
mistero, sono una diagnosi.

### 5.3 Cosa sta in Git e cosa no

| In Git | Fuori da Git |
|---|---|
| `config/**` (tutti gli input) | `runs/**` (tutti gli output) |
| `src/**`, `tests/**`, `scripts/**` | mesh, campi, STEP/STL generati |
| `environment.yml`, `conda-lock` | il file del meccanismo cinetico (solo lo SHA256 sta in Git) |
| `docs/**`, `ROADMAP.md` | il database SQLite |

Il principio: **in Git ciò che è scritto da un umano, fuori ciò che è
riproducibile da una macchina.** Se un output non è riproducibile, il problema è
l'output, non Git.

### 5.4 Database

SQLite (`runs.db`) come indice, Parquet per i campi. Tabella principale `runs`:
colonne scalari (una per ogni voce di `x`, ogni voce di `L0Result`, ogni voce di
`Objectives`) + `run_id` PK + `fidelity` + `env_hash` + timestamp.

Motivo: le query che vuoi fare sono
`SELECT run_id FROM runs WHERE plug_trunc BETWEEN 0.2 AND 0.4 AND margin_yield > 2`,
cioè SQL su scalari. Una colonna JSON renderebbe questa query lenta e non
indicizzabile. I campi spaziali (grandi, mai filtrati per valore) stanno in
Parquet, referenziati per path.

---

## 6. Convenzioni di naming

### 6.1 Codice
- moduli e funzioni `snake_case`; classi `PascalCase`; costanti `UPPER_SNAKE`.
- **Le variabili fisiche portano il nome del simbolo, non una parafrasi**:
  `p_c`, `T_ad`, `mdot_air`, `c_star`, `A_t`. Serve a rendere il codice
  confrontabile con la carta su cui hai fatto il conto.
- Suffisso di unità **solo** quando non è SI, e questo può accadere solo al
  confine I/O: `p_c_bar`, `Isp_s`. Dentro il codice il suffisso non c'è perché
  l'unità è SI per contratto.

### 6.2 Patch/boundary — nomi canonici
Vincolati qui perché sono l'interfaccia fra Gmsh, OpenFOAM e CalculiX, e un
refuso lì produce una condizione al contorno silenziosamente sbagliata:

```
inlet_air, inlet_fuel_core, inlet_fuel_film,
wall_chamber, wall_convergent, wall_throat, wall_plug, wall_cowl,
outlet_far, axis, periodic_a, periodic_b
```

`l1.mesh` valida che l'insieme dei nomi prodotti da Gmsh sia **esattamente**
questo insieme, e fallisce altrimenti.

### 6.3 File
`<run_id>_<contenuto>.<ext>`, tutto minuscolo, mai spazi.
Il meccanismo cinetico si chiama `<nome>_<sha256[:8]>.yaml`: il nome del file
porta il suo hash, così è impossibile usarne uno diverso senza accorgersene.

---

## 7. Dati fisici: provenienza e assenze

Regola del progetto: **nessun dato inventato**. Ogni numero ha una provenienza
tracciabile o è un `None` che fa fallire il codice.

### 7.1 Cosa è calcolato (non assunto)
- `AFR_stoich`, `T_ad`, composizione all'equilibrio, `γ`, `MW`, `c*`, `C_F`:
  calcolati da Cantera a partire dai polinomi NASA del meccanismo dichiarato.
- Contorno dell'aerospike: derivato analiticamente (§4.3).
- Caso di verifica del modulo L0: propano/aria stechiometrico, 300 K, 1 atm →
  `AFR = 15.671`, `T_ad = 2267.4 K`. Entrambi verificabili a mano
  (§ `tests/test_l0_reference.py`, che contiene il conto in commento).

### 7.2 Sorgenti termodinamiche e cinetiche
- **Equilibrio (default): `gri30.yaml`**, distribuito con Cantera. Contiene
  `C3H8` e un set di prodotti C/H/O/N completo. Nota importante: GRI-Mech 3.0 è
  un meccanismo *cinetico* tarato su metano, ma qui se ne usano **solo i dati
  termodinamici e l'elenco di specie**, e l'equilibrio non usa alcuna costante
  di velocità. È quindi un uso legittimo — a due condizioni, entrambe
  verificate dal codice: che il combustibile sia presente fra le specie, e che
  la miscela non sia troppo ricca (per `φ > 1.5` le specie C2+ e i precursori di
  fuliggine mancanti falsano l'equilibrio → **warning esplicito**).
- **Butano: assente da `gri30.yaml` e dal database NASA distribuito con
  Cantera.** Se `FuelSpec.composition` contiene `C4H10`, il codice **solleva
  `MissingThermoData`** e indica di eseguire `scripts/fetch_mechanism.py`. Non
  esiste alcun percorso in cui il butano viene sostituito con propano di
  nascosto.
- **Cinetica (τ_ign, S_L): San Diego mech**, non vendorizzato.
  `scripts/fetch_mechanism.py` lo scarica dalla fonte ufficiale UCSD, lo
  converte con `ck2yaml`, ne registra lo SHA256 in `config/mechanisms.lock` e
  rinomina il file col proprio hash. Finché il file non c'è, `l0.kinetics`
  solleva `MissingThermoData` e `L0Result.tau_chem` resta `None`.

---

## 8. L'impianto e i dati che ancora mancano

### 8.1 Che cosa l'impianto permette

Dati noti: **compressore 3 cavalli, serbatoio 100 L, max 10 bar**; **bombola GPL
da barbecue, indicata a 8 bar**. Da questi si ricava molto per via puramente
termodinamica (`zefiro.feed`, `scripts/plant_report.py`), senza stimare nulla.

**Portata continua.** Il lavoro minimo di compressione da 1 a 10 bar a 20 °C è
quello isotermo, `w = R T ln(10) = 192.7 kJ/kg`. Con 3 HP all'albero
(2237 W) questo dà un **limite superiore assoluto di 11.6 g/s**, irraggiungibile
perché richiederebbe raffreddamento perfetto. Una compressione adiabatica
monostadio ideale dà 8.2 g/s. Una macchina a pistoni reale sta sotto: attendersi
**4–7 g/s**. Il numero vero è sulla targhetta come *aria resa* in l/min
(300 l/min ⇒ 6.0 g/s).

**Portata a raffica — è qui che il serbatoio cambia il progetto.** 100 L a
10 bar contengono 1188 g d'aria. Scaricando fino a 6 bar se ne estraggono
**359 g** (limite adiabatico; 475 g nel limite isotermo, ottimistico perché una
raffica di pochi secondi non ha tempo di scambiare calore con le pareti).
Quindi:

| ṁ_aria | durata raffica | spinta a p_c = 4 bar | D_gola | Ø foro GPL | Ø foro film |
|---|---|---|---|---|---|
| 6 g/s (continuo) | illimitata | 8.2 N | 5.1 mm | 0.20 mm | 0.15 mm |
| 20 g/s | 26 s | 27 N | 9.2 mm | 0.36 mm | 0.27 mm |
| 50 g/s | 8.2 s | 68 N | 14.6 mm | 0.57 mm | 0.42 mm |
| 100 g/s | 3.8 s | 137 N | 20.7 mm | 0.80 mm | 0.60 mm |

La colonna che decide non è la spinta, è **il diametro dei fori**. In
funzionamento continuo iniettori e fori di film cooling scendono sotto i due
decimi di millimetro, cioè sotto qualunque limite SLM ragionevole: *il motore
non sarebbe stampabile*. A partire da ~50 g/s le dimensioni entrano in un
intervallo fabbricabile. **È il serbatoio, non il compressore, a rendere Zefiro
costruibile.**

### 8.2 La bombola: la pressione misura la composizione

La pressione di una bombola **non dipende da quanto liquido è rimasto**: finché
c'è liquido, dipende solo da temperatura e composizione. Per legge di Raoult su
una miscela binaria propano/butano:

```
p = x_C3H8 · p_sat,C3H8(T) + (1 − x_C3H8) · p_sat,C4H10(T)
```

e la relazione si inverte (`feed.composition_from_pressure`). Con p_sat da EOS
di riferimento (propano: Lemmon 2009; n-butano: Bücker & Wagner 2006), 8 bar
assoluti implicano:

| T bombola | x_propano dedotto |
|---|---|
| 10 °C | 1.34 → **impossibile** |
| 15 °C | 1.12 → **impossibile** |
| 20 °C | 0.94 |
| 25 °C | 0.79 |
| 30 °C | 0.65 |

Cioè: **un manometro e un termometro letti insieme sostituiscono un'analisi di
laboratorio.** Finché non hai quella coppia, la composizione resta un TODO.

### 8.3 L'autorefrigerazione non è il problema

Prelevare vapore richiede calore latente, che viene dalla bombola stessa:
`(m_liq c_liq + m_guscio c_guscio) dT/dt = −ṁ h_fg(T)`. Con 7.5 g/s di propano
(cioè ṁ_aria = 100 g/s) da una bombola con 6 kg di liquido residuo, la potenza
sottratta è 2.6 kW e la bombola perde **1.3 K in 10 secondi**: la pressione
scende da 8.36 a 8.08 bar. Trascurabile.

Il collo di bottiglia è invece **il riduttore**: a 100 g/s d'aria servono
27 kg/h di GPL. Un riduttore da barbecue standard eroga 30 mbar e circa 1 kg/h —
non è utilizzabile qui, dove servono 5–8 bar. Serve un riduttore di alta
pressione, e la sua portata massima è un vincolo hard sul punto operativo.

### 8.4 TODO ancora aperti

Stanno come `null` in `config/`; `zefiro.config.load_operating_point()`
fallisce con `MissingDatum` elencandoli tutti insieme.

| # | Dato | Perché blocca | Come ottenerlo |
|---|---|---|---|
| A | Aria resa del compressore (l/min) | Fissa la portata continua e quindi la ricarica durante la raffica | targhetta. Attenzione: alcuni dichiarano l'aria *aspirata*, che è 1.3–1.6× l'aria resa |
| B | Pressione minima utile a valle del serbatoio | Fissa la massa utilizzabile e quindi la durata della raffica | dipende dal riduttore che monti |
| C | Pressione **e** temperatura della bombola, misurate insieme | Identifica la composizione (§8.2). Se c'è butano serve un meccanismo con C4H10 | manometro + termometro, bombola in equilibrio termico |
| D | Prelievo gassoso o liquido | Se liquido serve un modello di flash a monte di L0 | com'è fatta la presa |
| E | Tara e peso totale della bombola | Dà la massa di liquido, che entra nel calcolo di autorefrigerazione | la tara è stampigliata sul collare, il peso con una bilancia |
| F | Portata massima del riduttore GPL (kg/h) | È il vincolo hard sul punto operativo a raffica (§8.3) | datasheet del riduttore |
| G | `T_air_in`, `T_fuel_in` a valle dei riduttori | L'espansione raffredda: non è temperatura ambiente | misura |
| H | Punto di progetto: continuo o raffica, e a che ṁ_aria | Dimensiona tutto il motore | **decisione tua** |
| I | `C_d` degli iniettori | Determina Δp e se p_c è raggiungibile | geometria + taratura a freddo |
| J | AISI 316L da SLM: σ_y(T), E(T), α(T), k(T), ρ, σ_ammissibile | È il FEM. I dati SLM sono anisotropi e diversi dal 316L laminato | prove sui tuoi provini, o datasheet per i tuoi parametri di processo |
| K | Spessore minimo e angolo di overhang della tua macchina SLM | Vincoli geometrici hard: §8.1 mostra che sono attivi | chi gestisce la stampante |
| L | Sistema di accensione | Determina se modellare l'ignizione o solo la combustione stabilizzata | decisione da fare |

## 9. Storia delle versioni di schema

| Versione | Data | Cambiamento |
|---|---|---|
| `zefiro-schema-0.1.0` | 2026-08-26 | Prima definizione. |
| — | 2026-08-26 | Aggiunto `zefiro.feed` (impianto di alimentazione). Nessun contratto scambiato modificato, quindi `SCHEMA_VERSION` invariata. |
