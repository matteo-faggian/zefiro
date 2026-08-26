# Zefiro

Software di progettazione e ottimizzazione di un combustore a GPL con aria
compressa come ossidante, alimentato a pressione (nessuna turbopompa), con
ugello **aerospike** e camera raffreddata a **film cooling**. Realizzazione
prevista in **SLM, AISI 316L**.

Il software è un **loop chiuso**: la geometria è input dei solver, non output.

```
parametri x → geometria CAD → mesh → CFD reattivo → FEM termo-strutturale
           → obiettivi/vincoli → ottimizzatore → (torna a x)
           → export STEP/STL solo a convergenza
```

- **Architettura e contratti dati**: [`docs/architettura.md`](docs/architettura.md)
- **Fasi e criteri di chiusura**: [`ROADMAP.md`](ROADMAP.md)

## Stato

| Blocco | Stato |
|---|---|
| Contratti dati (`schemas.py`) | implementato |
| Geometria parametrica → STEP/STL | implementato, testato |
| Aerospike (metodo di Angelino) | implementato, testato |
| L0 termochimica (Cantera) | implementato, testato su caso verificabile a mano |
| Impianto di alimentazione (`feed.py`) | implementato, testato |
| Identità e riproducibilità delle run | implementato, verificato byte per byte |
| Database delle run (SQLite + Parquet) | implementato, testato |
| DOE Latin Hypercube deterministico | implementato, testato |
| Obiettivi e vincoli L0 | implementato, testato |
| Carico termico e parete transitoria | implementato, verificato su soluzioni analitiche |
| Stime strutturali analitiche | implementato, testato |
| L1 (OpenFOAM / CalculiX) | stub con contratti fissati |
| Ottimizzazione e surrogato | stub con contratti fissati |

## Installazione

Target unico: **WSL2 / Ubuntu**. La motivazione è in
`docs/architettura.md` §2.2 — in breve: OpenFOAM e CalculiX girano solo lì, e
attraversare il confine Windows/WSL a ogni valutazione costa tempo e rompe la
riproducibilità binaria.

```bash
conda env create -f environment.yml
conda activate zefiro
pip install -e ".[dev]"
pytest -q
```

## Uso

```bash
zefiro-l0                    # valutazione termochimica
zefiro-geometry              # genera STEP + STL
zefiro-geometry --sector     # settore periodico 1/N per la CFD

python scripts/plant_report.py --fad 300 --bottle-T 20   # che motore permette il banco
python scripts/sweep_l0.py --n 500 --seed 0 --mdot-air 0.05   # DOE su L0 -> runs/runs.db
```

Il database si interroga in SQL:

```bash
sqlite3 runs/runs.db "SELECT run_id, thrust, Isp_s, p_c FROM runs
                      WHERE feasible = 1 ORDER BY thrust DESC LIMIT 10"
```

Entrambi falliscono con un elenco dei dati mancanti finché
`config/operating_point.yaml` contiene dei `null`. **È voluto**: vedi la regola
qui sotto.

## Regola del progetto

> **Nessun dato fisico inventato.** Ogni numero ha una provenienza tracciabile
> oppure è un `null` che fa fallire il codice con un messaggio che dice quale
> dato manca e chi lo sa.

I dati ancora mancanti sono elencati in `docs/architettura.md` §8.
