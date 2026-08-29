# Avviare Zefiro

Testato su un'installazione pulita: **installazione 1 min 15 s, test 1 min 3 s.**

---

## 1. Dove metti le cose

Due posti diversi, e la distinzione conta.

| | dove | perché |
|---|---|---|
| **il repo** | `I:\AA_ENGINE`, cioè `/mnt/i/AA_ENGINE` da WSL | è codice, è piccolo, e ti serve vederlo da Windows |
| **gli artefatti** (`runs/`) | nel filesystem di WSL, es. `~/zefiro-runs` | su `/mnt/...` l'I/O su file piccoli e numerosi è 5–10 volte più lento, ed è esattamente quello che produce una mesh |

La seconda si imposta con la variabile d'ambiente `ZEFIRO_RUNS`. Il punto 3 la mette
nel `.bashrc` una volta per tutte.

---

## 2. Prerequisiti

Se WSL2 non c'è ancora, da **PowerShell come amministratore**:

```powershell
wsl --install -d Ubuntu
```

poi riavvia e apri Ubuntu. Da qui in avanti tutti i comandi sono **dentro WSL**.

```bash
sudo apt update
sudo apt install -y python3 python3-venv python3-pip git build-essential
```

> **Perché WSL e non Windows.** OpenFOAM e CalculiX girano solo lì, e far
> attraversare a ogni valutazione il confine Windows/WSL con mesh da centinaia di
> MB è lento e introduce tre sorgenti di errore silenzioso (traduzione dei path,
> permessi POSIX inventati da drvfs, line ending). In più Cantera e OCCT hanno
> binari diversi sui due sistemi, e con un solo target di compilazione
> l'affermazione «stesso input, stesso output» resta verificabile.
> Il ragionamento completo è in `docs/architettura.md` §2.2.

---

## 3. Installazione

```bash
cd /mnt/i/AA_ENGINE

# ambiente isolato: non installare mai nel Python di sistema
python3 -m venv .venv
source .venv/bin/activate

pip install --upgrade pip
pip install -e ".[dev]"
```

L'ultimo comando scarica circa 500 MB (OCCT, il kernel CAD) e ci mette un paio di
minuti la prima volta.

Poi, una volta sola:

```bash
mkdir -p ~/zefiro-runs
echo 'export ZEFIRO_RUNS=$HOME/zefiro-runs' >> ~/.bashrc
source ~/.bashrc

git init && git add -A && git commit -m "stato iniziale"
```

> **Il `git init` non è opzionale.** L'identità di una run (`run_id`) include la
> revisione del codice. Senza un repository git, cambiare un modulo non cambia il
> `run_id`, e ti ritroveresti a confrontare risultati prodotti da versioni diverse
> credendoli uguali. `scripts/doctor.py` lo segnala come bloccante apposta.

---

## 4. Verifica

```bash
python scripts/doctor.py
```

Elenca in un colpo solo che cosa funziona e che cosa manca, ed esce con codice 1 se
c'è qualcosa di bloccante. Le voci marcate `nota` non sono problemi: Gmsh, OpenFOAM
e CalculiX non servono fino alla fase 3.

```bash
pytest -q
```

Devono passare tutti. Ci mette circa un minuto: la maggior parte del tempo se ne va
in Cantera e nelle booleane di OCCT.

---

## 5. Primo giro

```bash
# che motore permette il tuo banco
python scripts/plant_report.py --fad 300 --bottle-T 20

# valutazione termochimica del punto di progetto
zefiro-l0

# geometria: STEP + STL
zefiro-geometry
zefiro-geometry --sector          # settore periodico 1/N per la CFD

# piano sperimentale: 500 punti nel database
python scripts/sweep_l0.py --n 500 --seed 0
sqlite3 $ZEFIRO_RUNS/runs.db "SELECT run_id, thrust, Isp_s FROM runs
                              WHERE feasible=1 ORDER BY thrust DESC LIMIT 10"
```

### Se `zefiro-l0` si ferma dicendo che mancano dei dati, sta funzionando

```
MissingDatum: Dati operativi mancanti (vedi docs/architettura.md sezione 8):
T_air_in, T_fuel_in
```

I valori che non sono stati misurati stanno come `null` in `config/` e il codice si
ferma elencandoli tutti insieme, invece di riempirli con numeri plausibili. Per
provare il software prima di avere le misure, copia il file e metti dei valori
**dichiaratamente provvisori**:

```bash
cp config/operating_point.yaml /tmp/prova.yaml
# apri /tmp/prova.yaml e sostituisci i null che ti servono
zefiro-l0 --operating /tmp/prova.yaml
```

Così i valori finti restano fuori dal repo e non finiscono in una run che poi
crederai vera.

---

## 5b. Interfaccia web

```bash
pip install -e ".[web]"
python scripts/serve.py
```

Apre il browser su `http://127.0.0.1:8000`. Se i dati operativi in `config/`
sono ancora `null`, la pagina non lancia una valutazione destinata a fallire:
evidenzia in ambra i campi da misurare e ti lascia inserire valori provvisori,
che restano in memoria e **non toccano i file**.

I lavori lunghi (geometria, piano sperimentale) vanno in coda con barra di
avanzamento invece di bloccare la richiesta. Un solo lavoro alla volta, di
proposito: OCCT e i solutori numerici non sono thread-safe e sono limitati
dalla memoria, non dalla CPU.

## 6. Se qualcosa non va

| sintomo | causa | rimedio |
|---|---|---|
| `pip install` fallisce su `cadquery-ocp` | rete lenta, il wheel è grosso | `pip install --timeout 300 --retries 5 -e ".[dev]"` |
| `ls /mnt/i` è vuoto | il disco `I:` non è montato in WSL | `sudo mkdir -p /mnt/i && sudo mount -t drvfs I: /mnt/i` |
| `zefiro-l0: command not found` | ambiente non attivo | `source .venv/bin/activate` |
| `ImportError: OCP` | manca una libreria di sistema | `sudo apt install -y libgl1 libglu1-mesa libxrender1 libxi6` |
| il `run_id` finisce in `-dirty` | modifiche non committate | normale mentre sviluppi. Quelle run vanno in `runs_dirty` e non compaiono nelle query di sintesi: servono a esplorare, non a concludere |
| i test sono lentissimi | `ZEFIRO_RUNS` punta su `/mnt/...` | spostalo nel filesystem di WSL |
| l'interfaccia web dice `INTERFACCIA_ASSENTE` | manca la cartella `web/` | verifica di aver estratto tutto il repo |
| il 3D resta vuoto | il browser non ha WebGL | aggiorna il browser; il resto della pagina funziona lo stesso |
| da Windows non apri `localhost:8000` | WSL2 di norma inoltra la porta da solo; se non lo fa, lancia con `--host 0.0.0.0` e usa l'IP di WSL (`hostname -I`) | |

---

## 7. Più avanti: i solutori esterni

Servono dalla **fase 3**, non adesso. Quando ci arriviamo:

```bash
sudo apt install -y gmsh calculix-ccx
# OpenFOAM: repository dedicato, vedi openfoam.org
```

A quel punto conviene passare da `venv` a `conda` con `environment.yml`, perché
gmsh e pymoo si gestiscono meglio da conda-forge e serve un lockfile che descriva
l'intero ambiente numerico. Finché il layer è solo Python, `venv` è più semplice e
fa esattamente lo stesso.

---

## Nota su Windows

Il layer Python **funziona anche su Windows** — Cantera, build123d e CoolProp hanno
i wheel — quindi se ti fa comodo aprire il repo in un IDE su Windows e sperimentare,
fallo. Regola sola: **i risultati prodotti da Windows non entrano nel database
delle run.** Il campo `env_fingerprint` in `env.json` serve proprio a rendere questa
regola verificabile invece che sperata.

---

## 7. Far girare l'ottimizzazione (fase 5)

Tre comandi, in ordine. Tutti da dentro WSL, con l'ambiente attivo.

### 7.1 Prima: verificare che gli obiettivi abbiano ancora senso

```bash
python scripts/objective_screening.py -n 400 --seed 20260829
```

Campiona 400 punti del box e stampa la matrice di correlazione di rango fra le
metriche candidate. **Va rifatto ogni volta che cambi i bound o il punto
operativo**: due obiettivi che oggi sono indipendenti possono diventare
degeneri con un box diverso, e un fronte degenere non lo vedi guardandolo — lo
vedi solo qui. Costa ~2 minuti (la parte lenta è tabulare il tempo chimico).

### 7.2 Poi: il fronte di Pareto

```bash
python scripts/optimize_l0.py --pop 200 --gen 250 --seed 101 --seeds 3 --jobs 23
```

`--seeds 3` ripete la run con tre semi diversi e confronta gli ipervolumi. Non
è un lusso: NSGA-II è stocastico, e **un fronte da un solo seme è un campione,
non un risultato**. Se la dispersione supera il 5 % lo script lo dice, e la
risposta giusta è alzare `--gen`, non fidarsi lo stesso.

Costo sulla tua macchina: ~13 ms per valutazione su 23 processi, cioè **~30 s
per seme** a 50 000 valutazioni. Puoi permetterti `--pop 400 --gen 600 --seeds
10` (2.4 milioni di valutazioni) in una ventina di minuti. Il budget non è più
il vincolo: usalo per i semi, non solo per le generazioni.

### 7.3 Infine: leggere il fronte

```bash
python scripts/pareto_report.py runs/pareto/pareto_seed10{1,2,3}.json \
       --plot runs/pareto/fronte.png
```

Un fronte è un elenco di compromessi, non una risposta. Il report separa:

* **gli estremi** e cosa costa passare dall'uno all'altro;
* il **ginocchio**, cioè il punto più vicino all'ideale in norma normalizzata —
  un punto di partenza dichiarato, non «la» risposta;
* i parametri di **compromesso** contro quelli **decisi** (stesso valore su
  tutto il fronte: quelli non sono compromessi, sono conclusioni). Se un
  parametro deciso sta al bordo del box, **il bound è stretto e il vero ottimo
  è fuori**: va rimesso in discussione, non subìto;
* il **prezzo** di ciascun obiettivo. È la sezione che conta di più: se
  spostare un obiettivo di un ordine di grandezza costa il 2 % su un altro,
  quello non è un compromesso ma una decisione già presa.

### 7.4 E verificare un punto, in modo indipendente

```bash
python scripts/verify_pareto_point.py runs/pareto/pareto_seed101.json \
       --cad runs/pareto/cad
```

Riparte dal vettore `x` e **rifà tutto da capo** — termochimica, geometria, CAD
vero con OCCT, vincoli — senza riusare niente di ciò che l'ottimizzatore aveva
in memoria. Se l'ottimizzatore avesse un vincolo col segno sbagliato o leggesse
la colonna sbagliata, questo controllo lo vedrebbe; rileggere i suoi numeri no.
Esce con codice 1 se qualcosa non torna, quindi lo puoi mettere in uno script.
