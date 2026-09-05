# Far girare la CFD del mescolamento sul tuo PC

Il caso lo genera lo stesso modello che progetta il motore: le quote NON sono
scritte a mano nel caso, si chiedono a `progetta()`. Se il motore cambia, il
caso cambia con lui.

## 1. Installare OpenFOAM in WSL (una volta sola, ~5 minuti)

Serve la tua password perche' `apt` vuole i diritti di amministratore. Apri un
terminale WSL e incolla:

```bash
sudo apt-get update
sudo apt-get install -y --no-install-recommends software-properties-common curl
curl -s https://dl.openfoam.com/add-debian-repo.sh | sudo bash
sudo apt-get install -y openfoam2312-default
```

Se il repository di openfoam.com non risponde, l'alternativa e' il pacchetto
Debian, piu' vecchio ma sufficiente:

```bash
sudo apt-get install -y openfoam
```

Verifica:

```bash
source /usr/lib/openfoam/openfoam2312/etc/bashrc   # oppure /usr/share/openfoam/etc/bashrc
simpleFoam -help | head -3
```

## 2. Generare il caso

Dal repo (Windows o WSL, e' python puro):

```bash
cd /mnt/i/AA_ENGINE
python3 scripts/cfd/caso_mescolamento.py \
    --caso runs/cfd/mescolamento \
    --fine 4.0e-5 --grossa 3.0e-4 \
    --tempo 2.0e-3 --scrittura 5.0e-5 \
    --proc 30
```

`--fine 4.0e-5` sono 40 micron al getto: **18 celle sul diametro del foro**,
che e' la risoluzione a cui la coppia di vortici controrotanti si vede davvero.
Dovrebbero uscire circa 2-3 milioni di celle. Sulla mia macchina non ci
starebbero in memoria; sui tuoi 23 GB si'.

## 3. Lanciare

```bash
source /usr/share/openfoam/etc/bashrc
cd /mnt/i/AA_ENGINE/runs/cfd
gmshToFoam iniettore.msh -case mescolamento
python3 ../../scripts/cfd/patch_bordi.py mescolamento    # simmetrie e pareti
checkMesh -case mescolamento | tail -20                  # deve dire "Mesh OK"
cd mescolamento
decomposePar
mpirun -np 30 reactingFoam -parallel | tee log.foam
reconstructPar
```

**Stima:** con 2.5 milioni di celle su 30 processi, circa **un'ora e mezza**
per i 2 ms simulati. Sulla mia macchina (2 core) sarebbero tre giorni: e' per
questo che il tuo PC serve davvero.

Controlla ogni tanto che il numero di Courant resti sotto 1 e che i residui
scendano; se il passo temporale crolla, il caso sta divergendo e conviene
fermarlo invece di lasciarlo macinare.

## 4. Guardare i risultati

```bash
foamToVTK -case mescolamento
python3 /mnt/i/AA_ENGINE/scripts/cfd/analizza_mescolamento.py \
    --caso runs/cfd/mescolamento
```

Lo script misura la cosa per cui il caso esiste: **quanto e' ancora
disuniforme la frazione di massa di GPL, piano per piano lungo l'asse**, e la
confronta con i 3.31 mm che la correlazione di Holdeman promette. Produce anche
le immagini del piano meridiano, che sono quelle da guardare.

## Che cosa questo caso NON dice

- **Non brucia.** La chimica e' spenta di proposito: se mescolamento e
  combustione girassero insieme e il risultato fosse brutto, non si saprebbe
  quale dei due modelli incolpare.
- **E' RANS, non LES.** k-omega SST da' i valori medi; le strutture
  istantanee del getto in flusso trasversale sono piu' ricche di cosi'.
- **E' un settore di 90 gradi con simmetria sui fianchi.** Esatta se i quattro
  getti sono identici; se uno si otturasse, questo caso non lo vedrebbe.
