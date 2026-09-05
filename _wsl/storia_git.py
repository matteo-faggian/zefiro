#!/usr/bin/env python3
"""Ricostruisce la storia di Zefiro dai quattordici zip datati.

PERCHE' SI PUO' FARE ONESTAMENTE. Gli zip sono ISTANTANEE COMPLETE, non
consegne parziali: il conteggio dei file cresce in modo monotono da 48 a 118 e
`src/zefiro/l0/cycle.py` compare in tutti e quattordici. Se fossero stati
parziali, sostituire l'albero a ogni passo avrebbe fatto SPARIRE dei file a
ogni commit, e la storia avrebbe raccontato cancellazioni mai avvenute.

L'ORDINE E' QUELLO DELLE DATE, NON DEI NOMI, e non e' un dettaglio:
`zefiro_fase5.zip` e' del 29 agosto alle 18:06 e `zefiro_fase3.zip` delle
22:43 dello stesso giorno. Ordinare per nome metterebbe la fase 3 prima della
5 e la storia risulterebbe rovesciata proprio dove il progetto cambiava di
piu'.

DUE COSE CHE QUESTO SCRIPT NON PRETENDE DI ESSERE. Non e' la storia reale dei
salvataggi: e' la sequenza degli stati che sono stati archiviati, che e' tutto
quello che resta. E le date sono quelle di modifica dei file zip, cioe' quando
l'archivio e' stato scritto - non quando il codice e' stato pensato.
"""
from __future__ import annotations

import datetime
import os
import pathlib
import shutil
import subprocess
import zipfile

ARCHIVIO = pathlib.Path("/mnt/i/AA_ENGINE/(eliminare)/_da_cancellare")
MASTER = pathlib.Path("/mnt/i/AA_ENGINE")
REPO = pathlib.Path("/root/zefiro_git")

#: non entrano nel repo: rigenerabili, enormi, o non codice
ESCLUSI_DIR = {".git", "__pycache__", "runs", "(eliminare)", ".pytest_cache",
               ".ruff_cache", "_da_cancellare"}
ESCLUSI_EXT = {".stl", ".msh", ".step", ".parquet", ".db"}

NOME = "mat"
EMAIL = "gptmatandrew@gmail.com"


def git(*args, quando: datetime.datetime | None = None) -> str:
    amb = dict(os.environ)
    if quando is not None:
        stamp = quando.strftime("%Y-%m-%d %H:%M:%S")
        amb["GIT_AUTHOR_DATE"] = stamp
        amb["GIT_COMMITTER_DATE"] = stamp
    r = subprocess.run(["git", "-C", str(REPO), *args], env=amb,
                       capture_output=True, text=True)
    if r.returncode and "nothing to commit" not in r.stdout:
        raise RuntimeError(f"git {' '.join(args)}:\n{r.stdout}\n{r.stderr}")
    return r.stdout


def svuota_albero() -> None:
    """Toglie tutto tranne .git: ogni commit e' uno STATO, non una toppa."""
    for p in REPO.iterdir():
        if p.name == ".git":
            continue
        shutil.rmtree(p) if p.is_dir() else p.unlink()


def estrai(z: pathlib.Path) -> int:
    """Scompatta normalizzando la radice.

    Alcuni zip hanno tutto sotto `zefiro/`, altri hanno i file al primo
    livello. Senza normalizzare, la storia alternerebbe fra due alberi diversi
    e ogni commit risulterebbe una riscrittura totale del progetto.
    """
    with zipfile.ZipFile(z) as f:
        nomi = [n for n in f.namelist() if not n.endswith("/")]
        radici = {n.split("/")[0] for n in nomi}
        taglia = len("zefiro/") if radici == {"zefiro"} else 0
        n = 0
        for nome in nomi:
            rel = nome[taglia:]
            if not rel:
                continue
            parti = pathlib.PurePosixPath(rel).parts
            if any(p in ESCLUSI_DIR for p in parti):
                continue
            if pathlib.PurePosixPath(rel).suffix in ESCLUSI_EXT:
                continue
            dest = REPO / rel
            dest.parent.mkdir(parents=True, exist_ok=True)
            with f.open(nome) as src, open(dest, "wb") as out:
                shutil.copyfileobj(src, out)
            n += 1
    return n


def copia_stato_attuale() -> int:
    n = 0
    for src in MASTER.rglob("*"):
        rel = src.relative_to(MASTER)
        if any(p in ESCLUSI_DIR for p in rel.parts):
            continue
        if not src.is_file() or src.suffix in ESCLUSI_EXT:
            continue
        dest = REPO / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dest)
        n += 1
    return n


def main() -> int:
    if REPO.exists():
        shutil.rmtree(REPO)
    REPO.mkdir(parents=True)
    git("init", "-q", "-b", "main")
    git("config", "user.name", NOME)
    git("config", "user.email", EMAIL)

    gitignore = (MASTER / ".gitignore").read_text(encoding="utf-8", errors="replace")

    zips = sorted(ARCHIVIO.glob("*.zip"), key=lambda p: p.stat().st_mtime)
    print(f"{len(zips)} archivi da ricostruire\n")
    for z in zips:
        quando = datetime.datetime.fromtimestamp(z.stat().st_mtime)
        svuota_albero()
        n = estrai(z)
        #: il .gitignore e' SEMPRE quello di oggi, in tutti i commit. Cosi' le
        #: regole di esclusione non cambiano sotto i piedi lungo la storia e
        #: un file grosso non entra nel repo perche' in agosto non era ancora
        #: ignorato. E' l'unica cosa anacronistica, ed e' voluta.
        (REPO / ".gitignore").write_text(gitignore, encoding="utf-8")
        git("add", "-A")
        messaggio = (
            f"{z.stem}\n\n"
            f"Stato del progetto archiviato il "
            f"{quando:%d/%m/%Y alle %H:%M} ({n} file).\n\n"
            f"Commit RICOSTRUITO da {z.name}: il progetto non era sotto\n"
            f"controllo di versione e questi archivi datati sono l'unica\n"
            f"traccia rimasta della sua evoluzione. La data e' quella di\n"
            f"scrittura dell'archivio, non quella in cui il codice e' stato\n"
            f"pensato.")
        git("commit", "-q", "--allow-empty", "-m", messaggio, quando=quando)
        breve = git("log", "-1", "--format=%h").strip()
        print(f"  {quando:%Y-%m-%d %H:%M}  {breve}  {z.stem:22} {n:4d} file")

    #: --- lo stato di oggi ------------------------------------------------- #
    svuota_albero()
    n = copia_stato_attuale()
    git("add", "-A")
    git("commit", "-q", "-m", MESSAGGIO_FINALE.format(n=n))
    print(f"\n  {git('log', '-1', '--format=%h').strip()}  stato attuale  {n} file")
    print()
    print(git("log", "--oneline"))
    print(git("count-objects", "-vH").strip())
    return 0


MESSAGGIO_FINALE = """Stato attuale: mesh corretta, sonde di bordo, correzione dello sbocco

Da qui in avanti il progetto e' sotto controllo di versione. Questo commit
raccoglie il lavoro fatto dopo l'ultimo archivio (30/08) e contiene {n} file.

Le correzioni piu' importanti rispetto all'ultimo zip:

- mesh dell'iniettore: tolta la tasca cieca di 0.731 x 2.000 mm davanti a ogni
  getto (uno sfondamento di +2 mm pensato per un taglio booleano e finito
  dentro una fusione), e corretta la correzione: lo sconfinamento e' da una
  parte sola, perche' metterlo anche all'estremita' interna allungava il foro
  del 60 % e faceva sparire la patch d'ingresso del GPL - la mesh usciva con
  cinque bordi invece di sei e checkMesh diceva "Mesh OK".

- un gruppo di bordo vuoto ora solleva un errore invece di essere saltato in
  silenzio.

- sonde di bordo a passi invece che a tempo: con `adjustableRunTime` la corsa 4
  aveva prodotto DUE campioni in 250 microsecondi, e la domanda "gira con le
  pressioni giuste?" era rimasta senza risposta misurata.

- condizione allo sbocco: da `waveTransmissive` con lInf 2.7 volte il dominio
  (che non vincolava la pressione media, e ha fatto andare la corsa 5 alla
  deriva fino a strozzare il getto) a `fixedMean` ancorato a p_c, reso
  possibile dalla rampa di portata che toglie l'onda d'avvio alla radice.

- riscritti i test che fissavano l'IMPLEMENTAZIONE invece del REQUISITO: uno
  imponeva `waveTransmissive` per nome e avrebbe impedito proprio la
  correzione che serviva.

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01PqPTkjhc82NmVTJozaYyCp"""


if __name__ == "__main__":
    raise SystemExit(main())
