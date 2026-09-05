#!/usr/bin/env python3
"""Controllo dell'ambiente: dice che cosa funziona e che cosa manca.

Serve a trasformare "non parte" in un elenco di cose specifiche da sistemare.
Esce con codice 1 se qualcosa di BLOCCANTE non va, 0 altrimenti: cosi' si puo'
usare anche in uno script di avvio.

    python scripts/doctor.py
"""
from __future__ import annotations

import importlib
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path

OK, WARN, FAIL = "  ok  ", " nota ", "MANCA "
_blocking = 0
_warnings = 0


def line(status: str, what: str, detail: str = "") -> None:
    global _blocking, _warnings
    if status == FAIL:
        _blocking += 1
    elif status == WARN:
        _warnings += 1
    print(f"[{status}] {what:<34} {detail}")


def section(title: str) -> None:
    print(f"\n\033[1m{title}\033[0m" if sys.stdout.isatty() else f"\n{title}")
    print("-" * 72)


def check_import(name: str, label: str | None = None, blocking: bool = True,
                 hint: str = "") -> None:
    label = label or name
    try:
        mod = importlib.import_module(name)
        v = getattr(mod, "__version__", None)
        if v is None and name == "cantera":
            v = mod.__version__
        line(OK, label, str(v) if v else "")
    except ImportError:
        line(FAIL if blocking else WARN, label, hint or "non installato")


def main() -> int:
    section("Interprete")
    v = sys.version_info
    line(OK if v >= (3, 10) else FAIL, "Python", f"{v.major}.{v.minor}.{v.micro}")
    line(OK, "sistema", platform.platform())
    in_venv = sys.prefix != sys.base_prefix
    line(OK if in_venv else WARN, "ambiente isolato",
         sys.prefix if in_venv else "stai usando il Python di sistema")
    if "microsoft" in platform.release().lower() or "WSL_DISTRO_NAME" in os.environ:
        line(OK, "WSL", os.environ.get("WSL_DISTRO_NAME", "si'"))
    else:
        line(WARN, "WSL", "non rilevato: OpenFOAM e CalculiX gireranno solo in WSL")

    section("Dipendenze del layer Python")
    check_import("numpy")
    check_import("yaml", "pyyaml")
    check_import("cantera")
    check_import("build123d")
    check_import("CoolProp", "CoolProp")
    check_import("pyarrow")
    check_import("pytest", blocking=False, hint="serve solo per far girare i test")

    section("Il pacchetto zefiro")
    try:
        import zefiro
        line(OK, "import zefiro", zefiro.__version__)
        installed_editable = "site-packages" not in str(Path(zefiro.__file__).resolve())
        line(OK if installed_editable else WARN, "installato in modifica",
             "si'" if installed_editable else 'rifai  pip install -e ".[dev]"')
        from zefiro.geometry import build_and_export  # noqa: F401
        from zefiro.l0 import evaluate_l0  # noqa: F401
        line(OK, "moduli principali", "geometry, l0, feed, thermal, cooling, store")
    except Exception as exc:
        line(FAIL, "import zefiro", f"{type(exc).__name__}: {exc}")
        print("\n  -> dalla radice del repo:  pip install -e \".[dev]\"")
        return 1

    section("Riproducibilita'")
    git = shutil.which("git")
    line(OK if git else FAIL, "git", git or "installalo: senza, il run_id non "
         "distingue le versioni del codice")
    if git:
        try:
            r = subprocess.run(["git", "rev-parse", "--is-inside-work-tree"],
                               capture_output=True, text=True, timeout=10)
            inside = r.returncode == 0 and r.stdout.strip() == "true"
        except Exception:
            inside = False
        if inside:
            from zefiro.store.runid import git_revision
            rev = git_revision()
            line(OK, "repository git", rev)
            if rev.endswith("-dirty"):
                line(WARN, "working tree", "sporco: le run andranno in runs_dirty")
        else:
            line(FAIL, "repository git", "non inizializzato: fai  git init && git add -A "
                 "&& git commit -m 'stato iniziale'")

    runs = Path(os.environ.get("ZEFIRO_RUNS", "runs"))
    on_drvfs = str(runs.resolve()).startswith("/mnt/")
    line(WARN if on_drvfs else OK, "ZEFIRO_RUNS", str(runs.resolve())
         + ("  <- su drvfs: I/O lento, spostalo nel filesystem di WSL" if on_drvfs else ""))

    section("Dati di configurazione")
    try:
        from zefiro.config import load_material, load_operating_point
        op = load_operating_point()
        line(OK, "operating_point.yaml", "caricato")
        mancanti = [n for n in ("T_air_in", "T_fuel_in", "mdot_air_max",
                                "cd_injector_ox", "cd_injector_fuel", "burn_time")
                    if getattr(op, n) is None]
        if mancanti:
            line(WARN, "dati ancora da misurare", ", ".join(mancanti))
        else:
            line(OK, "dati operativi", "completi")
        mat = load_material()
        if mat["density_kg_m3"] is None:
            line(WARN, "AISI 316L (SLM)", "proprieta' non fornite: il FEM non partira'")
    except Exception as exc:
        line(FAIL, "configurazione", f"{type(exc).__name__}: {exc}")

    section("Prova di funzionamento")
    try:
        from zefiro.feed import compressor_bounds
        b = compressor_bounds(2237.0, 101325.0, 10e5, 293.15)
        line(OK, "termodinamica (feed)", f"limite compressore {b.mdot_isothermal*1e3:.1f} g/s")
    except Exception as exc:
        line(FAIL, "termodinamica (feed)", str(exc))
    try:
        from zefiro.l0.mixture import MixtureModel
        from zefiro.schemas import FuelSpec
        m = MixtureModel.from_fuel(FuelSpec(composition={"C3H8": 1.0},
                                            phase_at_injection="gas",
                                            thermo_source="gri30.yaml"))
        line(OK, "Cantera (L0)", f"AFR stechiometrico {m.afr_stoichiometric():.4f}")
    except Exception as exc:
        line(FAIL, "Cantera (L0)", str(exc))
    try:
        from build123d import Axis, Polyline, make_face, revolve
        s = revolve(make_face(Polyline((0, 0, 0), (0, 0, 5), (10, 0, 5), (10, 0, 0),
                                       close=True)), axis=Axis.X, revolution_arc=360)
        line(OK, "OCCT (geometria)", f"cilindro di prova, volume {s.volume:.1f} mm3")
    except Exception as exc:
        line(FAIL, "OCCT (geometria)", str(exc))

    section("Solutori esterni (servono dalla fase 3)")
    for exe, label in (("gmsh", "Gmsh"), ("simpleFoam", "OpenFOAM"), ("ccx", "CalculiX")):
        p = shutil.which(exe)
        line(OK if p else WARN, label, p or "non installato (non serve ancora)")

    print("\n" + "=" * 72)
    if _blocking:
        print(f"{_blocking} problemi BLOCCANTI. Risolvili e rilancia.")
        return 1
    print(f"Tutto a posto. {_warnings} note, nessun problema bloccante.")
    print("Prossimo passo:  pytest -q   poi   zefiro-l0")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
