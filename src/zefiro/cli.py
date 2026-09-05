"""Interfacce a riga di comando."""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

from zefiro.config import load_design_vector, load_material, load_operating_point
from zefiro.geometry import BuildOptions, build_and_export
from zefiro.geometry.parameters import check_manufacturability, derive
from zefiro.schemas import ZefiroError
from zefiro.store import env_fingerprint, make_run_id


def runs_root() -> Path:
    """Cartella degli artefatti di run.

    Impostabile con la variabile d'ambiente ZEFIRO_RUNS. Serve perche' sotto
    WSL il repo sta su /mnt/... (drvfs), dove l'I/O su file piccoli e numerosi
    e' 5-10 volte piu' lento del filesystem nativo — ed e' esattamente cio' che
    produce una mesh. Il repo puo' restare su /mnt perche' e' piccolo e vuoi
    vederlo da Windows; gli artefatti no.
    """
    return Path(os.environ.get("ZEFIRO_RUNS", "runs"))


def _common(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--operating", type=Path, default=None)
    parser.add_argument("--design", type=Path, default=None)
    parser.add_argument("--out", type=Path, default=None,
                        help="default: $ZEFIRO_RUNS, oppure ./runs")
    parser.add_argument("--mdot-air", type=float, default=None,
                        help="kg/s: sovrascrive mdot_air_max (utile finche' il TODO n.1 e' aperto)")


def _guard(fn):
    """Trasforma un errore di dominio in un messaggio leggibile e codice 2.

    Un TODO non compilato non e' un crash del programma: e' il programma che
    fa il suo lavoro. Mostrarlo come traceback lo farebbe sembrare un bug.
    """
    def wrapped(argv: list[str] | None = None) -> int:
        try:
            return fn(argv)
        except ZefiroError as exc:
            print(f"\n{type(exc).__name__}: {exc}\n", file=sys.stderr)
            return 2
    return wrapped


def main_l0(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Valutazione L0 (termochimica 0D/1D).")
    _common(ap)
    args = ap.parse_args(argv)

    op = load_operating_point(args.operating)
    x = load_design_vector(args.design)
    params, l0 = derive(x, op, mdot_air=args.mdot_air)

    run_id = make_run_id(x, op)
    out = (args.out or runs_root()) / run_id
    out.mkdir(parents=True, exist_ok=True)
    (out / "l0.json").write_text(l0.to_json(), encoding="utf-8")
    (out / "design.json").write_text(x.to_json(), encoding="utf-8")
    (out / "env.json").write_text(json.dumps(env_fingerprint(), indent=2), encoding="utf-8")

    print(f"run_id      {run_id}")
    print(f"T_ad        {l0.T_ad:9.1f} K")
    print(f"c*          {l0.c_star:9.1f} m/s")
    print(f"eps         {l0.epsilon:9.3f}     M_e {l0.M_e:.3f}   C_F {l0.C_F:.4f}")
    print(f"spinta      {l0.thrust:9.2f} N")
    print(f"Isp totale  {l0.Isp_s:9.1f} s   (su aria + GPL)")
    print(f"Isp GPL     {l0.Isp_fuel_s:9.1f} s   (solo GPL)")
    print(f"L*          {params.derived['L_star']:9.3f} m")
    for w in l0.warnings:
        print(f"AVVISO: {w}")
    return 0


def main_geometry(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Genera la geometria (STEP + STL).")
    _common(ap)
    ap.add_argument("--sector", action="store_true", help="ritaglia il settore periodico 1/N")
    args = ap.parse_args(argv)

    op = load_operating_point(args.operating)
    x = load_design_vector(args.design)
    params, l0 = derive(x, op, mdot_air=args.mdot_air)
    run_id = make_run_id(x, op)
    out = (args.out or runs_root()) / run_id

    art = build_and_export(params, out, run_id, BuildOptions(sector=args.sector))
    (out / "geometry_params.json").write_text(params.to_json(), encoding="utf-8")

    print(f"run_id       {run_id}")
    print(f"STEP         {art.step_path}")
    print(f"STL          {art.stl_path}  ({art.n_triangles} triangoli)")
    print(f"B-Rep valido {art.is_valid_brep}")
    print(f"water-tight  {art.is_watertight_mesh}")
    print(f"volume       {art.volume * 1e9:.1f} mm^3")
    print(f"area bagnata {art.wetted_area * 1e6:.1f} mm^2")

    mat = load_material()
    for msg in check_manufacturability(params, mat["process"]["min_feature_size_m"]):
        print(f"FABBRICABILITA': {msg}")
    return 0


main_l0 = _guard(main_l0)
main_geometry = _guard(main_geometry)
