#!/usr/bin/env python3
"""DOE deterministico su L0, con i risultati nel database delle run.

Migliaia di valutazioni L0 in pochi minuti: e' il livello dove si fa il 90 %
del lavoro di design. La geometria (--with-geometry) e' molto piu' lenta
(~2 s a candidato) e serve solo se vuoi i vincoli di fabbricabilita' e il
volume di parete.

Uso:
    python scripts/sweep_l0.py --n 500 --seed 0 --mdot-air 0.05
    python scripts/sweep_l0.py --n 40 --seed 0 --mdot-air 0.05 --with-geometry
"""
from __future__ import annotations

import argparse
import sys
import tempfile
from datetime import datetime, timezone
from hashlib import blake2b
from json import dumps
from pathlib import Path

from zefiro.cli import runs_root
from zefiro.config import load_design_vector, load_material, load_operating_point
from zefiro.doe import latin_hypercube
from zefiro.geometry import build_and_export
from zefiro.geometry.parameters import DESIGN_BOUNDS, INTEGER_PARAMETERS, derive
from zefiro.opt.objectives import objectives_l0
from zefiro.schemas import DesignVector, ZefiroError
from zefiro.store import RunStore, env_fingerprint, make_run_id
from zefiro.store.runid import git_revision


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--seed", type=int, required=True,
                    help="obbligatorio: un DOE non riproducibile non entra nel database")
    ap.add_argument("--db", type=Path, default=None,
                    help="default: $ZEFIRO_RUNS/runs.db, oppure ./runs/runs.db")
    ap.add_argument("--operating", type=Path, default=None)
    ap.add_argument("--mdot-air", type=float, default=None)
    ap.add_argument("--with-geometry", action="store_true")
    args = ap.parse_args(argv)
    db_path = args.db or (runs_root() / "runs.db")

    op = load_operating_point(args.operating)
    material = load_material()
    min_feature = material["process"]["min_feature_size_m"]

    samples = latin_hypercube(args.n, DESIGN_BOUNDS, args.seed, INTEGER_PARAMETERS)
    rev = git_revision()
    env = env_fingerprint()
    env_hash = blake2b(dumps(env, sort_keys=True).encode(), digest_size=8).hexdigest()
    now = datetime.now(timezone.utc).isoformat()

    ok = skipped = 0
    tmp = Path(tempfile.mkdtemp(prefix="zefiro-sweep-"))
    with RunStore(db_path) as store:
        for i, values in enumerate(samples):
            try:
                x = DesignVector(values=values, bounds=DESIGN_BOUNDS)
                params, l0 = derive(x, op, mdot_air=args.mdot_air)
                geom = None
                if args.with_geometry:
                    geom = build_and_export(params, tmp, f"sweep-{i:05d}")
                obj = objectives_l0(
                    make_run_id(x, op, rev), l0, op, values["p_c"],
                    geometry=geom, min_feature_size=min_feature,
                    derived=params.derived,
                )
                store.insert(obj.run_id, x, l0, now, rev, env_hash,
                             geometry=geom, objectives=obj, replace=True)
                ok += 1
            except (ZefiroError, ValueError, NotImplementedError) as exc:
                # Un candidato geometricamente o fisicamente impossibile NON e'
                # un errore del DOE: e' informazione. Si conta e si va avanti.
                skipped += 1
                if skipped <= 5:
                    print(f"  [skip {i}] {type(exc).__name__}: {exc}", file=sys.stderr)
            if (i + 1) % 50 == 0:
                print(f"  {i+1}/{len(samples)}  ok={ok} skip={skipped}", file=sys.stderr)

        print(f"\nvalutati {ok}, scartati {skipped}")
        print(f"database: {db_path}  ({store.count('runs')} run pulite, "
              f"{store.count('runs_dirty')} sporche)")
        feas = store.query("SELECT COUNT(*) c FROM runs WHERE feasible = 1")[0]["c"]
        print(f"fattibili: {feas}")
        best = store.pareto_candidates(limit=5)
        if best:
            print("\n  migliori 5 per spinta:")
            print(f"  {'run_id':<12}{'spinta':>9}{'Isp':>8}{'T_ad':>8}{'p_c':>8}{'phi':>7}{'f_film':>8}")
            for r in best:
                print(f"  {r['run_id'][:10]:<12}{r['thrust']:8.2f}N{r['Isp_s']:7.1f}s"
                      f"{r['T_ad']:7.0f}K{r['p_c']/1e5:7.2f}b{r['phi_core']:7.3f}{r['f_film']:8.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
