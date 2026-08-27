#!/usr/bin/env python3
"""Avvia l'interfaccia web locale.

    python scripts/serve.py
    python scripts/serve.py --port 8080 --reload

Il server ascolta su 127.0.0.1 e NON deve essere esposto in rete: non ha
autenticazione, esegue calcoli pesanti su richiesta e serve file dal disco.
E' uno strumento da scrivania, non un servizio.
"""
from __future__ import annotations

import argparse
import logging
import webbrowser
from threading import Timer


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8000)
    ap.add_argument("--reload", action="store_true", help="ricarica a ogni modifica (sviluppo)")
    ap.add_argument("--no-browser", action="store_true")
    ap.add_argument("--log-level", default="info")
    args = ap.parse_args()

    logging.basicConfig(
        level=args.log_level.upper(),
        format="%(asctime)s  %(levelname)-7s %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )

    import uvicorn

    if args.host not in ("127.0.0.1", "localhost", "::1"):
        logging.getLogger("zefiro.web").warning(
            "stai ascoltando su %s: il server non ha autenticazione, "
            "non esporlo in rete", args.host)

    url = f"http://{args.host if args.host != '0.0.0.0' else '127.0.0.1'}:{args.port}/"
    print(f"\n  Zefiro  ->  {url}\n  documentazione API  ->  {url}api/docs\n")
    if not args.no_browser:
        Timer(1.2, lambda: webbrowser.open(url)).start()

    uvicorn.run(
        "zefiro.web.app:create_app" if args.reload else create_target(),
        host=args.host, port=args.port, reload=args.reload,
        factory=args.reload, log_level=args.log_level,
    )
    return 0


def create_target():
    from zefiro.web import create_app
    return create_app()


if __name__ == "__main__":
    raise SystemExit(main())
