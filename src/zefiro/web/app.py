"""Applicazione FastAPI: solo rotte, nessuna fisica e nessuna logica.

Le rotte fanno tre cose e basta: validare l'ingresso, chiamare `service`, e
tradurre gli errori di dominio in risposte strutturate. Tutto cio' che calcola
sta in `service.py`, e tutto cio' che calcola *sul serio* sta nei moduli
`zefiro.*` che usa anche la riga di comando.

L'API e' versionata sotto `/api/v1`. Non e' burocrazia: il frontend e' un file
statico che l'utente puo' avere in cache, quindi cambiare la forma di una
risposta senza cambiare versione significa una pagina rotta senza spiegazione.
"""
from __future__ import annotations

import logging
import time
import uuid
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Any

from fastapi import APIRouter, FastAPI, Request
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

from zefiro import __version__
from zefiro.cli import runs_root
from zefiro.web import errors, service
from zefiro.web.errors import ApiError, to_api_error
from zefiro.web.jobs import JobRegistry

log = logging.getLogger("zefiro.web")
WEB_DIR = Path(__file__).resolve().parents[3] / "web"


# --------------------------------------------------------------------------- #
# Modelli di richiesta
# --------------------------------------------------------------------------- #
class DesignPayload(BaseModel):
    design: dict[str, float] = Field(..., description="valori in unita' di visualizzazione")
    operating: dict[str, float | None] = Field(default_factory=dict)


class GeometryPayload(DesignPayload):
    sector: bool = Field(False, description="ritaglia il settore periodico 1/N")


class ThermalPayload(DesignPayload):
    # Proprieta' del metallo. Sono il TODO n.J: l'interfaccia le chiede
    # esplicitamente e marca il risultato come basato su segnaposto.
    k_wall: float = Field(16.0, gt=0, description="W/(m K)")
    rho_wall: float = Field(7980.0, gt=0, description="kg/m3")
    cp_wall: float = Field(500.0, gt=0, description="J/(kg K)")
    water_cooled: bool = True
    h_coolant: float = Field(51676.0, ge=0, description="W/(m2 K) lato acqua")
    T_coolant: float = Field(288.15, gt=0, description="K")


class SweepPayload(BaseModel):
    n: int = Field(200, ge=1, le=5000)
    seed: int = Field(..., description="obbligatorio: un piano non riproducibile non entra nel database")
    operating: dict[str, float | None] = Field(default_factory=dict)


# --------------------------------------------------------------------------- #
# Applicazione
# --------------------------------------------------------------------------- #
def create_app(registry: JobRegistry | None = None) -> FastAPI:
    jobs = registry or JobRegistry()

    @asynccontextmanager
    async def lifespan(_: FastAPI):
        log.info("Zefiro %s — artefatti in %s", __version__, runs_root())
        yield
        jobs.shutdown()

    app = FastAPI(
        title="Zefiro",
        version=__version__,
        description=(
            "Interfaccia locale per la progettazione del combustore. "
            "Il backend chiama le stesse funzioni della riga di comando: ogni "
            "valutazione produce un run_id vero e finisce nello stesso database."
        ),
        lifespan=lifespan,
        docs_url="/api/docs",
        openapi_url="/api/openapi.json",
    )
    app.state.jobs = jobs
    errors.install(app)

    @app.middleware("http")
    async def correlate(request: Request, call_next):
        rid = uuid.uuid4().hex[:8]
        t0 = time.perf_counter()
        request.state.request_id = rid
        try:
            response = await call_next(request)
        except ApiError:
            raise
        except Exception as exc:  # noqa: BLE001
            try:
                api = to_api_error(exc)
            except Exception:
                log.exception("[%s] errore non gestito su %s", rid, request.url.path)
                return JSONResponse(
                    status_code=500,
                    content={"error": {"code": "ERRORE_INTERNO",
                                       "message": "errore interno del server",
                                       "details": {"request_id": rid}}},
                )
            api.details.setdefault("request_id", rid)
            log.warning("[%s] %s su %s: %s", rid, api.code, request.url.path, api.message)
            return JSONResponse(status_code=api.status, content=api.payload())
        dt = (time.perf_counter() - t0) * 1e3
        response.headers["X-Request-Id"] = rid
        response.headers["Server-Timing"] = f"app;dur={dt:.1f}"
        if dt > 500.0:
            log.info("[%s] %s %s in %.0f ms", rid, request.method, request.url.path, dt)
        return response

    api = APIRouter(prefix="/api/v1")

    # -- lettura ----------------------------------------------------------- #
    @api.get("/schema", summary="Parametri, bounds e dati mancanti")
    def get_schema() -> dict[str, Any]:
        return service.schema()

    @api.get("/health", summary="Stato del servizio")
    def health() -> dict[str, Any]:
        s = service.schema()
        return {
            "status": "ok",
            "version": __version__,
            "code_rev": s["code_rev"],
            "runs_root": s["runs_root"],
            "dati_mancanti": [f["name"] for f in s["operating"] if f["missing"]],
            "materiale_segnaposto": s["material_is_placeholder"],
        }

    @api.get("/runs", summary="Storico delle valutazioni")
    def get_runs(limit: int = 500) -> dict[str, Any]:
        return service.runs(limit=min(max(limit, 1), 5000))

    @api.get("/plant", summary="Che motore permette l'impianto")
    def get_plant(power_hp: float = 3.0, tank_l: float = 100.0, p_max: float = 10.0,
                  p_min: float = 6.0, burn: float = 5.0, T_amb: float = 20.0) -> dict[str, Any]:
        if p_min >= p_max:
            raise ApiError("PARAMETRO_NON_VALIDO",
                           "la pressione minima deve essere sotto la massima", 422,
                           {"fields": ["p_min", "p_max"]})
        return service.plant(power_hp, tank_l, p_max, p_min, burn, T_amb)

    # -- calcolo sincrono (~50 ms) ------------------------------------------ #
    @api.post("/evaluate", summary="Valutazione L0 (immediata)")
    def post_evaluate(payload: DesignPayload) -> dict[str, Any]:
        return service.evaluate(payload.design, payload.operating)

    @api.post("/thermal", summary="Carico termico e risposta della parete")
    def post_thermal(payload: ThermalPayload) -> dict[str, Any]:
        return service.thermal(
            payload.design, payload.operating, payload.k_wall, payload.rho_wall,
            payload.cp_wall, payload.water_cooled, payload.h_coolant, payload.T_coolant)

    # -- lavori in coda ------------------------------------------------------ #
    @api.post("/jobs/geometry", status_code=202, summary="Genera STEP e STL (in coda)")
    def post_geometry(payload: GeometryPayload) -> dict[str, Any]:
        # Si valida SUBITO, in modo che un progetto sbagliato dia 422 adesso e
        # non un lavoro fallito fra dieci secondi.
        service.design_to_si(payload.design)
        job = jobs.submit("geometry", _geometry_job, design=payload.design,
                          operating=payload.operating, sector=payload.sector)
        return job.public()

    @api.post("/jobs/sweep", status_code=202, summary="Piano sperimentale su L0 (in coda)")
    def post_sweep(payload: SweepPayload) -> dict[str, Any]:
        job = jobs.submit("sweep", _sweep_job, n=payload.n, seed=payload.seed,
                          operating=payload.operating)
        return job.public()

    @api.get("/jobs/{job_id}", summary="Stato di un lavoro")
    def get_job(job_id: str) -> dict[str, Any]:
        job = jobs.get(job_id)
        if job is None:
            raise ApiError("LAVORO_SCONOSCIUTO", f"nessun lavoro con id {job_id}", 404)
        return job.public()

    @api.delete("/jobs/{job_id}", summary="Annulla un lavoro")
    def delete_job(job_id: str) -> dict[str, Any]:
        if not jobs.cancel(job_id):
            raise ApiError("LAVORO_NON_ANNULLABILE",
                           "il lavoro non esiste o e' gia' concluso", 409)
        return {"cancelled": job_id}

    @api.get("/jobs", summary="Lavori recenti")
    def list_jobs(limit: int = 30) -> dict[str, Any]:
        return {"jobs": jobs.recent(limit)}

    # -- artefatti ----------------------------------------------------------- #
    @api.get("/geometry/{run_id}/{kind}", summary="Scarica STEP o STL")
    def geometry_file(run_id: str, kind: str):
        if kind not in ("stl", "step"):
            raise ApiError("TIPO_NON_VALIDO", "i tipi ammessi sono stl e step", 404)
        if not run_id.replace("-", "").isalnum():
            raise ApiError("RUN_ID_NON_VALIDO", "identificativo non valido", 400)
        f = runs_root() / run_id / f"{run_id}_geometry.{kind}"
        if not f.exists():
            raise ApiError("ARTEFATTO_ASSENTE",
                           "geometria non ancora generata per questa run", 404,
                           {"run_id": run_id})
        return FileResponse(f, media_type="model/stl" if kind == "stl" else "application/step",
                            filename=f.name)

    app.include_router(api)

    assets = WEB_DIR / "assets"
    if assets.is_dir():
        app.mount("/assets", StaticFiles(directory=str(assets)), name="assets")

    @app.get("/", response_class=HTMLResponse, include_in_schema=False)
    def index() -> HTMLResponse:
        f = WEB_DIR / "index.html"
        if not f.exists():
            raise ApiError("INTERFACCIA_ASSENTE", f"manca {f}", 500)
        # niente cache sull'HTML: l'API e' versionata, ma la pagina deve poter
        # cambiare senza che l'utente debba svuotare la cache del browser.
        return HTMLResponse(f.read_text(encoding="utf-8"),
                            headers={"Cache-Control": "no-cache"})

    return app


# I lavori girano in un thread separato e ricevono `progress` dal registro.
def _geometry_job(progress, design, operating, sector):
    return service.build_geometry(design, operating, sector, progress=progress)


def _sweep_job(progress, n, seed, operating):
    return service.sweep(n, seed, operating, progress=progress)
