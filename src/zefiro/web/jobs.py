"""Coda dei lavori lunghi.

Perche' serve: la valutazione L0 costa ~50 ms e si puo' servire nella richiesta,
ma la geometria costa ~2 s e uno sweep costa minuti. Bloccare il thread della
richiesta per quel tempo significa un browser che gira a vuoto senza sapere se
sta succedendo qualcosa, nessun modo di annullare, e un timeout del proxy che
uccide il lavoro a meta'.

Scelte, e i loro motivi:

* **Un solo worker.** Non e' una limitazione da rimuovere: OCCT e i solutori
  numerici non sono thread-safe e sono limitati dalla memoria, non dalla CPU.
  Su una macchina sola, serializzare e' la scelta corretta — evita che due
  build CAD concorrenti si corrompano a vicenda o saturino la RAM.
* **Registro in memoria, niente Redis.** Il server gira su localhost per un
  utente. Introdurre un broker esterno aggiungerebbe un pezzo da installare,
  configurare e tenere acceso, in cambio di nulla.
* **I lavori sopravvivono alla richiesta, non al processo.** Se riavvii il
  server i lavori in corso si perdono: gli artefatti pero' no, perche' ogni
  lavoro scrive su disco sotto il proprio `run_id`. Riprodurre un risultato non
  dipende mai dallo stato del server.
"""
from __future__ import annotations

import logging
import threading
import traceback
import uuid
from concurrent.futures import Future, ThreadPoolExecutor
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Callable

from zefiro.web.errors import to_api_error

log = logging.getLogger("zefiro.web.jobs")

PENDING, RUNNING, DONE, FAILED, CANCELLED = "pending", "running", "done", "failed", "cancelled"


@dataclass
class Job:
    id: str
    kind: str
    status: str = PENDING
    progress: float = 0.0
    message: str = ""
    result: Any = None
    error: dict[str, Any] | None = None
    created_utc: str = ""
    finished_utc: str = ""
    _future: Future | None = field(default=None, repr=False, compare=False)
    _cancel: threading.Event = field(default_factory=threading.Event, repr=False, compare=False)

    def public(self) -> dict[str, Any]:
        return {
            "id": self.id, "kind": self.kind, "status": self.status,
            "progress": round(self.progress, 4), "message": self.message,
            "created_utc": self.created_utc, "finished_utc": self.finished_utc,
            "result": self.result, "error": self.error,
        }


class Cancelled(Exception):
    """Sollevata dentro un lavoro quando l'utente lo annulla."""


class JobRegistry:
    def __init__(self, max_jobs: int = 200) -> None:
        self._jobs: dict[str, Job] = {}
        self._order: list[str] = []
        self._lock = threading.Lock()
        self._pool = ThreadPoolExecutor(max_workers=1, thread_name_prefix="zefiro-job")
        self._max = max_jobs

    def submit(self, kind: str, fn: Callable[..., Any], **kwargs: Any) -> Job:
        job = Job(id=uuid.uuid4().hex[:12], kind=kind,
                  created_utc=datetime.now(timezone.utc).isoformat())
        with self._lock:
            self._jobs[job.id] = job
            self._order.append(job.id)
            while len(self._order) > self._max:
                self._jobs.pop(self._order.pop(0), None)
        job._future = self._pool.submit(self._run, job, fn, kwargs)
        log.info("lavoro %s accodato (%s)", job.id, kind)
        return job

    def _run(self, job: Job, fn: Callable[..., Any], kwargs: dict[str, Any]) -> None:
        if job._cancel.is_set():
            job.status = CANCELLED
            return
        job.status = RUNNING

        def progress(frac: float, msg: str = "") -> None:
            if job._cancel.is_set():
                raise Cancelled
            job.progress = max(0.0, min(1.0, frac))
            if msg:
                job.message = msg

        try:
            job.result = fn(progress=progress, **kwargs)
            job.status = DONE
            job.progress = 1.0
            log.info("lavoro %s completato", job.id)
        except Cancelled:
            job.status = CANCELLED
            job.message = "annullato"
            log.info("lavoro %s annullato", job.id)
        except Exception as exc:  # noqa: BLE001
            job.status = FAILED
            try:
                api = to_api_error(exc)
                job.error = api.payload()["error"]
            except Exception:
                job.error = {"code": "ERRORE_INTERNO", "message": str(exc),
                             "details": {"type": type(exc).__name__}}
                log.error("lavoro %s fallito:\n%s", job.id, traceback.format_exc())
        finally:
            job.finished_utc = datetime.now(timezone.utc).isoformat()

    def get(self, job_id: str) -> Job | None:
        return self._jobs.get(job_id)

    def cancel(self, job_id: str) -> bool:
        job = self._jobs.get(job_id)
        if job is None or job.status in (DONE, FAILED, CANCELLED):
            return False
        job._cancel.set()
        if job._future is not None and job._future.cancel():
            job.status = CANCELLED
            job.finished_utc = datetime.now(timezone.utc).isoformat()
        return True

    def recent(self, limit: int = 30) -> list[dict[str, Any]]:
        with self._lock:
            ids = list(reversed(self._order))[:limit]
        return [self._jobs[i].public() for i in ids if i in self._jobs]

    def shutdown(self) -> None:
        self._pool.shutdown(wait=False, cancel_futures=True)
