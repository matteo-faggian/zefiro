"""Errori dell'API in forma strutturata.

Un errore che arriva al client come stringa libera e' inutilizzabile: il
frontend non puo' distinguere "ti manca un dato" da "hai sforato un vincolo" da
"il server e' rotto", e finisce per mostrare all'utente un messaggio generico.
Qui ogni errore ha un CODICE stabile, un messaggio leggibile e un campo
`details` che il frontend puo' usare per evidenziare i campi giusti.
"""
from __future__ import annotations

from typing import Any

from fastapi import Request
from fastapi.responses import JSONResponse

from zefiro.schemas import (
    ContractViolation,
    MissingDatum,
    MissingThermoData,
    ZefiroError,
)


class ApiError(Exception):
    """Errore con codice stabile. Il codice fa parte del contratto dell'API."""

    def __init__(self, code: str, message: str, status: int = 422,
                 details: dict[str, Any] | None = None) -> None:
        super().__init__(message)
        self.code = code
        self.message = message
        self.status = status
        self.details = details or {}

    def payload(self) -> dict[str, Any]:
        return {"error": {"code": self.code, "message": self.message,
                          "details": self.details}}


#: Mappa dagli errori di dominio ai codici dell'API. Aggiungere un errore di
#: dominio senza aggiungerlo qui lo farebbe arrivare al client come 500, cioe'
#: come se fosse colpa del server: il test `test_web.py` lo verifica.
DOMAIN_CODES: dict[type[ZefiroError], tuple[str, int]] = {
    MissingDatum: ("DATO_MANCANTE", 422),
    MissingThermoData: ("TERMODINAMICA_MANCANTE", 422),
    ContractViolation: ("CONTRATTO_VIOLATO", 422),
}


def to_api_error(exc: Exception) -> ApiError:
    for cls, (code, status) in DOMAIN_CODES.items():
        if isinstance(exc, cls):
            details: dict[str, Any] = {}
            if isinstance(exc, MissingDatum):
                testo = str(exc)
                if ":" in testo:
                    details["fields"] = [s.strip() for s in testo.rsplit(":", 1)[1].split(",")]
            return ApiError(code, str(exc), status, details)
    if isinstance(exc, ZefiroError):
        return ApiError("ERRORE_ZEFIRO", str(exc), 422)
    if isinstance(exc, NotImplementedError):
        return ApiError("NON_IMPLEMENTATO", str(exc) or "funzione non ancora implementata", 501)
    if isinstance(exc, ValueError):
        return ApiError("PARAMETRO_NON_VALIDO", str(exc), 422)
    raise exc


def install(app) -> None:
    @app.exception_handler(ApiError)
    async def _handler(_: Request, exc: ApiError) -> JSONResponse:  # noqa: RUF029
        return JSONResponse(status_code=exc.status, content=exc.payload())
