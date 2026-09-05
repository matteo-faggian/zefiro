"""Geometria parametrica: dal vettore di progetto a STEP/STL.

IMPORT PIGRO DI build123d, e non e' un dettaglio di stile.

`build.py` importa build123d, che tira dentro OCCT: centinaia di megabyte di
kernel CAD. Importandolo qui in cima, **ogni** cosa che tocca
`zefiro.geometry.parameters` - cioe' il sintetizzatore, l'ottimizzatore, il
DoE, i test - richiedeva OCCT installato, anche solo per calcolare un'area di
gola. Su una macchina senza build123d il sintetizzatore non partiva, e il
messaggio d'errore parlava di un modulo che non c'entrava niente con quello che
si stava facendo.

La geometria implicita (SDF) non usa OCCT affatto: serve solo per l'export STEP,
che e' un passo finale e opzionale. Quindi si importa quando lo si chiede.
"""
from zefiro.geometry.parameters import DESIGN_BOUNDS, derive  # noqa: F401

_ESPORTATI_DA_BUILD = frozenset({
    "BuildOptions", "build_and_export", "build_solid", "normalize_step_header",
    "stl_is_watertight", "wetted_area",
})


def __getattr__(nome: str):
    """PEP 562: `from zefiro.geometry import build_solid` continua a funzionare,
    ma OCCT viene caricato solo in quel momento."""
    if nome in _ESPORTATI_DA_BUILD:
        from zefiro.geometry import build
        return getattr(build, nome)
    raise AttributeError(f"module {__name__!r} has no attribute {nome!r}")


def __dir__():
    return sorted(set(globals()) | _ESPORTATI_DA_BUILD)
