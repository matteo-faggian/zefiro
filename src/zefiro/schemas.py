"""IL CONTRATTO DATI.

Ogni cosa che passa da un modulo all'altro e' definita qui e solo qui.
Regole (vedi docs/architettura.md sezione 4):
  * unita' SI stretto;
  * dataclass frozen, con `schema_version`;
  * niente default fisici (stanno nei YAML di config/);
  * serializzazione JSON con repr() sui float -> round-trip esatto.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Mapping

SCHEMA_VERSION = "zefiro-schema-0.3.0"


# --------------------------------------------------------------------------- #
# Errori del dominio
# --------------------------------------------------------------------------- #
class ZefiroError(Exception):
    """Base."""


class MissingDatum(ZefiroError):
    """Un dato fisico richiesto non e' stato fornito (e non va inventato)."""


class MissingThermoData(ZefiroError):
    """Il meccanismo dichiarato non contiene una specie/proprieta' necessaria."""


class ContractViolation(ZefiroError):
    """Un contratto dati di questo modulo e' stato violato."""


# --------------------------------------------------------------------------- #
# Helper di serializzazione
# --------------------------------------------------------------------------- #
def _jsonable(obj: Any) -> Any:
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, tuple):
        return list(obj)
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list,)):
        return [_jsonable(v) for v in obj]
    return obj


class _Serializable:
    def to_dict(self) -> dict[str, Any]:
        return _jsonable(asdict(self))  # type: ignore[call-overload]

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent, sort_keys=True)


# --------------------------------------------------------------------------- #
# Ingresso
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class FuelSpec(_Serializable):
    """Combustibile. `composition` in frazioni MOLARI.

    Molari e non massiche perche' le schede tecniche del GPL danno percentuali
    volumetriche, che in fase gassosa coincidono con le molari: cosi' il dato
    entra nel codice senza passaggi non tracciabili.
    """
    composition: Mapping[str, float]
    phase_at_injection: str            # "gas" | "liquid"
    thermo_source: str                 # nome file meccanismo per Cantera
    schema_version: str = SCHEMA_VERSION

    def __post_init__(self) -> None:
        if self.phase_at_injection not in ("gas", "liquid"):
            raise ContractViolation(f"phase_at_injection non valida: {self.phase_at_injection!r}")
        tot = sum(self.composition.values())
        if not self.composition or abs(tot - 1.0) > 1e-6:
            raise ContractViolation(f"composition deve sommare a 1.0, somma = {tot!r}")


@dataclass(frozen=True)
class OperatingPoint(_Serializable):
    """Condizioni al contorno dell'impianto. I `None` sono TODO bloccanti."""
    p_amb: float                       # Pa
    p_air_supply: float                # Pa, monte regolatore
    p_fuel_supply: float               # Pa, monte regolatore
    fuel: FuelSpec
    T_air_in: float | None = None      # K, valle regolatore
    T_fuel_in: float | None = None     # K
    mdot_air_max: float | None = None  # kg/s   <-- TODO n.1 (compressore)
    cd_injector_ox: float | None = None
    cd_injector_fuel: float | None = None
    burn_time: float | None = None     # s
    schema_version: str = SCHEMA_VERSION

    def require(self, *names: str) -> None:
        """Fallisce elencando TUTTI i dati mancanti fra quelli richiesti.

        Elencarli tutti insieme e' voluto: scoprirli uno per volta a ogni
        esecuzione fallita e' il modo piu' lento possibile di raccogliere dati.
        """
        missing = [n for n in names if getattr(self, n) is None]
        if missing:
            raise MissingDatum(
                "Dati operativi mancanti (vedi docs/architettura.md sezione 8): "
                + ", ".join(missing)
            )


@dataclass(frozen=True)
class DesignVector(_Serializable):
    """L'unica cosa che l'ottimizzatore muove."""
    values: Mapping[str, float]
    bounds: Mapping[str, tuple[float, float]]
    schema_version: str = SCHEMA_VERSION

    def __post_init__(self) -> None:
        if set(self.values) != set(self.bounds):
            raise ContractViolation(
                f"values e bounds hanno chiavi diverse: "
                f"{sorted(set(self.values) ^ set(self.bounds))}"
            )
        for k, v in self.values.items():
            lo, hi = self.bounds[k]
            if not (lo <= v <= hi):
                # Volutamente un errore e non un clamp: un ottimizzatore che
                # esce dai bounds ha un bug, e il clamp lo renderebbe invisibile.
                raise ContractViolation(f"{k} = {v!r} fuori dai bounds [{lo!r}, {hi!r}]")


# --------------------------------------------------------------------------- #
# Geometria
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class GeometryParams(_Serializable):
    """Record completo: parametri liberi + derivati. Va salvato accanto al CAD.

    E' cio' che rende un file STEP interpretabile fra sei mesi.
    """
    free: Mapping[str, float]
    derived: Mapping[str, float]
    plug_contour_x: tuple[float, ...]   # m
    plug_contour_r: tuple[float, ...]   # m
    schema_version: str = SCHEMA_VERSION


@dataclass(frozen=True)
class GeometryArtifact(_Serializable):
    run_id: str
    step_path: Path
    stl_path: Path
    params: GeometryParams
    volume: float                # m^3, materiale di parete
    wetted_area: float           # m^2
    is_valid_brep: bool
    is_watertight_mesh: bool
    n_triangles: int
    mesh_tolerance: float        # m
    sha256_step: str
    sha256_stl: str
    schema_version: str = SCHEMA_VERSION


# --------------------------------------------------------------------------- #
# L0
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class L0Result(_Serializable):
    # miscela
    phi_core: float
    phi_global: float
    AFR_stoich: float
    mdot_air: float
    mdot_fuel_core: float
    mdot_fuel_film: float
    # termochimica di camera
    T_ad: float
    X_eq: Mapping[str, float]
    gamma_c: float
    MW_c: float
    cp_c: float
    # prestazione 1D
    c_star: float
    M_e: float
    epsilon: float
    C_F: float
    Isp_s: float          # su massa TOTALE (aria + GPL) -> metrica di sistema
    Isp_fuel_s: float     # sul solo GPL              -> metrica di combustione
    thrust: float
    A_t: float
    # severita' termica in gola: viene da Bartz (empirica, +-30 %). Sta qui
    # perche' e' un OBIETTIVO di progetto, non un post-processing: a parita' di
    # spinta si sceglie il candidato che scalda meno.
    q_throat: float | None = None          # W/m^2
    T_wall_adiabatic: float | None = None  # K
    # proxy di massa: volume di parete stimato per Pappo-Guldino, senza CAD.
    # Serve a poter mettere la massa fra gli obiettivi in un ciclo veloce.
    wall_volume: float | None = None  # m^3, volume ESATTO del solido di
    #  rivoluzione (geometry/profile.py). Non e' una stima: coincide con OCCT
    #  a meno dei fori d'iniezione (test_objectives lo misura: ~1e-15).
    # tempi caratteristici
    tau_res: float | None = None
    tau_chem: float | None = None
    damkohler: float | None = None
    # diagnostica: il risultato porta con se' il proprio dominio di validita'
    warnings: tuple[str, ...] = ()
    assumptions: tuple[str, ...] = ()
    schema_version: str = SCHEMA_VERSION


# --------------------------------------------------------------------------- #
# L1 (contratti definiti ora, implementazione in fase 3-4)
# --------------------------------------------------------------------------- #
CANONICAL_BOUNDARIES: tuple[str, ...] = (
    "inlet_air", "inlet_fuel_core", "inlet_fuel_film",
    "wall_faceplate",
    "wall_chamber", "wall_convergent", "wall_throat", "wall_plug", "wall_cowl",
    "outlet_far", "axis", "periodic_a", "periodic_b",
)
# `wall_faceplate` e' la piastra d'iniezione al netto dei fori. Mancava nella
# prima stesura, ed era una svista con conseguenze: e' la parete su cui si e'
# scoperto che il film cooling in testa protegge una zona che non ne ha bisogno
# (docs/architettura.md 8.5). Una parete senza nome e' una parete di cui non si
# guarda il carico termico.
#
# `wall_throat` e `wall_cowl` restano nel registro ma NON sono prodotte da un
# aerospike a espansione esterna con labbro di spessore nullo: li' la gola e'
# delimitata dal plug e da uno spigolo, non da due superfici. Il registro
# elenca i nomi ammessi, non quelli obbligatori: vedi l1.mesh.REQUIRED_BOUNDARIES.


@dataclass(frozen=True)
class MeshArtifact(_Serializable):
    run_id: str
    msh_path: Path
    n_cells: int
    boundary_names: tuple[str, ...]
    is_periodic: bool
    sector_angle: float
    max_non_orthogonality: float
    max_skewness: float
    sha256: str
    #: Frazione di facce interne oltre i 70 gradi di non-ortogonalita'.
    #: Il MASSIMO da solo e' un pessimo riassunto di 47000 facce: puo' essere
    #: alto per una scheggia sola in un angolo che non interessa. La frazione
    #: dice se il problema e' locale o diffuso, e sono due situazioni diverse.
    frac_non_orthogonal: float = 0.0
    #: Scarto massimo fra le due facce periodiche dopo la rotazione, in
    #: multipli del raggio del labbro. Registrato e non solo confrontato con
    #: una soglia: un numero si puo' rivedere, un booleano no.
    periodic_mismatch: float = 0.0
    warnings: tuple[str, ...] = ()
    schema_version: str = SCHEMA_VERSION


@dataclass(frozen=True)
class CFDResult(_Serializable):
    run_id: str
    mesh_id: str
    converged: bool
    residual_final: Mapping[str, float]
    eta_combustion: float
    p_wall: Path
    q_wall: Path
    thrust_cfd: float
    solver_version: str
    schema_version: str = SCHEMA_VERSION


@dataclass(frozen=True)
class FEMResult(_Serializable):
    run_id: str
    cfd_id: str
    T_wall_max: float
    von_mises_max: float
    margin_yield: float
    margin_temp: float
    solver_version: str
    schema_version: str = SCHEMA_VERSION


# --------------------------------------------------------------------------- #
# Obiettivi
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Objectives(_Serializable):
    """Convenzione RIGIDA: minimizzare `f`, fattibile se `g <= 0`.

    Un obiettivo da massimizzare entra come -valore, e l'inversione la fa chi
    produce l'obiettivo, non l'ottimizzatore: cosi' cambiare ottimizzatore non
    cambia il significato dei numeri.
    """
    run_id: str
    fidelity: str                    # "L0" | "L1" | "L2"
    f: Mapping[str, float]
    g: Mapping[str, float]
    source: Mapping[str, str] = field(default_factory=dict)
    schema_version: str = SCHEMA_VERSION

    @property
    def feasible(self) -> bool:
        return all(v <= 0.0 for v in self.g.values())
