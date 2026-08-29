"""Assemblaggio del solido, export STEP/STL e verifiche di tenuta.

NOTA SULLE UNITA'. Il resto del codice e' in SI (metri). Dentro questo modulo
si lavora in MILLIMETRI, per due motivi concreti:
  * le tolleranze di OCCT (1e-7 di default) sono tarate su modelli in mm; con
    un modello in metri una gola da 6 mm diventa 6e-3 e le operazioni booleane
    iniziano a considerare coincidenti facce distinte;
  * lo standard STEP AP214 usa il millimetro, quindi non c'e' conversione
    all'export e nessun lettore CAD deve indovinare la scala.
La conversione avviene solo alle due frontiere di questo file (`_MM`).
"""
from __future__ import annotations

import hashlib
import math
import re
import struct
from dataclasses import dataclass
from pathlib import Path

from build123d import (
    Axis,
    Cylinder,
    Location,
    Polyline,
    Rotation,
    export_step,
    export_stl,
    make_face,
    revolve,
)

from zefiro.geometry.profile import meridian_polygon, wetted_area_from_contour
from zefiro.schemas import GeometryArtifact, GeometryParams

_MM = 1000.0          # m -> mm
#: Diametro minimo sotto il quale un foro non viene modellato. E' una soglia
#: NUMERICA (robustezza delle booleane OCCT in mm), non un'affermazione di
#: fabbricabilita': il limite di processo SLM e' il TODO n.7 e va verificato
#: separatamente con `check_manufacturability`.
_NUMERICAL_SLOT_FLOOR_MM = 0.05


@dataclass(frozen=True)
class BuildOptions:
    sector: bool = False              # True -> ritaglia il settore periodico 1/N
    stl_tolerance: float = 5.0e-6     # m, tolleranza lineare di tassellazione
    stl_angular_tolerance: float = 0.2  # rad
    include_injection: bool = True


def wall_profile_mm(p: GeometryParams) -> list[tuple[float, float, float]]:
    """Il poligono meridiano di `geometry.profile`, in mm e in forma 3D.

    Il poligono e' definito UNA VOLTA SOLA, in `profile.meridian_polygon`, che
    non importa OCCT. Averne due copie sarebbe il modo piu' rapido di far
    divergere il volume che l'ottimizzatore minimizza da quello che la
    stampante stampa.
    """
    poly = meridian_polygon(p.derived, p.plug_contour_x, p.plug_contour_r)
    return [(x * _MM, 0.0, r * _MM) for x, r in poly]


def wetted_area(p: GeometryParams) -> float:
    """Superficie bagnata dai gas [m^2]. Delega a `geometry.profile`, che non
    importa OCCT ed e' quindi usabile anche sul percorso veloce."""
    d = p.derived
    return wetted_area_from_contour(
        d["R_c"], d["R_lip"], d["L_c"], d["L_conv"], d["r_centerbody"],
        p.plug_contour_x, p.plug_contour_r,
    )


def build_solid(p: GeometryParams, opts: BuildOptions = BuildOptions()):
    """Costruisce il solido (mm). Ritorna l'oggetto build123d."""
    d = p.derived
    solid = revolve(
        make_face(Polyline(*wall_profile_mm(p), close=True)),
        axis=Axis.X,
        revolution_arc=360.0,
    )

    if opts.include_injection:
        solid = _cut_injection(solid, p)

    if opts.sector:
        N = int(d["N_inj"])
        span = 360.0 / N
        # cuneo abbondantemente piu' grande del pezzo, poi intersezione
        R_big = 4.0 * (d["R_c"] + d["t_wall"]) * _MM
        L_big = 4.0 * (d["L_c"] + d["L_conv"] + abs(d["x_tip_full"])) * _MM
        wedge = _wedge(R_big, L_big, span, x0=-2.0 * d["t_face"] * _MM)
        solid = solid & wedge
    return solid


def _wedge(radius: float, length: float, angle_deg: float, x0: float):
    """Cuneo angolare [0, angle_deg] attorno all'asse X, come solido di rivoluzione."""
    pts = [(x0, 0.0, 0.0), (x0, 0.0, radius), (x0 + length, 0.0, radius), (x0 + length, 0.0, 0.0)]
    return revolve(make_face(Polyline(*pts, close=True)), axis=Axis.X, revolution_arc=angle_deg)


def _cut_injection(solid, p: GeometryParams):
    """Fori ossidante (con swirl), fori combustibile e fessura di film cooling."""
    d = p.derived
    N = int(d["N_inj"])
    tf = d["t_face"] * _MM
    R_inj = d["R_inj"] * _MM
    d_ox = d["d_ox"] * _MM
    d_f = d["d_fuel"] * _MM
    theta_s = p.free["theta_swirl"]
    span = 2.0 * math.pi / N
    depth = 3.0 * tf / max(math.cos(theta_s), 0.2)

    for k in range(N):
        # ox a 1/4 del passo, fuel a 3/4: entrambi interni al settore [0, span)
        phi_ox = k * span + 0.25 * span
        phi_f = k * span + 0.75 * span
        solid = solid - _hole(R_inj, phi_ox, d_ox, depth, tf, theta_s)
        solid = solid - _hole(R_inj, phi_f, d_f, depth, tf, 0.0)

    d_film = d["d_film"] * _MM
    if d_film > _NUMERICAL_SLOT_FLOOR_MM:
        R_film = d["R_film"] * _MM
        for k in range(N):
            phi = k * span + 0.5 * span      # a meta' passo: interno al settore
            solid = solid - _hole(R_film, phi, d_film, depth, tf, 0.0)
    return solid


def _hole(R_inj: float, phi: float, diameter: float, depth: float, tf: float, swirl: float):
    """Foro cilindrico, eventualmente inclinato tangenzialmente di `swirl`.

    L'inclinazione e' attorno alla direzione RADIALE locale: produce una
    componente tangenziale di velocita' (swirl) senza componente radiale, che
    e' cio' che si vuole da un iniettore a swirl assiale.
    """
    cyl = Cylinder(radius=0.5 * diameter, height=depth, rotation=(0.0, 90.0, 0.0))
    cyl = Rotation(0.0, 0.0, 0.0) * cyl
    if swirl != 0.0:
        cyl = Rotation(0.0, 0.0, math.degrees(swirl)) * cyl
    cyl = Rotation(math.degrees(phi), 0.0, 0.0) * cyl
    pos = (
        -tf,
        R_inj * math.cos(phi),
        R_inj * math.sin(phi),
    )
    return Location(pos) * cyl


# --------------------------------------------------------------------------- #
# Verifica di water-tightness della MESH, indipendente da OCCT
# --------------------------------------------------------------------------- #
def stl_is_watertight(path: Path, quantum: float = 1.0e-6) -> tuple[bool, int]:
    """Vero se ogni spigolo della tassellazione e' condiviso da 2 triangoli.

    Verifica volutamente INDIPENDENTE da `Shape.is_valid` di OCCT: un B-Rep
    valido puo' tassellare in modo non chiuso se la tolleranza e' troppo lasca,
    ed e' esattamente quel caso a far fallire il mesher a valle.

    `quantum` [mm] e' la griglia su cui si arrotondano i vertici prima del
    confronto: serve perche' due triangoli adiacenti possono avere lo stesso
    vertice con ultimo bit diverso.
    """
    data = path.read_bytes()
    if data[:5].lstrip().lower().startswith(b"solid") and b"facet normal" in data[:2048]:
        tris = _parse_ascii_stl(data)
    else:
        tris = _parse_binary_stl(data)

    edges: dict[tuple, int] = {}
    for tri in tris:
        q = [tuple(round(c / quantum) for c in v) for v in tri]
        for a, b in ((q[0], q[1]), (q[1], q[2]), (q[2], q[0])):
            key = (a, b) if a <= b else (b, a)
            edges[key] = edges.get(key, 0) + 1
    return all(c == 2 for c in edges.values()), len(tris)


def _parse_binary_stl(data: bytes) -> list[tuple[tuple[float, ...], ...]]:
    (n,) = struct.unpack("<I", data[80:84])
    out = []
    off = 84
    for _ in range(n):
        vals = struct.unpack("<12fH", data[off:off + 50])
        out.append((vals[3:6], vals[6:9], vals[9:12]))
        off += 50
    return out


def _parse_ascii_stl(data: bytes) -> list[tuple[tuple[float, ...], ...]]:
    out, cur = [], []
    for line in data.decode("ascii", "replace").splitlines():
        s = line.strip()
        if s.startswith("vertex"):
            cur.append(tuple(float(v) for v in s.split()[1:4]))
            if len(cur) == 3:
                out.append(tuple(cur))
                cur = []
    return out


#: Timestamp fisso scritto nell'header STEP. OCCT ci mette l'ora di creazione,
#: che rende due esportazioni della STESSA geometria byte-diverse e fa fallire
#: il criterio di riproducibilita' (docs/architettura.md sezione 5). L'ora non
#: va persa: sta in `created_utc` nel database, dove e' un dato e non rumore.
_STEP_EPOCH = "1970-01-01T00:00:00"


def normalize_step_header(path: Path, run_id: str) -> None:
    """Rende l'header STEP deterministico e autoidentificante.

    Sostituisce il nome del modello col `run_id` e il timestamp con un
    sentinella fisso. Effetto: due export della stessa geometria danno file
    byte-identici, e aprendo il file si sa a quale run appartiene senza
    consultare nulla.

    Tocca SOLO l'header: la sezione DATA, cioe' la geometria, non viene
    sfiorata.
    """
    text = path.read_text(encoding="utf-8", errors="surrogateescape")
    head, sep, rest = text.partition("ENDSEC;")
    if not sep:
        raise ValueError(f"{path} non sembra un file STEP: manca ENDSEC nell'header.")
    head = re.sub(
        r"FILE_NAME\('[^']*','[^']*'",
        f"FILE_NAME('zefiro {run_id}','{_STEP_EPOCH}'",
        head,
        count=1,
    )
    path.write_text(head + sep + rest, encoding="utf-8", errors="surrogateescape")


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def build_and_export(
    p: GeometryParams,
    out_dir: Path,
    run_id: str,
    opts: BuildOptions = BuildOptions(),
) -> GeometryArtifact:
    """Costruisce, esporta STEP + STL e verifica la tenuta due volte."""
    out_dir.mkdir(parents=True, exist_ok=True)
    solid = build_solid(p, opts)

    step_path = out_dir / f"{run_id}_geometry.step"
    stl_path = out_dir / f"{run_id}_geometry.stl"
    export_step(solid, str(step_path))
    normalize_step_header(step_path, run_id)
    export_stl(
        solid,
        str(stl_path),
        tolerance=opts.stl_tolerance * _MM,
        angular_tolerance=opts.stl_angular_tolerance,
    )

    watertight, n_tri = stl_is_watertight(stl_path)
    return GeometryArtifact(
        run_id=run_id,
        step_path=step_path,
        stl_path=stl_path,
        params=p,
        volume=float(solid.volume) / _MM**3,      # mm^3 -> m^3
        wetted_area=wetted_area(p),
        is_valid_brep=bool(solid.is_valid),
        is_watertight_mesh=watertight,
        n_triangles=n_tri,
        mesh_tolerance=opts.stl_tolerance,
        sha256_step=_sha256(step_path),
        sha256_stl=_sha256(stl_path),
    )
