"""Geometria parametrica: dal vettore di progetto a STEP/STL."""
from zefiro.geometry.parameters import DESIGN_BOUNDS, derive  # noqa: F401
from zefiro.geometry.build import (  # noqa: F401
    BuildOptions,
    build_and_export,
    build_solid,
    stl_is_watertight,
    wetted_area,
)
