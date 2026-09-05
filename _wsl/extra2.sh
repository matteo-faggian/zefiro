#!/bin/bash
export MAMBA_ROOT_PREFIX=/root/mm
/root/mm/bin/micromamba install -y -p /root/zef -c conda-forge vtk ffmpeg 2>&1 | tail -5
/root/zef/bin/python -c "import vtk; print('vtk', vtk.vtkVersion.GetVTKVersion())"
