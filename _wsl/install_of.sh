#!/bin/bash
export MAMBA_ROOT_PREFIX=/root/mm
cd /root/mm
./bin/micromamba create -y -p /root/of -c conda-forge openfoam=2412
echo "RC_OPENFOAM=$?"
./bin/micromamba create -y -p /root/zef -c conda-forge python=3.12 numpy scipy pyyaml cantera coolprop gmsh python-gmsh scikit-image matplotlib
echo "RC_ZEF=$?"
echo FATTO
