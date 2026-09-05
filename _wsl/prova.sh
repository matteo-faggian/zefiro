#!/bin/bash
source /root/of/etc/bashrc
echo "WM_PROJECT_VERSION=$WM_PROJECT_VERSION"
echo "WM_PROJECT_DIR=$WM_PROJECT_DIR"
which reactingFoam gmshToFoam foamToVTK decomposePar
reactingFoam -help 2>&1 | head -3
echo "--- funzioni ---"
ls $FOAM_LIBBIN/libutilityFunctionObjects.so 2>/dev/null && echo "functionObjects presenti"
