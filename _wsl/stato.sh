#!/bin/bash
tail -c 300 /root/mm/install.log
echo
echo "--- zef ---"
ls /root/zef/bin/python 2>/dev/null && /root/zef/bin/python -c "import numpy,scipy,yaml,cantera,gmsh;print('py ok', cantera.__version__, gmsh.__file__)"
echo "--- libbin ---"
source /root/of/etc/bashrc >/dev/null 2>&1
echo "FOAM_LIBBIN=$FOAM_LIBBIN"
ls "$FOAM_LIBBIN" 2>/dev/null | grep -i functionobj
