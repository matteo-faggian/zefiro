#!/bin/bash
/root/zef/bin/python -m pip install --quiet meshio 2>&1 | tail -3
/root/zef/bin/python -c "import meshio; print('meshio', meshio.__version__)"
