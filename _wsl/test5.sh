#!/bin/bash
cd /root/zefiro
export PYTHONPATH=/root/zefiro/src:/root/zefiro/scripts/cfd
# build123d (geometria solida) non e' installato qui: quei tre moduli si
# provano nel container. Il resto, CFD compresa, gira in questo ambiente.
/root/zef/bin/python -m pytest -q -p no:cacheprovider \
    --ignore=tests/test_geometry_watertight.py \
    --ignore=tests/test_integration_l0_geometry.py \
    --ignore=tests/test_reproducibility.py \
    > /root/test5.out 2>&1
echo "EXIT=$?"
tail -25 /root/test5.out
