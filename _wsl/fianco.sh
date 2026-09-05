#!/bin/bash
# IL PIANO "A META' FRA DUE GETTI" ERA VUOTO, E NON PER CASO. Il dominio e' un
# settore di 90 gradi con il getto sulla BISETTRICE: i due fianchi del settore
# sono piani di simmetria, e uno di essi E' il piano z = 0. Chiedere un
# `cuttingPlane` con normale (0 0 1) passante per l'origine significa tagliare
# esattamente sulla faccia di bordo: non ci sono celle da intersecare e restano
# quattro triangoli degeneri. Il piano giusto e' lo stesso spostato dentro il
# dominio di una frazione di cella: a r = 4 mm uno scarto di 0.2 mm vale 2.9
# gradi contro i 45 che lo separano dal getto, quindi resta "meta' strada".
# Si ricampionano ENTRAMBI i piani allo stesso istante: confrontare due piani
# presi a istanti diversi sarebbe un confronto truccato.
set +e
source /mnt/i/AA_ENGINE/_wsl/ambiente.sh
C=/root/zefiro/runs/cfd/run4/caso
T=0.00025
cd "$C"
#: la normale del piano del getto si RILEGGE dal controlDict del caso, non si
#: riscrive a mano: due copie dello stesso numero divergono sempre.
NG=$(grep -A8 "meridiano_getto" system/controlDict | sed -n 's/.*normal *(\([^)]*\)).*/\1/p' | head -1)
[ -z "$NG" ] && { echo "NORMALE NON TROVATA nel controlDict"; exit 2; }
echo "normale del piano del getto: ($NG)"
cat > system/rifai_superfici <<DICT
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      rifai_superfici;
}
functions
{
superfici2
{
    type            surfaces;
    libs            (sampling);
    surfaceFormat   vtk;
    writeFormat     binary;
    interpolationScheme cellPoint;
    fields          (C3H8 O2 T U);
    surfaces
    {
        meridiano_getto
        {
            type            cuttingPlane;
            planeType       pointAndNormal;
            pointAndNormalDict { point (0 0 0); normal ($NG); }
            interpolate     true;
        }
        meridiano_fianco
        {
            type            cuttingPlane;
            planeType       pointAndNormal;
            pointAndNormalDict { point (0 0 0.0002); normal (0 0 1); }
            interpolate     true;
        }
    }
}
}
DICT
mpirun -np 24 postProcess -parallel -time "$T" -fields "(C3H8 O2 T U)" -dict system/rifai_superfici > log.fianco 2>&1
echo "rc=$?"
grep -iE "error|cannot|--> FOAM" log.fianco | head -5
ls -la postProcessing/superfici2/$T 2>/dev/null
