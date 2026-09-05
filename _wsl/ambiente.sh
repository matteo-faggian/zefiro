# ambiente OpenFOAM: si attiva l'ambiente conda, non si sorgente il bashrc a
# mano. Gli script di attivazione impostano WM_MPLIB=MPICH e il FOAM_MPI
# corrispondente: facendolo a mano si finisce con la libPstream "dummy", che
# in parallelo si rifiuta di partire (ed e' giusto che lo faccia).
export MAMBA_ROOT_PREFIX=/root/mm
eval "$(/root/mm/bin/micromamba shell hook -s bash)"
micromamba activate /root/of
