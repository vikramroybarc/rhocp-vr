
NPROC="$1"
RHOCPOPT=~/projects/rhocp-vr/rhocp-opt

if [ "$2" = "1" ]; then
    mpirun -n "$NPROC" "$RHOCPOPT" -i bcc_pxtal.i --recover
else
    mpirun -n "$NPROC" "$RHOCPOPT" -i bcc_pxtal.i
fi
