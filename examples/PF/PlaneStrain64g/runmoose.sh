
NPROC="$1"
RHOCPOPT=~/projects/rhocp-vr/rhocp-opt

if [ "$2" = "1" ]; then
    mpirun -n "$NPROC" "$RHOCPOPT" -i 2dplanestain.i --recover
else
    mpirun -n "$NPROC" "$RHOCPOPT" -i 2dplanestain.i
fi
