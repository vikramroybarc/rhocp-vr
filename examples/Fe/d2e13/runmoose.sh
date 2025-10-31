export RHOCPOPT=~/projects/rhocp-vr-250922/rhocp-opt

# cd /home/vikramroy/projects/rhocp-vr-250922/examples/Fe/400K_SR100_km_0p002
# mpirun -n 24 $RHOCPOPT -i bcc_pxtal.i --recover

# cd /home/vikramroy/projects/rhocp-vr-250922/examples/Fe/400K_SR10
# mpirun -n 24 $RHOCPOPT -i bcc_pxtal.i --recover

# cd /home/vikramroy/projects/rhocp-vr-250922/examples/Fe/400K_SR20
# mpirun -n 24 $RHOCPOPT -i bcc_pxtal.i --recover

cd /home/vikramroy/projects/rhocp-vr-250922/examples/Fe/400K_SR40
mpirun -n 24 $RHOCPOPT -i bcc_pxtal.i --recover