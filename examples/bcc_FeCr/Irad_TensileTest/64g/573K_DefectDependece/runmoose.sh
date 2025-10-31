# mpirun -n 24 ~/projects/rhocp-vr-250926/rhocp-opt -i bcc_pxtal.i \
# Materials/CPStressUpdate/propsFile=bcc_props_1p5.in \
# Outputs/file_base=out_573K_1p5dpa 


# mpirun -n 24 ~/projects/rhocp-vr-250926/rhocp-opt -i bcc_pxtal.i \
# Materials/CPStressUpdate/propsFile=bcc_props_1p5Loops.in \
# Outputs/file_base=out_573K_1p5dpaLoops



mpirun -n 24 ~/projects/rhocp-vr-250926/rhocp-opt -i bcc_pxtal.i \
Materials/CPStressUpdate/propsFile=bcc_props_1p5preci.in \
Outputs/file_base=out_573K_1p5dpapreci


# mpirun -n 24 ~/projects/rhocp-vr-250926/rhocp-opt -i bcc_pxtal.i \
# Materials/CPStressUpdate/propsFile=bcc_props_1p5voids.in \
# Outputs/file_base=out_573K_1p5dpavoids