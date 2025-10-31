mpirun -n 24 ~/projects/rhocp-vr-250420/rhocp-opt -i bcc_pxtal.i \
Materials/CPStressUpdate/propsFile=bcc_props.in \
Outputs/file_base=out_473K_3dpa

mpirun -n 24 ~/projects/rhocp-vr-250420/rhocp-opt -i bcc_pxtal.i \
Materials/CPStressUpdate/propsFile=bcc_propsLoops.in \
Outputs/file_base=out_473K_3dpa_Loops

mpirun -n 24 ~/projects/rhocp-vr-250420/rhocp-opt -i bcc_pxtal.i \
Materials/CPStressUpdate/propsFile=bcc_propsPreci.in \
Outputs/file_base=out_473K_3dpa_preci
