mpirun -n 20 ~/projects/rhocp-vr-250505/rhocp-opt -i bcc_pxtal.i \
Materials/CPStressUpdate/propsFile=bcc_props_0p6.in \
Outputs/file_base=out_573K_0p6dpa

mpirun -n 20 ~/projects/rhocp-vr-250505/rhocp-opt -i bcc_pxtal.i \
Materials/CPStressUpdate/propsFile=bcc_props_1p5.in \
Outputs/file_base=out_573K_1p5dpa 

mpirun -n 20 ~/projects/rhocp-vr-250505/rhocp-opt -i bcc_pxtal.i \
Materials/CPStressUpdate/propsFile=bcc_props_0p06.in \
Outputs/file_base=out_573K_0p06dpa

cd ~/projects/rhocp-vr-250505/examples/bcc_FeCr/Irad_TensileTest/512g/573K_0dpa

mpirun -n 20 ~/projects/rhocp-vr-250505/rhocp-opt -i bcc_pxtal.i \
Materials/CPStressUpdate/propsFile=bcc_props.in \
Outputs/file_base=out_573K_0dpa --recover
