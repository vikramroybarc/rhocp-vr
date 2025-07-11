# mpirun -n 24 ~/projects/rhocp-vr-250401/rhocp-opt -i bcc_pxtal.i \
# Functions/top_pull/expression='0.8*3e-3' \
# Executioner/end_time=17 \
# Materials/CPStressUpdate/temp=473 \
# Materials/CPStressUpdate/propsFile=bcc_props_473.in \
# Outputs/file_base=out_473K_3e-3

# mpirun -n 24 ~/projects/rhocp-vr-250401/rhocp-opt -i bcc_pxtal.i \
# Functions/top_pull/expression='0.8*3e-3' \
# Executioner/end_time=17 \
# Materials/CPStressUpdate/temp=673 \
# Materials/CPStressUpdate/propsFile=bcc_props_673.in \
# Outputs/file_base=out_673K_3e-3

# mpirun -n 24 ~/projects/rhocp-vr-250401/rhocp-opt -i bcc_pxtal.i \
# Functions/top_pull/expression='0.8*3e-3' \
# Executioner/end_time=17 \
# Materials/CPStressUpdate/temp=573 \
# Materials/CPStressUpdate/propsFile=bcc_props_573.in \
# Outputs/file_base=out_573K_3e-3


# mpirun -n 24 ~/projects/rhocp-vr-250401/rhocp-opt -i bcc_pxtal.i \
# Functions/top_pull/expression='0.8*1.26e-3' \
# Functions/dts/y='0.0001    0.04' \
# Executioner/end_time=10 \
# Materials/CPStressUpdate/temp=473 \
# Materials/CPStressUpdate/propsFile=bcc_props_473.in \
# Outputs/file_base=out_473K_1.26e-3

# mpirun -n 24 ~/projects/rhocp-vr-250401/rhocp-opt -i bcc_pxtal.i \
# Functions/top_pull/expression='0.8*1.26e-3' \
# Functions/dts/y='0.0001    0.04' \
# Executioner/end_time=10 \
# Materials/CPStressUpdate/temp=673 \
# Materials/CPStressUpdate/propsFile=bcc_props_673.in \
# Outputs/file_base=out_673K_1.26e-3

# mpirun -n 24 ~/projects/rhocp-vr-250401/rhocp-opt -i bcc_pxtal.i \
# Functions/top_pull/expression='0.8*1.26e-3' \
# Functions/dts/y='0.0001    0.04' \
# Executioner/end_time=10 \
# Materials/CPStressUpdate/temp=573 \
# Materials/CPStressUpdate/propsFile=bcc_props_573.in \
# Outputs/file_base=out_573K_1.26e-3

# mpirun -n 24 ~/projects/rhocp-vr-250401/rhocp-opt -i bcc_pxtal.i \
# Functions/top_pull/expression='0.8*1.e-4' \
# Functions/dts/y='0.0001    0.2' \
# Executioner/end_time=200 \
# Materials/CPStressUpdate/temp=573 \
# Materials/CPStressUpdate/propsFile=bcc_props_573.in \
# Outputs/file_base=out_573K_1.e-4

# mpirun -n 24 ~/projects/rhocp-vr-250401/rhocp-opt -i bcc_pxtal.i \
# Functions/top_pull/expression='0.8*1.e-4' \
# Functions/dts/y='0.0001    0.2' \
# Executioner/end_time=200 \
# Materials/CPStressUpdate/temp=473 \
# Materials/CPStressUpdate/propsFile=bcc_props_473.in \
# Outputs/file_base=out_473K_1.e-4


mpirun -n 24 ~/projects/rhocp-vr-250401/rhocp-opt -i bcc_pxtal.i \
Functions/top_pull/expression='0.8*1.26e-3' \
Functions/dts/y='0.0001    0.04' \
Executioner/end_time=10 \
Materials/CPStressUpdate/temp=723 \
Materials/CPStressUpdate/propsFile=bcc_props.in \
Outputs/file_base=out_723K_1.26e-3

mpirun -n 24 ~/projects/rhocp-vr-250401/rhocp-opt -i bcc_pxtal.i \
Functions/top_pull/expression='0.8*1.26e-3' \
Functions/dts/y='0.0001    0.04' \
Executioner/end_time=10 \
Materials/CPStressUpdate/temp=773 \
Materials/CPStressUpdate/propsFile=bcc_props.in \
Outputs/file_base=out_773K_1.26e-3

mpirun -n 24 ~/projects/rhocp-vr-250401/rhocp-opt -i bcc_pxtal.i \
Functions/top_pull/expression='0.8*1.26e-3' \
Functions/dts/y='0.0001    0.04' \
Executioner/end_time=10 \
Materials/CPStressUpdate/temp=823 \
Materials/CPStressUpdate/propsFile=bcc_props.in \
Outputs/file_base=out_823K_1.26e-3

mpirun -n 24 ~/projects/rhocp-vr-250401/rhocp-opt -i bcc_pxtal.i \
Functions/top_pull/expression='0.8*1.26e-3' \
Functions/dts/y='0.0001    0.04' \
Executioner/end_time=10 \
Materials/CPStressUpdate/temp=873 \
Materials/CPStressUpdate/propsFile=bcc_props.in \
Outputs/file_base=out_873K_1.26e-3

