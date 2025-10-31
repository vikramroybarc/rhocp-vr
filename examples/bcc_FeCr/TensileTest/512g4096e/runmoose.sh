NPROC=20

# for temp in {723..923..50}; do
#     mpirun -n "$NPROC" ~/projects/rhocp-vr-250902/rhocp-opt -i bcc_pxtal.i \
#         Functions/top_pull/expression='0.8*3e-3' \
#         Executioner/end_time=17 \
#         Materials/CPStressUpdate/temp=$temp \
#         Materials/CPStressUpdate/deltaH_eV=false\
#         Materials/CPStressUpdate/propsFile=bcc_props.in \
#         Outputs/file_base=out_${temp}K_3e-3
# done

# for temp in {573..573..100}; do
    temp=473
    mpirun -n "$NPROC" ~/projects/rhocp-vr-250902/rhocp-opt -i bcc_pxtal.i \
        Functions/top_pull/expression='0.8*3e-3' \
        Executioner/end_time=18 \
        Materials/CPStressUpdate/temp=$temp \
        Materials/CPStressUpdate/deltaH_eV=true\
        Materials/CPStressUpdate/propsFile=bcc_props_473.in \
        Functions/dts/y='0.0001    0.01' \
        Outputs/file_base=out_${temp}K_3e-3 --recover

    temp=573
    mpirun -n "$NPROC" ~/projects/rhocp-vr-250902/rhocp-opt -i bcc_pxtal.i \
        Functions/top_pull/expression='0.8*3e-3' \
        Executioner/end_time=18 \
        Materials/CPStressUpdate/temp=$temp \
        Materials/CPStressUpdate/deltaH_eV=true\
        Materials/CPStressUpdate/propsFile=bcc_props_473.in \
        Functions/dts/y='0.0001    0.01' \
        Outputs/file_base=out_${temp}K_3e-3 
# done




# for temp in {873..923..50}; do
#     mpirun -n "$NPROC" ~/projects/rhocp-vr-250902/rhocp-opt -i bcc_pxtal.i \
#         Functions/top_pull/expression='0.8*3e-3' \
#         Executioner/end_time=12 \
#         Materials/CPStressUpdate/temp=$temp \
#         Materials/CPStressUpdate/deltaH_eV=false\
#         Materials/CPStressUpdate/propsFile=bcc_props.in \
#         Outputs/file_base=out_${temp}K_3e-3
# done        

# cd /home/vikramroy/projects/rhocp-vr-250420/examples/bcc_FeCr/Irad_TensileTest/512g/473K_0dpa
# mpirun -n "$NPROC" ~/projects/rhocp-vr-250420/rhocp-opt -i bcc_pxtal.i 

# cd /home/vikramroy/projects/rhocp-vr-250420/examples/bcc_FeCr/Irad_TensileTest/512g/473K_3dpa
# mpirun -n "$NPROC" ~/projects/rhocp-vr-250420/rhocp-opt -i bcc_pxtal.i 

# cd /home/vikramroy/projects/rhocp-vr-250420/examples/bcc_FeCr/Irad_TensileTest/512g/573K_0dpa
# mpirun -n "$NPROC" ~/projects/rhocp-vr-250420/rhocp-opt -i bcc_pxtal.i 
    
    

# for temp in {823..923..50}; do
#     mpirun -n "$NPROC" ~/projects/rhocp-vr-250902/rhocp-opt -i bcc_pxtal.i \
#         Functions/top_pull/expression='0.8*1.26e-3' \
#         Executioner/end_time=10 \
#         Materials/CPStressUpdate/temp=$temp \
#         Materials/CPStressUpdate/deltaH_eV=false\
#         Materials/CPStressUpdate/propsFile=bcc_props.in \
#         Functions/dts/y='0.0001    0.04' \
#         Outputs/file_base=out_${temp}K_1.26e-3
# done

# for temp in {473..673..100}; do
#     mpirun -n "$NPROC" ~/projects/rhocp-vr-250902/rhocp-opt -i bcc_pxtal.i \
#         Functions/top_pull/expression='0.8*1.26e-3' \
#         Executioner/end_time=10 \
#         Materials/CPStressUpdate/temp=$temp \
#         Materials/CPStressUpdate/deltaH_eV=true\
#         Materials/CPStressUpdate/propsFile=bcc_props_473.in \
#         Functions/dts/y='0.0001    0.02' \
#         Outputs/file_base=out_${temp}K_1.26e-3
# done



# for temp in {473..673..100}; do
#     mpirun -n "$NPROC" ~/projects/rhocp-vr-250902/rhocp-opt -i bcc_pxtal.i \
#         Functions/top_pull/expression='0.4*5.e-5' \
#         Executioner/end_time=500 \
#         Materials/CPStressUpdate/temp=$temp \
#         Materials/CPStressUpdate/deltaH_eV=true\
#         Materials/CPStressUpdate/propsFile=bcc_props_473.in \
#         Functions/dts/y='0.0001    0.4' \
#         Outputs/file_base=out_${temp}K_5e-5
# done


# for temp in {723..973..50}; do
#     mpirun -n "$NPROC" ~/projects/rhocp-vr-250902/rhocp-opt -i bcc_pxtal.i \
#         Functions/top_pull/expression='0.4*5.e-5' \
#         Executioner/end_time=500 \
#         Materials/CPStressUpdate/temp=$temp \
#         Materials/CPStressUpdate/deltaH_eV=false\
#         Materials/CPStressUpdate/propsFile=bcc_props.in \
#         Functions/dts/y='0.0001    0.4' \
#         Outputs/file_base=out_${temp}K_5.e-5
# done




# Iradiation Hardening Simulations

# cd ~/projects/rhocp-vr-250420/examples/bcc_FeCr/Irad_TensileTest/512g/473K_0dpa/
# mpirun -n 24 ~/projects/rhocp-vr-250902/rhocp-opt -i bcc_pxtal.i \
# Materials/CPStressUpdate/propsFile=bcc_props.in \
# Outputs/file_base=out_473K_0dpa




# cd ~/projects/rhocp-vr-250420/examples/bcc_FeCr/Irad_TensileTest/512g/473K_3dpa/
# mpirun -n 24 ~/projects/rhocp-vr-250902/rhocp-opt -i bcc_pxtal.i \
# Materials/CPStressUpdate/propsFile=bcc_props.in \
# Outputs/file_base=out_473K_3dpa