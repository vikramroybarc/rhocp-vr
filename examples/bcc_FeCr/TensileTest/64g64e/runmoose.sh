NPROC=24
RHOCPOPT=~/projects/rhocp-vr-250926/rhocp-opt

# for temp in {773..923..50}; do
#     mpirun -n "$NPROC" $RHOCPOPT -i bcc_pxtal.i \
#         Functions/top_pull/expression='0.4*3e-3' \
#         Executioner/end_time=17 \
#         Materials/CPStressUpdate/temp=$temp \
#         Materials/CPStressUpdate/deltaH_eV=false\
#         Materials/CPStressUpdate/propsFile=bcc_props.in \
#         Outputs/file_base=out_${temp}K_3e-3
# done

# for temp in {473..673..100}; do
#     mpirun -n "$NPROC" $RHOCPOPT -i bcc_pxtal.i \
#         Functions/top_pull/expression='0.4*3e-3' \
#         Executioner/end_time=17 \
#         Materials/CPStressUpdate/temp=$temp \
#         Materials/CPStressUpdate/deltaH_eV=true\
#         Materials/CPStressUpdate/propsFile=bcc_props_473.in \
#         Functions/dts/y='0.0001    0.02' \
#         Outputs/file_base=out_${temp}K_3e-3
# done


for temp in {673..923..50}; do
    mpirun -n "$NPROC" $RHOCPOPT -i bcc_pxtal.i \
        Functions/top_pull/expression='0.4*3e-3' \
        Executioner/end_time=17 \
        Materials/CPStressUpdate/temp=$temp \
        Materials/CPStressUpdate/deltaH_eV=false\
        Materials/CPStressUpdate/propsFile=bcc_props.in \
        Outputs/file_base=out_${temp}K_3e-3
done

# for temp in {723..923..50}; do
#     mpirun -n "$NPROC" ~/projects/rhocp-vr-250401/rhocp-opt -i bcc_pxtal.i \
#         Functions/top_pull/expression='0.4*1.26e-3' \
#         Executioner/end_time=17 \
#         Materials/CPStressUpdate/temp=$temp \
#         Materials/CPStressUpdate/deltaH_eV=false\
#         Materials/CPStressUpdate/propsFile=bcc_props.in \
#         Functions/dts/y='0.0001    0.04' \
#         Outputs/file_base=out_${temp}K_1.26e-3
# done

# for temp in {473..673..100}; do
#     mpirun -n "$NPROC" ~/projects/rhocp-vr-250401/rhocp-opt -i bcc_pxtal.i \
#         Functions/top_pull/expression='0.4*1.26e-3' \
#         Executioner/end_time=17 \
#         Materials/CPStressUpdate/temp=$temp \
#         Materials/CPStressUpdate/deltaH_eV=true\
#         Materials/CPStressUpdate/propsFile=bcc_props_473.in \
#         Functions/dts/y='0.0001    0.04' \
#         Outputs/file_base=out_${temp}K_1.26e-3
# done



# for temp in {473..673..100}; do
#     mpirun -n "$NPROC" ~/projects/rhocp-vr-250401/rhocp-opt -i bcc_pxtal.i \
#         Functions/top_pull/expression='0.4*5.e-5' \
#         Executioner/end_time=500 \
#         Materials/CPStressUpdate/temp=$temp \
#         Materials/CPStressUpdate/deltaH_eV=true\
#         Materials/CPStressUpdate/propsFile=bcc_props_473.in \
#         Functions/dts/y='0.0001    0.4' \
#         Outputs/file_base=out_${temp}K_5e-5
# done


# for temp in {723..973..50}; do
#     mpirun -n "$NPROC" ~/projects/rhocp-vr-250401/rhocp-opt -i bcc_pxtal.i \
#         Functions/top_pull/expression='0.4*5.e-5' \
#         Executioner/end_time=500 \
#         Materials/CPStressUpdate/temp=$temp \
#         Materials/CPStressUpdate/deltaH_eV=false\
#         Materials/CPStressUpdate/propsFile=bcc_props.in \
#         Functions/dts/y='0.0001    0.4' \
#         Outputs/file_base=out_${temp}K_5.e-5
# done
