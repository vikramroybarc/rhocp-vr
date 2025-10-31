#!/bin/bash

# Define Simulation Parameters for Thermal Creep Simulation

NPROCS=24
STRESS=(9 26 39 53 66)  # Stress values in MPa
TIME1=(5e5 5e5 5e5 5e5 5e5)     # End times for initial simulations
TIME2=(7e7 7e7 7e7 7e7 7e7)     # End times for recovered simulations
TEMP=873

# # Iterate over the array index to pair STRESS, TIME1, and TIME2 correctly
# for i in "${!STRESS[@]}"
# do
#     stress=${STRESS[$i]}
#     time1=${TIME1[$i]}
#     time2=${TIME2[$i]}

#     echo "Running initial simulation for STRESS=${stress} MPa, TEMP=${TEMP} K, end_time=${time1}"

#     # Run the initial simulation
#     mpirun -n "$NPROCS" ~/projects/rhocp-vr-250401/rhocp-opt -i bcc_pxtal.i \
#         sigmaV="$stress" \
#         Materials/CPStressUpdate/temp="$TEMP" \
#         Materials/CPStressUpdate/climbmodel=true \
#         Executioner/end_time="$time1" \
#         Executioner/TimeStepper/log_dt=0.050 \
#         Outputs/file_base="out_Creep_${TEMP}K_${stress}MPa" > "log${stress}.run"

#     echo "Initial run complete. Starting recovered simulation with end_time=${time2}"

#     # Run the recovered simulation
#     mpirun -n "$NPROCS" ~/projects/rhocp-vr-250401/rhocp-opt -i bcc_pxtal.i \
#         sigmaV="$stress" \
#         Materials/CPStressUpdate/temp="$TEMP" \
#         Materials/CPStressUpdate/climbmodel=true \
#         Executioner/end_time="$time2" \
#         Executioner/TimeStepper/log_dt=0.01 \
#         Outputs/file_base="out_Creep_${TEMP}K_${stress}MPa" --recover > "log${stress}r1.run"
# done


# cd /home/vikramroy/projects/rhocp-vr-250420/examples/bcc_FeCr/Creep/512g/873K/UTS_7/
# ./runmoose.sh


    echo "Running initial simulation for STRESS=${stress} MPa, TEMP=${TEMP} K, end_time=${time1}"

    stress=${STRESS[0]}
    time1=${TIME1[0]}
    time2=${TIME2[0]}
    # Run the initial simulation
    mpirun -n "$NPROCS" ~/projects/rhocp-vr-250401/rhocp-opt -i bcc_pxtal.i \
        sigmaV="$stress" \
        Materials/CPStressUpdate/temp="$TEMP" \
        Materials/CPStressUpdate/climbmodel=true \
        Executioner/end_time="$time1" \
        Executioner/TimeStepper/log_dt=0.050 \
        Outputs/file_base="out_Creep_${TEMP}K_${stress}MPa" > "log${stress}.run"

    echo "Initial run complete. Starting recovered simulation with end_time=${time2}"

    # Run the recovered simulation
    mpirun -n "$NPROCS" ~/projects/rhocp-vr-250401/rhocp-opt -i bcc_pxtal.i \
        sigmaV="$stress" \
        Materials/CPStressUpdate/temp="$TEMP" \
        Materials/CPStressUpdate/climbmodel=true \
        Executioner/end_time="$time2" \
        Executioner/TimeStepper/log_dt=0.01 \
        Outputs/file_base="out_Creep_${TEMP}K_${stress}MPa" --recover > "log${stress}r1.run"