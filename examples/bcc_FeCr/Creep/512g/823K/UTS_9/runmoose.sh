#!/bin/bash

# Define Simulation Parameters for Thermal Creep Simulation

NPROCS=24
STRESS=(215)
TIME1=(5e4)     # End times for initial simulations
TIME2=(1e6)     # End times for recovered simulations
TEMP=823

# Iterate over the array index to pair STRESS, TIME1, and TIME2 correctly
for i in "${!STRESS[@]}"
do
    stress=${STRESS[$i]}
    time1=${TIME1[$i]}
    time2=${TIME2[$i]}

    echo "Running initial simulation for STRESS=${stress} MPa, TEMP=${TEMP} K, end_time=${time1}"

    # Run the initial simulation
    mpirun -n "$NPROCS" ~/projects/rhocp-vr-250401/rhocp-opt -i bcc_pxtal.i \
        sigmaV="$stress" \
        Materials/CPStressUpdate/temp="$TEMP" \
        Materials/CPStressUpdate/climbmodel=true \
        Executioner/end_time="$time1" \
        Executioner/TimeStepper/log_dt=0.010 \
        Outputs/file_base="out_Creep_${TEMP}K_${stress}MPa" > "log${stress}.run"

    echo "Initial run complete. Starting recovered simulation with end_time=${time2}"

    # Run the recovered simulation
    mpirun -n "$NPROCS" ~/projects/rhocp-vr-250401/rhocp-opt -i bcc_pxtal.i \
        sigmaV="$stress" \
        Materials/CPStressUpdate/temp="$TEMP" \
        Materials/CPStressUpdate/climbmodel=true \
        Executioner/end_time="$time2" \
        Executioner/TimeStepper/log_dt=0.001 \
        Outputs/file_base="out_Creep_${TEMP}K_${stress}MPa" --recover > "log${stress}r1.run"

    echo "Recovered run complete. Starting second recovered simulation with end_time=${time2}"
done



    # # Run the recovered simulation
    # mpirun -n 24 ~/projects/rhocp-vr-250420/rhocp-opt -i bcc_pxtal.i \
    #     sigmaV=200 \
    #     Materials/CPStressUpdate/temp="$TEMP" \
    #     Materials/CPStressUpdate/climbmodel=true \
    #     Executioner/end_time=3e6 \
    #     Executioner/TimeStepper/log_dt=0.001 \
    #     Outputs/file_base="out_Creep_823K_200MPa" --recover > "log200r1.run"