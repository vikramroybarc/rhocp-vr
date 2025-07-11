#!/bin/bash

# Define Simulation Parameters for Thermal Creep Simulation

NPROCS=20

TEMP=873

stress=150
time2=2e6
echo "Initial run complete. Starting recovered simulation with end_time=${time2}"
mpirun -n "$NPROCS" ~/projects/rhocp-vr-250505/rhocp-opt -i bcc_pxtal.i \
    sigmaV="$stress" \
    Materials/CPStressUpdate/temp="$TEMP" \
    Materials/CPStressUpdate/climbmodel=true \
    Executioner/end_time="$time2" \
    Executioner/TimeStepper/log_dt=0.001 \
    Outputs/file_base="out_Creep_${TEMP}K_${stress}MPa" --recover > "log${stress}r1.run"


STRESS=(175 200)
TIME1=(5e4 5e3 )     # End times for initial simulations
TIME2=(5e5 1e5 )     # End times for recovered simulations    

Iterate over the array index to pair STRESS, TIME1, and TIME2 correctly
for i in "${!STRESS[@]}"
do
    stress=${STRESS[$i]}
    time1=${TIME1[$i]}
    time2=${TIME2[$i]}

    echo "Running initial simulation for STRESS=${stress} MPa, TEMP=${TEMP} K, end_time=${time1}"

    # Run the initial simulation
    mpirun -n "$NPROCS" ~/projects/rhocp-vr-250505/rhocp-opt -i bcc_pxtal.i \
        sigmaV="$stress" \
        Materials/CPStressUpdate/temp="$TEMP" \
        Materials/CPStressUpdate/climbmodel=true \
        Executioner/end_time="$time1" \
        Executioner/TimeStepper/log_dt=0.10 \
        Outputs/file_base="out_Creep_${TEMP}K_${stress}MPa" > "log${stress}.run"

    echo "Initial run complete. Starting recovered simulation with end_time=${time2}"

    # Run the recovered simulation
    mpirun -n "$NPROCS" ~/projects/rhocp-vr-250505/rhocp-opt -i bcc_pxtal.i \
        sigmaV="$stress" \
        Materials/CPStressUpdate/temp="$TEMP" \
        Materials/CPStressUpdate/climbmodel=true \
        Executioner/end_time="$time2" \
        Executioner/TimeStepper/log_dt=0.01 \
        Outputs/file_base="out_Creep_${TEMP}K_${stress}MPa" --recover > "log${stress}r1.run"
done


