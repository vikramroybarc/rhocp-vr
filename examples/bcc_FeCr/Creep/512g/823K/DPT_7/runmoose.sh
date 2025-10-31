#!/bin/bash

# Define Simulation Parameters for Thermal Creep Simulation

NPROCS=24
STRESS=(121 87 52)
TEMP=823

for stress in "${STRESS[@]}"
do
    # Run the simulation with the specified parameters
    mpirun -n "$NPROCS" ~/projects/rhocp-vr-250401/rhocp-opt -i bcc_pxtal.i \
        sigmaV="$stress" \
        Materials/CPStressUpdate/temp="$TEMP" \
        Materials/CPStressUpdate/climbmodel=true \
        Executioner/end_time=5e5 \
        Executioner/TimeStepper/log_dt=0.10 \
        Outputs/file_base="out_Creep_${TEMP}K_${stress}MPa" > "log${stress}.run"

    # Run the simulation with a longer end time and smaller log_dt
    mpirun -n "$NPROCS" ~/projects/rhocp-vr-250401/rhocp-opt -i bcc_pxtal.i \
        sigmaV="$stress" \
        Materials/CPStressUpdate/temp="$TEMP" \
        Materials/CPStressUpdate/climbmodel=true \
        Executioner/end_time=5e6 \
        Executioner/TimeStepper/log_dt=0.05 \
        Outputs/file_base="out_Creep_${TEMP}K_${stress}MPa" --recover > "log${stress}r1.run"

        # Run the simulation with a longer end time and smaller log_dt
    mpirun -n "$NPROCS" ~/projects/rhocp-vr-250401/rhocp-opt -i bcc_pxtal.i \
        sigmaV="$stress" \
        Materials/CPStressUpdate/temp="$TEMP" \
        Materials/CPStressUpdate/climbmodel=true \
        Executioner/end_time=5e7 \
        Executioner/TimeStepper/log_dt=0.005 \
        Outputs/file_base="out_Creep_${TEMP}K_${stress}MPa" --recover > "log${stress}r2.run"
done


cd /home/vikramroy/projects/rhocp-vr-250420/examples/bcc_FeCr/Creep/512g/823K/UTS_7/
./runmoose.sh