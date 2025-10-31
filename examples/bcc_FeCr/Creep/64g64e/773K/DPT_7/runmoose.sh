#!/bin/bash

# Define Simulation Parameters for Thermal Creep Simulation

NPROCS=24
STRESS=(52 87 121 173)
TEMP=763

for stress in "${STRESS[@]}"
do
    # Run the simulation with the specified parameters
    mpirun -n "$NPROCS" ~/projects/rhocp-vr-250420/rhocp-opt -i bcc_pxtal.i \
        sigmaV="$stress" \
        Materials/CPStressUpdate/temp="$TEMP" \
        Materials/CPStressUpdate/climbmodel=true \
        Executioner/end_time=5e5 \
        Executioner/TimeStepper/log_dt=0.10 \
        Outputs/file_base="out_Creep_773K_${stress}MPa" > "log${stress}.run"

    # Run the simulation with a longer end time and smaller log_dt
    mpirun -n "$NPROCS" ~/projects/rhocp-vr-250420/rhocp-opt -i bcc_pxtal.i \
        sigmaV="$stress" \
        Materials/CPStressUpdate/temp="$TEMP" \
        Materials/CPStressUpdate/climbmodel=true \
        Executioner/end_time=5e7 \
        Executioner/TimeStepper/log_dt=0.005 \
        Outputs/file_base="out_Creep_773K_${stress}MPa" --recover > "log${stress}r1.run"
done


cd /home/vikramr/projects/rhocp-vr-250420/examples/bcc_FeCr/Creep/64g64e/823K/DPT_7
./runmoose.sh

cd /home/vikramr/projects/rhocp-vr-250420/examples/bcc_FeCr/Creep/64g64e/823K/UTS_7
./runmoose.sh

cd /home/vikramr/projects/rhocp-vr-250420/examples/bcc_FeCr/Creep/64g64e/873K/PTS_7
./runmoose.sh

cd /home/vikramr/projects/rhocp-vr-250420/examples/bcc_FeCr/Creep/64g64e/873K/UTS_7
./runmoose.sh

cd /home/vikramr/projects/rhocp-vr-250420/examples/bcc_FeCr/Creep/64g64e/923K/UTS_7
./runmoose.sh

cd /home/vikramr/projects/rhocp-vr-250420/examples/bcc_FeCr/Creep/64g64e/973K/UTS_7
./runmoose.sh