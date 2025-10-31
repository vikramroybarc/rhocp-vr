# !/bin/bash

# Define Simulation Parameters for Thermal Creep Simulation

# NPROCS=24
# STRESS=(140 150 175 200 220)
# TIME1=(1e5 5e4 5e4 5e3 5e2)     # End times for initial simulations
# TIME2=(2e6 1e6 5e5 5e4 5e3)     # End times for recovered simulations
# TEMP=873


    # stress=120
    # time2=5e7
    # mpirun -n "$NPROCS" ~/projects/rhocp-vr-250401/rhocp-opt -i bcc_pxtal.i \
    #     sigmaV="$stress" \
    #     Materials/CPStressUpdate/temp="$TEMP" \
    #     Materials/CPStressUpdate/climbmodel=true \
    #     Executioner/end_time="$time2" \
    #     Executioner/TimeStepper/log_dt=0.0005 \
    #     Outputs/file_base="out_Creep_${TEMP}K_${stress}MPa" --recover > "log${stress}r1.run"

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
#         Executioner/TimeStepper/log_dt=0.005 \
#         Outputs/file_base="out_Creep_${TEMP}K_${stress}MPa" --recover > "log${stress}r1.run"
# done



# NPROCS=24
# STRESS=(70 130 160)
# TIME1=(5e5 5e5 5e4)     # End times for initial simulations
# TIME2=(5e7 1e7 1e6)     # End times for recovered simulations
# TEMP=873

# Iterate over the array index to pair STRESS, TIME1, and TIME2 correctly
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
#         Executioner/TimeStepper/log_dt=0.005 \
#         Outputs/file_base="out_Creep_${TEMP}K_${stress}MPa" --recover > "log${stress}r1.run"
# done

NPROCS=24
TEMP=873
stress=150
time2=2e6
    mpirun -n "$NPROCS" ~/projects/rhocp-vr-250401/rhocp-opt -i bcc_pxtal.i \
        sigmaV="$stress" \
        Materials/CPStressUpdate/temp="$TEMP" \
        Materials/CPStressUpdate/climbmodel=true \
        Executioner/end_time="$time2" \
        Executioner/TimeStepper/log_dt=0.0005 \
        Outputs/file_base="out_Creep_${TEMP}K_${stress}MPa" --recover > "log${stress}r2.run"