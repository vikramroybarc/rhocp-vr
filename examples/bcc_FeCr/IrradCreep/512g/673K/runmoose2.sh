
cd 350MPa
mpirun -n 24 ~/projects/rhocp-vr/rhocp-opt -i bcc_pxtal.i sigmaV=350 Executioner/TimeStepper/log_dt=0.010 Executioner/end_time=1.0e3 Outputs/file_base=out_ICreep_673K_350MPa
mpirun -n 24 ~/projects/rhocp-vr/rhocp-opt -i bcc_pxtal.i sigmaV=350 Executioner/TimeStepper/log_dt=0.001 Executioner/end_time=1.0e5 Outputs/file_base=out_ICreep_673K_350MPa --recover
mpirun -n 24 ~/projects/rhocp-vr/rhocp-opt -i bcc_pxtal.i sigmaV=350 Executioner/TimeStepper/log_dt=0.001 Executioner/end_time=2.5e7 Outputs/file_base=out_ICreep_673K_350MPa --recover
