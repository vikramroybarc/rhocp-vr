[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  [./fmg]
    type = FileMeshGenerator
    file = 64grains_64elements.inp
    # type = GeneratedMeshGenerator
    # elem_type = HEX8
    # dim = 3
    # nx = 1
    # ny = 1
    # nz = 1
    # xmin = 0.0
    # xmax = 1.0
    # ymin = 0.0
    # ymax = 1.0
    # zmin = 0.0
    # zmax = 1.0    
  [../]
  [./bot_corner]
    type = ExtraNodesetGenerator
    new_boundary = bot_corner
    input = fmg
    coord = '0 0.0 0.0'
  [../]
  [./add_side_sets]
    type = SideSetsFromNormalsGenerator
    normals = '1  0  0
               0  1  0
               0  0  1
              -1  0  0
               0 -1  0
               0  0 -1'
    fixed_normal = false
    new_boundary = 'xp_face yp_face zp_face xn_face yn_face zn_face'
    input=bot_corner
  [../]  
[]

[Variables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
    scaling = 1e-4
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
    scaling = 1e-4
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
    scaling = 1e-4
  [../]
[]

[AuxVariables]
  [./strain_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strain_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vonmises]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./rho_m]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./rho_i]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./Ep_eff]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./eeq]
    order = FIRST
    family = MONOMIAL
  [../]  
  [./w]
    order = FIRST
    family = MONOMIAL
  [../]  
  [./tau_by_tauf]
    order = FIRST
    family = MONOMIAL
  [../]  
  [./taup]
    order = FIRST
    family = MONOMIAL
  [../]  
  [./k_r]
    order = FIRST
    family = MONOMIAL
  [../]  

[]

[Physics/SolidMechanics/QuasiStatic]
  [./all]
    strain = FINITE
    incremental = true
    use_finite_deform_jacobian = true
    volumetric_locking_correction = false
  [../]
[]

[AuxKernels]
  [./eeq]
    type = StateVariable
    variable = eeq
    sdv_id = 46
    execute_on = timestep_end
  [../]  
  [./strain_zz]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = strain_zz
    execute_on = timestep_end
    index_i = 2
    index_j = 2
  [../]
  [./strain_yy]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = strain_yy
    execute_on = timestep_end
    index_i = 1
    index_j = 1
  [../]
  [./strain_xx]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = strain_xx
    execute_on = timestep_end
    index_i = 0
    index_j = 0
  [../]
  [./stress_zz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    execute_on = timestep_end
    index_i = 2
    index_j = 2
  [../]
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    execute_on = timestep_end
    index_i = 1
    index_j = 1
  [../]
  [./stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    execute_on = timestep_end
    index_i = 0
    index_j = 0
  [../]
  [./vonmises]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = vonmises
    execute_on = timestep_end
    scalar_type = VonMisesStress
  [../]
  [./rho_m]
    type = StateVariable
    variable = rho_m
    sdv_id = 48
    execute_on = timestep_end
  [../]
  [./rho_i]
    type = StateVariable
    variable = rho_i
    sdv_id = 49
    execute_on = timestep_end
  [../]
  [./Ep_eff]
    type = StateVariable
    variable = Ep_eff
    sdv_id = 47
    execute_on = timestep_end
  [../]
  [./w]
    type = StateVariable
    variable = w
    sdv_id = 147
    execute_on = timestep_end
  [../]
  [./tau_by_tauf]
    type = StateVariable
    variable = tau_by_tauf
    sdv_id = 148
    execute_on = timestep_end
  [../]    
    [./taup]
    type = StateVariable
    variable = taup
    sdv_id = 149
    execute_on = timestep_end
  [../]   
  [./k_r]
    type = StateVariable
    variable = k_r
    sdv_id = 150
    execute_on = timestep_end
  [../]  
[]

[UserObjects]
  [./euler_angle]
    type = EulerAngleReader
    # file_name = orientations.in
    file_name = orientations_roltex1.in
    execute_on = 'initial'
  [../]
[]

[Functions]
  [./top_pull]
    type = ParsedFunction
    expression = '0.4*1' # 1 is the sample dimension, 10/s is the strain rate
  [../]

  [./dts]
    type = PiecewiseLinear
    x = '0      2e-5'
    y = '1e-6   5e-4'
  [../]
[]

[BCs]
  # roller BCs on lateral faces
  [./bottom_roller]
    type = DirichletBC
    preset = true
    variable = disp_y
    boundary = yn_face
    value = 0.0
  [../]
  [./left_roller]
    type = DirichletBC
    preset = true
    variable = disp_x
    boundary = xn_face
    value = 0.0
  [../]
  [./back_roller]
    type = DirichletBC
    preset = true
    variable = disp_z
    boundary = zn_face
    value = 0.0
  [../]

  # fixed BCs
  # corner node fixed in all DOFs
  [./bottom_nodes_x]
    type = DirichletBC
    preset = true
    variable = disp_x
    boundary = bot_corner
    value = 0.0
  [../]
  [./bottom_nodes_y]
    type = DirichletBC
    preset = true
    variable = disp_y
    boundary = bot_corner
    value = 0.0
  [../]
  [./bottom_nodes_z]
    type = DirichletBC
    preset = true
    variable = disp_z
    boundary = bot_corner
    value = 0.0
  [../]

  # tensile loading along x-direction
  [./y_pull_function]
    type = PresetVelocity
    variable = disp_x
    boundary = xp_face
    function = top_pull
  [../]
[]

[Materials]
  [./CPStressUpdate]
    type = DDCPStressUpdate2
    propsFile = bcc_props.in
    slipSysFile = bcc_slip_sys1.in
    num_slip_sys = 24
    num_state_vars = 150 # 50 + 3*num_slip_sys +4
    num_props = 37
    temp = 300 # K
    tol = 1e-3
    EulerAngFileReader = euler_angle
    #intMatFile = FeHintmat.in
    sliplaw = 2
    crsslaw = 2
    deltaH_eV = true
    
  [../]
  [./elasticity_tensor]
    type = ComputeCPElasticityTensor
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full=true
  [../]
[]

[Executioner]
  type = Transient

  solve_type = 'NEWTON'

  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu superlu_dist'
  # line_search = 'none'

  l_tol = 1e-8
  nl_abs_tol = 1e-7
  nl_rel_tol = 1e-6
  nl_max_its = 20
  nl_forced_its = 1
  l_max_its = 10

  start_time = 0.0
  end_time = 1e-1

  [./TimeStepper]
    type = FunctionDT
    function = dts
    min_dt = 1e-12
    cutback_factor_at_failure = 0.2
    growth_factor = 1.2
  [../]

  [./Predictor]
    type = SimplePredictor
    scale = 1
  [../]
[]

[Postprocessors]
  # [./strain_zz]
  #   type = ElementAverageValue
  #   variable = strain_zz
  # [../]
  # [./strain_yy]
  #   type = ElementAverageValue
  #   variable = strain_yy
  # [../]
  [./strain_xx]
    type = ElementAverageValue
    variable = strain_xx
  [../]
  # [./stress_zz]
  #   type = ElementAverageValue
  #   variable = stress_zz
  # [../]
  # [./stress_yy]
  #   type = ElementAverageValue
  #   variable = stress_yy
  # [../]
  [./stress_xx]
    type = ElementAverageValue
    variable = stress_xx
  [../]
  # [./vonmises]
  #   type = ElementAverageValue
  #   variable = vonmises
  # [../]
  [./rho_m]
    type = ElementAverageValue
    variable = rho_m
  [../]
  [./rho_i]
    type = ElementAverageValue
    variable = rho_i
  [../]
  [./Ep_eff]
    type = ElementAverageValue
    variable = Ep_eff
  [../]
  [./eeq]
    type = ElementAverageValue
    variable = eeq
  [../]  
  [./w]
    type = ElementAverageValue
    variable = w
  [../]  
  [./tau_by_tauf]
    type = ElementAverageValue
    variable = tau_by_tauf
  [../]   
  [./taup]
    type = ElementAverageValue
    variable = taup
  [../]
  [./k_r]
    type = ElementAverageValue
    variable = k_r
  [../]

[]

[Outputs]
  file_base = out_sr1
  csv = true
  print_linear_residuals = true
  perf_graph = true
  time_step_interval = 5
  [./exodus]
    type = Exodus
    time_step_interval = 10
  [../]
  [./cp]
    type = Checkpoint
    time_step_interval = 50
    num_files = 2
  [../]
[]

# !include'sdv.i'