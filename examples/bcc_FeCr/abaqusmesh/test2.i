[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]


[Mesh]
  [fmg]
  # type = FileMeshGenerator
  # file = Cantilever.inp
  type = GeneratedMeshGenerator
  elem_type = HEX8
  dim = 3
  nx = 1
  ny = 1
  nz = 1
  xmin = 0.0
  xmax = 1.0
  ymin = 0.0
  ymax = 1.0
  zmin = 0.0
  zmax = 1.0  
  []
  [bot_corner]
    type = ExtraNodesetGenerator
    new_boundary = bot_corner
    input = fmg
    coord = '0.0 0.0 0.0'
  []
  [add_side_sets]
    type = SideSetsFromNormalsGenerator
    normals = '0 0 -1
               0 0  1'
    fixed_normal = false
    new_boundary = 'bc_face load_face'
    input = bot_corner
  []
[]


[Variables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
  [../]
[]



[Physics/SolidMechanics/QuasiStatic]
  [all]
    strain = FINITE
    add_variables = true
  [../]
[]

[Materials]
  [./fplastic]
    type = FiniteStrainPlasticMaterial
    yield_stress='0. 445. 0.05 610. 0.1 680. 0.38 810. 0.95 920. 2. 950.'
  [../]
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    C_ijkl = '2.827e5 1.21e5 1.21e5 2.827e5 1.21e5 2.827e5 0.808e5 0.808e5 0.808e5'
    fill_method = symmetric9
  [../]
  # [./strain]
  #   type = ComputeFiniteStrain
  #   block = 1
  #   displacements = 'disp_x disp_y disp_z'
  # [../]
[]


[Functions]
  [./topfunc]
    type = ParsedFunction
    expression = 't'
  [../]
  [top_pull]
    type = ParsedFunction
    expression =  '1*1e-1'
  []
  [./dts]
    type = PiecewiseLinear
    x = '0       1e-2'
    y = '1e-4    5e-1'
  [../] 
  [force]
    type = PiecewiseLinear
    x = '0 1 2'
    y = '0 27.78 55.56' 
  []
  [zpresspull]
    type = PiecewiseLinear
    x = '0 10'
    y = '0 1'
  []
[]

[BCs]
  [x_fix_corner_node]
    type = DirichletBC
    preset = true
    variable = disp_x
    boundary = bot_corner
    value = 0.0
  []
  [y_fix_corner_node]
    type = DirichletBC
    preset = true
    variable = disp_y
    boundary = bot_corner
    value = 0.0
  []  
  [z_fix_corner_node]
    type = DirichletBC
    preset = true
    variable = disp_z
    boundary = bot_corner
    value = 0.0
  []
  [y_bcface]
    type = DirichletBC
    preset = true
    variable = disp_y
    boundary = bc_face
    value = 0.0
  []
  [x_bcface]
    type = DirichletBC
    preset = true
    variable = disp_x
    boundary = bc_face
    value = 0.0
  []
  [z_bcface]
    type = DirichletBC
    preset = true
    variable = disp_z
    boundary = bc_face
    value = 0.0
  []    
  # Z Pull by Displacement
  # [z_pull]
  #   type = PresetVelocity
  #   variable = disp_z
  #   boundary = load_face
  #   function = top_pull
  # []  
  # [z_pull_press]
  #   type = Pressure
  #   variable = disp_z
  #   boundary = load_face
  #   factor = -480
  #   function = zpresspull
  # [] 
[]

[NodalKernels]
  [force]
    type = UserForcingFunctionNodalKernel
    variable = disp_y
    boundary = load_face
    function = force
  []
[]



[AuxVariables]
  [./stress_zz]
    order = FIRST
    family = MONOMIAL
  [../]
    [./stress_yy]
      order = FIRST
      family = MONOMIAL
    [../]

  [./peeq]
    order = FIRST
    family = MONOMIAL
  [../]
  [./pe11]
    order = FIRST
    family = MONOMIAL
  [../]
  [./pe22]
    order = FIRST
    family = MONOMIAL
  [../]
  [./pe33]
    order = FIRST
    family = MONOMIAL
  [../]
[]


[AuxKernels]
  [./stress_zz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    index_i = 2
    index_j = 2
  [../]
    [./stress_yy]
      type = RankTwoAux
      rank_two_tensor = stress
      variable = stress_yy
      index_i = 1
      index_j = 1
    [../]    
  [./pe11]
    type = RankTwoAux
    rank_two_tensor = plastic_strain
    variable = pe11
    index_i = 0
    index_j = 0
    execute_on = TIMESTEP_END
  [../]
    [./pe22]
    type = RankTwoAux
    rank_two_tensor = plastic_strain
    variable = pe22
    index_i = 1
    index_j = 1
    execute_on = TIMESTEP_END
  [../]
  [./pe33]
    type = RankTwoAux
    rank_two_tensor = plastic_strain
    variable = pe33
    index_i = 2
    index_j = 2
    execute_on = TIMESTEP_END
  [../]
  [./eqv_plastic_strain]
    type = MaterialRealAux
    property = eqv_plastic_strain
    variable = peeq
    execute_on = TIMESTEP_END
  [../]
[]

[Postprocessors]
  [stress_zz]
    type = ElementAverageValue
    variable = stress_zz
  []
  [stress_yy]
    type = ElementAverageValue
    variable = stress_yy
  []  
  [peeq]
    type = ElementAverageValue
    variable = peeq
  []
  [reaction_z]
    type = SidesetReaction
    direction = '0 0 1'
    stress_tensor = stress
    boundary = bc_face
  []
[]

[Preconditioning]
  [./SMP]
   type = SMP
   full=true
  [../]
[]





[Executioner]
  type = Transient

  dt=1e-4
  dtmax=1
  dtmin=1e-8
  end_time=2

  solve_type = 'NEWTON'

  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu superlu_dist'
  line_search = 'none'

  l_tol = 1e-8
  nl_abs_tol = 1e-7
  nl_rel_tol = 1e-6
  nl_max_its = 20
  nl_forced_its = 1
  l_max_its = 10
  
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




[Outputs]
  file_base = test2
  exodus = true
[]

