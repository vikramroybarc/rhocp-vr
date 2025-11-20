#include "ComputeStressBase.h"
#include "EulerAngleReader.h"
#include "GrainAreaSize.h"

class DDJ2StressUpdate2 : public ComputeStressBase
{
public:
  static InputParameters validParams();

  DDJ2StressUpdate2(const InputParameters & parameters);
  virtual ~DDJ2StressUpdate2();

protected:
  FileName _propsFile;
  FileName _slipSysFile;

  unsigned int _num_props;
  unsigned int _num_slip_sys;  
  unsigned int _num_state_vars;

  int _grainid;

  const Real _tol;
  const VariableValue & _temp;


  const EulerAngleReader * _EulerAngFileReader;
  const GrainAreaSize * _GrainAreaSize;

  Real _grain_size;
  int _isEulerRadian;
  int _isEulerBunge;



  virtual void initQpStatefulProperties();
  virtual void computeQpStress();

  Real max_val(Real a,Real b);


  MaterialProperty<Point> & _euler_ang;

  // The eigenstrains
  std::vector<MaterialPropertyName> _eigenstrain_names;
  std::vector<const MaterialProperty<RankTwoTensor> *> _eigenstrains;
  std::vector<const MaterialProperty<RankTwoTensor> *> _eigenstrains_old;  



  MaterialProperty<std::vector<Real> > & _state_var;
  const MaterialProperty<std::vector<Real> > & _state_var_old;
  MaterialProperty<std::vector<Real> > & _properties;
  const MaterialProperty<std::vector<Real> > & _properties_old;

  // Miller indices of slip plane normals
  MaterialProperty<std::vector<std::vector<Real> > > & _y;
  const MaterialProperty<std::vector<std::vector<Real> > > & _y_old;
  // Miller indices of slip directions
  MaterialProperty<std::vector<std::vector<Real> > > & _z;
  const MaterialProperty<std::vector<std::vector<Real> > > & _z_old;  

  const MaterialProperty<RankTwoTensor> & _deformation_gradient;
  const MaterialProperty<RankTwoTensor> & _deformation_gradient_old;
  const MaterialProperty<RankTwoTensor> & _strain_increment;
  const MaterialProperty<RankTwoTensor> & _rotation_increment;
  const MaterialProperty<RankTwoTensor> & _stress_old;
  MaterialProperty<RankFourTensor> & _Cel_cp;

  // Name of the Phase Field Variable
  const VariableName _d_name;

  /// Decomposition types
  const enum class Decomposition {none, voldev} _decomposition;

  // @ {Strain Energy and its Derivatives with respect to damage
  const MaterialPropertyName _psie_name;
  ADMaterialProperty<Real> & _psie;
  ADMaterialProperty<Real> & _psie_active;
  ADMaterialProperty<Real> & _dpsie_dd;
  //}


  //@ {Plastic Strain Energy and its Derivatives with respect to damage
  const MaterialPropertyName _psip_name;
  ADMaterialProperty<Real> & _psip;
  ADMaterialProperty<Real> & _psip_active;
  ADMaterialProperty<Real> & _dpsip_dd;  
  // }


  // @ {degradation function of elastic energy and its derivatives with respect to damage
  const MaterialPropertyName _g_name;
  const ADMaterialProperty<Real> & _g;
  const ADMaterialProperty<Real> & _dg_dd;
  // }


  // @{ Plastic Heat Generation
  const ADMaterialProperty<Real> & _heat;
  //}  



  

  // parameters/variables used for calculations
  void readPropsFile();
  void assignProperties();
  void normalize_vector(Real*, Real*, Real*);
  void rotate_4th(Real a[3][3], Real b[3][3][3][3], Real (&c)[3][3][3][3]);
  void forth_to_Voigt(Real a[3][3][3][3], Real (&b)[6][6]);
  void Voigt_to_forth(Real b[6][6], Real (&a)[3][3][3][3]);
  void aaaa_dot_dot_bbbb(Real a[3][3][3][3], Real b[3][3][3][3], Real (&product)[3][3][3][3]);
  void aaaa_dot_dot_bb(Real a[3][3][3][3], Real b[3][3], Real (&product)[3][3]);
  void aaaa_dot_dot_bb2(Real a[3][3][3][3], RankTwoTensor bT, RankTwoTensor (&productT));
  void aa_dot_dot_bbbb(Real a[3][3], Real b[3][3][3][3], Real (&product)[3][3]);
  void aa_dot_bb(Real a[3][3], Real b[3][3], Real (&product)[3][3]);
  Real aa_dot_dot_bb(Real a[3][3], Real b[3][3]);
  Real power(Real x, Real y);
  Real sgn(Real x);

  void NR_residual_J2 (
    unsigned int num_slip_sys, 
    std::vector<std::vector<Real>> &xs0, 
    std::vector<std::vector<Real>> &xm0,     
    Real temp, 
    Real dt, 
    RankTwoTensor F1, 
    RankTwoTensor &F_p_inv, 
    RankTwoTensor F_p_inv_0, 
    Real C[3][3][3][3], 
    RankTwoTensor Np_star, 
    Real &rho_m0, 
    Real &rho_m, 
    Real &rho_i0, 
    Real &rho_i, 
    RankTwoTensor &bstress0, 
    RankTwoTensor &bstress, 
    RankTwoTensor &sig, 
    Real &sig_vm, 
    Real &sig_vm_star, 
    Real &s_a, 
    Real &s_t, 
    Real &gamma_dot_trial, 
    Real &residual);

  Real tolerance;

  Real act_vol;
  Real delF0;

  // @ Stress Decomposition Methods
  RankTwoTensor computeCauchyStress(
    const RankTwoTensor & F_el, 
    Real  C[3][3][3][3], 
    const bool converged);
  RankTwoTensor computeCauchyStressNoDecomposition(
    const RankTwoTensor & F_el, 
    Real  C[3][3][3][3], 
    const bool converged);
  RankTwoTensor computeCauchyStressVolDevDecomposition(
    const RankTwoTensor & F_el, 
    Real  C[3][3][3][3], 
    const bool converged);

  //

  // Variables in Props File
  Real 
  YM,                                                           // 1
  YM_perK,                                                      // 2
  poisson_ratio,                                                // 3  
  poisson_ratio_perK,                                           // 4
  G,                                                            // 5
  G_perK,                                                       // 6
  b_mag,                                                        // 7
  gammadot0g,                                                   // 8
  enthalpy_const,                                               // 9  
  p,                                                            // 10
  q,                                                            // 11
  tau0,                                                         // 12
  hp_coeff,                                                     // 13
  grain_size,                                                   // 14
  frictional_stress,                                            // 15  
  p0,                                                           // 16    
  k_rho,                                                        // 17  
  k_I,                                                          // 18
  rho_m_zero,                                                   // 19 
  rho_i_zero,                                                   // 20
  d_disl_zero,                                                  // 21
  k_M,                                                          // 22
  R_c,                                                          // 23
  k_ann,                                                        // 24
  k_D,                                                          // 25
  k_bs1,                                                        // 26
  k_bs2,                                                        // 27   
  TQfact;                                                       // 28 Taylor Quiney Factor
  
  Real B_k, freq;
};
