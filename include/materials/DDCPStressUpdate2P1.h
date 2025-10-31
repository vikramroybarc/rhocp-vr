#include "ComputeStressBase.h"
#include "EulerAngleReader.h"
#include "EBSDMeshReader.h"
#include "GrainAreaSize.h"
//#include "IntMatFileReader.h"

class DDCPStressUpdate2P1 : public ComputeStressBase
{
public:
  static InputParameters validParams();

  DDCPStressUpdate2P1(const InputParameters & parameters);
  virtual ~DDCPStressUpdate2P1();

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
  const EBSDMeshReader * _EBSDFileReader;
  const GrainAreaSize * _GrainAreaSize;

  Real _grain_size;
  int _isEulerRadian;
  int _isEulerBunge;

  
  unsigned int _sliplaw;
  unsigned int _crsslaw;
  bool _deltaH_eV;

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

  // parameters/variables used for calculations
  static const int max_loops = 20;
  void readPropsFile();
  void assignProperties();
  void normalize_vector(Real*, Real*, Real*);
  void rotate_4th(Real a[3][3], Real b[3][3][3][3], Real (&c)[3][3][3][3]);
  void forth_to_Voigt(Real a[3][3][3][3], Real (&b)[6][6]);
  void Voigt_to_forth(Real b[6][6], Real (&a)[3][3][3][3]);
  void aaaa_dot_dot_bbbb(Real a[3][3][3][3], Real b[3][3][3][3], Real (&product)[3][3][3][3]);
  void aaaa_dot_dot_bb(Real a[3][3][3][3], Real b[3][3], Real (&product)[3][3]);
  void aa_dot_bb(Real a[3][3], Real b[3][3], Real (&product)[3][3]);
  Real aa_dot_dot_bb(Real a[3][3], Real b[3][3]);
  void bunge_angles(Real (&array1)[3][3], Real (&psi0)[3]);
  Real power(Real x, Real y);
  Real sgn(Real x);

  void NR_residual (
    unsigned int num_slip_sys, 
    std::vector<std::vector<Real>> &xs0, 
    std::vector<std::vector<Real>> &xm0, 
    Real temp, 
    Real dt, 
    std::vector<Real> gamma_dot, 
    std::vector<Real> gammadot_climb,
    RankTwoTensor F1, 
    RankTwoTensor &F_el, 
    RankTwoTensor &F_p_inv, 
    RankTwoTensor F_p_inv0, 
    Real C[3][3][3][3], 
    std::vector<Real> &rho_m0, 
    std::vector<Real> &rho_m, 
    std::vector<Real> &rho_i0, 
    std::vector<Real> &rho_i, 
    std::vector<Real> &bstress0, 
    std::vector<Real> &bstress, 
    RankTwoTensor &sig, 
    std::vector<Real> &tau, 
    std::vector<Real> &tau_eff, 
    std::vector<Real> &s_a, 
    std::vector<Real> &s_t, 
    std::vector<std::vector<Real>> A, 
    std::vector<std::vector<Real>> H, 
    Real w,
    Real tau_by_tauf,
    std::vector<Real> &tau_p,
    Real &max_taup,
    std::vector<Real> &residual, 
    Real &sse);

    // Function for Calculation of CRSS
    void calc_crss(
      unsigned int num_slip_sys,
      std::vector<Real> &rho_m,
      std::vector<Real> &rho_i,
      std::vector<Real> &s_a,
      std::vector<std::vector<Real>> A);

    // Function for Calculation of Derivative of CRSS
void calc_dsadgb(
  unsigned int num_slip_sys,
  std::vector<Real> &rho_m,
  std::vector<Real> &rho_i,
  std::vector<Real> &s_a,
  std::vector<std::vector<Real>> &A,
  std::vector<std::vector<Real>> &drhomdgb,
  std::vector<std::vector<Real>> &drhoidgb,
  std::vector<std::vector<Real>> &dsadgb);
    
    // Function for Calculation of Slip Rate
    void calc_sliprate(
      unsigned int num_slip_sys,
      std::vector<Real> &tau_eff,
      std::vector<Real> &s_a,
      std::vector<Real> &s_t,
      std::vector<Real> &gamma_dot_g,
      std::vector<Real> &gamma_dot);

    // Function for Calculation of Climb Rate
    void calc_climbrate(
      unsigned int num_slip_sys,
      Real Di,
      Real Dv,
      Real vacdensity_th,
      Real &interdensity,
      Real &vacdensity,
      std::vector<Real> &climbstress,
      std::vector<Real> &zi,
      std::vector<Real> &zv,
      std::vector<Real> &Kis,
      std::vector<Real> &Kvs,
      std::vector<Real> &vacdensity_0,
      std::vector<Real> &vel_climb,      
      std::vector<Real> &rho_m,
      std::vector<Real> &rho_i,
      std::vector<Real> &d_disl2,
      std::vector<Real> &gammadot_c);

  Real tolerance;

  Real act_vol;
  Real delF0;

  // Constants
  Real
  B_k,              // Boltzman Constant
  freq;             // Debyde Freq

  // Variables Inside Props File
  Real                  // Line no in Props File 
  C11,                  //1 Elastic Constant
  C11_perK,             //2 Temp Variation of C11
  C12,                  //3 Elastic Constant
  C12_perK,             //4 Temperature Variation of C12
  C44,                  //5 Elastic Constant
  C44_perK,             //6 Temperature Variation of C44
  G,                    //7 Shear Modulus
  G_perK,               //8 Shear Modulus Temp variation
  b_mag,                //9 Burgers Vector Magnitude
  gammadot0g,           //10 reference strain rate
  enthalpy_const,       //11 Multiplication Const for Enthalpy Term
  p,                    //12 Inner Exponent of Slip Law
  q,                    //13 Outer Exponent of Slip Law
  tau0,                 //14 Athermal Slip Resistance Constant
  hp_coeff,             //15 Hall Petch Coefficient
  grain_size,           //16 Grain Size
  frictional_stress,    //17 Lattice Frictional Resistance
  p0,                   //18 Dislocation Barrier Strength
  k_rho,               //19 Taylor Hardening Parameter - Immobile Dislocations
  k_rhom,               //20 Taylor Hardening Parameter - Mobile Dislocations
  k_I0,                  //21 Mean Free Path Constant
  Alatent,              //22 Latent Hardening Coefficient
  rho_m_zero,           //23 Initial Mobile Dislocation Density
  rho_i_zero,           //24 Initial Immobile Dislocation Density
  d_disl_zero,          //25 Initial Dislocation Line Length
  k_M,                  //26 Dislocation Line Generation Constant
  R_c,                  //27 Critical Capture Radius for Dislocation Annhilation
  k_ann,                //28 Dislocation Evlution annhilation constant
  k_r0,                 //29 Immobile Dislocation Evolution Release Constant
  k_bs1,                //30 Backstress Evolution Constant 1
  k_bs2,                //31 Backstress Evolution Constant 2
  phi,                  //32 Steady State Order parameter related to average link length
  k_r1,                //33 Inner Exponent of Dislocation Relase Term
  corr_factor;         //34 Correlation Factor for mm to m conversion




  std::vector<std::vector<Real>> sintmat; // Interaction Matrix

  int numrows(std::string);
  int numcols(std::string);
  void readfile(std::string);

  Real calc_w(
  Real phi,
  Real &tau_by_tauf,
  std::vector<Real> &tau_eff,
  std::vector<Real> &s_a,
  std::vector<Real> &rho_i,
  std::vector<Real> &rho_m);

  void calc_taup(
  unsigned int num_slip_sys,
  std::vector<Real> &s_t,
  std::vector<Real> &gamma_dot,
  std::vector<Real> &tau_p,
  Real &max_taup);

  


  

  Real sse;
};
