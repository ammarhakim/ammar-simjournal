#pragma once 
// Auto-generated code. Do not edit by hand 

#include <math.h>


static inline double euler_sq(double x) { return x*x; } 


struct euler_eqn { 
  double gas_gamma; 
}; 

static void euler_mulby_phi_prime(void *ctx, const double *q, const double *vin, double *vout); 
static void euler_mulby_phi_prime_inv(void *ctx, const double *q, const double *vin, double *vout); 
static void euler_flux_0(void *ctx, const double *q, double *fout); 
static void euler_flux_1(void *ctx, const double *q, double *fout); 
static void euler_flux_2(void *ctx, const double *q, double *fout); 
static void euler_flux_dot_vec(void *ctx, const double *vec, const double *q, double *fout); 
static void euler_projon_left_ev_0(void *ctx, const double *q, const double *vin, double *vout); 
static void euler_projon_left_ev_1(void *ctx, const double *q, const double *vin, double *vout); 
static void euler_projon_left_ev_2(void *ctx, const double *q, const double *vin, double *vout); 
static void euler_recwith_right_ev_0(void *ctx, const double *q, const double *vin, double *vout); 
static void euler_recwith_right_ev_1(void *ctx, const double *q, const double *vin, double *vout); 
static void euler_recwith_right_ev_2(void *ctx, const double *q, const double *vin, double *vout); 
static void euler_rescale_right_ev_0(void *ctx, const double *q, const double *w, double *revout); 
static void euler_rescale_right_ev_1(void *ctx, const double *q, const double *w, double *revout); 
static void euler_rescale_right_ev_2(void *ctx, const double *q, const double *w, double *revout); 
static double euler_ev_0(void *ctx, const double *q, double *ev); 
static double euler_ev_1(void *ctx, const double *q, double *ev); 
static double euler_ev_2(void *ctx, const double *q, double *ev); 
static void euler_ravg(void *ctx, const double *ql, const double *qr, double *qavg); 
static double euler_roe_fluct_0(void *ctx, const double *ql, const double *qr, double *apdq, double *amdq); 
static double euler_roe_fluct_1(void *ctx, const double *ql, const double *qr, double *apdq, double *amdq); 
static double euler_roe_fluct_2(void *ctx, const double *ql, const double *qr, double *apdq, double *amdq); 

static inline void 
euler_mulby_phi_prime(void *ctx, const double *q, const double *vin, double *vout) 
{ 
  struct euler_eqn *eqn = ctx; 
  double gas_gamma = eqn->gas_gamma; 
  double rho = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double uz = q[3]/q[0]; 
  double u2 = ux*ux + uy*uy + uz*uz; 
  vout[0] = vin[0]; 
  vout[1] = vin[1]*rho+vin[0]*ux; 
  vout[2] = vin[2]*rho+vin[0]*uy; 
  vout[3] = vin[3]*rho+vin[0]*uz; 
  vout[4] = vin[4]/((-1.0)+gas_gamma)+0.5*vin[0]*u2+vin[1]*rho*ux+vin[2]*rho*uy+vin[3]*rho*uz; 
} 

static inline void 
euler_mulby_phi_prime_inv(void *ctx, const double *q, const double *vin, double *vout) 
{ 
  struct euler_eqn *eqn = ctx; 
  double gas_gamma = eqn->gas_gamma; 
  double rho = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double uz = q[3]/q[0]; 
  double u2 = ux*ux + uy*uy + uz*uz; 
  vout[0] = vin[0]; 
  vout[1] = -(1.0*((-1.0*vin[1])+vin[0]*ux))/rho; 
  vout[2] = -(1.0*((-1.0*vin[2])+vin[0]*uy))/rho; 
  vout[3] = -(1.0*((-1.0*vin[3])+vin[0]*uz))/rho; 
  vout[4] = -0.5*((-1.0)+gas_gamma)*((-2.0*vin[4])-1.0*vin[0]*u2+2.0*vin[1]*ux+2.0*vin[2]*uy+2.0*vin[3]*uz); 
} 

static inline void 
euler_flux_0(void *ctx, const double *q, double *fout) 
{ 
  struct euler_eqn *eqn = ctx; 
  double gas_gamma = eqn->gas_gamma; 
  double rho = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double uz = q[3]/q[0]; 
  double u2 = ux*ux + uy*uy + uz*uz; 
  double E = q[4]; 
  double p = (gas_gamma-1)*(E-0.5*rho*u2); 
  fout[0] = rho*ux; 
  fout[1] = p+rho*euler_sq(ux); 
  fout[2] = rho*ux*uy; 
  fout[3] = rho*ux*uz; 
  fout[4] = p*ux+E*ux; 
} 

static inline void 
euler_flux_1(void *ctx, const double *q, double *fout) 
{ 
  struct euler_eqn *eqn = ctx; 
  double gas_gamma = eqn->gas_gamma; 
  double rho = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double uz = q[3]/q[0]; 
  double u2 = ux*ux + uy*uy + uz*uz; 
  double E = q[4]; 
  double p = (gas_gamma-1)*(E-0.5*rho*u2); 
  fout[0] = rho*uy; 
  fout[1] = rho*ux*uy; 
  fout[2] = p+rho*euler_sq(uy); 
  fout[3] = rho*uy*uz; 
  fout[4] = p*uy+E*uy; 
} 

static inline void 
euler_flux_2(void *ctx, const double *q, double *fout) 
{ 
  struct euler_eqn *eqn = ctx; 
  double gas_gamma = eqn->gas_gamma; 
  double rho = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double uz = q[3]/q[0]; 
  double u2 = ux*ux + uy*uy + uz*uz; 
  double E = q[4]; 
  double p = (gas_gamma-1)*(E-0.5*rho*u2); 
  fout[0] = rho*uz; 
  fout[1] = rho*ux*uz; 
  fout[2] = rho*uy*uz; 
  fout[3] = p+rho*euler_sq(uz); 
  fout[4] = p*uz+E*uz; 
} 

static inline void 
euler_flux_dot_vec(void *ctx, const double *vec, const double *q, double *fout) 
{ 
  struct euler_eqn *eqn = ctx; 
  double gas_gamma = eqn->gas_gamma; 
  double rho = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double uz = q[3]/q[0]; 
  double u2 = ux*ux + uy*uy + uz*uz; 
  double E = q[4]; 
  double p = (gas_gamma-1)*(E-0.5*rho*u2); 
  fout[0] = vec[0]*rho*ux+vec[1]*rho*uy+vec[2]*rho*uz; 
  fout[1] = vec[0]*p+rho*ux*(vec[0]*ux+vec[1]*uy+vec[2]*uz); 
  fout[2] = vec[1]*p+rho*uy*(vec[0]*ux+vec[1]*uy+vec[2]*uz); 
  fout[3] = vec[2]*p+rho*uz*(vec[0]*ux+vec[1]*uy+vec[2]*uz); 
  fout[4] = vec[0]*p*ux+vec[0]*E*ux+vec[1]*p*uy+vec[1]*E*uy+vec[2]*p*uz+vec[2]*E*uz; 
} 

static void 
euler_projon_left_ev_0(void *ctx, const double *q, const double *vin, double *vout) 
{ 
  struct euler_eqn *eqn = ctx; 
  double gas_gamma = eqn->gas_gamma; 
  double rho = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double uz = q[3]/q[0]; 
  double u2 = ux*ux + uy*uy + uz*uz; 
  double E = q[4]; 
  double p = (gas_gamma-1)*(E-0.5*rho*u2); 
  double cs = sqrt(gas_gamma*p/rho); 
  double v1[5]; 
  euler_mulby_phi_prime_inv(ctx, q, vin, v1); 
  vout[0] = (-(0.5*v1[1]*rho)/cs)+(0.5*v1[4])/euler_sq(cs); 
  vout[1] = (0.5*v1[1]*rho)/cs+(0.5*v1[4])/euler_sq(cs); 
  vout[2] = v1[0]-(1.0*v1[4])/euler_sq(cs); 
  vout[3] = v1[2]; 
  vout[4] = v1[3]; 
} 

static void 
euler_projon_left_ev_1(void *ctx, const double *q, const double *vin, double *vout) 
{ 
  struct euler_eqn *eqn = ctx; 
  double gas_gamma = eqn->gas_gamma; 
  double rho = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double uz = q[3]/q[0]; 
  double u2 = ux*ux + uy*uy + uz*uz; 
  double E = q[4]; 
  double p = (gas_gamma-1)*(E-0.5*rho*u2); 
  double cs = sqrt(gas_gamma*p/rho); 
  double v1[5]; 
  euler_mulby_phi_prime_inv(ctx, q, vin, v1); 
  vout[0] = (-(0.5*v1[2]*rho)/cs)+(0.5*v1[4])/euler_sq(cs); 
  vout[1] = (0.5*v1[2]*rho)/cs+(0.5*v1[4])/euler_sq(cs); 
  vout[2] = v1[0]-(1.0*v1[4])/euler_sq(cs); 
  vout[3] = v1[1]; 
  vout[4] = v1[3]; 
} 

static void 
euler_projon_left_ev_2(void *ctx, const double *q, const double *vin, double *vout) 
{ 
  struct euler_eqn *eqn = ctx; 
  double gas_gamma = eqn->gas_gamma; 
  double rho = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double uz = q[3]/q[0]; 
  double u2 = ux*ux + uy*uy + uz*uz; 
  double E = q[4]; 
  double p = (gas_gamma-1)*(E-0.5*rho*u2); 
  double cs = sqrt(gas_gamma*p/rho); 
  double v1[5]; 
  euler_mulby_phi_prime_inv(ctx, q, vin, v1); 
  vout[0] = (-(0.5*v1[3]*rho)/cs)+(0.5*v1[4])/euler_sq(cs); 
  vout[1] = (0.5*v1[3]*rho)/cs+(0.5*v1[4])/euler_sq(cs); 
  vout[2] = v1[0]-(1.0*v1[4])/euler_sq(cs); 
  vout[3] = v1[1]; 
  vout[4] = v1[2]; 
} 

static inline void 
euler_recwith_right_ev_0(void *ctx, const double *q, const double *vin, double *vout) 
{ 
  struct euler_eqn *eqn = ctx; 
  double gas_gamma = eqn->gas_gamma; 
  double rho = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double uz = q[3]/q[0]; 
  double u2 = ux*ux + uy*uy + uz*uz; 
  double E = q[4]; 
  double p = (gas_gamma-1)*(E-0.5*rho*u2); 
  double cs = sqrt(gas_gamma*p/rho); 
  double v1[5]; 
  v1[0] = vin[0]+vin[1]+vin[2]; 
  v1[1] = (-(1.0*vin[0]*cs)/rho)+(vin[1]*cs)/rho; 
  v1[2] = vin[3]; 
  v1[3] = vin[4]; 
  v1[4] = vin[0]*euler_sq(cs)+vin[1]*euler_sq(cs); 
  euler_mulby_phi_prime(ctx, q, v1, vout); 
} 

static inline void 
euler_recwith_right_ev_1(void *ctx, const double *q, const double *vin, double *vout) 
{ 
  struct euler_eqn *eqn = ctx; 
  double gas_gamma = eqn->gas_gamma; 
  double rho = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double uz = q[3]/q[0]; 
  double u2 = ux*ux + uy*uy + uz*uz; 
  double E = q[4]; 
  double p = (gas_gamma-1)*(E-0.5*rho*u2); 
  double cs = sqrt(gas_gamma*p/rho); 
  double v1[5]; 
  v1[0] = vin[0]+vin[1]+vin[2]; 
  v1[1] = vin[3]; 
  v1[2] = (-(1.0*vin[0]*cs)/rho)+(vin[1]*cs)/rho; 
  v1[3] = vin[4]; 
  v1[4] = vin[0]*euler_sq(cs)+vin[1]*euler_sq(cs); 
  euler_mulby_phi_prime(ctx, q, v1, vout); 
} 

static inline void 
euler_recwith_right_ev_2(void *ctx, const double *q, const double *vin, double *vout) 
{ 
  struct euler_eqn *eqn = ctx; 
  double gas_gamma = eqn->gas_gamma; 
  double rho = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double uz = q[3]/q[0]; 
  double u2 = ux*ux + uy*uy + uz*uz; 
  double E = q[4]; 
  double p = (gas_gamma-1)*(E-0.5*rho*u2); 
  double cs = sqrt(gas_gamma*p/rho); 
  double v1[5]; 
  v1[0] = vin[0]+vin[1]+vin[2]; 
  v1[1] = vin[3]; 
  v1[2] = vin[4]; 
  v1[3] = (-(1.0*vin[0]*cs)/rho)+(vin[1]*cs)/rho; 
  v1[4] = vin[0]*euler_sq(cs)+vin[1]*euler_sq(cs); 
  euler_mulby_phi_prime(ctx, q, v1, vout); 
} 

static void 
euler_rescale_right_ev_0(void *ctx, const double *q, const double *w, double *revout) 
{ 
  struct euler_eqn *eqn = ctx; 
  double gas_gamma = eqn->gas_gamma; 
  double rho = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double uz = q[3]/q[0]; 
  double u2 = ux*ux + uy*uy + uz*uz; 
  double E = q[4]; 
  double p = (gas_gamma-1)*(E-0.5*rho*u2); 
  double cs = sqrt(gas_gamma*p/rho); 
  double v1[5]; 
  double *vout; 
  v1[0] = w[0]; 
  v1[1] = -(1.0*w[0]*cs)/rho; 
  v1[2] = 0.0; 
  v1[3] = 0.0; 
  v1[4] = w[0]*euler_sq(cs); 
  vout = &revout[5*0]; 
  euler_mulby_phi_prime(ctx, q, v1, vout); 

  v1[0] = w[1]; 
  v1[1] = (w[1]*cs)/rho; 
  v1[2] = 0.0; 
  v1[3] = 0.0; 
  v1[4] = w[1]*euler_sq(cs); 
  vout = &revout[5*1]; 
  euler_mulby_phi_prime(ctx, q, v1, vout); 

  v1[0] = w[2]; 
  v1[1] = 0.0; 
  v1[2] = 0.0; 
  v1[3] = 0.0; 
  v1[4] = 0.0; 
  vout = &revout[5*2]; 
  euler_mulby_phi_prime(ctx, q, v1, vout); 

  v1[0] = 0.0; 
  v1[1] = 0.0; 
  v1[2] = w[3]; 
  v1[3] = 0.0; 
  v1[4] = 0.0; 
  vout = &revout[5*3]; 
  euler_mulby_phi_prime(ctx, q, v1, vout); 

  v1[0] = 0.0; 
  v1[1] = 0.0; 
  v1[2] = 0.0; 
  v1[3] = w[4]; 
  v1[4] = 0.0; 
  vout = &revout[5*4]; 
  euler_mulby_phi_prime(ctx, q, v1, vout); 

} 

static void 
euler_rescale_right_ev_1(void *ctx, const double *q, const double *w, double *revout) 
{ 
  struct euler_eqn *eqn = ctx; 
  double gas_gamma = eqn->gas_gamma; 
  double rho = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double uz = q[3]/q[0]; 
  double u2 = ux*ux + uy*uy + uz*uz; 
  double E = q[4]; 
  double p = (gas_gamma-1)*(E-0.5*rho*u2); 
  double cs = sqrt(gas_gamma*p/rho); 
  double v1[5]; 
  double *vout; 
  v1[0] = w[0]; 
  v1[1] = 0.0; 
  v1[2] = -(1.0*w[0]*cs)/rho; 
  v1[3] = 0.0; 
  v1[4] = w[0]*euler_sq(cs); 
  vout = &revout[5*0]; 
  euler_mulby_phi_prime(ctx, q, v1, vout); 

  v1[0] = w[1]; 
  v1[1] = 0.0; 
  v1[2] = (w[1]*cs)/rho; 
  v1[3] = 0.0; 
  v1[4] = w[1]*euler_sq(cs); 
  vout = &revout[5*1]; 
  euler_mulby_phi_prime(ctx, q, v1, vout); 

  v1[0] = w[2]; 
  v1[1] = 0.0; 
  v1[2] = 0.0; 
  v1[3] = 0.0; 
  v1[4] = 0.0; 
  vout = &revout[5*2]; 
  euler_mulby_phi_prime(ctx, q, v1, vout); 

  v1[0] = 0.0; 
  v1[1] = w[3]; 
  v1[2] = 0.0; 
  v1[3] = 0.0; 
  v1[4] = 0.0; 
  vout = &revout[5*3]; 
  euler_mulby_phi_prime(ctx, q, v1, vout); 

  v1[0] = 0.0; 
  v1[1] = 0.0; 
  v1[2] = 0.0; 
  v1[3] = w[4]; 
  v1[4] = 0.0; 
  vout = &revout[5*4]; 
  euler_mulby_phi_prime(ctx, q, v1, vout); 

} 

static void 
euler_rescale_right_ev_2(void *ctx, const double *q, const double *w, double *revout) 
{ 
  struct euler_eqn *eqn = ctx; 
  double gas_gamma = eqn->gas_gamma; 
  double rho = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double uz = q[3]/q[0]; 
  double u2 = ux*ux + uy*uy + uz*uz; 
  double E = q[4]; 
  double p = (gas_gamma-1)*(E-0.5*rho*u2); 
  double cs = sqrt(gas_gamma*p/rho); 
  double v1[5]; 
  double *vout; 
  v1[0] = w[0]; 
  v1[1] = 0.0; 
  v1[2] = 0.0; 
  v1[3] = -(1.0*w[0]*cs)/rho; 
  v1[4] = w[0]*euler_sq(cs); 
  vout = &revout[5*0]; 
  euler_mulby_phi_prime(ctx, q, v1, vout); 

  v1[0] = w[1]; 
  v1[1] = 0.0; 
  v1[2] = 0.0; 
  v1[3] = (w[1]*cs)/rho; 
  v1[4] = w[1]*euler_sq(cs); 
  vout = &revout[5*1]; 
  euler_mulby_phi_prime(ctx, q, v1, vout); 

  v1[0] = w[2]; 
  v1[1] = 0.0; 
  v1[2] = 0.0; 
  v1[3] = 0.0; 
  v1[4] = 0.0; 
  vout = &revout[5*2]; 
  euler_mulby_phi_prime(ctx, q, v1, vout); 

  v1[0] = 0.0; 
  v1[1] = w[3]; 
  v1[2] = 0.0; 
  v1[3] = 0.0; 
  v1[4] = 0.0; 
  vout = &revout[5*3]; 
  euler_mulby_phi_prime(ctx, q, v1, vout); 

  v1[0] = 0.0; 
  v1[1] = 0.0; 
  v1[2] = w[4]; 
  v1[3] = 0.0; 
  v1[4] = 0.0; 
  vout = &revout[5*4]; 
  euler_mulby_phi_prime(ctx, q, v1, vout); 

} 

static inline double 
euler_ev_0(void *ctx, const double *q, double *ev) 
{ 
  struct euler_eqn *eqn = ctx; 
  double gas_gamma = eqn->gas_gamma; 
  double rho = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double uz = q[3]/q[0]; 
  double u2 = ux*ux + uy*uy + uz*uz; 
  double E = q[4]; 
  double p = (gas_gamma-1)*(E-0.5*rho*u2); 
  double cs = sqrt(gas_gamma*p/rho); 
  ev[0] = (-cs)+ux; 
  ev[1] = cs+ux; 
  ev[2] = ux; 
  ev[3] = ux; 
  ev[4] = ux; 
  return fabs(q[1]/q[0]) + cs; 
} 

static inline double 
euler_ev_1(void *ctx, const double *q, double *ev) 
{ 
  struct euler_eqn *eqn = ctx; 
  double gas_gamma = eqn->gas_gamma; 
  double rho = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double uz = q[3]/q[0]; 
  double u2 = ux*ux + uy*uy + uz*uz; 
  double E = q[4]; 
  double p = (gas_gamma-1)*(E-0.5*rho*u2); 
  double cs = sqrt(gas_gamma*p/rho); 
  ev[0] = (-cs)+uy; 
  ev[1] = cs+uy; 
  ev[2] = uy; 
  ev[3] = uy; 
  ev[4] = uy; 
  return fabs(q[2]/q[0]) + cs; 
} 

static inline double 
euler_ev_2(void *ctx, const double *q, double *ev) 
{ 
  struct euler_eqn *eqn = ctx; 
  double gas_gamma = eqn->gas_gamma; 
  double rho = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double uz = q[3]/q[0]; 
  double u2 = ux*ux + uy*uy + uz*uz; 
  double E = q[4]; 
  double p = (gas_gamma-1)*(E-0.5*rho*u2); 
  double cs = sqrt(gas_gamma*p/rho); 
  ev[0] = (-cs)+uz; 
  ev[1] = cs+uz; 
  ev[2] = uz; 
  ev[3] = uz; 
  ev[4] = uz; 
  return fabs(q[3]/q[0]) + cs; 
} 

static inline void 
euler_ravg(void *ctx, const double *ql, const double *qr, double *qavg) 
{ 
  struct euler_eqn *eqn = ctx; 
  double gas_gamma = eqn->gas_gamma; 
  double rhol = ql[0]; 
  double uxl = ql[1]/ql[0]; 
  double uyl = ql[2]/ql[0]; 
  double uzl = ql[3]/ql[0]; 
  double u2l = uxl*uxl + uyl*uyl + uzl*uzl; 
  double El = ql[4]; 
  double pl = (gas_gamma-1)*(El-0.5*rhol*u2l); 
  double hl = (El+pl)/rhol; 
  double rhor = qr[0]; 
  double uxr = qr[1]/qr[0]; 
  double uyr = qr[2]/qr[0]; 
  double uzr = qr[3]/qr[0]; 
  double u2r = uxr*uxr + uyr*uyr + uzr*uzr; 
  double Er = qr[4]; 
  double pr = (gas_gamma-1)*(Er-0.5*rhor*u2r); 
  double hr = (Er+pr)/rhor; 
  double rho = 0.5*(sqrt(rhol)+sqrt(rhor)); 
  double ux = (uxl*sqrt(rhol)+uxr*sqrt(rhor))/(sqrt(rhol)+sqrt(rhor)); 
  double uy = (uyl*sqrt(rhol)+uyr*sqrt(rhor))/(sqrt(rhol)+sqrt(rhor)); 
  double uz = (uzl*sqrt(rhol)+uzr*sqrt(rhor))/(sqrt(rhol)+sqrt(rhor)); 
  double h = (hl*sqrt(rhol)+hr*sqrt(rhor))/(sqrt(rhol)+sqrt(rhor)); 
  double u2 = ux*ux + uy*uy + uz*uz; 
  double p = (gas_gamma-1)/gas_gamma*(rho*h-0.5*rho*u2); 
  qavg[0] = rho; 
  qavg[1] = rho*ux; 
  qavg[2] = rho*uy; 
  qavg[3] = rho*uz; 
  qavg[4] = 0.5*rho*u2 + p/(gas_gamma-1); 
} 

static inline double 
euler_roe_fluct_0(void *ctx, const double *ql, const double *qr, double *apdq, double *amdq) 
{ 
  struct euler_eqn *eqn = ctx; 
  double gas_gamma = eqn->gas_gamma; 
  double dq[5]; 

  for (int i=0; i<5; ++i) dq[i] = qr[i]-ql[i];

  double qavg[5];
  euler_ravg(eqn, ql, qr, qavg);

  double w[5];
  euler_projon_left_ev_0(eqn, qavg, dq, w);
  double wrev[5*5];
  euler_rescale_right_ev_0(eqn, qavg, w, wrev);
    
  double ev[5];
  euler_ev_0(eqn, qavg, ev);
  double amax = 0.0;

  for (int i=0; i<5; ++i) {
    amdq[i] = apdq[i] = 0.0;
      
    for (int j=0; j<5; ++j){
      apdq[i] += fmax(0.0, ev[j])*wrev[j*5+i];
      amdq[i] += fmin(0.0, ev[j])*wrev[j*5+i];
    }
    amax = fmax(amax, fabs(ev[i]));
  }
      return amax; 
} 

static inline double 
euler_roe_fluct_1(void *ctx, const double *ql, const double *qr, double *apdq, double *amdq) 
{ 
  struct euler_eqn *eqn = ctx; 
  double gas_gamma = eqn->gas_gamma; 
  double dq[5]; 

  for (int i=0; i<5; ++i) dq[i] = qr[i]-ql[i];

  double qavg[5];
  euler_ravg(eqn, ql, qr, qavg);

  double w[5];
  euler_projon_left_ev_1(eqn, qavg, dq, w);
  double wrev[5*5];
  euler_rescale_right_ev_1(eqn, qavg, w, wrev);
    
  double ev[5];
  euler_ev_1(eqn, qavg, ev);
  double amax = 0.0;

  for (int i=0; i<5; ++i) {
    amdq[i] = apdq[i] = 0.0;
      
    for (int j=0; j<5; ++j){
      apdq[i] += fmax(0.0, ev[j])*wrev[j*5+i];
      amdq[i] += fmin(0.0, ev[j])*wrev[j*5+i];
    }
    amax = fmax(amax, fabs(ev[i]));
  }
      return amax; 
} 

static inline double 
euler_roe_fluct_2(void *ctx, const double *ql, const double *qr, double *apdq, double *amdq) 
{ 
  struct euler_eqn *eqn = ctx; 
  double gas_gamma = eqn->gas_gamma; 
  double dq[5]; 

  for (int i=0; i<5; ++i) dq[i] = qr[i]-ql[i];

  double qavg[5];
  euler_ravg(eqn, ql, qr, qavg);

  double w[5];
  euler_projon_left_ev_2(eqn, qavg, dq, w);
  double wrev[5*5];
  euler_rescale_right_ev_2(eqn, qavg, w, wrev);
    
  double ev[5];
  euler_ev_2(eqn, qavg, ev);
  double amax = 0.0;

  for (int i=0; i<5; ++i) {
    amdq[i] = apdq[i] = 0.0;
      
    for (int j=0; j<5; ++j){
      apdq[i] += fmax(0.0, ev[j])*wrev[j*5+i];
      amdq[i] += fmin(0.0, ev[j])*wrev[j*5+i];
    }
    amax = fmax(amax, fabs(ev[i]));
  }
      return amax; 
} 

