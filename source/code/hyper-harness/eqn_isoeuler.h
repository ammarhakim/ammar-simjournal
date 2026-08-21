#pragma once 
// Auto-generated code. Do not edit by hand 

#include <math.h>


static inline double isoeuler_sq(double x) { return x*x; } 


struct isoeuler_eqn { 
  double cs; 
}; 

static void isoeuler_mulby_phi_prime(void *ctx, const double *q, const double *vin, double *vout); 
static void isoeuler_mulby_phi_prime_inv(void *ctx, const double *q, const double *vin, double *vout); 
static void isoeuler_flux_0(void *ctx, const double *q, double *fout); 
static void isoeuler_flux_1(void *ctx, const double *q, double *fout); 
static void isoeuler_flux_2(void *ctx, const double *q, double *fout); 
static void isoeuler_flux_dot_vec(void *ctx, const double *vec, const double *q, double *fout); 
static void isoeuler_projon_left_ev_0(void *ctx, const double *q, const double *vin, double *vout); 
static void isoeuler_projon_left_ev_1(void *ctx, const double *q, const double *vin, double *vout); 
static void isoeuler_projon_left_ev_2(void *ctx, const double *q, const double *vin, double *vout); 
static void isoeuler_recwith_right_ev_0(void *ctx, const double *q, const double *vin, double *vout); 
static void isoeuler_recwith_right_ev_1(void *ctx, const double *q, const double *vin, double *vout); 
static void isoeuler_recwith_right_ev_2(void *ctx, const double *q, const double *vin, double *vout); 
static void isoeuler_rescale_right_ev_0(void *ctx, const double *q, const double *w, double *revout); 
static void isoeuler_rescale_right_ev_1(void *ctx, const double *q, const double *w, double *revout); 
static void isoeuler_rescale_right_ev_2(void *ctx, const double *q, const double *w, double *revout); 
static double isoeuler_ev_0(void *ctx, const double *q, double *ev); 
static double isoeuler_ev_1(void *ctx, const double *q, double *ev); 
static double isoeuler_ev_2(void *ctx, const double *q, double *ev); 
static void isoeuler_ravg(void *ctx, const double *ql, const double *qr, double *qavg); 
static double isoeuler_roe_fluct_0(void *ctx, const double *ql, const double *qr, double *apdq, double *amdq); 
static double isoeuler_roe_fluct_1(void *ctx, const double *ql, const double *qr, double *apdq, double *amdq); 
static double isoeuler_roe_fluct_2(void *ctx, const double *ql, const double *qr, double *apdq, double *amdq); 

static inline void 
isoeuler_mulby_phi_prime(void *ctx, const double *q, const double *vin, double *vout) 
{ 
  struct isoeuler_eqn *eqn = ctx; 
  double cs = eqn->cs; 
  double rho = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double uz = q[3]/q[0]; 
  double u2 = ux*ux + uy*uy + uz*uz; 
  vout[0] = vin[0]; 
  vout[1] = vin[1]*rho+vin[0]*ux; 
  vout[2] = vin[2]*rho+vin[0]*uy; 
  vout[3] = vin[3]*rho+vin[0]*uz; 
} 

static inline void 
isoeuler_mulby_phi_prime_inv(void *ctx, const double *q, const double *vin, double *vout) 
{ 
  struct isoeuler_eqn *eqn = ctx; 
  double cs = eqn->cs; 
  double rho = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double uz = q[3]/q[0]; 
  double u2 = ux*ux + uy*uy + uz*uz; 
  vout[0] = vin[0]; 
  vout[1] = -(1.0*((-1.0*vin[1])+vin[0]*ux))/rho; 
  vout[2] = -(1.0*((-1.0*vin[2])+vin[0]*uy))/rho; 
  vout[3] = -(1.0*((-1.0*vin[3])+vin[0]*uz))/rho; 
} 

static inline void 
isoeuler_flux_0(void *ctx, const double *q, double *fout) 
{ 
  struct isoeuler_eqn *eqn = ctx; 
  double cs = eqn->cs; 
  double rho = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double uz = q[3]/q[0]; 
  double u2 = ux*ux + uy*uy + uz*uz; 
  fout[0] = rho*ux; 
  fout[1] = isoeuler_sq(cs)*rho+rho*isoeuler_sq(ux); 
  fout[2] = rho*ux*uy; 
  fout[3] = rho*ux*uz; 
} 

static inline void 
isoeuler_flux_1(void *ctx, const double *q, double *fout) 
{ 
  struct isoeuler_eqn *eqn = ctx; 
  double cs = eqn->cs; 
  double rho = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double uz = q[3]/q[0]; 
  double u2 = ux*ux + uy*uy + uz*uz; 
  fout[0] = rho*uy; 
  fout[1] = rho*ux*uy; 
  fout[2] = isoeuler_sq(cs)*rho+rho*isoeuler_sq(uy); 
  fout[3] = rho*uy*uz; 
} 

static inline void 
isoeuler_flux_2(void *ctx, const double *q, double *fout) 
{ 
  struct isoeuler_eqn *eqn = ctx; 
  double cs = eqn->cs; 
  double rho = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double uz = q[3]/q[0]; 
  double u2 = ux*ux + uy*uy + uz*uz; 
  fout[0] = rho*uz; 
  fout[1] = rho*ux*uz; 
  fout[2] = rho*uy*uz; 
  fout[3] = isoeuler_sq(cs)*rho+rho*isoeuler_sq(uz); 
} 

static inline void 
isoeuler_flux_dot_vec(void *ctx, const double *vec, const double *q, double *fout) 
{ 
  struct isoeuler_eqn *eqn = ctx; 
  double cs = eqn->cs; 
  double rho = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double uz = q[3]/q[0]; 
  double u2 = ux*ux + uy*uy + uz*uz; 
  fout[0] = vec[0]*rho*ux+vec[1]*rho*uy+vec[2]*rho*uz; 
  fout[1] = vec[0]*isoeuler_sq(cs)*rho+rho*ux*(vec[0]*ux+vec[1]*uy+vec[2]*uz); 
  fout[2] = vec[1]*isoeuler_sq(cs)*rho+rho*uy*(vec[0]*ux+vec[1]*uy+vec[2]*uz); 
  fout[3] = vec[2]*isoeuler_sq(cs)*rho+rho*uz*(vec[0]*ux+vec[1]*uy+vec[2]*uz); 
} 

static void 
isoeuler_projon_left_ev_0(void *ctx, const double *q, const double *vin, double *vout) 
{ 
  struct isoeuler_eqn *eqn = ctx; 
  double cs = eqn->cs; 
  double rho = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double uz = q[3]/q[0]; 
  double u2 = ux*ux + uy*uy + uz*uz; 
   
  double v1[4]; 
  isoeuler_mulby_phi_prime_inv(ctx, q, vin, v1); 
  vout[0] = 0.5*v1[0]-(0.5*v1[1]*rho)/cs; 
  vout[1] = 0.5*v1[0]+(0.5*v1[1]*rho)/cs; 
  vout[2] = v1[2]; 
  vout[3] = v1[3]; 
} 

static void 
isoeuler_projon_left_ev_1(void *ctx, const double *q, const double *vin, double *vout) 
{ 
  struct isoeuler_eqn *eqn = ctx; 
  double cs = eqn->cs; 
  double rho = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double uz = q[3]/q[0]; 
  double u2 = ux*ux + uy*uy + uz*uz; 
   
  double v1[4]; 
  isoeuler_mulby_phi_prime_inv(ctx, q, vin, v1); 
  vout[0] = 0.5*v1[0]-(0.5*v1[2]*rho)/cs; 
  vout[1] = 0.5*v1[0]+(0.5*v1[2]*rho)/cs; 
  vout[2] = v1[1]; 
  vout[3] = v1[3]; 
} 

static void 
isoeuler_projon_left_ev_2(void *ctx, const double *q, const double *vin, double *vout) 
{ 
  struct isoeuler_eqn *eqn = ctx; 
  double cs = eqn->cs; 
  double rho = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double uz = q[3]/q[0]; 
  double u2 = ux*ux + uy*uy + uz*uz; 
   
  double v1[4]; 
  isoeuler_mulby_phi_prime_inv(ctx, q, vin, v1); 
  vout[0] = 0.5*v1[0]-(0.5*v1[3]*rho)/cs; 
  vout[1] = 0.5*v1[0]+(0.5*v1[3]*rho)/cs; 
  vout[2] = v1[1]; 
  vout[3] = v1[2]; 
} 

static inline void 
isoeuler_recwith_right_ev_0(void *ctx, const double *q, const double *vin, double *vout) 
{ 
  struct isoeuler_eqn *eqn = ctx; 
  double cs = eqn->cs; 
  double rho = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double uz = q[3]/q[0]; 
  double u2 = ux*ux + uy*uy + uz*uz; 
  double v1[4]; 
  v1[0] = vin[0]+vin[1]; 
  v1[1] = (-(1.0*vin[0]*cs)/rho)+(vin[1]*cs)/rho; 
  v1[2] = vin[2]; 
  v1[3] = vin[3]; 
  isoeuler_mulby_phi_prime(ctx, q, v1, vout); 
} 

static inline void 
isoeuler_recwith_right_ev_1(void *ctx, const double *q, const double *vin, double *vout) 
{ 
  struct isoeuler_eqn *eqn = ctx; 
  double cs = eqn->cs; 
  double rho = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double uz = q[3]/q[0]; 
  double u2 = ux*ux + uy*uy + uz*uz; 
  double v1[4]; 
  v1[0] = vin[0]+vin[1]; 
  v1[1] = vin[2]; 
  v1[2] = (-(1.0*vin[0]*cs)/rho)+(vin[1]*cs)/rho; 
  v1[3] = vin[3]; 
  isoeuler_mulby_phi_prime(ctx, q, v1, vout); 
} 

static inline void 
isoeuler_recwith_right_ev_2(void *ctx, const double *q, const double *vin, double *vout) 
{ 
  struct isoeuler_eqn *eqn = ctx; 
  double cs = eqn->cs; 
  double rho = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double uz = q[3]/q[0]; 
  double u2 = ux*ux + uy*uy + uz*uz; 
  double v1[4]; 
  v1[0] = vin[0]+vin[1]; 
  v1[1] = vin[2]; 
  v1[2] = vin[3]; 
  v1[3] = (-(1.0*vin[0]*cs)/rho)+(vin[1]*cs)/rho; 
  isoeuler_mulby_phi_prime(ctx, q, v1, vout); 
} 

static void 
isoeuler_rescale_right_ev_0(void *ctx, const double *q, const double *w, double *revout) 
{ 
  struct isoeuler_eqn *eqn = ctx; 
  double cs = eqn->cs; 
  double rho = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double uz = q[3]/q[0]; 
  double u2 = ux*ux + uy*uy + uz*uz; 
  double v1[4]; 
  double *vout; 
  v1[0] = w[0]; 
  v1[1] = -(1.0*w[0]*cs)/rho; 
  v1[2] = 0.0; 
  v1[3] = 0.0; 
  vout = &revout[4*0]; 
  isoeuler_mulby_phi_prime(ctx, q, v1, vout); 

  v1[0] = w[1]; 
  v1[1] = (w[1]*cs)/rho; 
  v1[2] = 0.0; 
  v1[3] = 0.0; 
  vout = &revout[4*1]; 
  isoeuler_mulby_phi_prime(ctx, q, v1, vout); 

  v1[0] = 0.0; 
  v1[1] = 0.0; 
  v1[2] = w[2]; 
  v1[3] = 0.0; 
  vout = &revout[4*2]; 
  isoeuler_mulby_phi_prime(ctx, q, v1, vout); 

  v1[0] = 0.0; 
  v1[1] = 0.0; 
  v1[2] = 0.0; 
  v1[3] = w[3]; 
  vout = &revout[4*3]; 
  isoeuler_mulby_phi_prime(ctx, q, v1, vout); 

} 

static void 
isoeuler_rescale_right_ev_1(void *ctx, const double *q, const double *w, double *revout) 
{ 
  struct isoeuler_eqn *eqn = ctx; 
  double cs = eqn->cs; 
  double rho = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double uz = q[3]/q[0]; 
  double u2 = ux*ux + uy*uy + uz*uz; 
  double v1[4]; 
  double *vout; 
  v1[0] = w[0]; 
  v1[1] = 0.0; 
  v1[2] = -(1.0*w[0]*cs)/rho; 
  v1[3] = 0.0; 
  vout = &revout[4*0]; 
  isoeuler_mulby_phi_prime(ctx, q, v1, vout); 

  v1[0] = w[1]; 
  v1[1] = 0.0; 
  v1[2] = (w[1]*cs)/rho; 
  v1[3] = 0.0; 
  vout = &revout[4*1]; 
  isoeuler_mulby_phi_prime(ctx, q, v1, vout); 

  v1[0] = 0.0; 
  v1[1] = w[2]; 
  v1[2] = 0.0; 
  v1[3] = 0.0; 
  vout = &revout[4*2]; 
  isoeuler_mulby_phi_prime(ctx, q, v1, vout); 

  v1[0] = 0.0; 
  v1[1] = 0.0; 
  v1[2] = 0.0; 
  v1[3] = w[3]; 
  vout = &revout[4*3]; 
  isoeuler_mulby_phi_prime(ctx, q, v1, vout); 

} 

static void 
isoeuler_rescale_right_ev_2(void *ctx, const double *q, const double *w, double *revout) 
{ 
  struct isoeuler_eqn *eqn = ctx; 
  double cs = eqn->cs; 
  double rho = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double uz = q[3]/q[0]; 
  double u2 = ux*ux + uy*uy + uz*uz; 
  double v1[4]; 
  double *vout; 
  v1[0] = w[0]; 
  v1[1] = 0.0; 
  v1[2] = 0.0; 
  v1[3] = -(1.0*w[0]*cs)/rho; 
  vout = &revout[4*0]; 
  isoeuler_mulby_phi_prime(ctx, q, v1, vout); 

  v1[0] = w[1]; 
  v1[1] = 0.0; 
  v1[2] = 0.0; 
  v1[3] = (w[1]*cs)/rho; 
  vout = &revout[4*1]; 
  isoeuler_mulby_phi_prime(ctx, q, v1, vout); 

  v1[0] = 0.0; 
  v1[1] = w[2]; 
  v1[2] = 0.0; 
  v1[3] = 0.0; 
  vout = &revout[4*2]; 
  isoeuler_mulby_phi_prime(ctx, q, v1, vout); 

  v1[0] = 0.0; 
  v1[1] = 0.0; 
  v1[2] = w[3]; 
  v1[3] = 0.0; 
  vout = &revout[4*3]; 
  isoeuler_mulby_phi_prime(ctx, q, v1, vout); 

} 

static inline double 
isoeuler_ev_0(void *ctx, const double *q, double *ev) 
{ 
  struct isoeuler_eqn *eqn = ctx; 
  double cs = eqn->cs; 
  double rho = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double uz = q[3]/q[0]; 
  double u2 = ux*ux + uy*uy + uz*uz; 
  ev[0] = (-cs)+ux; 
  ev[1] = cs+ux; 
  ev[2] = ux; 
  ev[3] = ux; 
  return fabs(q[1]/q[0]) + cs; 
} 

static inline double 
isoeuler_ev_1(void *ctx, const double *q, double *ev) 
{ 
  struct isoeuler_eqn *eqn = ctx; 
  double cs = eqn->cs; 
  double rho = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double uz = q[3]/q[0]; 
  double u2 = ux*ux + uy*uy + uz*uz; 
  ev[0] = (-cs)+uy; 
  ev[1] = cs+uy; 
  ev[2] = uy; 
  ev[3] = uy; 
  return fabs(q[2]/q[0]) + cs; 
} 

static inline double 
isoeuler_ev_2(void *ctx, const double *q, double *ev) 
{ 
  struct isoeuler_eqn *eqn = ctx; 
  double cs = eqn->cs; 
  double rho = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double uz = q[3]/q[0]; 
  double u2 = ux*ux + uy*uy + uz*uz; 
  ev[0] = (-cs)+uz; 
  ev[1] = cs+uz; 
  ev[2] = uz; 
  ev[3] = uz; 
  return fabs(q[3]/q[0]) + cs; 
} 

static inline void 
isoeuler_ravg(void *ctx, const double *ql, const double *qr, double *qavg) 
{ 
  struct isoeuler_eqn *eqn = ctx; 
  double cs = eqn->cs; 
  double rhol = ql[0]; 
  double uxl = ql[1]/ql[0]; 
  double uyl = ql[2]/ql[0]; 
  double uzl = ql[3]/ql[0]; 
  double u2l = uxl*uxl + uyl*uyl + uzl*uzl; 
  double rhor = qr[0]; 
  double uxr = qr[1]/qr[0]; 
  double uyr = qr[2]/qr[0]; 
  double uzr = qr[3]/qr[0]; 
  double u2r = uxr*uxr + uyr*uyr + uzr*uzr; 
  double rho = 0.5*(sqrt(rhol)+sqrt(rhor)); 
  double ux = (uxl*sqrt(rhol)+uxr*sqrt(rhor))/(sqrt(rhol)+sqrt(rhor)); 
  double uy = (uyl*sqrt(rhol)+uyr*sqrt(rhor))/(sqrt(rhol)+sqrt(rhor)); 
  double uz = (uzl*sqrt(rhol)+uzr*sqrt(rhor))/(sqrt(rhol)+sqrt(rhor)); 
  qavg[0] = rho; 
  qavg[1] = rho*ux; 
  qavg[2] = rho*uy; 
  qavg[3] = rho*uz; 
} 

static inline double 
isoeuler_roe_fluct_0(void *ctx, const double *ql, const double *qr, double *apdq, double *amdq) 
{ 
  struct isoeuler_eqn *eqn = ctx; 
  double dq[4]; 

  for (int i=0; i<4; ++i) dq[i] = qr[i]-ql[i];

  double qavg[4];
  isoeuler_ravg(eqn, ql, qr, qavg);

  double w[4];
  isoeuler_projon_left_ev_0(eqn, qavg, dq, w);
  double wrev[4*4];
  isoeuler_rescale_right_ev_0(eqn, qavg, w, wrev);
    
  double ev[4];
  isoeuler_ev_0(eqn, qavg, ev);
  double amax = 0.0;

  for (int i=0; i<4; ++i) {
    amdq[i] = apdq[i] = 0.0;
      
    for (int j=0; j<4; ++j){
      apdq[i] += fmax(0.0, ev[j])*wrev[j*4+i];
      amdq[i] += fmin(0.0, ev[j])*wrev[j*4+i];
    }
    amax = fmax(amax, fabs(ev[i]));
  }
      return amax; 
} 

static inline double 
isoeuler_roe_fluct_1(void *ctx, const double *ql, const double *qr, double *apdq, double *amdq) 
{ 
  struct isoeuler_eqn *eqn = ctx; 
  double dq[4]; 

  for (int i=0; i<4; ++i) dq[i] = qr[i]-ql[i];

  double qavg[4];
  isoeuler_ravg(eqn, ql, qr, qavg);

  double w[4];
  isoeuler_projon_left_ev_1(eqn, qavg, dq, w);
  double wrev[4*4];
  isoeuler_rescale_right_ev_1(eqn, qavg, w, wrev);
    
  double ev[4];
  isoeuler_ev_1(eqn, qavg, ev);
  double amax = 0.0;

  for (int i=0; i<4; ++i) {
    amdq[i] = apdq[i] = 0.0;
      
    for (int j=0; j<4; ++j){
      apdq[i] += fmax(0.0, ev[j])*wrev[j*4+i];
      amdq[i] += fmin(0.0, ev[j])*wrev[j*4+i];
    }
    amax = fmax(amax, fabs(ev[i]));
  }
      return amax; 
} 

static inline double 
isoeuler_roe_fluct_2(void *ctx, const double *ql, const double *qr, double *apdq, double *amdq) 
{ 
  struct isoeuler_eqn *eqn = ctx; 
  double dq[4]; 

  for (int i=0; i<4; ++i) dq[i] = qr[i]-ql[i];

  double qavg[4];
  isoeuler_ravg(eqn, ql, qr, qavg);

  double w[4];
  isoeuler_projon_left_ev_2(eqn, qavg, dq, w);
  double wrev[4*4];
  isoeuler_rescale_right_ev_2(eqn, qavg, w, wrev);
    
  double ev[4];
  isoeuler_ev_2(eqn, qavg, ev);
  double amax = 0.0;

  for (int i=0; i<4; ++i) {
    amdq[i] = apdq[i] = 0.0;
      
    for (int j=0; j<4; ++j){
      apdq[i] += fmax(0.0, ev[j])*wrev[j*4+i];
      amdq[i] += fmin(0.0, ev[j])*wrev[j*4+i];
    }
    amax = fmax(amax, fabs(ev[i]));
  }
      return amax; 
} 

