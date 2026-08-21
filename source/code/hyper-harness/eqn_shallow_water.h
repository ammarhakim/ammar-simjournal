#pragma once 
// Auto-generated code. Do not edit by hand 

#include <math.h>


static inline double shallow_water_sq(double x) { return x*x; } 


struct shallow_water_eqn { 
  double grav; 
}; 

static void shallow_water_mulby_phi_prime(void *ctx, const double *q, const double *vin, double *vout); 
static void shallow_water_mulby_phi_prime_inv(void *ctx, const double *q, const double *vin, double *vout); 
static void shallow_water_flux_0(void *ctx, const double *q, double *fout); 
static void shallow_water_flux_1(void *ctx, const double *q, double *fout); 
static void shallow_water_flux_dot_vec(void *ctx, const double *vec, const double *q, double *fout); 
static void shallow_water_projon_left_ev_0(void *ctx, const double *q, const double *vin, double *vout); 
static void shallow_water_projon_left_ev_1(void *ctx, const double *q, const double *vin, double *vout); 
static void shallow_water_recwith_right_ev_0(void *ctx, const double *q, const double *vin, double *vout); 
static void shallow_water_recwith_right_ev_1(void *ctx, const double *q, const double *vin, double *vout); 
static void shallow_water_rescale_right_ev_0(void *ctx, const double *q, const double *w, double *revout); 
static void shallow_water_rescale_right_ev_1(void *ctx, const double *q, const double *w, double *revout); 
static double shallow_water_ev_0(void *ctx, const double *q, double *ev); 
static double shallow_water_ev_1(void *ctx, const double *q, double *ev); 
static void shallow_water_ravg(void *ctx, const double *ql, const double *qr, double *qavg); 
static double shallow_water_roe_fluct_0(void *ctx, const double *ql, const double *qr, double *apdq, double *amdq); 
static double shallow_water_roe_fluct_1(void *ctx, const double *ql, const double *qr, double *apdq, double *amdq); 

static inline void 
shallow_water_mulby_phi_prime(void *ctx, const double *q, const double *vin, double *vout) 
{ 
  struct shallow_water_eqn *eqn = ctx; 
  double g = eqn->grav; 
  double h = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  double u2 = ux*ux + uy*uy; 
  vout[0] = vin[0]; 
  vout[1] = vin[1]*h+vin[0]*ux; 
  vout[2] = vin[2]*h+vin[0]*uy; 
} 

static inline void 
shallow_water_mulby_phi_prime_inv(void *ctx, const double *q, const double *vin, double *vout) 
{ 
  struct shallow_water_eqn *eqn = ctx; 
  double g = eqn->grav; 
  double h = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  
  double u2 = ux*ux + uy*uy; 
  vout[0] = vin[0]; 
  vout[1] = -(1.0*((-1.0*vin[1])+vin[0]*ux))/h; 
  vout[2] = -(1.0*((-1.0*vin[2])+vin[0]*uy))/h; 
} 

static inline void 
shallow_water_flux_0(void *ctx, const double *q, double *fout) 
{ 
  struct shallow_water_eqn *eqn = ctx; 
  double g = eqn->grav; 
  double h = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  
  double u2 = ux*ux + uy*uy; 
  fout[0] = h*ux; 
  fout[1] = 0.5*g*shallow_water_sq(h)+h*shallow_water_sq(ux); 
  fout[2] = h*ux*uy; 
} 

static inline void 
shallow_water_flux_1(void *ctx, const double *q, double *fout) 
{ 
  struct shallow_water_eqn *eqn = ctx; 
  double g = eqn->grav; 
  double h = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  
  double u2 = ux*ux + uy*uy; 
  fout[0] = h*uy; 
  fout[1] = h*ux*uy; 
  fout[2] = 0.5*g*shallow_water_sq(h)+h*shallow_water_sq(uy); 
} 

static void 
shallow_water_projon_left_ev_0(void *ctx, const double *q, const double *vin, double *vout) 
{ 
  struct shallow_water_eqn *eqn = ctx; 
  double g = eqn->grav; 
  double h = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  
  double u2 = ux*ux + uy*uy; 
  double v1[3]; 
  shallow_water_mulby_phi_prime_inv(ctx, q, vin, v1); 
  vout[0] = 0.5*v1[0]-(0.5*v1[1]*h)/sqrt(g*h); 
  vout[1] = 0.5*v1[0]+(0.5*v1[1]*h)/sqrt(g*h); 
  vout[2] = v1[2]; 
} 

static void 
shallow_water_projon_left_ev_1(void *ctx, const double *q, const double *vin, double *vout) 
{ 
  struct shallow_water_eqn *eqn = ctx; 
  double g = eqn->grav; 
  double h = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  
  double u2 = ux*ux + uy*uy; 
  double v1[3]; 
  shallow_water_mulby_phi_prime_inv(ctx, q, vin, v1); 
  vout[0] = 0.5*v1[0]-(0.5*v1[2]*h)/sqrt(g*h); 
  vout[1] = 0.5*v1[0]+(0.5*v1[2]*h)/sqrt(g*h); 
  vout[2] = v1[1]; 
} 

static inline void 
shallow_water_recwith_right_ev_0(void *ctx, const double *q, const double *vin, double *vout) 
{ 
  struct shallow_water_eqn *eqn = ctx; 
  double g = eqn->grav; 
  double h = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  
  double u2 = ux*ux + uy*uy; 
  double v1[3]; 
  v1[0] = vin[0]+vin[1]; 
  v1[1] = (-(1.0*vin[0]*sqrt(g*h))/h)+(vin[1]*sqrt(g*h))/h; 
  v1[2] = vin[2]; 
  shallow_water_mulby_phi_prime(ctx, q, v1, vout); 
} 

static inline void 
shallow_water_recwith_right_ev_1(void *ctx, const double *q, const double *vin, double *vout) 
{ 
  struct shallow_water_eqn *eqn = ctx; 
  double g = eqn->grav; 
  double h = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  
  double u2 = ux*ux + uy*uy; 
  double v1[3]; 
  v1[0] = vin[0]+vin[1]; 
  v1[1] = vin[2]; 
  v1[2] = (-(1.0*vin[0]*sqrt(g*h))/h)+(vin[1]*sqrt(g*h))/h; 
  shallow_water_mulby_phi_prime(ctx, q, v1, vout); 
} 

static void 
shallow_water_rescale_right_ev_0(void *ctx, const double *q, const double *w, double *revout) 
{ 
  struct shallow_water_eqn *eqn = ctx; 
  double g = eqn->grav; 
  double h = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  
  double u2 = ux*ux + uy*uy; 
  double v1[3]; 
  double *vout; 
  v1[0] = w[0]; 
  v1[1] = -(1.0*w[0]*sqrt(g*h))/h; 
  v1[2] = 0.0; 
  vout = &revout[3*0]; 
  shallow_water_mulby_phi_prime(ctx, q, v1, vout); 

  v1[0] = w[1]; 
  v1[1] = (w[1]*sqrt(g*h))/h; 
  v1[2] = 0.0; 
  vout = &revout[3*1]; 
  shallow_water_mulby_phi_prime(ctx, q, v1, vout); 

  v1[0] = 0.0; 
  v1[1] = 0.0; 
  v1[2] = w[2]; 
  vout = &revout[3*2]; 
  shallow_water_mulby_phi_prime(ctx, q, v1, vout); 

} 

static void 
shallow_water_rescale_right_ev_1(void *ctx, const double *q, const double *w, double *revout) 
{ 
  struct shallow_water_eqn *eqn = ctx; 
  double g = eqn->grav; 
  double h = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  
  double u2 = ux*ux + uy*uy; 
  double v1[3]; 
  double *vout; 
  v1[0] = w[0]; 
  v1[1] = 0.0; 
  v1[2] = -(1.0*w[0]*sqrt(g*h))/h; 
  vout = &revout[3*0]; 
  shallow_water_mulby_phi_prime(ctx, q, v1, vout); 

  v1[0] = w[1]; 
  v1[1] = 0.0; 
  v1[2] = (w[1]*sqrt(g*h))/h; 
  vout = &revout[3*1]; 
  shallow_water_mulby_phi_prime(ctx, q, v1, vout); 

  v1[0] = 0.0; 
  v1[1] = w[2]; 
  v1[2] = 0.0; 
  vout = &revout[3*2]; 
  shallow_water_mulby_phi_prime(ctx, q, v1, vout); 

} 

static inline double 
shallow_water_ev_0(void *ctx, const double *q, double *ev) 
{ 
  struct shallow_water_eqn *eqn = ctx; 
  double g = eqn->grav; 
  double h = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  
  double u2 = ux*ux + uy*uy; 
  ev[0] = (-sqrt(g*h))+ux; 
  ev[1] = sqrt(g*h)+ux; 
  ev[2] = ux; 
  return fabs(q[1]/q[0]) + sqrt(g*h); 
} 

static inline double 
shallow_water_ev_1(void *ctx, const double *q, double *ev) 
{ 
  struct shallow_water_eqn *eqn = ctx; 
  double g = eqn->grav; 
  double h = q[0]; 
  double ux = q[1]/q[0]; 
  double uy = q[2]/q[0]; 
  
  double u2 = ux*ux + uy*uy; 
  ev[0] = (-sqrt(g*h))+uy; 
  ev[1] = sqrt(g*h)+uy; 
  ev[2] = uy; 
  return fabs(q[2]/q[0]) + sqrt(g*h); 
} 

static inline void 
shallow_water_ravg(void *ctx, const double *ql, const double *qr, double *qavg) 
{ 
  struct shallow_water_eqn *eqn = ctx; 
  double g = eqn->grav; 
  double hl = ql[0]; 
  double uxl = ql[1]/ql[0]; 
  double uyl = ql[2]/ql[0]; 
  double u2l = uxl*uxl + uyl*uyl; 
  double hr = qr[0]; 
  double uxr = qr[1]/qr[0]; 
  double uyr = qr[2]/qr[0]; 
  double u2r = uxr*uxr + uyr*uyr; 
  double h = 0.5*(hl+hr); 
  double ux = (uxl*sqrt(hl)+uxr*sqrt(hr))/(sqrt(hl)+sqrt(hr)); 
  double uy = (uyl*sqrt(hl)+uyr*sqrt(hr))/(sqrt(hl)+sqrt(hr)); 
  qavg[0] = h; 
  qavg[1] = h*ux; 
  qavg[2] = h*uy; 
} 

static inline double 
shallow_water_roe_fluct_0(void *ctx, const double *ql, const double *qr, double *apdq, double *amdq) 
{ 
  struct shallow_water_eqn *eqn = ctx; 
  double dq[3]; 

  for (int i=0; i<3; ++i) dq[i] = qr[i]-ql[i];

  double qavg[3];
  shallow_water_ravg(eqn, ql, qr, qavg);

  double w[3];
  shallow_water_projon_left_ev_0(eqn, qavg, dq, w);
  double wrev[3*3];
  shallow_water_rescale_right_ev_0(eqn, qavg, w, wrev);
    
  double ev[3];
  shallow_water_ev_0(eqn, qavg, ev);
  double amax = 0.0;

  for (int i=0; i<3; ++i) {
    amdq[i] = apdq[i] = 0.0;
      
    for (int j=0; j<3; ++j){
      apdq[i] += fmax(0.0, ev[j])*wrev[j*3+i];
      amdq[i] += fmin(0.0, ev[j])*wrev[j*3+i];
    }
    amax = fmax(amax, fabs(ev[i]));
  }
      return amax; 
} 

static inline double 
shallow_water_roe_fluct_1(void *ctx, const double *ql, const double *qr, double *apdq, double *amdq) 
{ 
  struct shallow_water_eqn *eqn = ctx; 
  double dq[3]; 

  for (int i=0; i<3; ++i) dq[i] = qr[i]-ql[i];

  double qavg[3];
  shallow_water_ravg(eqn, ql, qr, qavg);

  double w[3];
  shallow_water_projon_left_ev_1(eqn, qavg, dq, w);
  double wrev[3*3];
  shallow_water_rescale_right_ev_1(eqn, qavg, w, wrev);
    
  double ev[3];
  shallow_water_ev_1(eqn, qavg, ev);
  double amax = 0.0;

  for (int i=0; i<3; ++i) {
    amdq[i] = apdq[i] = 0.0;
      
    for (int j=0; j<3; ++j){
      apdq[i] += fmax(0.0, ev[j])*wrev[j*3+i];
      amdq[i] += fmin(0.0, ev[j])*wrev[j*3+i];
    }
    amax = fmax(amax, fabs(ev[i]));
  }
      return amax; 
} 

