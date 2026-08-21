#pragma once 
// Auto-generated code. Do not edit by hand 

#include <math.h>


static inline double burgers_sq(double x) { return x*x; } 


struct burgers_eqn { 
  
}; 

static void burgers_flux_0(void *ctx, const double *q, double *fout); 
static void burgers_flux_1(void *ctx, const double *q, double *fout); 
static void burgers_flux_2(void *ctx, const double *q, double *fout); 
static void burgers_flux_dot_vec(void *ctx, const double *vec, const double *q, double *fout); 
static void burgers_projon_left_ev_0(void *ctx, const double *q, const double *vin, double *vout); 
static void burgers_projon_left_ev_1(void *ctx, const double *q, const double *vin, double *vout); 
static void burgers_projon_left_ev_2(void *ctx, const double *q, const double *vin, double *vout); 
static void burgers_recwith_right_ev_0(void *ctx, const double *q, const double *vin, double *vout); 
static void burgers_recwith_right_ev_1(void *ctx, const double *q, const double *vin, double *vout); 
static void burgers_recwith_right_ev_2(void *ctx, const double *q, const double *vin, double *vout); 
static void burgers_rescale_right_ev_0(void *ctx, const double *q, const double *w, double *revout); 
static void burgers_rescale_right_ev_1(void *ctx, const double *q, const double *w, double *revout); 
static void burgers_rescale_right_ev_2(void *ctx, const double *q, const double *w, double *revout); 
static double burgers_ev_0(void *ctx, const double *q, double *ev); 
static double burgers_ev_1(void *ctx, const double *q, double *ev); 
static double burgers_ev_2(void *ctx, const double *q, double *ev); 
static void burgers_ravg(void *ctx, const double *ql, const double *qr, double *qavg); 
static double burgers_fluct_0(void *ctx, const double *ql, const double *qr, double *apdq, double *amdq); 
static double burgers_fluct_1(void *ctx, const double *ql, const double *qr, double *apdq, double *amdq); 
static double burgers_fluct_2(void *ctx, const double *ql, const double *qr, double *apdq, double *amdq); 

static inline void 
burgers_flux_0(void *ctx, const double *q, double *fout) 
{ 
  fout[0] = q[0]*q[0]/2; 
} 

static inline void 
burgers_flux_1(void *ctx, const double *q, double *fout) 
{ 
  fout[0] = q[0]*q[0]/2; 
} 

static inline void 
burgers_flux_2(void *ctx, const double *q, double *fout) 
{ 
  fout[0] = q[0]*q[0]/2; 
} 

static inline void 
burgers_flux_dot_vec(void *ctx, const double *vec, const double *q, double *fout) 
{ 
  fout[0] = q[0]*q[0]/2; 
} 

static void 
burgers_projon_left_ev_0(void *ctx, const double *q, const double *vin, double *vout) 
{ 
  vout[0] = vin[0]; 
} 

static void 
burgers_projon_left_ev_1(void *ctx, const double *q, const double *vin, double *vout) 
{ 
  vout[0] = vin[0]; 
} 

static void 
burgers_projon_left_ev_2(void *ctx, const double *q, const double *vin, double *vout) 
{ 
  vout[0] = vin[0]; 
} 

static inline void 
burgers_recwith_right_ev_0(void *ctx, const double *q, const double *vin, double *vout) 
{ 
  vout[0] = vin[0]; 
} 

static inline void 
burgers_recwith_right_ev_1(void *ctx, const double *q, const double *vin, double *vout) 
{ 
  vout[0] = vin[0]; 
} 

static inline void 
burgers_recwith_right_ev_2(void *ctx, const double *q, const double *vin, double *vout) 
{ 
  vout[0] = vin[0]; 
} 

static void 
burgers_rescale_right_ev_0(void *ctx, const double *q, const double *w, double *revout) 
{ 
  revout[0] = w[0]; 
} 

static void 
burgers_rescale_right_ev_1(void *ctx, const double *q, const double *w, double *revout) 
{ 
  revout[0] = w[0]; 
} 

static void 
burgers_rescale_right_ev_2(void *ctx, const double *q, const double *w, double *revout) 
{ 
  revout[0] = w[0]; 
} 

static inline double 
burgers_ev_0(void *ctx, const double *q, double *ev) 
{ 
  ev[0] = q[0]; 
  return fabs(ev[0]); 
} 

static inline double 
burgers_ev_1(void *ctx, const double *q, double *ev) 
{ 
  ev[0] = q[0]; 
  return fabs(ev[0]); 
} 

static inline double 
burgers_ev_2(void *ctx, const double *q, double *ev) 
{ 
  ev[0] = q[0]; 
  return fabs(ev[0]); 
} 

static inline void 
burgers_ravg(void *ctx, const double *ql, const double *qr, double *qavg) 
{ 
  qavg[0] = 0.5*(ql[0]+qr[0]); 
} 

static inline double 
burgers_fluct_0(void *ctx, const double *ql, const double *qr, double *apdq, double *amdq) 
{ 

      double c = 0.5*(ql[0]+qr[0]);
      apdq[0] = fmax(0.0, c)*(qr[0]-ql[0]);
      amdq[0] = fmin(0.0, c)*(qr[0]-ql[0]);      
      return fabs(c); 
} 

static inline double 
burgers_fluct_1(void *ctx, const double *ql, const double *qr, double *apdq, double *amdq) 
{ 

      double c = 0.5*(ql[0]+qr[0]);
      apdq[0] = fmax(0.0, c)*(qr[0]-ql[0]);
      amdq[0] = fmin(0.0, c)*(qr[0]-ql[0]);      
      return fabs(c); 
} 

static inline double 
burgers_fluct_2(void *ctx, const double *ql, const double *qr, double *apdq, double *amdq) 
{ 

      double c = 0.5*(ql[0]+qr[0]);
      apdq[0] = fmax(0.0, c)*(qr[0]-ql[0]);
      amdq[0] = fmin(0.0, c)*(qr[0]-ql[0]);      
      return fabs(c); 
} 

