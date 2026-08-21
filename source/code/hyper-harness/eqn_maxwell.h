#pragma once 
// Auto-generated code. Do not edit by hand 

#include <math.h>


static inline double maxwell_sq(double x) { return x*x; } 


struct maxwell_eqn { 
  double c, efact, bfact; 
}; 

static void maxwell_flux_0(void *ctx, const double *q, double *fout); 
static void maxwell_flux_1(void *ctx, const double *q, double *fout); 
static void maxwell_flux_2(void *ctx, const double *q, double *fout); 
static void maxwell_flux_dot_vec(void *ctx, const double *vec, const double *q, double *fout); 
static void maxwell_projon_left_ev_0(void *ctx, const double *q, const double *vin, double *vout); 
static void maxwell_projon_left_ev_1(void *ctx, const double *q, const double *vin, double *vout); 
static void maxwell_projon_left_ev_2(void *ctx, const double *q, const double *vin, double *vout); 
static void maxwell_recwith_right_ev_0(void *ctx, const double *q, const double *vin, double *vout); 
static void maxwell_recwith_right_ev_1(void *ctx, const double *q, const double *vin, double *vout); 
static void maxwell_recwith_right_ev_2(void *ctx, const double *q, const double *vin, double *vout); 
static void maxwell_rescale_right_ev_0(void *ctx, const double *q, const double *w, double *revout); 
static void maxwell_rescale_right_ev_1(void *ctx, const double *q, const double *w, double *revout); 
static void maxwell_rescale_right_ev_2(void *ctx, const double *q, const double *w, double *revout); 
static double maxwell_ev_0(void *ctx, const double *q, double *ev); 
static double maxwell_ev_1(void *ctx, const double *q, double *ev); 
static double maxwell_ev_2(void *ctx, const double *q, double *ev); 
static void maxwell_ravg(void *ctx, const double *ql, const double *qr, double *qavg); 
static double maxwell_fluct_0(void *ctx, const double *ql, const double *qr, double *apdq, double *amdq); 
static double maxwell_fluct_1(void *ctx, const double *ql, const double *qr, double *apdq, double *amdq); 
static double maxwell_fluct_2(void *ctx, const double *ql, const double *qr, double *apdq, double *amdq); 

static inline void 
maxwell_flux_0(void *ctx, const double *q, double *fout) 
{ 
  struct maxwell_eqn *eqn = ctx; 
  double c = eqn->c; 
  double c2 = eqn->c*eqn->c; 
  double efact = eqn->efact; 
  double bfact = eqn->bfact; 
  fout[0] = q[6]*c2*efact; 
  fout[1] = q[5]*c2; 
  fout[2] = -1.0*q[4]*c2; 
  fout[3] = q[7]*bfact; 
  fout[4] = -1.0*q[2]; 
  fout[5] = q[1]; 
  fout[6] = q[0]*efact; 
  fout[7] = q[3]*bfact*c2; 
} 

static inline void 
maxwell_flux_1(void *ctx, const double *q, double *fout) 
{ 
  struct maxwell_eqn *eqn = ctx; 
  double c = eqn->c; 
  double c2 = eqn->c*eqn->c; 
  double efact = eqn->efact; 
  double bfact = eqn->bfact; 
  fout[0] = -1.0*q[5]*c2; 
  fout[1] = q[6]*c2*efact; 
  fout[2] = q[3]*c2; 
  fout[3] = q[2]; 
  fout[4] = q[7]*bfact; 
  fout[5] = -1.0*q[0]; 
  fout[6] = q[1]*efact; 
  fout[7] = q[4]*bfact*c2; 
} 

static inline void 
maxwell_flux_2(void *ctx, const double *q, double *fout) 
{ 
  struct maxwell_eqn *eqn = ctx; 
  double c = eqn->c; 
  double c2 = eqn->c*eqn->c; 
  double efact = eqn->efact; 
  double bfact = eqn->bfact; 
  fout[0] = q[4]*c2; 
  fout[1] = -1.0*q[3]*c2; 
  fout[2] = q[6]*c2*efact; 
  fout[3] = -1.0*q[1]; 
  fout[4] = q[0]; 
  fout[5] = q[7]*bfact; 
  fout[6] = q[2]*efact; 
  fout[7] = q[5]*bfact*c2; 
} 


static void 
maxwell_projon_left_ev_0(void *ctx, const double *q, const double *vin, double *vout) 
{ 
  struct maxwell_eqn *eqn = ctx; 
  double c = eqn->c; 
  double c2 = eqn->c*eqn->c; 
  double efact = eqn->efact; 
  double bfact = eqn->bfact; 
  vout[0] = 0.5*vin[3]-(0.5*vin[7])/c; 
  vout[1] = 0.5*vin[3]+(0.5*vin[7])/c; 
  vout[2] = 0.5*vin[0]-0.5*vin[6]*c; 
  vout[3] = 0.5*vin[0]+0.5*vin[6]*c; 
  vout[4] = 0.5*vin[1]-0.5*vin[5]*c; 
  vout[5] = 0.5*vin[2]+0.5*vin[4]*c; 
  vout[6] = 0.5*vin[1]+0.5*vin[5]*c; 
  vout[7] = 0.5*vin[2]-0.5*vin[4]*c; 
} 

static void 
maxwell_projon_left_ev_1(void *ctx, const double *q, const double *vin, double *vout) 
{ 
  struct maxwell_eqn *eqn = ctx; 
  double c = eqn->c; 
  double c2 = eqn->c*eqn->c; 
  double efact = eqn->efact; 
  double bfact = eqn->bfact; 
  vout[0] = 0.5*vin[4]-(0.5*vin[7])/c; 
  vout[1] = 0.5*vin[4]+(0.5*vin[7])/c; 
  vout[2] = 0.5*vin[1]-0.5*vin[6]*c; 
  vout[3] = 0.5*vin[1]+0.5*vin[6]*c; 
  vout[4] = 0.5*vin[0]+0.5*vin[5]*c; 
  vout[5] = 0.5*vin[2]-0.5*vin[3]*c; 
  vout[6] = 0.5*vin[0]-0.5*vin[5]*c; 
  vout[7] = 0.5*vin[2]+0.5*vin[3]*c; 
} 

static void 
maxwell_projon_left_ev_2(void *ctx, const double *q, const double *vin, double *vout) 
{ 
  struct maxwell_eqn *eqn = ctx; 
  double c = eqn->c; 
  double c2 = eqn->c*eqn->c; 
  double efact = eqn->efact; 
  double bfact = eqn->bfact; 
  vout[0] = 0.5*vin[5]-(0.5*vin[7])/c; 
  vout[1] = 0.5*vin[5]+(0.5*vin[7])/c; 
  vout[2] = 0.5*vin[2]-0.5*vin[6]*c; 
  vout[3] = 0.5*vin[2]+0.5*vin[6]*c; 
  vout[4] = 0.5*vin[0]-0.5*vin[4]*c; 
  vout[5] = 0.5*vin[1]+0.5*vin[3]*c; 
  vout[6] = 0.5*vin[0]+0.5*vin[4]*c; 
  vout[7] = 0.5*vin[1]-0.5*vin[3]*c; 
} 

static inline void 
maxwell_recwith_right_ev_0(void *ctx, const double *q, const double *vin, double *vout) 
{ 
  struct maxwell_eqn *eqn = ctx; 
  double c = eqn->c; 
  double c2 = eqn->c*eqn->c; 
  double efact = eqn->efact; 
  double bfact = eqn->bfact; 
  vout[0] = vin[2]+vin[3]; 
  vout[1] = vin[4]+vin[6]; 
  vout[2] = vin[5]+vin[7]; 
  vout[3] = vin[0]+vin[1]; 
  vout[4] = vin[5]/c-(1.0*vin[7])/c; 
  vout[5] = (-(1.0*vin[4])/c)+vin[6]/c; 
  vout[6] = (-(1.0*vin[2])/c)+vin[3]/c; 
  vout[7] = (-1.0*vin[0]*c)+vin[1]*c; 
} 

static inline void 
maxwell_recwith_right_ev_1(void *ctx, const double *q, const double *vin, double *vout) 
{ 
  struct maxwell_eqn *eqn = ctx; 
  double c = eqn->c; 
  double c2 = eqn->c*eqn->c; 
  double efact = eqn->efact; 
  double bfact = eqn->bfact; 
  vout[0] = vin[4]+vin[6]; 
  vout[1] = vin[2]+vin[3]; 
  vout[2] = vin[5]+vin[7]; 
  vout[3] = (-(1.0*vin[5])/c)+vin[7]/c; 
  vout[4] = vin[0]+vin[1]; 
  vout[5] = vin[4]/c-(1.0*vin[6])/c; 
  vout[6] = (-(1.0*vin[2])/c)+vin[3]/c; 
  vout[7] = (-1.0*vin[0]*c)+vin[1]*c; 
} 

static inline void 
maxwell_recwith_right_ev_2(void *ctx, const double *q, const double *vin, double *vout) 
{ 
  struct maxwell_eqn *eqn = ctx; 
  double c = eqn->c; 
  double c2 = eqn->c*eqn->c; 
  double efact = eqn->efact; 
  double bfact = eqn->bfact; 
  vout[0] = vin[4]+vin[6]; 
  vout[1] = vin[5]+vin[7]; 
  vout[2] = vin[2]+vin[3]; 
  vout[3] = vin[5]/c-(1.0*vin[7])/c; 
  vout[4] = (-(1.0*vin[4])/c)+vin[6]/c; 
  vout[5] = vin[0]+vin[1]; 
  vout[6] = (-(1.0*vin[2])/c)+vin[3]/c; 
  vout[7] = (-1.0*vin[0]*c)+vin[1]*c; 
} 

static void 
maxwell_rescale_right_ev_0(void *ctx, const double *q, const double *w, double *revout) 
{ 
  struct maxwell_eqn *eqn = ctx; 
  double c = eqn->c; 
  double c2 = eqn->c*eqn->c; 
  double efact = eqn->efact; 
  double bfact = eqn->bfact; 
  double v1[8]; 
  double *vout; 
  vout = &revout[8*0]; 
  vout[0] = 0.0; 
  vout[1] = 0.0; 
  vout[2] = 0.0; 
  vout[3] = w[0]; 
  vout[4] = 0.0; 
  vout[5] = 0.0; 
  vout[6] = 0.0; 
  vout[7] = -1.0*w[0]*c; 

  vout = &revout[8*1]; 
  vout[0] = 0.0; 
  vout[1] = 0.0; 
  vout[2] = 0.0; 
  vout[3] = w[1]; 
  vout[4] = 0.0; 
  vout[5] = 0.0; 
  vout[6] = 0.0; 
  vout[7] = w[1]*c; 

  vout = &revout[8*2]; 
  vout[0] = w[2]; 
  vout[1] = 0.0; 
  vout[2] = 0.0; 
  vout[3] = 0.0; 
  vout[4] = 0.0; 
  vout[5] = 0.0; 
  vout[6] = -(1.0*w[2])/c; 
  vout[7] = 0.0; 

  vout = &revout[8*3]; 
  vout[0] = w[3]; 
  vout[1] = 0.0; 
  vout[2] = 0.0; 
  vout[3] = 0.0; 
  vout[4] = 0.0; 
  vout[5] = 0.0; 
  vout[6] = w[3]/c; 
  vout[7] = 0.0; 

  vout = &revout[8*4]; 
  vout[0] = 0.0; 
  vout[1] = w[4]; 
  vout[2] = 0.0; 
  vout[3] = 0.0; 
  vout[4] = 0.0; 
  vout[5] = -(1.0*w[4])/c; 
  vout[6] = 0.0; 
  vout[7] = 0.0; 

  vout = &revout[8*5]; 
  vout[0] = 0.0; 
  vout[1] = 0.0; 
  vout[2] = w[5]; 
  vout[3] = 0.0; 
  vout[4] = w[5]/c; 
  vout[5] = 0.0; 
  vout[6] = 0.0; 
  vout[7] = 0.0; 

  vout = &revout[8*6]; 
  vout[0] = 0.0; 
  vout[1] = w[6]; 
  vout[2] = 0.0; 
  vout[3] = 0.0; 
  vout[4] = 0.0; 
  vout[5] = w[6]/c; 
  vout[6] = 0.0; 
  vout[7] = 0.0; 

  vout = &revout[8*7]; 
  vout[0] = 0.0; 
  vout[1] = 0.0; 
  vout[2] = w[7]; 
  vout[3] = 0.0; 
  vout[4] = -(1.0*w[7])/c; 
  vout[5] = 0.0; 
  vout[6] = 0.0; 
  vout[7] = 0.0; 

} 

static void 
maxwell_rescale_right_ev_1(void *ctx, const double *q, const double *w, double *revout) 
{ 
  struct maxwell_eqn *eqn = ctx; 
  double c = eqn->c; 
  double c2 = eqn->c*eqn->c; 
  double efact = eqn->efact; 
  double bfact = eqn->bfact; 
  double v1[8]; 
  double *vout; 
  vout = &revout[8*0]; 
  vout[0] = 0.0; 
  vout[1] = 0.0; 
  vout[2] = 0.0; 
  vout[3] = 0.0; 
  vout[4] = w[0]; 
  vout[5] = 0.0; 
  vout[6] = 0.0; 
  vout[7] = -1.0*w[0]*c; 

  vout = &revout[8*1]; 
  vout[0] = 0.0; 
  vout[1] = 0.0; 
  vout[2] = 0.0; 
  vout[3] = 0.0; 
  vout[4] = w[1]; 
  vout[5] = 0.0; 
  vout[6] = 0.0; 
  vout[7] = w[1]*c; 

  vout = &revout[8*2]; 
  vout[0] = 0.0; 
  vout[1] = w[2]; 
  vout[2] = 0.0; 
  vout[3] = 0.0; 
  vout[4] = 0.0; 
  vout[5] = 0.0; 
  vout[6] = -(1.0*w[2])/c; 
  vout[7] = 0.0; 

  vout = &revout[8*3]; 
  vout[0] = 0.0; 
  vout[1] = w[3]; 
  vout[2] = 0.0; 
  vout[3] = 0.0; 
  vout[4] = 0.0; 
  vout[5] = 0.0; 
  vout[6] = w[3]/c; 
  vout[7] = 0.0; 

  vout = &revout[8*4]; 
  vout[0] = w[4]; 
  vout[1] = 0.0; 
  vout[2] = 0.0; 
  vout[3] = 0.0; 
  vout[4] = 0.0; 
  vout[5] = w[4]/c; 
  vout[6] = 0.0; 
  vout[7] = 0.0; 

  vout = &revout[8*5]; 
  vout[0] = 0.0; 
  vout[1] = 0.0; 
  vout[2] = w[5]; 
  vout[3] = -(1.0*w[5])/c; 
  vout[4] = 0.0; 
  vout[5] = 0.0; 
  vout[6] = 0.0; 
  vout[7] = 0.0; 

  vout = &revout[8*6]; 
  vout[0] = w[6]; 
  vout[1] = 0.0; 
  vout[2] = 0.0; 
  vout[3] = 0.0; 
  vout[4] = 0.0; 
  vout[5] = -(1.0*w[6])/c; 
  vout[6] = 0.0; 
  vout[7] = 0.0; 

  vout = &revout[8*7]; 
  vout[0] = 0.0; 
  vout[1] = 0.0; 
  vout[2] = w[7]; 
  vout[3] = w[7]/c; 
  vout[4] = 0.0; 
  vout[5] = 0.0; 
  vout[6] = 0.0; 
  vout[7] = 0.0; 

} 

static void 
maxwell_rescale_right_ev_2(void *ctx, const double *q, const double *w, double *revout) 
{ 
  struct maxwell_eqn *eqn = ctx; 
  double c = eqn->c; 
  double c2 = eqn->c*eqn->c; 
  double efact = eqn->efact; 
  double bfact = eqn->bfact; 
  double v1[8]; 
  double *vout; 
  vout = &revout[8*0]; 
  vout[0] = 0.0; 
  vout[1] = 0.0; 
  vout[2] = 0.0; 
  vout[3] = 0.0; 
  vout[4] = 0.0; 
  vout[5] = w[0]; 
  vout[6] = 0.0; 
  vout[7] = -1.0*w[0]*c; 

  vout = &revout[8*1]; 
  vout[0] = 0.0; 
  vout[1] = 0.0; 
  vout[2] = 0.0; 
  vout[3] = 0.0; 
  vout[4] = 0.0; 
  vout[5] = w[1]; 
  vout[6] = 0.0; 
  vout[7] = w[1]*c; 

  vout = &revout[8*2]; 
  vout[0] = 0.0; 
  vout[1] = 0.0; 
  vout[2] = w[2]; 
  vout[3] = 0.0; 
  vout[4] = 0.0; 
  vout[5] = 0.0; 
  vout[6] = -(1.0*w[2])/c; 
  vout[7] = 0.0; 

  vout = &revout[8*3]; 
  vout[0] = 0.0; 
  vout[1] = 0.0; 
  vout[2] = w[3]; 
  vout[3] = 0.0; 
  vout[4] = 0.0; 
  vout[5] = 0.0; 
  vout[6] = w[3]/c; 
  vout[7] = 0.0; 

  vout = &revout[8*4]; 
  vout[0] = w[4]; 
  vout[1] = 0.0; 
  vout[2] = 0.0; 
  vout[3] = 0.0; 
  vout[4] = -(1.0*w[4])/c; 
  vout[5] = 0.0; 
  vout[6] = 0.0; 
  vout[7] = 0.0; 

  vout = &revout[8*5]; 
  vout[0] = 0.0; 
  vout[1] = w[5]; 
  vout[2] = 0.0; 
  vout[3] = w[5]/c; 
  vout[4] = 0.0; 
  vout[5] = 0.0; 
  vout[6] = 0.0; 
  vout[7] = 0.0; 

  vout = &revout[8*6]; 
  vout[0] = w[6]; 
  vout[1] = 0.0; 
  vout[2] = 0.0; 
  vout[3] = 0.0; 
  vout[4] = w[6]/c; 
  vout[5] = 0.0; 
  vout[6] = 0.0; 
  vout[7] = 0.0; 

  vout = &revout[8*7]; 
  vout[0] = 0.0; 
  vout[1] = w[7]; 
  vout[2] = 0.0; 
  vout[3] = -(1.0*w[7])/c; 
  vout[4] = 0.0; 
  vout[5] = 0.0; 
  vout[6] = 0.0; 
  vout[7] = 0.0; 

} 

static inline double 
maxwell_ev_0(void *ctx, const double *q, double *ev) 
{ 
  struct maxwell_eqn *eqn = ctx; 
  double c = eqn->c; 
  double c2 = eqn->c*eqn->c; 
  double efact = eqn->efact; 
  double bfact = eqn->bfact; 
  ev[0] = -c*bfact; 
  ev[1] = c*bfact; 
  ev[2] = -c*efact; 
  ev[3] = c*efact; 
  ev[4] = -c; 
  ev[5] = -c; 
  ev[6] = c; 
  ev[7] = c; 
  return c; 
} 

static inline double 
maxwell_ev_1(void *ctx, const double *q, double *ev) 
{ 
  struct maxwell_eqn *eqn = ctx; 
  double c = eqn->c; 
  double c2 = eqn->c*eqn->c; 
  double efact = eqn->efact; 
  double bfact = eqn->bfact; 
  ev[0] = -c*bfact; 
  ev[1] = c*bfact; 
  ev[2] = -c*efact; 
  ev[3] = c*efact; 
  ev[4] = -c; 
  ev[5] = -c; 
  ev[6] = c; 
  ev[7] = c; 
  return c; 
} 

static inline double 
maxwell_ev_2(void *ctx, const double *q, double *ev) 
{ 
  struct maxwell_eqn *eqn = ctx; 
  double c = eqn->c; 
  double c2 = eqn->c*eqn->c; 
  double efact = eqn->efact; 
  double bfact = eqn->bfact; 
  ev[0] = -c*bfact; 
  ev[1] = c*bfact; 
  ev[2] = -c*efact; 
  ev[3] = c*efact; 
  ev[4] = -c; 
  ev[5] = -c; 
  ev[6] = c; 
  ev[7] = c; 
  return c; 
} 


static inline void 
maxwell_ravg(void *ctx, const double *ql, const double *qr, double *qavg) 
{ 
  qavg[0] = 0.5*(ql[0]+qr[0]); 
  qavg[1] = 0.5*(ql[1]+qr[1]); 
  qavg[2] = 0.5*(ql[2]+qr[2]); 
  qavg[3] = 0.5*(ql[3]+qr[3]); 
  qavg[4] = 0.5*(ql[4]+qr[4]); 
  qavg[5] = 0.5*(ql[5]+qr[5]); 
  qavg[6] = 0.5*(ql[6]+qr[6]); 
  qavg[7] = 0.5*(ql[7]+qr[7]); 
} 

static inline double 
maxwell_fluct_0(void *ctx, const double *ql, const double *qr, double *apdq, double *amdq) 
{ 
  struct maxwell_eqn *eqn = ctx; 
  double c = eqn->c; 
  double c2 = eqn->c*eqn->c; 
  double efact = eqn->efact; 
  double bfact = eqn->bfact; 
  double dq[8]; 
  dq[0] = qr[0]-ql[0]; 
  dq[1] = qr[1]-ql[1]; 
  dq[2] = qr[2]-ql[2]; 
  dq[3] = qr[3]-ql[3]; 
  dq[4] = qr[4]-ql[4]; 
  dq[5] = qr[5]-ql[5]; 
  dq[6] = qr[6]-ql[6]; 
  dq[7] = qr[7]-ql[7]; 
  apdq[0] = c*(0.5*dq[0]+0.5*dq[6]*c)*efact; 
  apdq[1] = c*(0.5*dq[1]+0.5*dq[5]*c); 
  apdq[2] = c*(0.5*dq[2]-0.5*dq[4]*c); 
  apdq[3] = (0.5*dq[3]+(0.5*dq[7])/c)*c*bfact; 
  apdq[4] = (0.5*dq[4]-(0.5*dq[2])/c)*c; 
  apdq[5] = (0.5*dq[5]+(0.5*dq[1])/c)*c; 
  apdq[6] = (0.5*dq[6]+(0.5*dq[0])/c)*c*efact; 
  apdq[7] = c*(0.5*dq[7]+0.5*dq[3]*c)*bfact; 
  amdq[0] = -1.0*c*(0.5*dq[0]-0.5*dq[6]*c)*efact; 
  amdq[1] = -1.0*c*(0.5*dq[1]-0.5*dq[5]*c); 
  amdq[2] = -1.0*c*(0.5*dq[2]+0.5*dq[4]*c); 
  amdq[3] = -1.0*(0.5*dq[3]-(0.5*dq[7])/c)*c*bfact; 
  amdq[4] = -1.0*(0.5*dq[4]+(0.5*dq[2])/c)*c; 
  amdq[5] = -1.0*(0.5*dq[5]-(0.5*dq[1])/c)*c; 
  amdq[6] = -1.0*(0.5*dq[6]-(0.5*dq[0])/c)*c*efact; 
  amdq[7] = -1.0*c*(0.5*dq[7]-0.5*dq[3]*c)*bfact; 
  return c; 
} 

static inline double 
maxwell_fluct_1(void *ctx, const double *ql, const double *qr, double *apdq, double *amdq) 
{ 
  struct maxwell_eqn *eqn = ctx; 
  double c = eqn->c; 
  double c2 = eqn->c*eqn->c; 
  double efact = eqn->efact; 
  double bfact = eqn->bfact; 
  double dq[8]; 
  dq[0] = qr[0]-ql[0]; 
  dq[1] = qr[1]-ql[1]; 
  dq[2] = qr[2]-ql[2]; 
  dq[3] = qr[3]-ql[3]; 
  dq[4] = qr[4]-ql[4]; 
  dq[5] = qr[5]-ql[5]; 
  dq[6] = qr[6]-ql[6]; 
  dq[7] = qr[7]-ql[7]; 
  apdq[0] = c*(0.5*dq[0]-0.5*dq[5]*c); 
  apdq[1] = c*(0.5*dq[1]+0.5*dq[6]*c)*efact; 
  apdq[2] = c*(0.5*dq[2]+0.5*dq[3]*c); 
  apdq[3] = (0.5*dq[3]+(0.5*dq[2])/c)*c; 
  apdq[4] = (0.5*dq[4]+(0.5*dq[7])/c)*c*bfact; 
  apdq[5] = (0.5*dq[5]-(0.5*dq[0])/c)*c; 
  apdq[6] = (0.5*dq[6]+(0.5*dq[1])/c)*c*efact; 
  apdq[7] = c*(0.5*dq[7]+0.5*dq[4]*c)*bfact; 
  amdq[0] = -1.0*c*(0.5*dq[0]+0.5*dq[5]*c); 
  amdq[1] = -1.0*c*(0.5*dq[1]-0.5*dq[6]*c)*efact; 
  amdq[2] = -1.0*c*(0.5*dq[2]-0.5*dq[3]*c); 
  amdq[3] = -1.0*(0.5*dq[3]-(0.5*dq[2])/c)*c; 
  amdq[4] = -1.0*(0.5*dq[4]-(0.5*dq[7])/c)*c*bfact; 
  amdq[5] = -1.0*(0.5*dq[5]+(0.5*dq[0])/c)*c; 
  amdq[6] = -1.0*(0.5*dq[6]-(0.5*dq[1])/c)*c*efact; 
  amdq[7] = -1.0*c*(0.5*dq[7]-0.5*dq[4]*c)*bfact; 
  return c; 
} 

static inline double 
maxwell_fluct_2(void *ctx, const double *ql, const double *qr, double *apdq, double *amdq) 
{ 
  struct maxwell_eqn *eqn = ctx; 
  double c = eqn->c; 
  double c2 = eqn->c*eqn->c; 
  double efact = eqn->efact; 
  double bfact = eqn->bfact; 
  double dq[8]; 
  dq[0] = qr[0]-ql[0]; 
  dq[1] = qr[1]-ql[1]; 
  dq[2] = qr[2]-ql[2]; 
  dq[3] = qr[3]-ql[3]; 
  dq[4] = qr[4]-ql[4]; 
  dq[5] = qr[5]-ql[5]; 
  dq[6] = qr[6]-ql[6]; 
  dq[7] = qr[7]-ql[7]; 
  apdq[0] = c*(0.5*dq[0]+0.5*dq[4]*c); 
  apdq[1] = c*(0.5*dq[1]-0.5*dq[3]*c); 
  apdq[2] = c*(0.5*dq[2]+0.5*dq[6]*c)*efact; 
  apdq[3] = (0.5*dq[3]-(0.5*dq[1])/c)*c; 
  apdq[4] = (0.5*dq[4]+(0.5*dq[0])/c)*c; 
  apdq[5] = (0.5*dq[5]+(0.5*dq[7])/c)*c*bfact; 
  apdq[6] = (0.5*dq[6]+(0.5*dq[2])/c)*c*efact; 
  apdq[7] = c*(0.5*dq[7]+0.5*dq[5]*c)*bfact; 
  amdq[0] = -1.0*c*(0.5*dq[0]-0.5*dq[4]*c); 
  amdq[1] = -1.0*c*(0.5*dq[1]+0.5*dq[3]*c); 
  amdq[2] = -1.0*c*(0.5*dq[2]-0.5*dq[6]*c)*efact; 
  amdq[3] = -1.0*(0.5*dq[3]+(0.5*dq[1])/c)*c; 
  amdq[4] = -1.0*(0.5*dq[4]-(0.5*dq[0])/c)*c; 
  amdq[5] = -1.0*(0.5*dq[5]-(0.5*dq[7])/c)*c*bfact; 
  amdq[6] = -1.0*(0.5*dq[6]-(0.5*dq[2])/c)*c*efact; 
  amdq[7] = -1.0*c*(0.5*dq[7]-0.5*dq[5]*c)*bfact; 
  return c; 
} 

