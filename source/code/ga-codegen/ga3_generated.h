#pragma once 
// Ga(3) auto-generated code. Do not edit by hand. 
// Include the ga.h header rather than this directly. 

// k-vector types 
typedef struct ga3_vec0 { double x[1]; } ga3_vec0; 
typedef struct ga3_vec1 { double x[3]; } ga3_vec1; 
typedef struct ga3_vec2 { double x[3]; } ga3_vec2; 
typedef struct ga3_vec3 { double x[1]; } ga3_vec3; 

typedef struct ga3_mv { 
  ga3_vec0 v0; 
  ga3_vec1 v1; 
  ga3_vec2 v2; 
  ga3_vec3 v3; 
  } ga3_mv; 

typedef struct ga3_even {
  ga3_vec0 v0;
  ga3_vec2 v2;
} ga3_even; 

static inline int idx3x3(int r, int c) { return c*3+r; }
static ga3_vec1 ga3_scale_vec1(double s, ga3_vec1 a); 
static ga3_vec2 ga3_scale_vec2(double s, ga3_vec2 a); 
static ga3_vec3 ga3_scale_vec3(double s, ga3_vec3 a); 
static ga3_even ga3_scale_even(double s, ga3_even a); 

static ga3_vec1 ga3_reverse_vec1(ga3_vec1 a); 
static ga3_vec2 ga3_reverse_vec2(ga3_vec2 a); 
static ga3_vec3 ga3_reverse_vec3(ga3_vec3 a); 
static ga3_even ga3_reverse_even(ga3_even a); 

static ga3_vec3 ga3_hodge_vec0(ga3_vec0 a); 
static ga3_vec2 ga3_hodge_vec1(ga3_vec1 a); 
static ga3_vec1 ga3_hodge_vec2(ga3_vec2 a); 
static ga3_vec0 ga3_hodge_vec3(ga3_vec3 a); 

static ga3_vec0 ga3_sum_vec0(double sa, ga3_vec0 a, double sb, ga3_vec0 b); 
static ga3_vec1 ga3_sum_vec1(double sa, ga3_vec1 a, double sb, ga3_vec1 b); 
static ga3_vec2 ga3_sum_vec2(double sa, ga3_vec2 a, double sb, ga3_vec2 b); 
static ga3_vec3 ga3_sum_vec3(double sa, ga3_vec3 a, double sb, ga3_vec3 b); 

static ga3_vec0 ga3_dot_v1_v1(ga3_vec1 a, ga3_vec1 b); 
static ga3_vec1 ga3_dot_v1_v2(ga3_vec1 a, ga3_vec2 b); 
static ga3_vec2 ga3_dot_v1_v3(ga3_vec1 a, ga3_vec3 b); 
static ga3_vec1 ga3_dot_v2_v1(ga3_vec2 a, ga3_vec1 b); 
static ga3_vec0 ga3_dot_v2_v2(ga3_vec2 a, ga3_vec2 b); 
static ga3_vec1 ga3_dot_v2_v3(ga3_vec2 a, ga3_vec3 b); 
static ga3_vec2 ga3_dot_v3_v1(ga3_vec3 a, ga3_vec1 b); 
static ga3_vec1 ga3_dot_v3_v2(ga3_vec3 a, ga3_vec2 b); 
static ga3_vec0 ga3_dot_v3_v3(ga3_vec3 a, ga3_vec3 b); 

static ga3_vec2 ga3_wedge_v1_v1(ga3_vec1 a, ga3_vec1 b); 
static ga3_vec3 ga3_wedge_v1_v2(ga3_vec1 a, ga3_vec2 b); 
static ga3_vec3 ga3_wedge_v2_v1(ga3_vec2 a, ga3_vec1 b); 

static ga3_vec0 ga3_sel_0(ga3_mv mv); 
static ga3_vec1 ga3_sel_1(ga3_mv mv); 
static ga3_vec2 ga3_sel_2(ga3_mv mv); 
static ga3_vec3 ga3_sel_3(ga3_mv mv); 

static ga3_even ga3_gp_even(ga3_even M, ga3_even N); 
static ga3_mv ga3_gp_mv_mv(ga3_mv M, ga3_mv N); 
static double ga3_scalar_mv_mv(ga3_mv M, ga3_mv N); 

// R is a rotor: R reverse(R) = 1 
static void ga3_rotmat_vec1(ga3_even R, double rmat[3*3]);
// Scale k-vectors 
static inline ga3_vec1 
ga3_scale_vec1(double s, ga3_vec1 a) 
{
  ga3_vec1 rv; 
  rv.x[0] = s*a.x[0]; 
  rv.x[1] = s*a.x[1]; 
  rv.x[2] = s*a.x[2]; 
  return rv;
}

static inline ga3_vec2 
ga3_scale_vec2(double s, ga3_vec2 a) 
{
  ga3_vec2 rv; 
  rv.x[0] = s*a.x[0]; 
  rv.x[1] = s*a.x[1]; 
  rv.x[2] = s*a.x[2]; 
  return rv;
}

static inline ga3_vec3 
ga3_scale_vec3(double s, ga3_vec3 a) 
{
  ga3_vec3 rv; 
  rv.x[0] = s*a.x[0]; 
  return rv;
}

static inline ga3_even 
ga3_scale_even(double s, ga3_even a) 
{
  ga3_even rv; 
  rv.v0.x[0] = s*a.v0.x[0]; 
  rv.v2 = ga3_scale_vec2(s, a.v2); 
  return rv;
}

// Reverse k-vectors 
static inline ga3_vec1 
ga3_reverse_vec1(ga3_vec1 a) 
{
  const double *ax = a.x; 
  ga3_vec1 rv; 
  rv.x[0] = ax[0]; 
  rv.x[1] = ax[1]; 
  rv.x[2] = ax[2]; 
  return rv;
}

static inline ga3_vec2 
ga3_reverse_vec2(ga3_vec2 a) 
{
  const double *ax = a.x; 
  ga3_vec2 rv; 
  rv.x[0] = -ax[0]; 
  rv.x[1] = -ax[1]; 
  rv.x[2] = -ax[2]; 
  return rv;
}

static inline ga3_vec3 
ga3_reverse_vec3(ga3_vec3 a) 
{
  const double *ax = a.x; 
  ga3_vec3 rv; 
  rv.x[0] = -ax[0]; 
  return rv;
}

static inline ga3_even 
ga3_reverse_even(ga3_even a) 
{
  ga3_even rv; 
  rv.v0.x[0] = a.v0.x[0]; 
  rv.v2 = ga3_reverse_vec2(a.v2); 
  return rv;
}

// Hodge-dual of k-vectors: A^dagger = A I_3 
static inline ga3_vec3 
ga3_hodge_vec0(ga3_vec0 a) 
{
  const double *ax = a.x; 
  ga3_vec3 rv; 
  rv.x[0] = ax[0]; 
  return rv;
}

static inline ga3_vec2 
ga3_hodge_vec1(ga3_vec1 a) 
{
  const double *ax = a.x; 
  ga3_vec2 rv; 
  rv.x[0] = ax[2]; 
  rv.x[1] = -ax[1]; 
  rv.x[2] = ax[0]; 
  return rv;
}

static inline ga3_vec1 
ga3_hodge_vec2(ga3_vec2 a) 
{
  const double *ax = a.x; 
  ga3_vec1 rv; 
  rv.x[0] = -ax[2]; 
  rv.x[1] = ax[1]; 
  rv.x[2] = -ax[0]; 
  return rv;
}

static inline ga3_vec0 
ga3_hodge_vec3(ga3_vec3 a) 
{
  const double *ax = a.x; 
  ga3_vec0 rv; 
  rv.x[0] = -ax[0]; 
  return rv;
}

// Linear combination of k-vectors 
static ga3_vec0 
ga3_sum_vec0(double sa, ga3_vec0 a, double sb, ga3_vec0 b) 
{
  ga3_vec0 rv; 
  rv.x[0] = sa*a.x[0] + sb*b.x[0]; 
  return rv;
}
static ga3_vec1 
ga3_sum_vec1(double sa, ga3_vec1 a, double sb, ga3_vec1 b) 
{
  ga3_vec1 rv; 
  rv.x[0] = sa*a.x[0] + sb*b.x[0]; 
  rv.x[1] = sa*a.x[1] + sb*b.x[1]; 
  rv.x[2] = sa*a.x[2] + sb*b.x[2]; 
  return rv;
}
static ga3_vec2 
ga3_sum_vec2(double sa, ga3_vec2 a, double sb, ga3_vec2 b) 
{
  ga3_vec2 rv; 
  rv.x[0] = sa*a.x[0] + sb*b.x[0]; 
  rv.x[1] = sa*a.x[1] + sb*b.x[1]; 
  rv.x[2] = sa*a.x[2] + sb*b.x[2]; 
  return rv;
}
static ga3_vec3 
ga3_sum_vec3(double sa, ga3_vec3 a, double sb, ga3_vec3 b) 
{
  ga3_vec3 rv; 
  rv.x[0] = sa*a.x[0] + sb*b.x[0]; 
  return rv;
}

// Dot products between k-vectors 
static inline ga3_vec0 
ga3_dot_v1_v1(ga3_vec1 a, ga3_vec1 b)
{
  const double *ax = a.x; 
  const double *bx = b.x; 
  ga3_vec0 rv; 
  rv.x[0] = ax[0]*bx[0]+ax[1]*bx[1]+ax[2]*bx[2]; 
  return rv;
}

static inline ga3_vec1 
ga3_dot_v1_v2(ga3_vec1 a, ga3_vec2 b)
{
  const double *ax = a.x; 
  const double *bx = b.x; 
  ga3_vec1 rv; 
  rv.x[0] = (-bx[0]*ax[1])-bx[1]*ax[2]; 
  rv.x[1] = ax[0]*bx[0]-ax[2]*bx[2]; 
  rv.x[2] = ax[0]*bx[1]+ax[1]*bx[2]; 
  return rv;
}

static inline ga3_vec2 
ga3_dot_v1_v3(ga3_vec1 a, ga3_vec3 b)
{
  const double *ax = a.x; 
  const double *bx = b.x; 
  ga3_vec2 rv; 
  rv.x[0] = bx[0]*ax[2]; 
  rv.x[1] = -bx[0]*ax[1]; 
  rv.x[2] = ax[0]*bx[0]; 
  return rv;
}


static inline ga3_vec1 
ga3_dot_v2_v1(ga3_vec2 a, ga3_vec1 b)
{
  const double *ax = a.x; 
  const double *bx = b.x; 
  ga3_vec1 rv; 
  rv.x[0] = ax[0]*bx[1]+ax[1]*bx[2]; 
  rv.x[1] = (-ax[0]*bx[0])+ax[2]*bx[2]; 
  rv.x[2] = (-bx[0]*ax[1])-bx[1]*ax[2]; 
  return rv;
}

static inline ga3_vec0 
ga3_dot_v2_v2(ga3_vec2 a, ga3_vec2 b)
{
  const double *ax = a.x; 
  const double *bx = b.x; 
  ga3_vec0 rv; 
  rv.x[0] = (-ax[0]*bx[0])-ax[1]*bx[1]-ax[2]*bx[2]; 
  return rv;
}

static inline ga3_vec1 
ga3_dot_v2_v3(ga3_vec2 a, ga3_vec3 b)
{
  const double *ax = a.x; 
  const double *bx = b.x; 
  ga3_vec1 rv; 
  rv.x[0] = -bx[0]*ax[2]; 
  rv.x[1] = bx[0]*ax[1]; 
  rv.x[2] = -ax[0]*bx[0]; 
  return rv;
}


static inline ga3_vec2 
ga3_dot_v3_v1(ga3_vec3 a, ga3_vec1 b)
{
  const double *ax = a.x; 
  const double *bx = b.x; 
  ga3_vec2 rv; 
  rv.x[0] = ax[0]*bx[2]; 
  rv.x[1] = -ax[0]*bx[1]; 
  rv.x[2] = ax[0]*bx[0]; 
  return rv;
}

static inline ga3_vec1 
ga3_dot_v3_v2(ga3_vec3 a, ga3_vec2 b)
{
  const double *ax = a.x; 
  const double *bx = b.x; 
  ga3_vec1 rv; 
  rv.x[0] = -ax[0]*bx[2]; 
  rv.x[1] = ax[0]*bx[1]; 
  rv.x[2] = -ax[0]*bx[0]; 
  return rv;
}

static inline ga3_vec0 
ga3_dot_v3_v3(ga3_vec3 a, ga3_vec3 b)
{
  const double *ax = a.x; 
  const double *bx = b.x; 
  ga3_vec0 rv; 
  rv.x[0] = -ax[0]*bx[0]; 
  return rv;
}



// Wedge products between k-vectors 
static inline ga3_vec2 
ga3_wedge_v1_v1(ga3_vec1 a, ga3_vec1 b)
{
  const double *ax = a.x; 
  const double *bx = b.x; 
  ga3_vec2 rv; 
  rv.x[0] = (-bx[0]*ax[1])+ax[0]*bx[1]; 
  rv.x[1] = (-bx[0]*ax[2])+ax[0]*bx[2]; 
  rv.x[2] = (-bx[1]*ax[2])+ax[1]*bx[2]; 
  return rv;
}

static inline ga3_vec3 
ga3_wedge_v1_v2(ga3_vec1 a, ga3_vec2 b)
{
  const double *ax = a.x; 
  const double *bx = b.x; 
  ga3_vec3 rv; 
  rv.x[0] = (-ax[1]*bx[1])+bx[0]*ax[2]+ax[0]*bx[2]; 
  return rv;
}


static inline ga3_vec3 
ga3_wedge_v2_v1(ga3_vec2 a, ga3_vec1 b)
{
  const double *ax = a.x; 
  const double *bx = b.x; 
  ga3_vec3 rv; 
  rv.x[0] = (-ax[1]*bx[1])+bx[0]*ax[2]+ax[0]*bx[2]; 
  return rv;
}



static inline ga3_vec0 
ga3_sel_0(ga3_mv mv) 
{
  return mv.v0; 
}
static inline ga3_vec1 
ga3_sel_1(ga3_mv mv) 
{
  return mv.v1; 
}
static inline ga3_vec2 
ga3_sel_2(ga3_mv mv) 
{
  return mv.v2; 
}
static inline ga3_vec3 
ga3_sel_3(ga3_mv mv) 
{
  return mv.v3; 
}

static inline ga3_even 
ga3_gp_even(ga3_even M, ga3_even N) 
{
  const double *mv0 = M.v0.x; 
  const double *mv2 = M.v2.x; 

  const double *nv0 = N.v0.x; 
  const double *nv2 = N.v2.x; 

  struct ga3_even rv; 
  double *rvv[4]; 
  rvv[0] = rv.v0.x; 
  rvv[2] = rv.v2.x; 

  rvv[0][0] = mv0[0]*nv0[0]-mv2[0]*nv2[0]-mv2[1]*nv2[1]-mv2[2]*nv2[2]; 
  rvv[2][0] = mv2[0]*nv0[0]+mv0[0]*nv2[0]+nv2[1]*mv2[2]-mv2[1]*nv2[2]; 
  rvv[2][1] = nv0[0]*mv2[1]+mv0[0]*nv2[1]-nv2[0]*mv2[2]+mv2[0]*nv2[2]; 
  rvv[2][2] = nv2[0]*mv2[1]-mv2[0]*nv2[1]+nv0[0]*mv2[2]+mv0[0]*nv2[2]; 

  return rv; 
}
static inline ga3_mv 
ga3_gp_mv_mv(ga3_mv M, ga3_mv N) 
{
  const double *mv0 = M.v0.x; 
  const double *mv1 = M.v1.x; 
  const double *mv2 = M.v2.x; 
  const double *mv3 = M.v3.x; 

  const double *nv0 = N.v0.x; 
  const double *nv1 = N.v1.x; 
  const double *nv2 = N.v2.x; 
  const double *nv3 = N.v3.x; 

  struct ga3_mv rv; 
  double *rvv[4]; 
  rvv[0] = rv.v0.x; 
  rvv[1] = rv.v1.x; 
  rvv[2] = rv.v2.x; 
  rvv[3] = rv.v3.x; 

  rvv[0][0] = mv0[0]*nv0[0]+mv1[0]*nv1[0]-mv2[0]*nv2[0]-mv3[0]*nv3[0]+mv1[1]*nv1[1]-mv2[1]*nv2[1]+mv1[2]*nv1[2]-mv2[2]*nv2[2]; 
  rvv[1][0] = mv1[0]*nv0[0]+mv0[0]*nv1[0]-nv2[0]*mv1[1]+mv2[0]*nv1[1]-nv2[1]*mv1[2]-nv3[0]*mv2[2]+mv2[1]*nv1[2]-mv3[0]*nv2[2]; 
  rvv[1][1] = (-mv2[0]*nv1[0])+mv1[0]*nv2[0]+nv0[0]*mv1[1]+nv3[0]*mv2[1]+mv0[0]*nv1[1]+mv3[0]*nv2[1]+mv2[2]*nv1[2]-mv1[2]*nv2[2]; 
  rvv[1][2] = (-mv3[0]*nv2[0])-mv2[0]*nv3[0]-nv1[0]*mv2[1]+mv1[0]*nv2[1]+nv0[0]*mv1[2]-nv1[1]*mv2[2]+mv0[0]*nv1[2]+mv1[1]*nv2[2]; 

  rvv[2][0] = mv2[0]*nv0[0]+mv0[0]*nv2[0]-nv1[0]*mv1[1]+mv1[0]*nv1[1]+nv3[0]*mv1[2]+nv2[1]*mv2[2]+mv3[0]*nv1[2]-mv2[1]*nv2[2]; 
  rvv[2][1] = (-nv3[0]*mv1[1])+nv0[0]*mv2[1]-mv3[0]*nv1[1]+mv0[0]*nv2[1]-nv1[0]*mv1[2]-nv2[0]*mv2[2]+mv1[0]*nv1[2]+mv2[0]*nv2[2]; 
  rvv[2][2] = mv3[0]*nv1[0]+mv1[0]*nv3[0]+nv2[0]*mv2[1]-mv2[0]*nv2[1]-nv1[1]*mv1[2]+nv0[0]*mv2[2]+mv1[1]*nv1[2]+mv0[0]*nv2[2]; 

  rvv[3][0] = mv3[0]*nv0[0]+mv0[0]*nv3[0]-mv2[1]*nv1[1]-mv1[1]*nv2[1]+nv2[0]*mv1[2]+nv1[0]*mv2[2]+mv2[0]*nv1[2]+mv1[0]*nv2[2]; 

  return rv; 
}

static inline double 
ga3_scalar_mv_mv(ga3_mv M, ga3_mv N) 
{
  const double *mv0 = M.v0.x; 
  const double *mv1 = M.v1.x; 
  const double *mv2 = M.v2.x; 
  const double *mv3 = M.v3.x; 

  const double *nv0 = N.v0.x; 
  const double *nv1 = N.v1.x; 
  const double *nv2 = N.v2.x; 
  const double *nv3 = N.v3.x; 

  return mv0[0]*nv0[0]+mv1[0]*nv1[0]-mv2[0]*nv2[0]-mv3[0]*nv3[0]+mv1[1]*nv1[1]-mv2[1]*nv2[1]+mv1[2]*nv1[2]-mv2[2]*nv2[2]; 
}


static void
ga3_rotmat_vec1(ga3_even M, double R[3*3]) 
{
  const double *mv0 = M.v0.x; 
  const double *mv2 = M.v2.x; 
  double mv0_0_2 = mv0[0]*mv0[0]; 
  double mv2_0_2 = mv2[0]*mv2[0]; 
  double mv2_1_2 = mv2[1]*mv2[1]; 
  double mv2_2_2 = mv2[2]*mv2[2]; 

  R[0] = mv0_0_2-mv2_0_2-mv2_1_2+mv2_2_2; 
  R[1] = (-2*mv0[0]*mv2[0])-2*mv2[1]*mv2[2]; 
  R[2] = (-2*mv0[0]*mv2[1])+2*mv2[0]*mv2[2]; 
  R[3] = 2*mv0[0]*mv2[0]-2*mv2[1]*mv2[2]; 
  R[4] = mv0_0_2-mv2_0_2+mv2_1_2-mv2_2_2; 
  R[5] = (-2*mv2[0]*mv2[1])-2*mv0[0]*mv2[2]; 
  R[6] = 2*mv0[0]*mv2[1]+2*mv2[0]*mv2[2]; 
  R[7] = (-2*mv2[0]*mv2[1])+2*mv0[0]*mv2[2]; 
  R[8] = mv0_0_2+mv2_0_2-mv2_1_2-mv2_2_2; 
}
