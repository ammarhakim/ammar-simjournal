#pragma once

#include <math.h>
#include <stdio.h>
#include <ga3_generated.h>

static inline double ga3_sq(double x) { return x*x; }

// some additional functions
static ga3_mv ga3_reverse_mv(ga3_mv mv);
static ga3_mv ga3_hodge_mv(ga3_mv mv);
static ga3_mv ga3_sum_mv(double sa, ga3_mv mva, double sb, ga3_mv mvb);
static void ga3_show_mv(const char *nm, FILE *fp, ga3_mv mv);
static void matvec_3x3(const double mat[9], const double v[3], double out[3]);

// functions to work with even-grade elements
static double ga3_norm_even(ga3_even M);

// functions to work with 3D vectors
static double ga3_dot(ga3_vec1 a, ga3_vec1 b);
static ga3_vec1 ga3_cross(ga3_vec1 a, ga3_vec1 b);
static double ga3_len(ga3_vec1 a);
static double ga3_triple(ga3_vec1 a, ga3_vec1 b, ga3_vec1 c);

// Implementation of static functions

static inline void
matvec_3x3(const double mat[9], const double v[3], double out[3])
{
  out[0] = mat[idx3x3(0,0)]*v[0] + mat[idx3x3(0,1)]*v[1] + mat[idx3x3(0,2)]*v[2];
  out[1] = mat[idx3x3(1,0)]*v[0] + mat[idx3x3(1,1)]*v[1] + mat[idx3x3(1,2)]*v[2];
  out[2] = mat[idx3x3(2,0)]*v[0] + mat[idx3x3(2,1)]*v[1] + mat[idx3x3(2,2)]*v[2];
}

static inline ga3_mv
ga3_reverse_mv(ga3_mv mv)
{
  return (ga3_mv) {
    .v0 = mv.v0,
    .v1 = mv.v1,
    .v2 = ga3_reverse_vec2(mv.v2),
    .v3 = ga3_reverse_vec3(mv.v3)
  };
}

static inline ga3_mv
ga3_hodge_mv(ga3_mv mv)
{
  return (ga3_mv) {
    .v0 = ga3_hodge_vec3(mv.v3),
    .v1 = ga3_hodge_vec2(mv.v2),
    .v2 = ga3_hodge_vec1(mv.v1),
    .v3 = ga3_hodge_vec0(mv.v0),
  };
}

static inline ga3_mv
ga3_sum_mv(double sa, ga3_mv mva, double sb, ga3_mv mvb)
{
  struct ga3_mv rv;
  rv.v0 = ga3_sum_vec0(sa, mva.v0, sb, mvb.v0);
  rv.v1 = ga3_sum_vec1(sa, mva.v1, sb, mvb.v1);
  rv.v2 = ga3_sum_vec2(sa, mva.v2, sb, mvb.v2);
  rv.v3 = ga3_sum_vec3(sa, mva.v3, sb, mvb.v3);
  return rv;
}

static void
ga3_show_mv(const char *nm, FILE *fp, ga3_mv mv)
{
  fprintf(fp, "-- %s --\n", nm);
  fprintf(fp, "v0: %lg\n", mv.v0.x[0]);
  fprintf(fp, "v1: [%lg, %lg, %lg]\n", mv.v1.x[0], mv.v1.x[1], mv.v1.x[2]);
  fprintf(fp, "v2: [%lg, %lg, %lg]\n", mv.v2.x[0], mv.v2.x[1], mv.v2.x[2]);
  fprintf(fp, "v3: %lg\n", mv.v3.x[0]);
}

static inline double
ga3_norm_even(ga3_even M)
{
  return ga3_sq(M.v0.x[0])+ga3_sq(M.v2.x[0])+ga3_sq(M.v2.x[1])+ga3_sq(M.v2.x[2]);
}

static inline double
ga3_dot(ga3_vec1 a, ga3_vec1 b)
{
  return ga3_dot_v1_v1(a, b).x[0];
}

static inline double
ga3_len(ga3_vec1 a)
{
  return sqrt( ga3_dot_v1_v1(a, a).x[0] );
}

static inline ga3_vec1
ga3_cross(ga3_vec1 a, ga3_vec1 b)
{
  // a x b = - (a wedge b ) I_3
  return ga3_scale_vec1(-1.0, ga3_hodge_vec2(ga3_wedge_v1_v1(a, b)));
}

// a \dot (b \times c)
static inline double
ga3_triple(ga3_vec1 a, ga3_vec1 b, ga3_vec1 c)
{
  return ga3_dot(a, ga3_cross(b, c));
}
