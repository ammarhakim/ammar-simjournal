#include <ga3.h>

#include <float.h>
#include <stdio.h>
#include <math.h>
#include <acutest.h>

int
compare_double(double a, double b, double eps)
{
  if (isnan(a) || isnan(b)) return 0;
  
  double absa = fabs(a), absb = fabs(b), diff = fabs(a-b);
  if (a == b) return 1;
  if (a == 0 || b == 0 || (absa+absb < DBL_MIN)) return diff < eps;
  if (absa < eps) return diff < eps;
  if (absb < eps) return diff < eps;
  return diff/fmin(absa+absb, DBL_MAX) < eps;  
}
    
void
test_v1(void)
{
  ga3_vec1 a = { .x = { 1.0, 2.0, 3.0 } };
  ga3_vec1 b = { .x = { 2.0, 3.0, 1.0 } };
  ga3_vec1 c = { .x = { 8.0, 9.0, 10.0 } };

  ga3_vec0 adotb = ga3_dot_v1_v1(a, b);
  TEST_CHECK( adotb.x[0] == (2+6+3) );

  ga3_vec2 awedgeb = ga3_wedge_v1_v1(a, b);
  TEST_CHECK( awedgeb.x[0] == -1.0 );
  TEST_CHECK( awedgeb.x[1] == -5.0 );
  TEST_CHECK( awedgeb.x[2] == -7.0 );

  // checks if a | ( b & c) = (a|b) c - (a|c) b
  ga3_vec1 t1 = ga3_dot_v1_v2(a, ga3_wedge_v1_v1(b, c));

  ga3_vec1 t2a = ga3_scale_vec1(ga3_dot_v1_v1(a, b).x[0], c);
  ga3_vec1 t2b = ga3_scale_vec1(ga3_dot_v1_v1(a, c).x[0], b);

  TEST_CHECK( t1.x[0] == (t2a.x[0] - t2b.x[0]) );
  TEST_CHECK( t1.x[1] == (t2a.x[1] - t2b.x[1]) );
  TEST_CHECK( t1.x[2] == (t2a.x[2] - t2b.x[2]) );

  // reverse
  ga3_vec1 ar = ga3_reverse_vec1(a);
  for (int i=0; i<3; ++i)
    TEST_CHECK( ar.x[i] == a.x[i] );

  // dual
  ga3_vec2 ad = ga3_hodge_vec1(a);
  TEST_CHECK( ad.x[0] == a.x[2] );
  TEST_CHECK( ad.x[1] == -a.x[1] );
  TEST_CHECK( ad.x[2] == a.x[0] );

  ga3_vec1 add = ga3_hodge_vec2(ga3_hodge_vec1(a));
  for (int i=0; i<3; ++i)
    TEST_CHECK( fabs(ar.x[i]) == fabs(a.x[i]) );

  ga3_mv agab = ga3_gp_mv_mv( (ga3_mv) { .v1 = a }, (ga3_mv) { .v1 = b });

  TEST_CHECK( agab.v0.x[0] == adotb.x[0] );

  for (int i=0; i<3; ++i)
    TEST_CHECK( agab.v2.x[i] == awedgeb.x[i] );
}

void
test_v2(void)
{
  ga3_vec2 A = { .x = { 1.0, 2.0, 3.0 } };
  ga3_vec2 B = { .x = { 2.0, 3.0, 1.0 } };
  ga3_vec1 c = { .x = { 8.0, 9.0, 10.0 } };

  ga3_vec0 AdotB = ga3_dot_v2_v2(A, B);
  TEST_CHECK( AdotB.x[0] == -A.x[0]*B.x[0]-A.x[1]*B.x[1]-A.x[2]*B.x[2] );

  ga3_vec3 cwedgeA = ga3_wedge_v1_v2(c, A);
  TEST_CHECK( cwedgeA.x[0] == -A.x[1]*c.x[1] + c.x[0]*A.x[2] + A.x[0]*c.x[2] );

  ga3_vec1 cdotA = ga3_dot_v1_v2(c, A);
  TEST_CHECK( cdotA.x[0] == -A.x[0]*c.x[1]-A.x[1]*c.x[2] );
  TEST_CHECK( cdotA.x[1] == A.x[0]*c.x[0]-A.x[2]*c.x[2] );
  TEST_CHECK( cdotA.x[2] == A.x[1]*c.x[0]+A.x[2]*c.x[1] );

  ga3_vec1 Adotc = ga3_dot_v2_v1(A, c);
  for (int i=0; i<3; ++i)
    TEST_CHECK( Adotc.x[i] == -cdotA.x[i] );
}

void
test_v3(void)
{
  ga3_vec3 A = { .x = { 3.0 } };
  ga3_vec1 c = { .x = { 8.0, 9.0, 10.0 } };

  ga3_vec2 cdotA = ga3_dot_v1_v3(c, A);
  TEST_CHECK( cdotA.x[0] == A.x[0]*c.x[2] );
  TEST_CHECK( cdotA.x[1] == -A.x[0]*c.x[1] );
  TEST_CHECK( cdotA.x[2] == A.x[0]*c.x[0] );
}
  
void
test_even()
{
  ga3_even e1 = {
    .v0 = { 1.0 },
    .v2 = { 5.0, 6.0, 7.0 }
  };

  ga3_even e1r = ga3_reverse_even(e1);
  ga3_even e1e1r = ga3_gp_even(e1, e1r);

  // no bivector part
  for (int i=0; i<3; ++i)
    TEST_CHECK( e1e1r.v2.x[0] == 0.0 );

  double e1norm = ga3_norm_even(e1);
  TEST_CHECK( e1norm == e1e1r.v0.x[0] );

  ga3_even R = ga3_scale_even(1/sqrt(e1norm), e1);
  TEST_CHECK( ga3_norm_even(R) == 1.0 );

  double rmat[9];
  ga3_rotmat_vec1(R, rmat);
      
  double isunit[9] = { 0.0 };
  for (int j=0; j<3; ++j) {
    for (int i=0; i<3; ++i) {
      isunit[idx3x3(i,j)] = 0.0;
      for (int k=0; k<3; ++k)
        isunit[idx3x3(i,j)] += rmat[idx3x3(i,k)]*rmat[idx3x3(j,k)];
    }
  }

  for (int k=0; k<3; ++k)
    TEST_CHECK( compare_double(1.0, isunit[idx3x3(k,k)], 1e-15) );

  ga3_vec1 a = { .x = { 1.0, 2.0, 3.0 } };
  ga3_vec1 ap;
  matvec_3x3(rmat, a.x, ap.x);
  TEST_CHECK( compare_double(ga3_len(a), ga3_len(ap), 1e-15) );
}

void
test_mv(void)
{
  ga3_vec1 E = { .x = { 1.0, 2.0, 3.0 } };
  ga3_vec1 B = { .x = { 2.0, 3.0, 1.0 } };

  ga3_mv F = { .v1 = E, .v2 = ga3_hodge_vec1(B) };
  ga3_mv Fr = ga3_reverse_mv(F);

  ga3_mv FFr = ga3_gp_mv_mv(F, Fr);

  ga3_vec0 Er = ga3_sel_0(FFr);
  ga3_vec1 Mom = ga3_sel_1(FFr);

  double Er_t = ga3_dot_v1_v1(E,E).x[0] + ga3_dot_v1_v1(B,B).x[0];

  TEST_CHECK( Er_t == Er.x[0] );

  for (int i=0; i<3; ++i)
    TEST_CHECK( FFr.v2.x[0] == 0.0 );

  ga3_vec1 ExB = ga3_cross(E, B);

  for (int i=0; i<3; ++i)
    TEST_CHECK( 0.5*FFr.v1.x[i] == ExB.x[i] );
}

TEST_LIST = {
  { "test_v1", test_v1 },
  { "test_v2", test_v2 },
  { "test_v3", test_v3 },
  { "test_even", test_even },
  { "test_mv", test_mv },  
  { NULL, NULL },
};
