#include <acutest.h>
#include <eqn_isoeuler.h>
#include <app.h>

struct eqn_sys
eqn_isoeuler_init(void *ctx)
{
  return (struct eqn_sys) {
    .meqn = 4,
    .mwave = 4,
    .eqn_ctx = &ctx,
    .flux = { isoeuler_flux_0, isoeuler_flux_1, isoeuler_flux_2  },
    .projon_left_ev = { isoeuler_projon_left_ev_0, isoeuler_projon_left_ev_1, isoeuler_projon_left_ev_2  },
    .recwith_right_ev = { isoeuler_recwith_right_ev_0, isoeuler_recwith_right_ev_1, isoeuler_recwith_right_ev_2 },
    .rescale_right_ev = { isoeuler_rescale_right_ev_0, isoeuler_rescale_right_ev_1, isoeuler_rescale_right_ev_2 },
    .ev = { isoeuler_ev_0, isoeuler_ev_1, isoeuler_ev_2 },
    .fluct = { isoeuler_roe_fluct_0, isoeuler_roe_fluct_1, isoeuler_roe_fluct_2 },
    .rotate_to_local = isoeuler_rotate_to_local,
    .rotate_to_global = isoeuler_rotate_to_global,    
  };
}

static void
prim_to_cons(double rho, double ux, double uy, double uz, double q[4])
{
  q[0] = rho;
  q[1] = rho*ux;
  q[2] = rho*uy;
  q[3] = rho*uz;
}

static void
show_vec(const char *nm, const double *q)
{
  printf("* %s\n", nm);
  for (int i=0; i<4; ++i)
    printf(" %.10lg\n", q[i]);
  printf("**\n");
}

void
test_phi_prime(void)
{
  struct isoeuler_eqn isoeuler = {.cs = 1.5};

  double v1[4], v2[4];
  double vec[4] = {1.0, 20.0, 30.0, 40.0};
  double q[4] = { 1.0, 1.5, 2.5, 3.5 };

  isoeuler_mulby_phi_prime(&isoeuler, q, vec, v1);
  isoeuler_mulby_phi_prime_inv(&isoeuler, q, v1, v2);

  for (int i=0; i<4; ++i)
    TEST_CHECK( gkyl_compare_double(v2[i], vec[i], 1e-14) );
}

void
test_ev(void)
{
  struct isoeuler_eqn isoeuler = {.cs = 1.5};

  double v1[4], v2[4];
  double vec[4] = {1.0, 20.0, 30.0, 40.0};
  double q[4] = { 1.0, 1.5, 2.5, 3.5 };

  struct eqn_sys eqn_isoeuler = eqn_isoeuler_init(&isoeuler);  
  
  for (int d=0; d<3; ++d) {
    eqn_isoeuler.projon_left_ev[d](&isoeuler, q, vec, v1);
    eqn_isoeuler.recwith_right_ev[d](&isoeuler, q, v1, v2);
    for (int i=0; i<4; ++i)
      TEST_CHECK( gkyl_compare_double(v2[i], vec[i], 1e-14) );
  }
}

void
test_rescale_ev(void)
{
  struct isoeuler_eqn isoeuler = {.cs = 1.5};

  double v1[4], v2[4], rev[4*4];
  double vec[4] = {1.0, 20.0, 30.0, 40.0 };
  double q[4] = { 1.0, 1.5, 2.5, 3.5 };

  struct eqn_sys eqn_isoeuler = eqn_isoeuler_init(&isoeuler);  
  
  for (int d=0; d<3; ++d) {
    eqn_isoeuler.projon_left_ev[d](&isoeuler, q, vec, v1);
    eqn_isoeuler.rescale_right_ev[d](&isoeuler, q, v1, rev);
    
    for (int i=0; i<4; ++i) {
      v2[i] = 0.0;
      for (int j=0; j<4; ++j)
        v2[i] += rev[4*j+i];
    }
    
    for (int i=0; i<4; ++i)
      TEST_CHECK( gkyl_compare_double(v2[i], vec[i], 1e-12) );
  }
}


void
test_flux_jump_roe(void)
{
  struct isoeuler_eqn isoeuler = { .cs = 1.5 };
  struct eqn_sys eqn_isoeuler = eqn_isoeuler_init(&isoeuler);

  double rhol = 1.0, uxl = 1.5, uyl = 2.5, uzl = 3.5, pl = 0.1;
  double rhor = 10.0, uxr = 10.5, uyr = 25.5, uzr = 35.5, pr = 100.0;
  
  double ql[4], qr[4];
  prim_to_cons(rhol, uxl, uyl, uzl, ql);
  prim_to_cons(rhor, uxr, uyr, uzr, qr);

  double dq[4];
  for (int i=0; i<4; ++i) dq[i] = qr[i]-ql[i];

  for (int d=0; d<3; ++d) {
    double apdq[4], amdq[4];    
    eqn_isoeuler.fluct[d](&isoeuler, ql, qr, apdq, amdq);

    double fr[4], fl[4];
    eqn_isoeuler.flux[d](&isoeuler, qr, fr);
    eqn_isoeuler.flux[d](&isoeuler, ql, fl);
    
    for (int i=0; i<4; ++i)
      TEST_CHECK( gkyl_compare_double(apdq[i]+amdq[i], fr[i]-fl[i], 1e-12) );
  }
}

void
test_rotate1(void)
{
  struct isoeuler_eqn isoeuler = { .cs = 1.5 };
  struct eqn_sys eqn_isoeuler = eqn_isoeuler_init(&isoeuler);
  
  double n[3] = { 0, 1.0, 0.0 };
  double tau1[3] = { -1.0, 0.0, 0.0 };
  double tau2[3] = { 0.0, 0.0, 1.0 };

  double rhol = 1.0, uxl = 1.5, uyl = 2.5, uzl = 3.5;
  
  double qglobal[4];
  prim_to_cons(rhol, uxl, uyl, uzl, qglobal);

  double qlocal[4];
  eqn_isoeuler.rotate_to_local(&isoeuler, n, tau1, tau2, qglobal, qlocal);
  double flocal[4];
  eqn_isoeuler.flux[0](&isoeuler, qlocal, flocal);
  double fglobal[4];
  eqn_isoeuler.rotate_to_global(&isoeuler, n, tau1, tau2, flocal, fglobal);

  double fy[4];
  eqn_isoeuler.flux[1](&isoeuler, qglobal, fy);

  for (int i=0; i<4; ++i)
    TEST_CHECK( gkyl_compare_double( fglobal[i], fy[i], 1e-14) );
}

TEST_LIST = {
  { "test_phi_prime", test_phi_prime },
  { "test_ev", test_ev },
  { "test_rescale_ev", test_rescale_ev },  
  { "test_flux_jump_roe", test_flux_jump_roe },
  { "test_rotate", test_rotate1 },    
  { NULL, NULL },
};    
