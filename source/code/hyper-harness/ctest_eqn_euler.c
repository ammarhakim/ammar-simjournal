#include <acutest.h>
#include <eqn_euler.h>
#include <app.h>

struct eqn_sys
eqn_euler_init(void *ctx)
{
  return (struct eqn_sys) {
    .meqn = 5,
    .mwave = 5,
    .eqn_ctx = &ctx,
    .flux = { euler_flux_0, euler_flux_1, euler_flux_2  },
    .projon_left_ev = { euler_projon_left_ev_0, euler_projon_left_ev_1, euler_projon_left_ev_2  },
    .recwith_right_ev = { euler_recwith_right_ev_0, euler_recwith_right_ev_1, euler_recwith_right_ev_2 },
    .rescale_right_ev = { euler_rescale_right_ev_0, euler_rescale_right_ev_1, euler_rescale_right_ev_2 },
    .ev = { euler_ev_0, euler_ev_1, euler_ev_2 },
    .fluct = { euler_roe_fluct_0, euler_roe_fluct_1, euler_roe_fluct_2 }
  };
}

static void
prim_to_cons(double gg, double rho, double ux, double uy, double uz, double p, double q[5])
{
  q[0] = rho;
  q[1] = rho*ux;
  q[2] = rho*uy;
  q[3] = rho*uz;
  double u2 = ux*ux + uy*uy + uz*uz;
  q[4] = 0.5*rho*u2 + p/(gg-1.0);
}

static void
show_vec(const char *nm, const double *q)
{
  printf("* %s\n", nm);
  for (int i=0; i<5; ++i)
    printf(" %.10lg\n", q[i]);
  printf("**\n");
}

void
test_phi_prime(void)
{
  struct euler_eqn euler = {.gas_gamma = 1.4};

  double v1[5], v2[5];
  double vec[5] = {1.0, 20.0, 30.0, 40.0, 5.0};
  double q[5] = { 1.0, 1.5, 2.5, 3.5, 10.5 };

  euler_mulby_phi_prime(&euler, q, vec, v1);
  euler_mulby_phi_prime_inv(&euler, q, v1, v2);

  for (int i=0; i<5; ++i)
    TEST_CHECK( gkyl_compare_double(v2[i], vec[i], 1e-14) );
}

void
test_ev(void)
{
  struct euler_eqn euler = {.gas_gamma = 1.4};

  double v1[5], v2[5];
  double vec[5] = {1.0, 20.0, 30.0, 40.0, 5.0};
  double q[5] = { 1.0, 1.5, 2.5, 3.5, 10.5 };

  struct eqn_sys eqn_euler = eqn_euler_init(&euler);  
  
  for (int d=0; d<3; ++d) {
    eqn_euler.projon_left_ev[d](&euler, q, vec, v1);
    eqn_euler.recwith_right_ev[d](&euler, q, v1, v2);
    for (int i=0; i<5; ++i)
      TEST_CHECK( gkyl_compare_double(v2[i], vec[i], 1e-14) );
  }
}

void
test_rescale_ev(void)
{
  struct euler_eqn euler = {.gas_gamma = 1.4};

  double v1[5], v2[5], rev[5*5];
  double vec[5] = {1.0, 20.0, 30.0, 40.0, 5.0};
  double q[5] = { 1.0, 1.5, 2.5, 3.5, 10.5 };

  struct eqn_sys eqn_euler = eqn_euler_init(&euler);  
  
  for (int d=0; d<3; ++d) {
    eqn_euler.projon_left_ev[d](&euler, q, vec, v1);
    eqn_euler.rescale_right_ev[d](&euler, q, v1, rev);
    
    for (int i=0; i<5; ++i) {
      v2[i] = 0.0;
      for (int j=0; j<5; ++j)
        v2[i] += rev[5*j+i];
    }
    
    for (int i=0; i<5; ++i)
      TEST_CHECK( gkyl_compare_double(v2[i], vec[i], 1e-12) );
  }
}


void
test_flux_jump_roe(void)
{
  double g = 1.4;
  struct euler_eqn euler = { .gas_gamma = g };
  struct eqn_sys eqn_euler = eqn_euler_init(&euler);

  double rhol = 1.0, uxl = 1.5, uyl = 2.5, uzl = 3.5, pl = 0.1;
  double rhor = 10.0, uxr = 10.5, uyr = 25.5, uzr = 35.5, pr = 100.0;
  
  double ql[5], qr[5];
  prim_to_cons(g, rhol, uxl, uyl, uzl, pl, ql);
  prim_to_cons(g, rhor, uxr, uyr, uzr, pr, qr);

  double dq[5];
  for (int i=0; i<5; ++i) dq[i] = qr[i]-ql[i];

  for (int d=0; d<3; ++d) {
    double apdq[5], amdq[5];    
    eqn_euler.fluct[d](&euler, ql, qr, apdq, amdq);

    double fr[5], fl[5];
    eqn_euler.flux[d](&euler, qr, fr);
    eqn_euler.flux[d](&euler, ql, fl);
    
    for (int i=0; i<5; ++i)
      TEST_CHECK( gkyl_compare_double(apdq[i]+amdq[i], fr[i]-fl[i], 1e-12) );
  }
}

TEST_LIST = {
  { "test_phi_prime", test_phi_prime },
  { "test_ev", test_ev },
  { "test_rescale_ev", test_rescale_ev },
  { "test_flux_jump_roe", test_flux_jump_roe },   
  { NULL, NULL },
};    
