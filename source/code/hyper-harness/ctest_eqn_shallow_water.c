#include <acutest.h>
#include <eqn_shallow_water.h>
#include <app.h>

struct eqn_sys
eqn_shallow_water_init(void *ctx)
{
  return (struct eqn_sys) {
    .meqn = 3, 
    .mtotal = 4, // bottom topo is stored in the q
    .eqn_ctx = &ctx,
    .flux = { shallow_water_flux_0, shallow_water_flux_1 },
    .projon_left_ev = { shallow_water_projon_left_ev_0, shallow_water_projon_left_ev_1 },
    .recwith_right_ev = { shallow_water_recwith_right_ev_0, shallow_water_recwith_right_ev_1 },
    .rescale_right_ev = { shallow_water_rescale_right_ev_0, shallow_water_rescale_right_ev_1 },
    .ev = { shallow_water_ev_0, shallow_water_ev_1 },
    .fluct = { shallow_water_roe_fluct_0, shallow_water_roe_fluct_1 }
  };
}

static void
prim_to_cons(double rho, double ux, double uy, double q[3])
{
  q[0] = rho;
  q[1] = rho*ux;
  q[2] = rho*uy;
}

static void
show_vec(const char *nm, const double *q)
{
  printf("* %s\n", nm);
  for (int i=0; i<3; ++i)
    printf(" %.10lg\n", q[i]);
  printf("**\n");
}

void
test_phi_prime(void)
{
  struct shallow_water_eqn shallow_water = {.grav = 10.0};

  double v1[3], v2[3];
  double vec[3] = {1.0, 20.0, 30.0};
  double q[3] = { 1.0, 1.5, 2.5 };

  shallow_water_mulby_phi_prime(&shallow_water, q, vec, v1);
  shallow_water_mulby_phi_prime_inv(&shallow_water, q, v1, v2);

  for (int i=0; i<3; ++i)
    TEST_CHECK( gkyl_compare_double(v2[i], vec[i], 1e-14) );
}

void
test_ev(void)
{
  struct shallow_water_eqn shallow_water = {.grav = 10.0};

  double v1[3], v2[3];
  double vec[3] = {1.0, 20.0, 30.0};
  double q[3] = { 1.0, 1.5, 2.5 };

  struct eqn_sys eqn_shallow_water = eqn_shallow_water_init(&shallow_water);  
  
  for (int d=0; d<2; ++d) {
    eqn_shallow_water.projon_left_ev[d](&shallow_water, q, vec, v1);
    eqn_shallow_water.recwith_right_ev[d](&shallow_water, q, v1, v2);
    for (int i=0; i<3; ++i)
      TEST_CHECK( gkyl_compare_double(v2[i], vec[i], 1e-14) );
  }
}

void
test_rescale_ev(void)
{
  struct shallow_water_eqn shallow_water = {.grav = 10.0};

  double v1[3], v2[3], rev[3*3];
  double vec[3] = {1.0, 20.0, 30.0 };
  double q[3] = { 1.0, 1.5, 2.5 };

  struct eqn_sys eqn_shallow_water = eqn_shallow_water_init(&shallow_water);  
  
  for (int d=0; d<2; ++d) {
    eqn_shallow_water.projon_left_ev[d](&shallow_water, q, vec, v1);
    eqn_shallow_water.rescale_right_ev[d](&shallow_water, q, v1, rev);
    
    for (int i=0; i<3; ++i) {
      v2[i] = 0.0;
      for (int j=0; j<3; ++j)
        v2[i] += rev[3*j+i];
    }
    
    for (int i=0; i<3; ++i)
      TEST_CHECK( gkyl_compare_double(v2[i], vec[i], 1e-12) );
  }
}


void
test_flux_jump_roe(void)
{
  struct shallow_water_eqn shallow_water = { .grav = 10.0 };
  struct eqn_sys eqn_shallow_water = eqn_shallow_water_init(&shallow_water);

  double rhol = 1.0, uxl = 1.5, uyl = 2.5;
  double rhor = 10.0, uxr = 10.5, uyr = 25.5;
  
  double ql[3], qr[3];
  prim_to_cons(rhol, uxl, uyl, ql);
  prim_to_cons(rhor, uxr, uyr, qr);

  double dq[3];
  for (int i=0; i<3; ++i) dq[i] = qr[i]-ql[i];

  for (int d=0; d<2; ++d) {
    double apdq[3], amdq[3];    
    eqn_shallow_water.fluct[d](&shallow_water, ql, qr, apdq, amdq);

    double fr[3], fl[3];
    eqn_shallow_water.flux[d](&shallow_water, qr, fr);
    eqn_shallow_water.flux[d](&shallow_water, ql, fl);
    
    for (int i=0; i<3; ++i)
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
