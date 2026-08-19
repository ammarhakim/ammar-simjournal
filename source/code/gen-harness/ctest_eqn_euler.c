#include <acutest.h>
#include <eqn_euler.h>

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

  euler_projon_left_ev_0(&euler, q, vec, v1);
  euler_recwith_right_ev_0(&euler, q, v1, v2);
  for (int i=0; i<5; ++i)
    TEST_CHECK( gkyl_compare_double(v2[i], vec[i], 1e-14) );

  euler_projon_left_ev_1(&euler, q, vec, v1);
  euler_recwith_right_ev_1(&euler, q, v1, v2);
  for (int i=0; i<5; ++i)
    TEST_CHECK( gkyl_compare_double(v2[i], vec[i], 1e-14) );

  euler_projon_left_ev_2(&euler, q, vec, v1);
  euler_recwith_right_ev_2(&euler, q, v1, v2);
  for (int i=0; i<5; ++i)
    TEST_CHECK( gkyl_compare_double(v2[i], vec[i], 1e-14) );
}

void
test_rescale_ev(void)
{
  struct euler_eqn euler = {.gas_gamma = 1.4};

  double v1[5], v2[5], rev[5*5];
  double vec[5] = {1.0, 20.0, 30.0, 40.0, 5.0};
  double q[5] = { 1.0, 1.5, 2.5, 3.5, 10.5 };

  do {
    euler_projon_left_ev_0(&euler, q, vec, v1);
    euler_rescale_right_ev_0(&euler, q, v1, rev);
    
    for (int i=0; i<5; ++i) {
      v2[i] = 0.0;
      for (int j=0; j<5; ++j)
        v2[i] += rev[5*j+i];
    }
    
    for (int i=0; i<5; ++i)
      TEST_CHECK( gkyl_compare_double(v2[i], vec[i], 1e-14) );
  } while (0);

  do {
    euler_projon_left_ev_1(&euler, q, vec, v1);
    euler_rescale_right_ev_1(&euler, q, v1, rev);
    
    for (int i=0; i<5; ++i) {
      v2[i] = 0.0;
      for (int j=0; j<5; ++j)
        v2[i] += rev[5*j+i];
    }
    
    for (int i=0; i<5; ++i)
      TEST_CHECK( gkyl_compare_double(v2[i], vec[i], 1e-10) );
  } while (0);

  do {
    euler_projon_left_ev_2(&euler, q, vec, v1);
    euler_rescale_right_ev_2(&euler, q, v1, rev);
    
    for (int i=0; i<5; ++i) {
      v2[i] = 0.0;
      for (int j=0; j<5; ++j)
        v2[i] += rev[5*j+i];
    }
    
    for (int i=0; i<5; ++i)
      TEST_CHECK( gkyl_compare_double(v2[i], vec[i], 1e-14) );
  } while (0);
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
test_flux_jump(void)
{
  double g = 1.4;
  struct euler_eqn euler = { .gas_gamma = g };

  double rhol = 1.0, uxl = 1.5, uyl = 2.5, uzl = 3.5, pl = 0.1;
  double rhor = 10.0, uxr = 10.5, uyr = 25.5, uzr = 35.5, pr = 100.0;
  
  double ql[5], qr[5];
  prim_to_cons(g, rhol, uxl, uyl, uzl, pl, ql);
  prim_to_cons(g, rhor, uxr, uyr, uzr, pr, qr);

  double dq[5];
  for (int i=0; i<5; ++i) dq[i] = qr[i]-ql[i];

  do {
    double fl[5], fr[5];
    euler_flux_0(&euler, ql, fl);
    euler_flux_0(&euler, qr, fr);

    double df[5];
    for (int i=0; i<5; ++i) df[i] = fr[i]-fl[i];

    double qavg[5];
    euler_ravg(&euler, ql, qr, qavg);

    double w[5];
    euler_projon_left_ev_0(&euler, qavg, dq, w);
    double wrev[5*5];
    euler_rescale_right_ev_0(&euler, qavg, w, wrev);
    
    double ev[5];
    euler_ev_0(&euler, qavg, ev);

    double apdq[5] = { 0.0, 0.0, 0.0, 0.0, 0.0};
    double amdq[5] = { 0.0, 0.0, 0.0, 0.0, 0.0};

    for (int i=0; i<5; ++i) {
      apdq[i] = 0.0;
      
      for (int j=0; j<5; ++j){
        apdq[i] += fmax(0.0, ev[j])*wrev[j*5+i];
        amdq[i] += fmin(0.0, ev[j])*wrev[j*5+i];
      }        
    }

    double adq[5];
    for (int i=0; i<5; ++i) adq[i] = apdq[i]+amdq[i];

    for (int i=0; i<5; ++i) 
      TEST_CHECK( gkyl_compare_double(adq[i], df[i], 1e-12) );
    
  } while (0);

  do {
    double fl[5], fr[5];
    euler_flux_1(&euler, ql, fl);
    euler_flux_1(&euler, qr, fr);

    double df[5];
    for (int i=0; i<5; ++i) df[i] = fr[i]-fl[i];

    double qavg[5];
    euler_ravg(&euler, ql, qr, qavg);

    double w[5];
    euler_projon_left_ev_1(&euler, qavg, dq, w);
    double wrev[5*5];
    euler_rescale_right_ev_1(&euler, qavg, w, wrev);
    
    double ev[5];
    euler_ev_1(&euler, qavg, ev);

    double apdq[5] = { 0.0, 0.0, 0.0, 0.0, 0.0};
    double amdq[5] = { 0.0, 0.0, 0.0, 0.0, 0.0};

    for (int i=0; i<5; ++i) {
      apdq[i] = 0.0;
      
      for (int j=0; j<5; ++j){
        apdq[i] += fmax(0.0, ev[j])*wrev[j*5+i];
        amdq[i] += fmin(0.0, ev[j])*wrev[j*5+i];
      }        
    }

    double adq[5];
    for (int i=0; i<5; ++i) adq[i] = apdq[i]+amdq[i];

    for (int i=0; i<5; ++i) 
      TEST_CHECK( gkyl_compare_double(adq[i], df[i], 1e-12) );
    
  } while (0);

  do {
    double fl[5], fr[5];
    euler_flux_2(&euler, ql, fl);
    euler_flux_2(&euler, qr, fr);

    double df[5];
    for (int i=0; i<5; ++i) df[i] = fr[i]-fl[i];

    double qavg[5];
    euler_ravg(&euler, ql, qr, qavg);

    double w[5];
    euler_projon_left_ev_2(&euler, qavg, dq, w);
    double wrev[5*5];
    euler_rescale_right_ev_2(&euler, qavg, w, wrev);
    
    double ev[5];
    euler_ev_2(&euler, qavg, ev);

    double apdq[5] = { 0.0, 0.0, 0.0, 0.0, 0.0};
    double amdq[5] = { 0.0, 0.0, 0.0, 0.0, 0.0};

    for (int i=0; i<5; ++i) {
      apdq[i] = 0.0;
      
      for (int j=0; j<5; ++j){
        apdq[i] += fmax(0.0, ev[j])*wrev[j*5+i];
        amdq[i] += fmin(0.0, ev[j])*wrev[j*5+i];
      }
    }

    double adq[5];
    for (int i=0; i<5; ++i) adq[i] = apdq[i]+amdq[i];

    for (int i=0; i<5; ++i) 
      TEST_CHECK( gkyl_compare_double(adq[i], df[i], 1e-12) );
    
  } while (0);  
  
}
  
TEST_LIST = {
  { "test_phi_prime", test_phi_prime },
  { "test_ev", test_ev },
  { "test_rescale_ev", test_rescale_ev },
  { "test_flux_jump", test_flux_jump },  
  { NULL, NULL },
};    
