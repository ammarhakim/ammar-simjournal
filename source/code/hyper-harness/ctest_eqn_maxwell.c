#include <acutest.h>
#include <eqn_maxwell.h>
#include <app.h>

struct eqn_sys
eqn_maxwell_init(void *ctx)
{
  return (struct eqn_sys) {
    .meqn = 8,
    .mtotal = 8,
    .eqn_ctx = &ctx,
    .flux = { maxwell_flux_0, maxwell_flux_1, maxwell_flux_2  },
    .projon_left_ev = { maxwell_projon_left_ev_0, maxwell_projon_left_ev_1, maxwell_projon_left_ev_2  },
    .recwith_right_ev = { maxwell_recwith_right_ev_0, maxwell_recwith_right_ev_1, maxwell_recwith_right_ev_2 },
    .rescale_right_ev = { maxwell_rescale_right_ev_0, maxwell_rescale_right_ev_1, maxwell_rescale_right_ev_2 },
    .ev = { maxwell_ev_0, maxwell_ev_1, maxwell_ev_2 },
    .fluct = { maxwell_fluct_0, maxwell_fluct_1, maxwell_fluct_2 },
    .rotate_to_local = maxwell_rotate_to_local,
    .rotate_to_global = maxwell_rotate_to_global,
  };
}

void
test_ev(void)
{
  struct maxwell_eqn maxwell = { .c = 1.0, .efact = 1.0, .bfact = 1.0 };
  int meqn = 8;

  double v1[8], v2[8];
  double vec[8] = {1.0, 20.0, 30.0, 40.0, 5.0, 1.0, 2.0, 3.0};
  double q[8] = { 1.0, 1.5, 2.5, 3.5, 10.5, 1.0, 1.5, 2.5 };

  struct eqn_sys eqn_maxwell = eqn_maxwell_init(&maxwell);  
  
  for (int d=0; d<3; ++d) {
    eqn_maxwell.projon_left_ev[d](&maxwell, q, vec, v1);
    eqn_maxwell.recwith_right_ev[d](&maxwell, q, v1, v2);
    for (int i=0; i<8; ++i)
      TEST_CHECK( gkyl_compare_double(v2[i], vec[i], 1e-14) );
  }
}

void
test_rescale_ev(void)
{
  struct maxwell_eqn maxwell = { .c = 1.0, .efact = 1.0, .bfact = 1.0 };
  int meqn = 8;  

  double v1[8], v2[8], rev[8*8];
  double vec[8] = {1.0, 20.0, 30.0, 40.0, 5.0, 1.0, 2.0, 3.0};
  double q[8] = { 1.0, 1.5, 2.5, 3.5, 10.5, 1.0, 1.5, 2.5 };

  struct eqn_sys eqn_maxwell = eqn_maxwell_init(&maxwell);  
  
  for (int d=0; d<3; ++d) {
    eqn_maxwell.projon_left_ev[d](&maxwell, q, vec, v1);
    eqn_maxwell.rescale_right_ev[d](&maxwell, q, v1, rev);
    
    for (int i=0; i<8; ++i) {
      v2[i] = 0.0;
      for (int j=0; j<8; ++j)
        v2[i] += rev[8*j+i];
    }
    
    for (int i=0; i<8; ++i)
      TEST_CHECK( gkyl_compare_double(v2[i], vec[i], 1e-12) );
  }
}

void
test_flux_jump_roe(void)
{
  struct maxwell_eqn maxwell = { .c = 1.0, .efact = 1.0, .bfact = 1.0 };
  struct eqn_sys eqn_maxwell = eqn_maxwell_init(&maxwell);

  int meqn = 8;  

  double rhol = 1.0, uxl = 1.5, uyl = 2.5, uzl = 3.5, pl = 0.1;
  double rhor = 10.0, uxr = 10.5, uyr = 25.5, uzr = 35.5, pr = 100.0;

  double ql[8] = { 1.0, 1.5, 2.5, 3.5, 10.5, 1.0, 1.5, 2.5 };
  double qr[8] = { 2.0, 3.5, -5.5, 30.5, 1.5, 10.0, 15.5, 32.5 };

  double dq[8];
  for (int i=0; i<8; ++i) dq[i] = qr[i]-ql[i];

  for (int d=0; d<3; ++d) {
    double apdq[8], amdq[8];    
    eqn_maxwell.fluct[d](&maxwell, ql, qr, apdq, amdq);

    double fr[8], fl[8];
    eqn_maxwell.flux[d](&maxwell, qr, fr);
    eqn_maxwell.flux[d](&maxwell, ql, fl);
    
    for (int i=0; i<8; ++i)
      TEST_CHECK( gkyl_compare_double(apdq[i]+amdq[i], fr[i]-fl[i], 1e-12) );
  }
}

void
test_rotate(void)
{
  struct maxwell_eqn maxwell = { .c = 1.0, .efact = 1.0, .bfact = 1.0 };
  struct eqn_sys eqn_maxwell = eqn_maxwell_init(&maxwell);
  
  double n[3] = { 0, 1.0, 0.0 };
  double tau1[3] = { -1.0, 0.0, 0.0 };
  double tau2[3] = { 0.0, 0.0, 1.0 };

  double qglobal[8] = { 1.0, 1.5, 2.5, 3.5, 10.5, 1.0, 1.5, 2.5 };  
 
  double qlocal[8];
  eqn_maxwell.rotate_to_local(&maxwell, n, tau1, tau2, qglobal, qlocal);
  double flocal[8];
  eqn_maxwell.flux[0](&maxwell, qlocal, flocal);
  double fglobal[8];
  eqn_maxwell.rotate_to_global(&maxwell, n, tau1, tau2, flocal, fglobal);

  double fy[8];
  eqn_maxwell.flux[1](&maxwell, qglobal, fy);

  for (int i=0; i<8; ++i)
    TEST_CHECK( gkyl_compare_double( fglobal[i], fy[i], 1e-14) );
}

TEST_LIST = {
  { "test_ev", test_ev },
  { "test_rescale_ev", test_rescale_ev },
  { "test_flux_jump_roe", test_flux_jump_roe },
  { "test_rotate", test_rotate },  
  { NULL, NULL },
};    
