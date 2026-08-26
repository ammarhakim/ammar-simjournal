#include <acutest.h>
#include <eqn_burgers.h>
#include <app.h>

struct eqn_sys
eqn_burgers_init(void *ctx)
{
  return (struct eqn_sys) {
    .meqn = 1,
    .mtotal = 1,
    .eqn_ctx = &ctx,
    .flux = { burgers_flux_0, burgers_flux_1, burgers_flux_2  },
    .projon_left_ev = { burgers_projon_left_ev_0, burgers_projon_left_ev_1, burgers_projon_left_ev_2  },
    .recwith_right_ev = { burgers_recwith_right_ev_0, burgers_recwith_right_ev_1, burgers_recwith_right_ev_2 },
    .rescale_right_ev = { burgers_rescale_right_ev_0, burgers_rescale_right_ev_1, burgers_rescale_right_ev_2 },
    .ev = { burgers_ev_0, burgers_ev_1, burgers_ev_2 },
    .fluct = { burgers_fluct_0, burgers_fluct_1, burgers_fluct_2 }
  };
}

void
test_ev(void)
{
  double v1[1], v2[1];
  double vec[1] = { 1.0 };
  double q[1] = { 2.5 };

  struct eqn_sys eqn_burgers = eqn_burgers_init(0);  
  
  for (int d=0; d<3; ++d) {
    eqn_burgers.projon_left_ev[d](0, q, vec, v1);
    eqn_burgers.recwith_right_ev[d](0, q, v1, v2);
    for (int i=0; i<1; ++i) {
      TEST_CHECK( gkyl_compare_double(v2[i], vec[i], 1e-14) );
    }
  }
}

void
test_rescale_ev(void)
{
  double v1[1], v2[1], rev[1*1];
  double vec[1] = { 1.0 };
  double q[1] = { -2.5 };

  struct eqn_sys eqn_burgers = eqn_burgers_init(0);  
  
  for (int d=0; d<3; ++d) {
    eqn_burgers.projon_left_ev[d](0, q, vec, v1);
    eqn_burgers.rescale_right_ev[d](0, q, v1, rev);
    
    for (int i=0; i<1; ++i) {
      v2[i] = 0.0;
      for (int j=0; j<1; ++j)
        v2[i] += rev[1*j+i];
    }
    
    for (int i=0; i<1; ++i)
      TEST_CHECK( gkyl_compare_double(v2[i], vec[i], 1e-12) );
  }
}

void
test_flux_jump_roe(void)
{
  struct eqn_sys eqn_burgers = eqn_burgers_init(0);

  double rhol = 1.0, uxl = 1.5, uyl = 2.5, uzl = 3.5, pl = 0.1;
  double rhor = 10.0, uxr = 10.5, uyr = 25.5, uzr = 35.5, pr = 100.0;

  double ql[1] = { -1.0 };
  double qr[1] = { 3.5 };

  double dq[1];
  for (int i=0; i<1; ++i) dq[i] = qr[i]-ql[i];

  for (int d=0; d<3; ++d) {
    double apdq[1], amdq[1];    
    eqn_burgers.fluct[d](0, ql, qr, apdq, amdq);

    double fr[1], fl[1];
    eqn_burgers.flux[d](0, qr, fr);
    eqn_burgers.flux[d](0, ql, fl);
    
    for (int i=0; i<1; ++i)
      TEST_CHECK( gkyl_compare_double(apdq[i]+amdq[i], fr[i]-fl[i], 1e-12) );
  }
}

TEST_LIST = {
  { "test_ev", test_ev },
  { "test_rescale_ev", test_rescale_ev },
  { "test_flux_jump_roe", test_flux_jump_roe },
  { NULL, NULL },
};    
