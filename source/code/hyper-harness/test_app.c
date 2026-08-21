#include <app.h>
#include <eqn_euler.h>

#include <math.h>

struct euler_ctx {
  double gas_gamma;
};

void
init_euler_sod(double t, const double * restrict xn, double* restrict fout, void *ctx)
{
  struct euler_ctx *ectx = ctx;

  double x = xn[0];  
  double xloc = 0.5;
  double gas_gamma = ectx->gas_gamma;  

  double rhol = 3.0, ul = 0.0, pl = 3.0;
  double rhor = 1.0, ur = 0.0, pr = 1.0;

  if (x < xloc) {
    fout[0] = rhol;
    fout[1] = rhol*ul;
    fout[2] = fout[3] = 0.0;
    fout[4] = 0.5*rhol*ul*ul + pl/(gas_gamma - 1);
  }
  else {
    fout[0] = rhor;
    fout[1] = rhor*ur;
    fout[2] = fout[3] = 0.0;    
    fout[4] = 0.5*rhor*ur*ur + pr/(gas_gamma - 1); 
  }
}

static inline void
noop_projon_left_ev(void *ctx, const double *q, const double *vin, double *vout)
{
  for (int i=0; i<5; ++i) vout[i] = vin[i];
}

static inline void
noop_recwith_right_ev_0(void *ctx, const double *q, const double *vin, double *vout)
{
  for (int i=0; i<5; ++i) vout[i] = vin[i];
}

struct eqn_sys
eqn_euler_init(void *ctx)
{
  return (struct eqn_sys) {
    .meqn = 5,
    .mwave = 5,
    .eqn_ctx = &ctx,
    .flux = { euler_flux_0, euler_flux_1, euler_flux_2  },
    .recwith_right_ev = { euler_recwith_right_ev_0, euler_recwith_right_ev_1, euler_recwith_right_ev_2 },
    .rescale_right_ev = { euler_rescale_right_ev_0, euler_rescale_right_ev_1, euler_rescale_right_ev_2 },
    .ev = { euler_ev_0, euler_ev_1, euler_ev_2 },
    .fluct = { euler_roe_fluct_0, euler_roe_fluct_1, euler_roe_fluct_2 }
  };
}

int
main(void)
{
  struct euler_ctx ctx = {
    .gas_gamma = 1.4    
  };

  struct eqn_sys euler = eqn_euler_init(&ctx);
  
  struct hyper_app_inp inp = {
    .name = "euler_sod",
    .ndim = 1,

    .tend = 0.2,
    .nframe = 1,
    
    .cells = { 100 },
    .lower = { 0.0 },
    .upper = { 1.0 },

    .scheme_type = SCHEME_FD,
    .hyper_eqn = euler,

    .has_diffusive_eqn = false,

    .init = init_euler_sod,
    .init_ctx = &ctx
  };

  hyper_app *app = hyper_app_new(&inp);

  hyper_app_release(app);
}
