#include <app.h>

#include <math.h>

void
init_f(double t, const double * restrict xn, double* restrict fout, void *ctx)
{
  fout[0] = sin(2*M_PI*xn[0])*cos(2*M_PI*xn[1]);
  fout[1] = sin(2*2*M_PI*xn[0])*cos(3*M_PI*xn[1]);
}

void
test_app_0(void)
{
  enum { FV_ARR_1, FV_ARR_2, DG_ARR_1, NARR };

  int ndim = 2;
  int poly_order = 2;
  struct gkyl_basis b1;
  gkyl_cart_modal_serendip(&b1, ndim, poly_order);  
  
  struct app_0_inp inp = {
    .ndim = ndim,
    .cells = { 100, 100 },

    .lower = { 0.0, 0.0 },
    .upper = { 1.0, 1.0 },
    
    .nghost = 2,

  };
  
  inp.narray = NARR;
  inp.arr_info[FV_ARR_1] = (struct app_0_array_inp) { .ncomp = 2, .ncoeff = 1 };
  inp.arr_info[FV_ARR_2] = (struct app_0_array_inp) { .ncomp = 2, .ncoeff = 1 };
  inp.arr_info[DG_ARR_1] = (struct app_0_array_inp) { .ncomp = 2, .ncoeff = b1.num_basis };

  app_0 *app = app_0_new(&inp);

  // FV
  app_0_fv_init(app, FV_ARR_1, 1.5, init_f, 0);
  app_0_write(app, &(struct app_0_write_inp) {
     .n = FV_ARR_1,
     .frame = 1,
     .tm = 0.0,
     .nm = "f_0.gkyl",
     .array_type = ARRAY_FV
    }
  );

  // DG
  app_0_dg_init(app, &b1, DG_ARR_1, 1.5, init_f, 0);
  app_0_write(app, &(struct app_0_write_inp) {
     .n = DG_ARR_1,
     .frame = 1,
     .tm = 0.0,
     .nm = "f_1.gkyl",
     .array_type = ARRAY_DG,
     .basis_type = b1.id,
     .poly_order = b1.poly_order
    }
  );

  app_0_release(app);
}

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

int
main(void)
{
  struct euler_ctx ctx = {
    .gas_gamma = 1.4
  };

  struct eqn_sys euler = {
    .meqn = 5,
    .mwave = 3,
    .eqn_ctx = &ctx
  };
  
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
