#include "gkyl_array.h"
#include "gkyl_array_ops.h"
#include "gkyl_array_rio.h"
#include "gkyl_elem_type.h"
#include "gkyl_fv_proj.h"
#include "gkyl_range.h"
#include "gkyl_rect_grid.h"

#include <math.h>
#include <stdio.h>
#include <gkylzero.h>

#define SQ(x) ((x) * (x))

// parameters for source initialization
struct src_params {
  double xc, yc, sig;
};

// parameters for diffusion tensor initialization
struct D_param {
  // bx^2 + by^2 must be 1.0
  double (*bx)(double x, double y);
  double (*by)(double x, double y);
  double kpar, kperp;
};

// App input
struct app_inp {
  int cells[2]; // number of cells in each direction
  double lower[2], upper[2]; // lower/upper bounds of domain
  double tend; // end time for simulation

  struct src_params src_inp; // source parameters
  struct D_param D_inp; // input to setup Dij
};

// Ansio diffusion "app"
struct ad_app {
  // grid and for cell-centered quantities
  struct gkyl_rect_grid grid;
  struct gkyl_range range;

  // grid and for nodal quantities
  struct gkyl_rect_grid nc_grid;
  struct gkyl_range nc_range;

  // temperature (nodal)
  struct gkyl_array *T, *Tnew;
  // heat source (nodal)
  struct gkyl_array *S;

  // heat-flux vector (cell-centered)
  struct gkyl_array *Q;
  // Diffusion tensor (cell-centered)
  struct gkyl_array *Dij;
};

////// source initialization

void
src_func(double t, const double * restrict xn, double* restrict fout, void *ctx)
{
  // source is a gaussian
  struct src_params *sp = ctx;
  
  double xc = sp->xc, yc = sp->yc;
  double sig = sp->sig;
  double x = xn[0], y = xn[1];
  double r2 = SQ(x-xc) + SQ(y-yc);
  fout[0] = exp(-r2/(2.0*SQ(sig)));
}

void
init_source(struct gkyl_rect_grid grid, struct gkyl_range range,
  struct gkyl_array *S, struct src_params sc)
{
  int num_quad = 1; // ensures cell-center value
  gkyl_fv_proj *fv_proj = gkyl_fv_proj_new(&grid, 1, 1, src_func, &sc);
  gkyl_fv_proj_advance(fv_proj, 0.0, &range, S);
  gkyl_fv_proj_release(fv_proj);
}

////// Diffusion tensor initialization

// O-point field
double
op_bx(double x, double y)
{
  double x0 = 0.5, y0 = 0.5;
  return -(y-y0)/sqrt(SQ(x-x0)+SQ(y-y0));
}
double
op_by(double x, double y)
{
  double x0 = 0.5, y0 = 0.5;
  return (x-x0)/sqrt(SQ(x-x0)+SQ(y-y0));
}

void
D_func(double t, const double * restrict xn, double* restrict fout, void *ctx)
{
  struct D_param *dp = ctx;
  
  double kpar = dp->kpar, kperp = dp->kperp;
  double x = xn[0], y = xn[1];
  double bx = dp->bx(x, y), by = dp->by(x, y);
  
  fout[0] = kpar*bx*bx + (1-bx*bx)*kperp; // Dxx
  fout[1] = kpar*bx*by + (0-bx*by)*kperp; // Dxy
  fout[2] = kpar*by*by + (1-by*by)*kperp; // Dyy
}

void
init_diff(struct gkyl_rect_grid grid, struct gkyl_range range,
  struct gkyl_array *Dij, struct D_param dp)
{
  int num_quad = 1; // ensures cell-center value
  gkyl_fv_proj *fv_proj = gkyl_fv_proj_new(&grid, 1, 3, D_func, &dp);
  gkyl_fv_proj_advance(fv_proj, 0.0, &range, Dij);
  gkyl_fv_proj_release(fv_proj);
}


// Return new anisotropic-diffusion app object
struct ad_app *
ad_app_new(struct app_inp inp)
{
  struct ad_app *ad = gkyl_malloc(sizeof(struct ad_app));

  // create grid and range to on cell-centered grid  
  gkyl_rect_grid_init(&ad->grid, 2, inp.lower, inp.upper, inp.cells);
  gkyl_range_init_from_shape(&ad->range, ad->grid.ndim, ad->grid.cells);

  double dx = ad->grid.dx[0], dy = ad->grid.dx[0];

  // on 'nc_grid' the nodes of 'grid' are cell-centers: this is needed
  // as G0 and pgkyl assume that we only work with cell-centered
  // fields
  gkyl_rect_grid_init(&ad->nc_grid, 2,
    (double[]) { ad->grid.lower[0]-0.5*dx, ad->grid.lower[1]-0.5*dy },
    (double[]) { ad->grid.upper[0]+0.5*dx, ad->grid.upper[1]+0.5*dy },
    (int[]) { inp.cells[0]+1, inp.cells[1]+1 }
  );

  // create range to represent nodes
  gkyl_range_init_from_shape(&ad->nc_range, ad->nc_grid.ndim, ad->nc_grid.cells);

  // allocate temperature and source on nodes
  ad->T = gkyl_array_new(GKYL_DOUBLE, 1, ad->nc_range.volume);
  ad->Tnew = gkyl_array_new(GKYL_DOUBLE, 1, ad->nc_range.volume);
  ad->S = gkyl_array_new(GKYL_DOUBLE, 1, ad->nc_range.volume);  
  
  gkyl_array_clear(ad->T, 0.0); gkyl_array_clear(ad->Tnew, 0.0);
  
  // allocate heat-flux vector and diffusion tensor on cells
  ad->Q = gkyl_array_new(GKYL_DOUBLE, ad->grid.ndim, ad->range.volume);
  ad->Dij = gkyl_array_new(GKYL_DOUBLE, 3, ad->range.volume);

  return ad;
}

// initialize sources and diffusion tensor
void
ad_app_init(struct ad_app *ad, struct src_params sp, struct D_param dp)
{
  // initialize source
  init_source(ad->nc_grid, ad->nc_range, ad->S, sp);
  gkyl_grid_sub_array_write(&ad->nc_grid, &ad->nc_range, ad->S, "source.gkyl");

  // initialize diffusion tensor
  init_diff(ad->grid, ad->range, ad->Dij, dp);
  gkyl_grid_sub_array_write(&ad->grid, &ad->range, ad->Dij, "Dij.gkyl");
}

// compute heat-flux vector from temperature (and diffusion tensor)
void
ad_app_calc_heat_flux(struct ad_app *ad)
{
  enum { I, R, T, TR }; // (i,j), (i+1,j), (i,j+1), (i+1,j+1)
  enum { XX, XY, YY }; // components of Dij

  // compute offsets
  long offsets[4];
  offsets[I] = 0; // (i,j)
  offsets[R] = gkyl_range_offset(&ad->nc_range, (int[]) { 1,0 } ); // (i+1,j)
  offsets[T] = gkyl_range_offset(&ad->nc_range, (int[]) { 0,1 } ); // (i,j+1)
  offsets[TR] = gkyl_range_offset(&ad->nc_range, (int[]) { 1,1 } ); // (i+1,j+1)

  double dx = ad->grid.dx[0], dy = ad->grid.dx[1];

  // iterator over cells
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &ad->range);

  while (gkyl_range_iter_next(&iter)) {

    // fetch temperature at the four corners of cell
    long tidx = gkyl_range_idx(&ad->nc_range, iter.idx);
    
    const double *T_I = gkyl_array_cfetch(ad->T, tidx+offsets[I]);
    const double *T_R = gkyl_array_cfetch(ad->T, tidx+offsets[R]);
    const double *T_T = gkyl_array_cfetch(ad->T, tidx+offsets[T]);
    const double *T_TR = gkyl_array_cfetch(ad->T, tidx+offsets[TR]);

    // compute gradients
    double dx_T = (0.5*(T_TR[0]+T_R[0]) - 0.5*(T_T[0]+T_I[0]))/dx;
    double dy_T = (0.5*(T_TR[0]+T_T[0]) - 0.5*(T_R[0]+T_I[0]))/dy;

    // compute heat-flux
    long qidx = gkyl_range_idx(&ad->range, iter.idx);

    const double *D = gkyl_array_cfetch(ad->Dij, qidx);
    double *Q = gkyl_array_fetch(ad->Q, qidx);

    Q[0] = D[XX]*dx_T + D[XY]*dy_T;
    Q[1] = D[XY]*dx_T + D[YY]*dy_T;
  }
}

void
ad_app_write(struct ad_app *ad, int frame)
{
  const char *fmt = "T_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, frame);  
  gkyl_grid_sub_array_write(&ad->nc_grid, &ad->nc_range, ad->T, fileNm);
}

// free memory for app
void
ad_app_release(struct ad_app *ad)
{
  gkyl_array_release(ad->T);
  gkyl_array_release(ad->Tnew);
  gkyl_array_release(ad->S);
  gkyl_array_release(ad->Q);
  gkyl_array_release(ad->Dij);

  gkyl_free(ad);
}

// read app input (eventually from input file)
struct app_inp
get_app_inp(int argc, char *argv[])
{
  return (struct app_inp) {
    .cells = { 64, 64 },
    .lower = { 0, 0 },
    .upper = { 1.0, 1.0 },

    .tend = 1.0,

    .src_inp = {
      .xc = 0.25,
      .yc = 0.5,
      .sig = 0.05,
    },

    .D_inp = {
      .bx = op_bx,
      .by = op_by,
      .kpar = 1.0,
      .kperp = 1.0e-9
    }
  };
}

int
main(int argc, char *argv[])
{
  struct app_inp inp = get_app_inp(argc, argv);

  struct ad_app *ad = ad_app_new(inp);
  ad_app_init(ad, inp.src_inp, inp.D_inp);

  ad_app_release(ad);

  return 0;
}


/*
  // for testing initialize T
  gkyl_array_copy(ad->T, ad->S);
  ad_app_write(ad, 0);
  
  ad_app_calc_heat_flux(ad);
  gkyl_grid_sub_array_write(&ad->grid, &ad->range, ad->Q, "Q.gkyl");

*/
