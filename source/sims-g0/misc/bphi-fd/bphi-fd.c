#include <math.h>
#include <stdio.h>
#include <gkylzero.h>
#include <rt_arg_parse.h>

#define SQ(x) ((x) * (x))

// App input
struct app_inp {
  int cells[2]; // number of cells in each direction
  double lower[2], upper[2]; // lower/upper bounds of domain

  int nframe; // number of data-frames to write
  double tend; // end time for simulation

  double eta; // diffusion constant
  double jz, jrad; // current, radius of current
};

struct init_ctx {
  double jz, jrad;
};

// SImple Bphi app
struct lm_app {
  // grid and ranges for cell-centered quantities
  struct gkyl_rect_grid grid;
  struct gkyl_range ext_range, range;

  // geometry
  struct gkyl_wave_geom *geom;
  // boundary conditions
  struct gkyl_wv_apply_bc *bc_lower[2], *bc_upper[2];

  // temperature and RHS of induction equation (nodal)
  struct gkyl_array *Bphi, *Bphi_new, *rhs;
  struct gkyl_array *J;
};

// boundary condition functions
void
bc_axis(double t, int ncomp, const double *skin, double *ghost, void *ctx)
{
  ghost[0] = -skin[0];
}
void
bc_copy(double t, int ncomp, const double *skin, double *ghost, void *ctx)
{
  ghost[0] = skin[0];
}

// initial condition on Bphi
void
init_Bphi(double t, const double *xn, double *fout, void *ctx)
{
  struct init_ctx *ic = ctx;
  double r = xn[0];  

  double jrad = ic->jrad, jz = ic->jz;
  double Bphi = r<jrad ? jz*r/2 : jz*jrad*jrad/2/r;
  fout[0] = Bphi;
}

// apply boundary condition
void
lm_app_apply_bc(struct lm_app *ad, struct gkyl_array *Bphi)
{
  gkyl_wv_apply_bc_advance(ad->bc_lower[0], 0.0, &ad->range, Bphi);
  gkyl_wv_apply_bc_advance(ad->bc_upper[0], 0.0, &ad->range, Bphi);
  gkyl_wv_apply_bc_advance(ad->bc_lower[1], 0.0, &ad->range, Bphi);
}

// Return new anisotropic-diffusion app object
struct lm_app *
lm_app_new(struct app_inp inp)
{
  struct lm_app *ad = gkyl_malloc(sizeof(struct lm_app));

  // create grid and range to on cell-centered grid  
  gkyl_rect_grid_init(&ad->grid, 2, inp.lower, inp.upper, inp.cells);

  int nghost[GKYL_MAX_DIM] = { 1, 1 };
  gkyl_create_grid_ranges(&ad->grid, nghost, &ad->ext_range, &ad->range);

  // compute geometry
  ad->geom = gkyl_wave_geom_new(&ad->grid, &ad->ext_range, 0, 0);

  // using brugers equation object just to get hold of rotations
  struct gkyl_wv_eqn *eqn = gkyl_wv_burgers_new();
  
  // left wall is axis
  ad->bc_lower[0] = gkyl_wv_apply_bc_new(&ad->grid, eqn, ad->geom,
    0, GKYL_LOWER_EDGE, nghost, bc_axis, 0);
  // bottom wall is copy
  ad->bc_lower[1] = gkyl_wv_apply_bc_new(&ad->grid, eqn, ad->geom,
    1, GKYL_LOWER_EDGE, nghost, bc_copy, 0);
  
  // right wall is copy
  ad->bc_upper[0] = gkyl_wv_apply_bc_new(&ad->grid, eqn, ad->geom,
    0, GKYL_UPPER_EDGE, nghost, bc_copy, 0);

  ad->Bphi = gkyl_array_new(GKYL_DOUBLE, 1, ad->ext_range.volume);
  ad->Bphi_new = gkyl_array_new(GKYL_DOUBLE, 1, ad->ext_range.volume);
  ad->rhs = gkyl_array_new(GKYL_DOUBLE, 1, ad->ext_range.volume);

  ad->J = gkyl_array_new(GKYL_DOUBLE, 2, ad->ext_range.volume);

  // we are setting the initial field to that from a wire, and then ...
  struct init_ctx ctx = { .jrad = inp.jrad, .jz = inp.jz };
  gkyl_fv_proj *ic = gkyl_fv_proj_new(&ad->grid, 2, 1, init_Bphi, &ctx);
  gkyl_fv_proj_advance(ic, 0.0, &ad->ext_range, ad->Bphi);

  gkyl_array_copy(ad->Bphi_new, ad->Bphi);

  // .. clearing out the interior so only ghost cells are set
  gkyl_array_clear_range(ad->Bphi, 0.0, ad->range);
  gkyl_array_clear_range(ad->Bphi_new, 0.0, ad->range);

  lm_app_apply_bc(ad, ad->Bphi);
  lm_app_apply_bc(ad, ad->Bphi_new);

  gkyl_fv_proj_release(ic);
  gkyl_wv_eqn_release(eqn);

  return ad;
}

// compute RHS of induction equation
void
lm_app_calc_rhs(struct lm_app *ad)
{
  enum { I, L, R, T, B }; // cells indexing

  // compute offsets
  long offsets[5];
  offsets[I] = 0; // (i,j)
  
  offsets[L] = gkyl_range_offset(&ad->range, (int[]) { -1,0 } ); // (i-1,j)
  offsets[R] = gkyl_range_offset(&ad->range, (int[]) { 1,0 } ); // (i+1,j)

  offsets[T] = gkyl_range_offset(&ad->range, (int[]) { 0,1 } ); // (i,j+1)
  offsets[B] = gkyl_range_offset(&ad->range, (int[]) { 0,-1 } ); // (i,j-1)

  double xc[GKYL_MAX_CDIM];
  double dr = ad->grid.dx[0], dz = ad->grid.dx[1];

  gkyl_array_clear(ad->rhs, 0.0);

  // iterator over cells
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &ad->range);

  while (gkyl_range_iter_next(&iter)) {

    gkyl_rect_grid_cell_center(&ad->grid, iter.idx, xc);
    double r = xc[0]; // radial coordinate of cell-center

    long bidx = gkyl_range_idx(&ad->range, iter.idx);

    const double *B_I = gkyl_array_cfetch(ad->Bphi, bidx+offsets[I]);
    const double *B_L = gkyl_array_cfetch(ad->Bphi, bidx+offsets[L]);
    const double *B_R = gkyl_array_cfetch(ad->Bphi, bidx+offsets[R]);
    const double *B_T = gkyl_array_cfetch(ad->Bphi, bidx+offsets[T]);
    const double *B_B = gkyl_array_cfetch(ad->Bphi, bidx+offsets[B]);

    double *rhs = gkyl_array_fetch(ad->rhs, bidx);

    rhs[0] = (B_R[0]-2*B_I[0]+B_L[0])/(dr*dr)
      + (B_R[0]-B_L[0])/(2*r*dr)
      + (B_T[0]-2*B_I[0]+B_B[0])/(dz*dz)
      - B_I[0]/(r*r)
      ;
  }
}

// compute RHS of induction equation
void
lm_app_calc_J(struct lm_app *ad)
{
  enum { I, L, R, T, B }; // cells indexing

  // compute offsets
  long offsets[5];
  offsets[I] = 0; // (i,j)
  
  offsets[L] = gkyl_range_offset(&ad->range, (int[]) { -1,0 } ); // (i-1,j)
  offsets[R] = gkyl_range_offset(&ad->range, (int[]) { 1,0 } ); // (i+1,j)

  offsets[T] = gkyl_range_offset(&ad->range, (int[]) { 0,1 } ); // (i,j+1)
  offsets[B] = gkyl_range_offset(&ad->range, (int[]) { 0,-1 } ); // (i,j-1)

  double xc[GKYL_MAX_CDIM];
  double dr = ad->grid.dx[0], dz = ad->grid.dx[1];

  gkyl_array_clear(ad->rhs, 0.0);

  // iterator over cells
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &ad->range);

  while (gkyl_range_iter_next(&iter)) {

    gkyl_rect_grid_cell_center(&ad->grid, iter.idx, xc);
    double r = xc[0]; // radial coordinate of cell-center

    long bidx = gkyl_range_idx(&ad->range, iter.idx);

    const double *B_I = gkyl_array_cfetch(ad->Bphi, bidx+offsets[I]);
    const double *B_L = gkyl_array_cfetch(ad->Bphi, bidx+offsets[L]);
    const double *B_R = gkyl_array_cfetch(ad->Bphi, bidx+offsets[R]);
    const double *B_T = gkyl_array_cfetch(ad->Bphi, bidx+offsets[T]);
    const double *B_B = gkyl_array_cfetch(ad->Bphi, bidx+offsets[B]);

    double *J = gkyl_array_fetch(ad->J, bidx);

    J[0] = -(B_T[0]-B_B[0])/(2*dz);
    J[1] = (B_R[0]-B_L[0])/(2*dr) + B_I[0]/r;
  }
}

void
lm_app_write(struct lm_app *ad, int frame)
{
  const char *fmt = "Bphi_%d.gkyl";
  int sz = gkyl_calc_strlen(fmt, frame);
  char fileNm[sz+1]; // ensures no buffer overflow
  snprintf(fileNm, sizeof fileNm, fmt, frame);  
  gkyl_grid_sub_array_write(&ad->grid, &ad->range, ad->Bphi, fileNm);

  // compute current
  lm_app_calc_J(ad);
  const char *jfmt = "J_%d.gkyl";
  sz = gkyl_calc_strlen(fmt, frame);
  char jfileNm[sz+1]; // ensures no buffer overflow
  snprintf(jfileNm, sizeof jfileNm, jfmt, frame);
  gkyl_grid_sub_array_write(&ad->grid, &ad->range, ad->J, jfileNm);
  
}

// free memory for app
void
lm_app_release(struct lm_app *ad)
{
  gkyl_array_release(ad->Bphi);
  gkyl_array_release(ad->Bphi_new);
  gkyl_array_release(ad->rhs);
  gkyl_array_release(ad->J);
  
  gkyl_wave_geom_release(ad->geom);
  
  gkyl_wv_apply_bc_release(ad->bc_lower[0]);
  gkyl_wv_apply_bc_release(ad->bc_lower[1]);
  gkyl_wv_apply_bc_release(ad->bc_upper[0]);

  gkyl_free(ad);
}

// read app input (eventually from input file)
struct app_inp
get_app_inp(void)
{
  return (struct app_inp) {
    .cells = { 40*5, 40 },
    .lower = { 0, 0 },
    .upper = { 5.0, 1.0 },

    .nframe = 50,
    .tend = 5.0,

    .jrad = 0.5,
    .jz = -1.0,
    .eta = 1.0,
  };
}

int
main(int argc, char *argv[])
{
  struct gkyl_app_args app_args = parse_app_args(argc, argv);

  if (app_args.trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }  
  
  struct app_inp inp = get_app_inp();

  struct lm_app *ad = lm_app_new(inp);

  double cfl = 0.25;
  double eta = inp.eta;

  // comput CFL frequency
  double cfl_freq = eta/SQ(ad->grid.dx[0]) + eta/SQ(ad->grid.dx[1]);
  double dt = cfl/cfl_freq;

  // create trigger for IO
  struct gkyl_tm_trigger io_trig = { .dt = inp.tend/(inp.nframe+1) };

  if (gkyl_tm_trigger_check_and_bump(&io_trig, 0.0))
    lm_app_write(ad, io_trig.curr-1);  

  // main loop
  double tcurr = 0.0, tend = inp.tend;
  long step = 1;
  while ((tcurr < tend) && (step <= app_args.num_steps)) {
    printf("Taking time-step %ld at t = %g with dt = %g ...\n", step, tcurr, dt);

    // compute rhs and do first-order Euler update
    lm_app_calc_rhs(ad);
    gkyl_array_accumulate_range(ad->Bphi, dt, ad->rhs, ad->range);
    lm_app_apply_bc(ad, ad->Bphi);

    tcurr += dt; 
    step += 1;

    // write data if needed
    if (gkyl_tm_trigger_check_and_bump(&io_trig, tcurr))
      lm_app_write(ad, io_trig.curr-1);

    dt = fmin(dt, (tend-tcurr)*(1+1e-6));
  }
  lm_app_write(ad, 1000);

  lm_app_release(ad);

  return 0;
}
