#include <adiff.h>

// std includes
#include <assert.h>
#include <float.h>
#include <limits.h>
#include <string.h>

// Different numerical flux functions for advection

// First-order upwind
static inline double
calc_flux_upwind_1o(double vel, double fll, double fl, double fr, double frr)
{
  return 0.5*vel*(fr+fl) - 0.5*fabs(vel)*(fr-fl);
}

// Second-order central
static inline double
calc_flux_central_2o(double vel, double fll, double fl, double fr, double frr)
{
  return 0.5*vel*(fr+fl);
}

// Third-order upwind
static inline double
calc_flux_upwind_3o(double vel, double fll, double fl, double fr, double frr)
{
  double f = 0;
  if (vel>0)
    f = vel*(3*fr + 6*fl - fll)/8.0;
  else
    f = vel*(3*fl + 6*fr - frr)/8.0;
  
  return f;
}

// Different stencils for second-order derivatives: the cell-spacing
// is not included

// Second-order central
static inline double
calc_diff2_central_2o(double dx, double fll, double fl, double f0, double fr, double frr)
{
  return (fl-2.0*f0+fr)/(dx*dx);
}

// Fourth-order central
static inline double
calc_diff2_central_4o(double dx, double fll, double fl, double f0, double fr, double frr)
{
  return 0.0;
}

// function pointer types for computing advective interface flux
typedef double (*advect_flux_t)(double vel, double fll, double fl, double fr, double frr);

// function pointer types for computing diffusive term
typedef double (*diff_stencil_t)(double dx, double fll, double fl, double f0, double fr, double frr);

// Update status
struct update_status {
  bool success; // status of update
  double dt_actual; // actual time-step taken
  double dt_suggested; // suggested stable time-step
};

struct adiff_app {
  char name[128];

  int max_steps;  
  int nframe;
  double tend;
  double cfl;

  int ndim;
  struct gkyl_rect_grid grid;
  struct gkyl_range local, local_ext;
  struct gkyl_range face_range;

  enum adiff_advection_scheme scheme;  
  
  void *init_ctx; // context object for ICs
  evalf_t init_func; // initial conditions

  void *vel_ctx; // context object for advection velocity
  evalf_t vel_func; // function for advection velocity
  
  double alpha; // diffusion coefficient

  struct gkyl_array *f, *f1, *f2, *fnew;
  struct gkyl_array *vel;
  struct gkyl_array *cfl_freq;

  gkyl_dynvec fint; // integrated f, and integrated f2

  struct gkyl_comm *comm; // parallel communicator  
};

adiff_app *
adiff_app_new(const struct adiff_app_inp *inp)
{
  adiff_app *app = gkyl_malloc(sizeof *app);
  strncpy(app->name, inp->name, sizeof(inp->name));
  
  app->nframe = inp->nframe;
  app->tend = inp->tend;
  int ndim = app->ndim = inp->ndim > 0 ? inp->ndim : 2;

  app->max_steps = inp->max_steps > 0 ? inp->max_steps : INT_MAX;
  
  double cfl_frac = inp->cfl_frac > 0 ? inp->cfl_frac : 1.0;
  app->cfl = cfl_frac/ndim;
  
  // create grid and range to on cell-centered grid
  gkyl_rect_grid_init(&app->grid, ndim, inp->lower, inp->upper, inp->cells);
  int nghost[GKYL_MAX_CDIM] = { 2, 2, 2 };
  gkyl_create_grid_ranges(&app->grid, nghost, &app->local_ext, &app->local);

  // construct range for velocity
  int elo[] = { 0, 0, 0 };
  int eup[] = { 1, 1, 1 };
  struct gkyl_range vel_range;
  gkyl_range_extend(&vel_range, &app->local, elo, eup);
  gkyl_sub_range_intersect(&app->face_range, &app->local_ext, &vel_range);

  app->scheme = inp->scheme;

  app->init_ctx = inp->init_ctx;
  app->init_func = inp->init;

  app->vel_ctx = inp->velocity_ctx;
  app->vel_func = inp->velocity;

  app->f = gkyl_array_new(GKYL_DOUBLE, 1, app->local_ext.volume);
  app->f1 = gkyl_array_new(GKYL_DOUBLE, 1, app->local_ext.volume);
  app->fnew = gkyl_array_new(GKYL_DOUBLE, 1, app->local_ext.volume);
  app->f2 = gkyl_array_new(GKYL_DOUBLE, 1, app->local_ext.volume);
  app->vel = gkyl_array_new(GKYL_DOUBLE, ndim, app->local_ext.volume);
  app->cfl_freq = gkyl_array_new(GKYL_DOUBLE, 1, app->local_ext.volume);

  app->fint = gkyl_dynvec_new(GKYL_DOUBLE, 2);

  app->alpha = inp->alpha;

  // create communicator for periodic BCs (works only in serial)
  int cuts[GKYL_MAX_CDIM] = { 1, 1, 1 };
  struct gkyl_rect_decomp *decomp = gkyl_rect_decomp_new_from_cuts(ndim,
    cuts, &app->local);
  app->comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .decomp = decomp
    }
  );
  gkyl_rect_decomp_release(decomp);  

  return app;
}

static void
apply_bc(adiff_app *app, struct gkyl_array *arr)
{
  int per_dirs[] = { 0, 1, 2 };
  gkyl_comm_array_per_sync(app->comm, &app->local, &app->local_ext,
    app->ndim, per_dirs, arr);
}

static void
apply_ic(adiff_app *app)
{
  struct gkyl_offset_descr offsets[] = {
    { 0.0, 0.0, 0.0 } // f is at cell-center
  };
  
  struct gkyl_eval_offset_fd_inp inp = {
    .grid = &app->grid,
    .num_ret_vals = 1,
    .offsets = offsets,
    .eval = app->init_func,
    .ctx = app->init_ctx
  };

  gkyl_eval_offset_fd *fev = gkyl_eval_offset_fd_new(&inp);
  gkyl_eval_offset_fd_advance(fev, 0.0, &app->local, app->f);
  gkyl_eval_offset_fd_release(fev);
  apply_bc(app, app->f);
}

static void
calc_diagnostics(adiff_app *app, double tm)
{
  // integrated quantity and L2 norm
  double fint[2] = { 0.0, 0.0 };

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&app->local, iter.idx);

    const double *f_d = gkyl_array_cfetch(app->f, loc);
    fint[0] += f_d[0];
    fint[1] += f_d[0]*f_d[0]/2.0;
  }

  double vol = app->grid.cellVolume;
  fint[0] *= vol; fint[1] *= vol;

  gkyl_dynvec_append(app->fint, tm, fint);
}

static void
calc_velocity_field(adiff_app *app, double tm)
{
  struct gkyl_offset_descr offsets[] = {
    { -0.5, 0.0, 0.0 }, // ux is at left face
    { 0.0, -0.5, 0.0 }, // uy is at bottom face
    { 0.0, 0.0, -0.5 }, // ux is at back face
  };
  
  struct gkyl_eval_offset_fd_inp inp = {
    .grid = &app->grid,
    .num_ret_vals = app->ndim,
    .offsets = offsets,
    .eval = app->vel_func,
    .ctx = app->vel_ctx
  };

  gkyl_eval_offset_fd *vel_ev = gkyl_eval_offset_fd_new(&inp);
  gkyl_eval_offset_fd_advance(vel_ev, 0.0, &app->face_range, app->vel);
  gkyl_eval_offset_fd_release(vel_ev);
}

static void
write_field(adiff_app *app, double tm, int frame, const struct gkyl_array *arr, const char *fname)
{
  struct gkyl_msgpack_map_elem elist[] = {
    { .key = "time", .elem_type = GKYL_MP_DOUBLE, .dval = tm },
    { .key = "frame", .elem_type = GKYL_MP_INT, .ival = frame },
  };
  
  struct gkyl_msgpack_data *meta = gkyl_msgpack_create(2, elist); 
  gkyl_grid_sub_array_write(&app->grid, &app->local, meta, arr, fname);
  gkyl_msgpack_data_release(meta);  
}

static void
write_f(adiff_app *app, double tm, int frame)
{
  char fileNm[1024];
  const char *fmt = "%s_f_%d.gkyl";
  snprintf(fileNm, sizeof fileNm, fmt, app->name, frame);
  write_field(app, tm, frame, app->f, fileNm);
}

static void
write_vel(adiff_app *app, double tm, int frame)
{
  // NOTE: this method is not completely correct. As the velocity is
  // defined at cell-faces the below output will miss the upper layer
  // of velocity values. Could be fixed, but for now this is good
  // enough.

  char fileNm[1024];
  const char *fmt = "%s_vel_%d.gkyl";
  snprintf(fileNm, sizeof fileNm, fmt, app->name, frame);
  write_field(app, tm, frame, app->vel, fileNm);
}

static void
write_data(struct gkyl_tm_trigger* iot, adiff_app* app, double t_curr)
{
  if (gkyl_tm_trigger_check_and_bump(iot, t_curr)) {
    int frame = iot->curr - 1;
    write_f(app, t_curr, frame);
    if (frame == 0)
      write_vel(app, t_curr, frame);
  }
}

static void
write_diagnostics(adiff_app* app)
{
  char fileNm[1024];
  const char *fmt = "%s_diag.gkyl";
  snprintf(fileNm, sizeof fileNm, fmt, app->name);
  gkyl_dynvec_write(app->fint, fileNm);
}

static inline long
get_offset(int dir, int loc, const struct gkyl_range *range)
{
  int idx[GKYL_MAX_CDIM] = { 0, 0, 0 };
  idx[dir] = loc;
  return gkyl_range_offset(range, idx);
}

static double
calc_max_dt(adiff_app* app)
{
  double dx[3];
  for (int d=0; d<app->ndim; ++d)
    dx[d] = app->grid.dx[d];

  // compute maximum time-step due to diffusion term
  double alpha = app->alpha;
  double cfl_freq = 0.0;
  for (int d=0; d<app->ndim; ++d)
    cfl_freq += fmax(alpha, DBL_MIN)/(dx[d]*dx[d]);
  double max_dt_diff = 0.5*app->cfl/cfl_freq;

  gkyl_array_clear(app->cfl_freq, 0.0);
  // compute maximum time-step due to advection term
  for (int d=0; d<app->ndim; ++d) {
  
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &app->face_range);
    while (gkyl_range_iter_next(&iter)) {
      long loc = gkyl_range_idx(&app->local, iter.idx);
      const double *vel = gkyl_array_cfetch(app->vel, loc);

      double *cfl_freq_p = gkyl_array_fetch(app->cfl_freq, loc);
      cfl_freq_p[0] += fabs(vel[d])/dx[d];      
    }
  }
  double global_cfl_freq[1];
  gkyl_array_reduce_range(global_cfl_freq, app->cfl_freq, GKYL_MAX, &app->face_range);
  double max_dt_adv = app->cfl/global_cfl_freq[0];
  
  return fmin(max_dt_adv, max_dt_diff);
}

// Computes contribution to RHS from advection term, returning maximum
// stable time-step
static void
advection_rhs(adiff_app* app, const struct gkyl_array *fin,
  struct gkyl_array *rhs)
{
  int ndim = app->ndim;

  advect_flux_t flux_func = 0;

  switch (app->scheme) {
    case SCHEME_C2_C2:
      flux_func = calc_flux_central_2o;
      break;

    case SCHEME_U1_C2:
      flux_func = calc_flux_upwind_1o;
      break;

    case SCHEME_U3_C4:
      flux_func = calc_flux_upwind_3o;
      break;
  };
  
  gkyl_array_clear(app->cfl_freq, 0.0);
  
  double dx[3];
  for (int d=0; d<ndim; ++d) dx[d] = app->grid.dx[d];

  // labels for four cells attached to face: two on left, two on right
  enum { ILL, IL, IR, IRR };
  
  // outer loop is over directions
  for (int d=0; d<ndim; ++d) {
    long offsets[4];

    offsets[ILL] = get_offset(d, -2, &app->local);
    offsets[IL] = get_offset(d, -1, &app->local);
    offsets[IR] = get_offset(d, 0, &app->local);
    offsets[IRR] = get_offset(d, 1, &app->local);

    // inner loop is over faces: we update right and left cells
    // attached to face
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &app->face_range);
    while (gkyl_range_iter_next(&iter)) {
      long loc = gkyl_range_idx(&app->local, iter.idx);

      const double *vel = gkyl_array_cfetch(app->vel, loc+offsets[IR]);

      const double *fll = gkyl_array_cfetch(fin, loc+offsets[ILL]);
      const double *fl = gkyl_array_cfetch(fin, loc+offsets[IL]);
      const double *fr = gkyl_array_cfetch(fin, loc+offsets[IR]);
      const double *frr = gkyl_array_cfetch(fin, loc+offsets[IRR]);
      
      // compute flux
      double flux = flux_func(vel[d], fll[0], fl[0], fr[0], frr[0]);

      // we update two cells
      double *rhsl_p = gkyl_array_fetch(rhs, loc+offsets[IL]);
      double *rhsr_p = gkyl_array_fetch(rhs, loc+offsets[IR]);

      rhsl_p[0] += -flux/dx[d];
      rhsr_p[0] += flux/dx[d];
    }
  }
 }

// Computes contribution to RHS from diffusion term
static void
diffusion_rhs(adiff_app* app, const struct gkyl_array *fin,
  struct gkyl_array *rhs)
{
  int ndim = app->ndim;

  diff_stencil_t diff_stencil = 0;

  switch (app->scheme) {
    case SCHEME_C2_C2:
      diff_stencil = calc_diff2_central_2o;
      break;

    case SCHEME_U1_C2:
      diff_stencil = calc_diff2_central_2o;
      break;

    case SCHEME_U3_C4:
      diff_stencil = calc_diff2_central_4o;
      break;
  };  
  
  double alpha = app->alpha;
  double dx[3];
  for (int d=0; d<ndim; ++d) dx[d] = app->grid.dx[d];

  // labels for fice cells: two to the left, itself, two to the right
  enum { ILL, IL, I0, IR, IRR };
  
  // outer loop is over directions: the RHS is updated direction by
  // direction
  for (int d=0; d<ndim; ++d) {

    long offsets[5];
    offsets[ILL] = get_offset(d, -2, &app->local);
    offsets[IL] = get_offset(d, -1, &app->local);
    offsets[I0] = get_offset(d, 0, &app->local);
    offsets[IR] = get_offset(d, 1, &app->local);
    offsets[IRR] = get_offset(d, 2, &app->local);

    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &app->local);
    while (gkyl_range_iter_next(&iter)) {
      long loc = gkyl_range_idx(&app->local, iter.idx);

      double *rhs_p = gkyl_array_fetch(rhs, loc+offsets[I0]);

      const double *fll = gkyl_array_cfetch(fin, loc+offsets[ILL]);
      const double *fl = gkyl_array_cfetch(fin, loc+offsets[IL]);
      const double *f0 = gkyl_array_cfetch(fin, loc+offsets[I0]);
      const double *fr = gkyl_array_cfetch(fin, loc+offsets[IR]);
      const double *frr = gkyl_array_cfetch(fin, loc+offsets[IRR]);

      rhs_p[0] += alpha*diff_stencil(dx[d], fll[0], fl[0], f0[0], fr[0], frr[0]);
    }
  }
}

static void
forward_euler(adiff_app* app, double tcurr, double dt,
  const struct gkyl_array *fin, struct gkyl_array *fout,
  struct update_status *st)
{
  gkyl_array_clear_range(fout, 0.0, &app->local);
  // fout will have the RHS of advection and diffusion semi-discrete
  // operators accumulated
  advection_rhs(app, fin, fout);
  diffusion_rhs(app, fin, fout);

  // complete forward Euler update:
  // fout = fin + dt*fout
  gkyl_array_accumulate_range(gkyl_array_scale_range(fout, dt, &app->local),
    1.0, fin, &app->local);
  apply_bc(app, fout);

  st->dt_actual = st->dt_suggested = dt;
}

static struct update_status
update_ssp_rk3(adiff_app* app, double tcurr, double dt0)
{
  struct update_status st = { .success = true };

  // Stage 1:
  // f1 = FE[fn]
  forward_euler(app, tcurr, dt0, app->f, app->f1, &st);
  apply_bc(app, app->f1);

  // Stage 2:
  // f2 = 3/4*fn + 1/4*FE[f1]
  forward_euler(app, tcurr+dt0, dt0, app->f1, app->f2, &st);
  gkyl_array_accumulate_range(gkyl_array_scale_range(app->f2, 1.0/4.0, &app->local),
    3.0/4.0, app->f, &app->local);
  apply_bc(app, app->f2);

  // Stage 3:
  // fnew = 1/3*fn + 2/3*FE[f2]
  forward_euler(app, tcurr+0.5*dt0, dt0, app->f2, app->fnew, &st);
  gkyl_array_accumulate_range(gkyl_array_scale_range(app->fnew, 2.0/3.0, &app->local),
    1.0/3.0, app->f, &app->local);
  apply_bc(app, app->fnew);

  // copy solution over for use in next time-step
  gkyl_array_copy_range(app->f, app->fnew, &app->local_ext);  

  return st;
}

void
adiff_app_run(adiff_app *app)
{
  int num_frames = app->nframe;
  struct gkyl_tm_trigger io_trig = { .dt = app->tend/num_frames, .tcurr = 0.0 };

  apply_ic(app);
  calc_velocity_field(app, 0.0);
  calc_diagnostics(app, 0.0);
  write_data(&io_trig, app, 0.0);

  double tcurr = 0.0, tend = app->tend;
  double dt = calc_max_dt(app);
  
  long step = 1;
  while ((tcurr < tend) && (step <= app->max_steps) ) {
    fprintf(stdout, "Taking step %6ld at t = %#11.8g with dt %#11.8g ....\n", step, tcurr, dt);
    struct update_status status = update_ssp_rk3(app, tcurr,  dt);

    tcurr += status.dt_actual;
    calc_diagnostics(app, tcurr);
    write_data(&io_trig, app, tcurr);
    
    dt = status.dt_suggested;
    // readjust dt to hit tend exactly
    if (tcurr+dt > tend)
      dt = (tend-tcurr)*(1+1e-6);
    
    step += 1;
  }

  write_diagnostics(app);
}

void
adiff_app_release(adiff_app *app)
{
  gkyl_array_release(app->f);
  gkyl_array_release(app->f1);
  gkyl_array_release(app->fnew);
  gkyl_array_release(app->f2);
  gkyl_array_release(app->vel);
  gkyl_array_release(app->cfl_freq);

  gkyl_dynvec_release(app->fint);
  
  gkyl_comm_release(app->comm);
  gkyl_free(app);
}
