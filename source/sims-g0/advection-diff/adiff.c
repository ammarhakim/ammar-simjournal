#include <adiff.h>

// std includes
#include <float.h>
#include <limits.h>
#include <string.h>

// Update status
struct update_status {
  bool success; // status of update
  double dt_actual; // actual time-step taken
  double dt_suggested; // suggested stable time-step
};

struct adiff_app {
  char name[128]; // name of simulation
  
  int nframe; // number of output frames to write
  double tend; // time to run simulation
  double cfl; // CFL condition

  int max_steps; // maxmimum steps to run (for debugging only)

  int ndim;
  struct gkyl_rect_grid grid;
  struct gkyl_range local, local_ext;
  struct gkyl_range face_range; // range for quantities faces

  enum adiff_advection_scheme scheme;  
  
  void *init_ctx; // context object for ICs
  evalf_t init_func; // initial conditions

  bool is_flow_dynamic; // set to true if velocity is time-dependent
  void *vel_ctx; // context object for advection velocity
  evalf_t vel_func; // function for advection velocity
  
  double alpha; // diffusion coefficient

  // arrays needed in the update
  struct gkyl_array *f, *f1, *fnew, *rhs;
  struct gkyl_array *vel;
  struct gkyl_array *cfl_freq;

  // integrated diagnostics
  gkyl_dynvec fint;

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

  app->is_flow_dynamic = inp->is_flow_dynamic;
  app->vel_ctx = inp->velocity_ctx;
  app->vel_func = inp->velocity;

  app->f = gkyl_array_new(GKYL_DOUBLE, 1, app->local_ext.volume);
  app->f1 = gkyl_array_new(GKYL_DOUBLE, 1, app->local_ext.volume);
  app->fnew = gkyl_array_new(GKYL_DOUBLE, 1, app->local_ext.volume);
  app->rhs = gkyl_array_new(GKYL_DOUBLE, 1, app->local_ext.volume);
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

static inline double
calc_flux_upwind(double vel, double fl, double fr)
{
  return 0.5*vel*(fr+fl) - 0.5*fabs(vel)*(fr-fl);
}

static inline double
calc_flux_central(double vel, double fl, double fr)
{
  return 0.5*vel*(fr+fl);
}

// Computes contribution to RHS from advection term, returning maximum
// stable time-step
static double
advection_rhs(adiff_app* app, const struct gkyl_array *fin,
  struct gkyl_array *rhs)
{
  int ndim = app->ndim;

  double (*flux_func)(double vel, double fl, double fr);
  if (app->scheme == ADIFF_CENTRAL_2)
    flux_func = calc_flux_central;
  else if (app->scheme == ADIFF_UPWIND_1)
    flux_func = calc_flux_upwind;

  gkyl_array_clear(app->cfl_freq, 0.0);
  
  double cfl_freq = 0.0;
  double dx[3];
  for (int d=0; d<ndim; ++d)
    dx[d] = app->grid.dx[d];

  // labels for two cells attached to face: left, right
  enum { IL, IR };
  
  // outer loop is over directions
  for (int d=0; d<ndim; ++d) {

    long offsets[3];
    offsets[IL] = get_offset(d, -1, &app->local);
    offsets[IR] = get_offset(d, 0, &app->local);

    // inner loop is over faces: we update right and left cells
    // attached to face
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &app->face_range);
    while (gkyl_range_iter_next(&iter)) {
      long loc = gkyl_range_idx(&app->local, iter.idx);

      const double *vel = gkyl_array_cfetch(app->vel, loc+offsets[IR]);
      const double *fl = gkyl_array_cfetch(fin, loc+offsets[IL]);
      const double *fr = gkyl_array_cfetch(fin, loc+offsets[IR]);
      // compute flux
      double flux = flux_func(vel[d], fl[0], fr[0]);

      // we update two cells
      double *rhsl_p = gkyl_array_fetch(rhs, loc+offsets[IL]);
      double *rhsr_p = gkyl_array_fetch(rhs, loc+offsets[IR]);

      rhsl_p[0] += -flux/dx[d];
      rhsr_p[0] += flux/dx[d];

      double *cfl_freq_p = gkyl_array_fetch(app->cfl_freq, loc+offsets[IR]);
      cfl_freq_p[0] += fabs(vel[d])/dx[d];
    }
  }

  double global_cfl_freq[1];
  gkyl_array_reduce_range(global_cfl_freq, app->cfl_freq, GKYL_MAX, &app->local);

  return app->cfl/global_cfl_freq[0];
  
  return DBL_MAX;
}

// Computes contribution to RHS from diffusion term, returning maximum
// stable time-step
static double
diffusion_rhs(adiff_app* app, const struct gkyl_array *fin,
  struct gkyl_array *rhs)
{
  int ndim = app->ndim;

  double alpha = app->alpha;
  double cfl_freq = 0.0;
  double dx[3];
  for (int d=0; d<ndim; ++d) {
    dx[d] = app->grid.dx[d];
    // not sure why there is a 2 here
    cfl_freq += 2.0*fmax(alpha, DBL_MIN)/(dx[d]*dx[d]);
  }

  // labels for three cells: left, center, right
  enum { IL, I0, IR };
  
  // outer loop is over directions: the RHS is updated direction by
  // direction
  for (int d=0; d<ndim; ++d) {

    long offsets[3];
    offsets[IL] = get_offset(d, -1, &app->local);
    offsets[I0] = get_offset(d, 0, &app->local);
    offsets[IR] = get_offset(d, 1, &app->local);

    double alpha_dx2 = alpha/(dx[d]*dx[d]);

    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &app->local);
    while (gkyl_range_iter_next(&iter)) {
      long loc = gkyl_range_idx(&app->local, iter.idx);

      double *rhs_p = gkyl_array_fetch(rhs, loc+offsets[I0]);
      
      const double *fl = gkyl_array_cfetch(fin, loc+offsets[IL]);
      const double *f0 = gkyl_array_cfetch(fin, loc+offsets[I0]);
      const double *fr = gkyl_array_cfetch(fin, loc+offsets[IR]);

      rhs_p[0] += alpha_dx2*(fl[0]-2*f0[0]+fr[0]);
    }
  }

  return app->cfl/cfl_freq;
}

// Take a forward Euler step with the suggested time-step dt. This may
// not be the actual time-step taken. However, the function will never
// take a time-step larger than dt even if it is allowed by
// stability. The actual time-step and dt_suggested are returned in
// the status object.
static void
forward_euler(adiff_app* app, double tcurr, double dt,
  const struct gkyl_array *fin, struct gkyl_array *fout,
  struct update_status *st)
{
  double dt1, dtmin = DBL_MAX;

  gkyl_array_clear_range(fout, 0.0, &app->local);
  
  dt1 = advection_rhs(app, fin, fout);
  dtmin = fmin(dt1, dtmin);

  dt1 = diffusion_rhs(app, fin, fout);
  dtmin = fmin(dt1, dtmin);

  double dt_max_rel_diff = 0.01;
  // check if dtmin is slightly smaller than dt. Use dt if it is
  // (avoids retaking steps if dt changes are very small).
  double dt_rel_diff = (dt-dtmin)/dt;
  if (dt_rel_diff > 0 && dt_rel_diff < dt_max_rel_diff)
    dtmin = dt;

  // don't take a time-step larger that input dt
  double dta = st->dt_actual = dt < dtmin ? dt : dtmin;
  st->dt_suggested = dtmin;

  // complete the forward Euler update
  gkyl_array_accumulate_range(gkyl_array_scale_range(fout, dta, &app->local),
    1.0, fin, &app->local);
  apply_bc(app, fout);
}

static struct update_status
update_ssp_rk3(adiff_app* app, double dt0)
{
  struct update_status st = { .success = true };
  
  return st;
}

static struct update_status
update_rk1(adiff_app* app, double tcurr, double dt0)
{
  struct update_status st = { .success = true };
  forward_euler(app, tcurr, dt0, app->f, app->fnew, &st);
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

  struct update_status (*stepper)(adiff_app* app, double tcurr, double dt0)
    = update_rk1;

  double tcurr = 0.0, tend = app->tend;
  double dt = tend-tcurr;

  long step = 1;
  while ((tcurr < tend) && (step <= app->max_steps) ) {
    fprintf(stdout, "Taking step %ld at t = %g ....\n", step, tcurr);
    struct update_status status = stepper(app, tcurr,  dt);

    tcurr += status.dt_actual;
    calc_diagnostics(app, tcurr);
    write_data(&io_trig, app, tcurr);
    
    dt = status.dt_suggested;
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
  gkyl_array_release(app->rhs);
  gkyl_array_release(app->vel);
  gkyl_array_release(app->cfl_freq);

  gkyl_dynvec_release(app->fint);
  
  gkyl_comm_release(app->comm);
  gkyl_free(app);
}
