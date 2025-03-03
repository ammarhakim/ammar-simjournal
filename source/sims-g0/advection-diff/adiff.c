#include <adiff.h>

// std includes
#include <string.h>

struct adiff_app {
  char name[128]; // name of simulation
  
  int nframe; // number of output frames to write
  double tend; // time to run simulation
  double cfl; // CFL condition

  int ndim;
  struct gkyl_rect_grid grid;
  struct gkyl_range local, local_ext;
  struct gkyl_range vel_range; // range for velocity (which are on nodes)

  void *init_ctx; // context object for ICs
  evalf_t init_func; // initial conditions

  bool is_flow_dynamic; // set to true if velocity is time-dependent
  void *vel_ctx; // context object for advection velocity
  evalf_t vel_func; // function for advection velocity
  
  double alpha; // diffusion coefficient

  // arrays needed in the update
  struct gkyl_array *f, *f1, *fnew, *rhs;
  struct gkyl_array *vel;

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
  gkyl_sub_range_intersect(&app->vel_range, &app->local_ext, &vel_range);

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
}

static void
calc_velocity_field(double tm, adiff_app *app)
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
  gkyl_eval_offset_fd_advance(vel_ev, 0.0, &app->vel_range, app->vel);
  gkyl_eval_offset_fd_release(vel_ev);
}

static void
write_f(adiff_app *app, double tm, int frame)
{
  struct gkyl_msgpack_map_elem elist[] = {
    { .key = "time", .elem_type = GKYL_MP_DOUBLE, .dval = tm },
    { .key = "frame", .elem_type = GKYL_MP_INT, .ival = frame },
  };
  
  struct gkyl_msgpack_data *meta = gkyl_msgpack_create(2, elist);

  char fileNm[1024];
  do {
    const char *fmt = "%s_f_%d.gkyl";
    snprintf(fileNm, sizeof fileNm, fmt, app->name, frame);
    gkyl_grid_sub_array_write(&app->grid, &app->local, meta, app->f, fileNm);
  } while (0);
  
  gkyl_msgpack_data_release(meta);
}

static void
write_vel(adiff_app *app, double tm, int frame)
{
  // NOTE: this method is not completely correct. As th velocity is
  // defined at cell-faces the below output will miss the upper layer
  // of velocity values. Could be fixed, but for now this is good
  // enough.
  
  struct gkyl_msgpack_map_elem elist[] = {
    { .key = "time", .elem_type = GKYL_MP_DOUBLE, .dval = tm },
    { .key = "frame", .elem_type = GKYL_MP_INT, .ival = frame },
  };
  
  struct gkyl_msgpack_data *meta = gkyl_msgpack_create(2, elist);

  char fileNm[1024];
  do {
    const char *fmt = "%s_vel_%d.gkyl";
    snprintf(fileNm, sizeof fileNm, fmt, app->name, frame);
    gkyl_grid_sub_array_write(&app->grid, &app->local, meta, app->vel, fileNm);
  } while (0);
  
  gkyl_msgpack_data_release(meta);
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

void
adiff_app_run(adiff_app *app)
{
  apply_ic(app);
  calc_velocity_field(0.0, app);
  
  int num_frames = app->nframe;
  struct gkyl_tm_trigger io_trig = { .dt = app->tend/num_frames, .tcurr = 0.0 };

  write_data(&io_trig, app, 0.0);
  
}

void
adiff_app_release(adiff_app *app)
{
  gkyl_array_release(app->f);
  gkyl_array_release(app->f1);
  gkyl_array_release(app->fnew);
  gkyl_array_release(app->rhs);
  gkyl_array_release(app->vel);
  
  gkyl_comm_release(app->comm);
  gkyl_free(app);
}
