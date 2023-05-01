#include <gkylzero.h>

// multivector in 3D GA
struct mv_array {
  struct gkyl_array *m0; // scalar
  struct gkyl_array *m1; // vector
  struct gkyl_array *m2; // bivectors
  struct gkyl_array *ps; // pseudoscalar
};

struct mv_array*
mv_array_new(size_t sz)
{
  struct mv_array *mva = gkyl_malloc(sizeof(*mva));

  mva->m0 = gkyl_array_new(GKYL_DOUBLE, 1, sz);
  gkyl_array_clear(mva->m0, 0.0);

  mva->m1 = gkyl_array_new(GKYL_DOUBLE, 3, sz);
  gkyl_array_clear(mva->m1, 0.0);

  mva->m2 = gkyl_array_new(GKYL_DOUBLE, 3, sz);
  gkyl_array_clear(mva->m2, 0.0);

  mva->ps = gkyl_array_new(GKYL_DOUBLE, 1, sz);
  gkyl_array_clear(mva->ps, 0.0);

  return mva;
}

void
mv_array_release(struct mv_array *mva)
{
  gkyl_array_release(mva->m0);
  gkyl_array_release(mva->m1);
  gkyl_array_release(mva->m2);
  gkyl_array_release(mva->ps);

  gkyl_free(mva);
}

// labels for cells involved in FD stencil
enum offset_labels {
  I, // current cell
  MII, PII, // left/right
  IMI, IPI, // bottom/top
  IIM, IIP, // back/front
  OFF_END
};

// offsets for cells involved in FD stencil
static int offset_idx[][3] = {
  [I] = { 0, 0, 0 },
  [MII] = { -1, 0, 0 },
  [PII] = { 1, 0, 0 },
  [IMI] = { 0, -1, 0 },
  [IPI] = { 0, 1, 0 },
  [IIM] = { 0, 0, -1 },
  [IIP] = { 0, 0, 1 }
};

// compute cell offsets
void
calc_offsets(const struct gkyl_range *range, long offsets[])
{
  for (int i=0; i<OFF_END; ++i)
    offsets[i] = gkyl_range_offset(range, offset_idx[i]);
}

struct dgc_inp {
  int ndim; // dimension

  // lower, upper bounds and cells
  double lower[GKYL_MAX_CDIM];
  double upper[GKYL_MAX_CDIM];
  int cells[GKYL_MAX_CDIM];

  // initial condition and evaluation context
  evalf_t init_E, init_B;
  void *ctx;
};

// Discrete GC Maxwell solver app
struct dgc_app {
  struct gkyl_rect_grid grid;
  struct gkyl_range local, local_ext;

  struct mv_array *Fhalf_new, *Fhalf;
  struct mv_array *Ffull_new, *Ffull;

  evalf_t init_E, init_B;
  void *ctx;
};

struct dgc_app*
dgc_app_new(const struct dgc_inp *inp)
{
  struct dgc_app *app = gkyl_malloc(sizeof(*app));

  gkyl_rect_grid_init(&app->grid, inp->ndim, inp->lower, inp->upper, inp->cells);

  int nghost[GKYL_MAX_CDIM] = { 1, 1, 1 };
  gkyl_create_grid_ranges(&app->grid, nghost, &app->local_ext, &app->local);

  app->Ffull = mv_array_new(app->local_ext.volume);
  app->Ffull_new = mv_array_new(app->local_ext.volume);

  app->Fhalf = mv_array_new(app->local_ext.volume);
  app->Fhalf_new = mv_array_new(app->local_ext.volume);

  app->init_E = inp->init_E;
  app->init_B = inp->init_B;
  app->ctx = inp->ctx;

  return app;
}

void
init_fields(struct dgc_app *app)
{
  // Initialize electric field
  struct gkyl_offset_descr offset_E[3] = {
    { 0.0, -0.5, -0.5 },
    { -0.5, 0.0, -0.5 },
    { -0.5, -0.5, 0.0 }
  };
  gkyl_eval_offset_fd *ev_E = gkyl_eval_offset_fd_new(
    &(struct gkyl_eval_offset_fd_inp) {
      .grid = &app->grid,
      .num_ret_vals = 3,
      .eval = app->init_E,
      .ctx = app->ctx,
      .offsets = offset_E
    }
  );
  gkyl_eval_offset_fd_advance(ev_E, 0.0, &app->local_ext, app->Ffull->m1);
  gkyl_eval_offset_fd_release(ev_E);

  // Initialize magnetic field
  struct gkyl_offset_descr offset_B[3] = {
    { -0.5, 0.0, 0.0 },
    { 0.0, -0.5, 0.0 },
    { 0.0, 0.0, -0.5 }
  };
  gkyl_eval_offset_fd *ev_B = gkyl_eval_offset_fd_new(
    &(struct gkyl_eval_offset_fd_inp) {
      .grid = &app->grid,
      .eval = app->init_B,
      .num_ret_vals = 3,
      .ctx = app->ctx,
      .offsets = offset_B
    }
  );
  gkyl_eval_offset_fd_advance(ev_B, 0.0, &app->local_ext, app->Ffull->m2);
  gkyl_eval_offset_fd_release(ev_B);  
}

void
write_fields(struct dgc_app *app, int nframe)
{
  do {
    const char *fmt = "m0_%d.gkyl";
    int sz = gkyl_calc_strlen(fmt, nframe);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, nframe);

    gkyl_grid_sub_array_write(&app->grid, &app->local, app->Ffull->m0, fileNm);
  } while (0);

  do {
    const char *fmt = "m1_%d.gkyl";
    int sz = gkyl_calc_strlen(fmt, nframe);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, nframe);

    gkyl_grid_sub_array_write(&app->grid, &app->local, app->Ffull->m1, fileNm);
  } while (0);

  do {
    const char *fmt = "m2_%d.gkyl";
    int sz = gkyl_calc_strlen(fmt, nframe);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, nframe);

    gkyl_grid_sub_array_write(&app->grid, &app->local, app->Ffull->m2, fileNm);
  } while (0);

  do {
    const char *fmt = "ps_%d.gkyl";
    int sz = gkyl_calc_strlen(fmt, nframe);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, nframe);

    gkyl_grid_sub_array_write(&app->grid, &app->local, app->Ffull->ps, fileNm);
  } while (0);    

}

void
copy_fields(struct dgc_app *app, struct mv_array *fout, const struct mv_array *finp)
{
  gkyl_array_copy(fout->m0, finp->m0);
  gkyl_array_copy(fout->m1, finp->m1);
  gkyl_array_copy(fout->m2, finp->m2);
  gkyl_array_copy(fout->ps, finp->ps);
}

// compute gradient of multivector f and store output in gradf
void
calc_grad_f(struct dgc_app *app, const struct mv_array *f, struct mv_array *gradf)
{
  int ndim = app->grid.ndim;
  enum { IX1, IX2, IX3 }; // indexing

  double dx[3] = {
    app->grid.dx[0], app->grid.dx[1], app->grid.dx[3]
  };
  
  long offsets[OFF_END];
  calc_offsets(&app->local, offsets);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &app->local);
  
  while ( gkyl_range_iter_next(&iter) ) {
    long lidx = gkyl_range_idx(&app->local, iter.idx);

    const double *alpha_I = gkyl_array_cfetch(f->m0, lidx+offsets[I]);
    const double *a_I = gkyl_array_cfetch(f->m1, lidx+offsets[I]);    
    const double *b_I = gkyl_array_cfetch(f->m2, lidx+offsets[I]);
    const double *beta_I = gkyl_array_cfetch(f->ps, lidx+offsets[I]);
    
    // grad(f) grade-0
    double *gf0 = gkyl_array_fetch(gradf->m0, lidx);
    gf0[0] = 0.0;

    if (ndim > 0) {
      const double *a_MII = gkyl_array_cfetch(f->m1, lidx+offsets[MII]);
      gf0[0] +=  (a_I[0]-a_MII[0])/dx[0];
    }
    if (ndim > 1) {
      const double *a_IMI = gkyl_array_cfetch(f->m1, lidx+offsets[IMI]);
      gf0[0] +=  (a_I[1]-a_IMI[1])/dx[1];
    }
    if (ndim > 2) {
      const double *a_IIM = gkyl_array_cfetch(f->m1, lidx+offsets[IIM]);
      gf0[0] +=  (a_I[2]-a_IIM[2])/dx[2];
    }

    // grad(f) grade-1
    double *gf1 = gkyl_array_fetch(gradf->m1, lidx);
    gf1[0] = gf1[1] = gf1[2] = 0.0;

    if (ndim > 0) {
      const double *alpha_PII = gkyl_array_cfetch(f->m0, lidx+offsets[PII]);
      const double *b_MII = gkyl_array_cfetch(f->m2, lidx+offsets[MII]);
      
      gf1[IX1] += (alpha_PII[0]-alpha_I[0])/dx[0];
      
      gf1[IX3] += -(b_I[IX2]-b_MII[IX2])/dx[0];
      gf1[IX2] += (b_I[IX3]-b_MII[IX3])/dx[0];
    }
    if (ndim > 1) {
      const double *alpha_IPI = gkyl_array_cfetch(f->m0, lidx+offsets[IPI]);
      const double *b_IMI = gkyl_array_cfetch(f->m2, lidx+offsets[IMI]);
      
      gf1[IX2] +=  (alpha_IPI[0]-alpha_I[0])/dx[1];
      
      gf1[IX3] += (b_I[IX1]-b_IMI[IX1])/dx[1];
      gf1[IX1] += -(b_I[IX3]-b_IMI[IX3])/dx[1];

    }
    if (ndim > 2) {
      const double *alpha_IIP = gkyl_array_cfetch(f->m0, lidx+offsets[IIP]);
      const double *b_IIM = gkyl_array_cfetch(f->m2, lidx+offsets[IIM]);
      
      gf1[IX3] +=  (alpha_IIP[0]-alpha_I[0])/dx[2];

      gf1[IX1] += (b_I[IX2]-b_IIM[IX2])/dx[2];
      gf1[IX2] += -(b_I[IX1]-b_IIM[IX1])/dx[2];
    }

    // grad(f) grade-2
    double *gf2 = gkyl_array_fetch(gradf->m2, lidx);
    gf2[0] = gf2[1] = gf2[2] = 0.0;

    if (ndim > 0) {
      const double *beta_MII = gkyl_array_cfetch(f->ps, lidx+offsets[MII]);
      const double *a_PII = gkyl_array_cfetch(f->m1, lidx+offsets[PII]);
      
      gf2[IX1] += (beta_I[0]-beta_MII[0])/dx[0];
      
      gf2[IX3] += (a_PII[IX2]-a_I[IX2])/dx[0];
      gf2[IX2] += -(a_PII[IX3]-a_I[IX3])/dx[0];
    }
    if (ndim > 1) {
      const double *beta_IMI = gkyl_array_cfetch(f->ps, lidx+offsets[IMI]);
      const double *a_IPI = gkyl_array_cfetch(f->m1, lidx+offsets[IPI]);
      
      gf2[IX2] += (beta_I[0]-beta_IMI[0])/dx[1];

      gf2[IX1] += (a_IPI[IX3]-a_I[IX3])/dx[1];
      gf2[IX3] += -(a_IPI[IX1]-a_I[IX1])/dx[1];
      
    }
    if (ndim > 2) {
      const double *beta_IIM = gkyl_array_cfetch(f->ps, lidx+offsets[IIM]);
      const double *a_IIP = gkyl_array_cfetch(f->m1, lidx+offsets[IIP]);
      
      gf2[IX3] += (beta_I[0]-beta_IIM[0])/dx[2];

      gf2[IX2] += (a_IIP[IX1]-a_I[IX1])/dx[2];
      gf2[IX1] += -(a_IIP[IX2]-a_I[IX2])/dx[2];
    }    

    // grad(f) pseudoscalar
    double *gfps = gkyl_array_fetch(gradf->ps, lidx);
    gfps[0] = 0.0;

    if (ndim > 0) {
      const double *b_PII = gkyl_array_cfetch(f->m2, lidx+offsets[PII]);
      gfps[0] = (b_PII[0]-b_I[0])/dx[0];
    }
    if (ndim > 1) {
      const double *b_IPI = gkyl_array_cfetch(f->m2, lidx+offsets[IPI]);
      gfps[0] +=  (b_IPI[1]-b_I[1])/dx[1];
    }
    if (ndim > 2) {
      const double *b_IIP = gkyl_array_cfetch(f->m2, lidx+offsets[IIP]);
      gfps[0] +=  (b_IIP[2]-b_I[2])/dx[2];
    }
  }
}

void
dgc_app_release(struct dgc_app *app)
{
  mv_array_release(app->Fhalf);
  mv_array_release(app->Fhalf_new);

  mv_array_release(app->Ffull);
  mv_array_release(app->Ffull_new);
  
  gkyl_free(app);
}

void
eval_Efld(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  fout[0] = 0.0;
  fout[1] = 0.0;
  fout[2] = sin(2*M_PI*x)*sin(2*M_PI*y);
}

void
eval_Bfld(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0], y = xn[1];
  fout[0] = 2*M_PI*sin(2*M_PI*x)*cos(2*M_PI*y);
  fout[1] = -2*M_PI*cos(2*M_PI*x)*sin(2*M_PI*y);
  fout[2] = sin(2*M_PI*x)*sin(2*M_PI*y);
}

int
main(void)
{
  struct dgc_inp inp = {
    .ndim = 2,
    
    .lower = { 0.0, 0.0 },
    .upper = { 1.0, 1.0 },
    .cells = { 16, 16 },

    .init_E = eval_Efld,
    .init_B = eval_Bfld,
  };

  struct dgc_app *app = dgc_app_new( &inp );

  init_fields(app);
  write_fields(app, 0);

  calc_grad_f(app, app->Ffull, app->Ffull_new);
  copy_fields(app, app->Ffull, app->Ffull_new);
  write_fields(app, 1);
  
  dgc_app_release(app);
  
  return 0;
}
