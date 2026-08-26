#include <app.h>

#include <gkyl_alloc.h>
#include <gkyl_array_rio.h>
#include <gkyl_fv_proj.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_rect_decomp.h>
#include <gkyl_array_ops.h>

#include <string.h>

//
// Base app
//

struct app_0 {
  int ndim;
  struct gkyl_rect_grid grid;
  struct gkyl_range local, local_ext;
  int narray;
  struct gkyl_array **f;
  struct app_0_array_inp *finfo;
  struct gkyl_ref_count ref_count;
};

static void
app_0_free(const struct gkyl_ref_count* rc)
{
  struct app_0 *app = container_of(rc, struct app_0, ref_count);
  
  for (int i=0; i<app->narray; ++i)
    gkyl_array_release(app->f[i]);
  gkyl_free(app->finfo);
  gkyl_free(app->f);
  gkyl_free(app);  
}

struct app_0 *
app_0_new(struct app_0_inp *inp)
{
  struct app_0 *app = gkyl_malloc(sizeof(*app));

  app->ndim = inp->ndim;
  gkyl_rect_grid_init(&app->grid, inp->ndim, inp->lower, inp->upper, inp->cells);

  int nghost[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<inp->ndim; ++d) nghost[d] = inp->nghost;
  gkyl_create_grid_ranges(&app->grid, nghost, &app->local_ext, &app->local);

  app->narray = inp->narray;
  app->f = gkyl_malloc(inp->narray*sizeof(struct gkyl_aray *));
  app->finfo = gkyl_malloc(sizeof(struct app_0_array_inp[inp->narray]));

  for (int i=0; i<inp->narray; ++i) {
    int nc = inp->arr_info[i].ncomp*inp->arr_info[i].ncoeff;
    app->f[i] = gkyl_array_new(GKYL_DOUBLE, nc, app->local_ext.volume);
    
    app->finfo[i].ncomp = inp->arr_info[i].ncomp;
    app->finfo[i].ncoeff = inp->arr_info[i].ncoeff;
  }

  app->ref_count = gkyl_ref_count_init(app_0_free);
  
  return app;
}

struct app_0 *
app_0_acquire(const struct app_0 *app)
{
  gkyl_ref_count_inc(&app->ref_count);
  return (struct app_0 *) app;
}

enum gkyl_array_rio_status
app_0_write(struct app_0 *app, struct app_0_write_inp *inp)
{
  int nelem = 0;
  struct gkyl_msgpack_map_elem elist[4];

  if (inp->array_type == ARRAY_FV) {
    elist[0] = GKYL_MSGPACK_MAP_ELEM("time", inp->tm);
    elist[1] = GKYL_MSGPACK_MAP_ELEM("frame", inp->frame);
    nelem = 2;
  }
  else if (inp->array_type == ARRAY_DG) {
    elist[0] = GKYL_MSGPACK_MAP_ELEM("time", inp->tm);
    elist[1] = GKYL_MSGPACK_MAP_ELEM("frame", inp->frame);
    elist[2] = GKYL_MSGPACK_MAP_ELEM("poly_order", inp->poly_order);
    elist[3] = GKYL_MSGPACK_MAP_ELEM("basis_type", inp->basis_type);
    nelem = 4;
  }

  struct gkyl_msgpack_data *meta = gkyl_msgpack_create(nelem, elist);
  enum gkyl_array_rio_status status =
    gkyl_grid_sub_array_write(&app->grid, &app->local, meta, app->f[inp->n], inp->nm);
  gkyl_msgpack_data_release(meta);

  return status;
}

void
app_0_fv_init(struct app_0 *app, int n, double tm, evalf_t init, void *ctx)
{
  int num_quad = 2;
  int nc = app->finfo[n].ncomp;
  gkyl_fv_proj *fv_proj = gkyl_fv_proj_new(&app->grid, num_quad, nc, init, ctx);
  gkyl_fv_proj_advance(fv_proj, tm, &app->local, app->f[n]);
  gkyl_fv_proj_release(fv_proj);
}

void
app_0_dg_init(struct app_0 *app, struct gkyl_basis *basis, int n, double tm, evalf_t init, void *ctx)
{
  int poly_order = basis->poly_order;
  int nc = app->finfo[n].ncomp;
  gkyl_proj_on_basis *pob = gkyl_proj_on_basis_new(&app->grid, basis,
    poly_order+1, nc, init, ctx);
  gkyl_proj_on_basis_advance(pob, tm, &app->local, app->f[n]);
  gkyl_proj_on_basis_release(pob);
}

void
app_0_release(struct app_0 *app)
{
  gkyl_ref_count_dec(&app->ref_count);
}

//
// Hyperbolic solver app
//

enum hyper_app_fields { HYPER_F0, HYPER_RHS, HYPER_LAST };

typedef void (*recovery_fn_t)(
  int meqn,
  const double *f3m, const double *f2m, const double *fm,
  const double *fp, const double *f2p, const double *f3p,
  double *outl, double *outr  
);

struct hyper_app {
  char name[128];
  struct app_0 *app0;
  enum hyper_scheme_type scheme_type;
  struct eqn_sys hyper_eqn;
  bool has_diffusive_eqn;
  struct eqn_sys diff_eqn;

  gkyl_mem_buff q_on_left_ev_buff, q_with_right_ev_buff;
  struct gkyl_array *ql, *qr;
  struct gkyl_array *apdq, *amdq;

  enum mp_recon recon_type;
  recovery_fn_t recon_fn;
  bool use_char_limiters;
};

static inline long
get_offset(int dir, int loc, const struct gkyl_range *range)
{
  int idx[GKYL_MAX_DIM] = { 0, 0, 0 };
  idx[dir] = loc;
  return gkyl_range_offset(range, idx);
}

static void
proj_on_left_ev(struct hyper_app *app, int dir, const double *qifc, const double *qin, double *qout)
{
  if (app->use_char_limiters)
    app->hyper_eqn.projon_left_ev[dir](app->hyper_eqn.eqn_ctx, qifc, qin, qout);
  else{ 
    for (int m=0; m<app->hyper_eqn.meqn; ++m)
      qout[m] = qin[m];
  }
}

static void
recwith_right_ev(struct hyper_app *app, int dir, const double *qifc, const double *qin, double *qout)
{
  if (app->use_char_limiters)
    app->hyper_eqn.recwith_right_ev[dir](app->hyper_eqn.eqn_ctx, qifc, qin, qout);
  else {
    for (int m=0; m<app->hyper_eqn.meqn; ++m)
      qout[m] = qin[m];
  }
}

// Each of the recovery methods below take 6 cells, three to the left
// and three to the right, and returns the values at the left/right of
// the interface. Note that depending on the scheme, some of the input
// values may be ignored.

static inline void
c2_recovery(int meqn,
  const double *f3m, const double *f2m, const double *fm,
  const double *fp, const double *f2p, const double *f3p,
  double *outl, double *outr)
{
  // c2 is symmetric 2nd order scheme, so outl and outr are same
  for (int m=0; m<meqn; ++m)
    outr[m] = outl[m] = fm[m]/2.0 + fp[m]/2.0;
}

static inline void
c4_recovery(int meqn,
  const double *f3m, const double *f2m, const double *fm,
  const double *fp, const double *f2p, const double *f3p,
  double *outl, double *outr)
{
  // c4 is symmetric 4th order scheme, so outl and outr are same
  for (int m=0; m<meqn; ++m)
    outr[m] = outl[m] = -f2m[m]/12.0 + 7.0*fm[m]/12.0 + 7.0*fp[m]/12.0 - f2p[m]/12.0;
}

static inline void
c6_recovery(int meqn,
  const double *f3m, const double *f2m, const double *fm,
  const double *fp, const double *f2p, const double *f3p,
  double *outl, double *outr)
{
  // c6 is symmetric 6th order scheme, so outl and outr are same
  for (int m=0; m<meqn; ++m)  
    outr[m] = outl[m] =
      37.0*fp[m]/60.0+37.0*fm[m]/60.0+f3p[m]/60.0+f3m[m]/60.0-2.0*f2p[m]/15.0-2.0*f2m[m]/15.0;
}

static inline void
u1_recovery(int meqn,
  const double *f3m, const double *f2m, const double *fm,
  const double *fp, const double *f2p, const double *f3p,
  double *outl, double *outr)
{
  // u1 is upwind-biased 1st order scheme
  for (int m=0; m<meqn; ++m) {
    outl[m] = fm[m];
    outr[m] = fp[m];
  }
}

static inline void
u3_recovery(int meqn,
  const double *f3m, const double *f2m, const double *fm,
  const double *fp, const double *f2p, const double *f3p,
  double *outl, double *outr)
{
  // u3 is upwind-biased 3rd order scheme
  for (int m=0; m<meqn; ++m) {
    outl[m] = -1.0/6.0*f2m[m] + 5.0/6.0*fm[m] + 1.0/3.0*fp[m];
    outr[m] = 1.0/3.0*fm[m] + 5.0/6.0*fp[m] - 1.0/6.0*f2p[m];
  }
}

static inline void
u5_recovery(int meqn,
  const double *f3m, const double *f2m, const double *fm,
  const double *fp, const double *f2p, const double *f3p,
  double *outl, double *outr)
{
  // u5 is upwind-biased 5th order scheme
  for (int m=0; m<meqn; ++m) {
    outl[m] = 1.0/30.0*f3m[m] - 13.0/60.0*f2m[m] + 47.0/60.0*fm[m] + 9.0/20.0*fp[m] - 1.0/20.0*f2p[m];
    outr[m] = -1.0/20.0*f2m[m] + 9.0/20.0*fm[m] + 47.0/60.0*fp[m] - 13.0/60.0*f2p[m] + 1.0/30.0*f3p[m];
  }
}

struct hyper_app *
hyper_app_new(struct hyper_app_inp *inp)
{
  struct hyper_app *app = gkyl_malloc(sizeof *app);
  strcpy(app->name, inp->name);

  // create base app
  struct app_0_inp inp0 = {
    .ndim = inp->ndim,
    .nghost = 3
  };
  
  for (int d=0; d<inp->ndim; ++d) {
    inp0.cells[d] = inp->cells[d];
    inp0.lower[d] = inp->lower[d];
    inp0.upper[d] = inp->upper[d];
  }

  inp0.narray = HYPER_LAST;
  for (int i=0; i<inp0.narray; ++i)
    inp0.arr_info[i] = (struct app_0_array_inp) { .ncomp = inp->hyper_eqn.mtotal, .ncoeff = 1 };

  app->app0 = app_0_new(&inp0);
  app_0_fv_init(app->app0, HYPER_F0, 0.0, inp->init, inp->init_ctx);

  do {
    const char *fmt = "%s-f_%d.gkyl";
    int sz = gkyl_calc_strlen(fmt, app->name, 0);
    char fileNm[sz+1]; // ensures no buffer overflow
    snprintf(fileNm, sizeof fileNm, fmt, app->name, 0);
    app_0_write(app->app0, &(struct app_0_write_inp) {
        .frame = 0,
        .tm = 0.0,
        .array_type = ARRAY_FV,
        .n = HYPER_F0,
        .nm = fileNm
      }
    );
    
  } while (0);

  app->hyper_eqn = inp->hyper_eqn;
  app->q_on_left_ev_buff = gkyl_mem_buff_new(sizeof(double[6*app->hyper_eqn.meqn]));
  app->q_with_right_ev_buff = gkyl_mem_buff_new(sizeof(double[6*app->hyper_eqn.meqn]));

  app->ql = gkyl_array_new(GKYL_DOUBLE, inp->hyper_eqn.meqn, app->app0->local_ext.volume);
  app->qr = gkyl_array_new(GKYL_DOUBLE, inp->hyper_eqn.meqn, app->app0->local_ext.volume);

  app->apdq = gkyl_array_new(GKYL_DOUBLE, inp->hyper_eqn.meqn, app->app0->local_ext.volume);
  app->amdq = gkyl_array_new(GKYL_DOUBLE, inp->hyper_eqn.meqn, app->app0->local_ext.volume);

  recovery_fn_t rfun[16];
  rfun[MP_C2] = c2_recovery;
  rfun[MP_C4] = c4_recovery;
  rfun[MP_C6] = c6_recovery;
  rfun[MP_U1] = u1_recovery;
  rfun[MP_U3] = u3_recovery;
  rfun[MP_U5] = u5_recovery;

  app->recon_type = inp->recon_type;
  app->recon_fn = rfun[inp->recon_type];
  app->use_char_limiters = inp->use_char_limiters;

  return app;
}

static inline double
minmod_2(double x, double y)
{
  if (x>0 && y>0)
    return fmin(x,y);
  if (x<0 && y<0)
    return fmax(x,y);
  return 0.0;
}

static inline double
minmod_4(double x, double y, double z, double w)
{
  if (x>0 && y>0 && z>0 && w>0)
    return fmin(fmin(x,y),fmin(z,w));
  if (x<0 && y<0 && z<0 && w<0)
    return fmax(fmax(x,y),fmax(z,w));
  return 0.0;
}

static inline double
median(double x, double y, double z)
{
  return x + minmod_2(y-x,z-x);
}

static inline double
min_3(double x, double y, double z)
{
  return fmin(x,fmin(y,z));
}

static inline double
max_3(double x, double y, double z)
{
  return fmax(x,fmax(y,z));
}

// MP limiter: See Eqns 3.44 - 3.57 of Peterson and Hammett SIAM
// J. Sci. Comput, vol 35, No 3 pp B576, 2013
static inline double
mp_limiter(double qe, double q2m, double q1m, double q0, double q1p, double q2p)
{
  double alpha = 4.0;
 // Suresh and Huynh recommend 1e-10, but that seems turns off the
 // limiter when the jumps are very tiny, leading to small-scale noise
 // in certain situations. Not sure what this should be. Greg Hammett
 // recommends eps = 0.0
  double eps = 0.0;
  // Eq 3.44
  double qmp = q0 + minmod_2(q1p-q0, alpha*(q0-q1m));

  //return median(qe, q0, qmp);

  // Eq 3.45
  if ((qe-q0)*(qe-qmp)<eps) return qe; // no need to apply limiter

  // Eq 3.46-3.48
  double d1m = q2m + q0 - 2*q1m;
  double d0 = q1m + q1p - 2*q0;
  double d1p = q0 + q2p - 2*q1p;
  // Eq 3.49 and 3.50
  double dp_m4 = minmod_4(4*d0-d1p, 4*d1p-d0, d0, d1p);
  double dm_m4 = minmod_4(4*d1m-d0, 4*d0-d1m, d1m, d0);

  // Eq 3.51-3.54
  double qul = q0 + alpha*(q0-q1m);
  double qavg = 0.5*(q0+q1p);
  double qmd = qavg - 0.5*dp_m4;
  double qlc = q0 + 0.5*(q0-q1m) + 4.0/3.0*dm_m4;

  // Eq 3.55 and 3.56
  double qmin = fmax(min_3(q0,q1p,qmd), min_3(q0,qul,qlc));
  double qmax = fmin(max_3(q0,q1p,qmd), max_3(q0,qul,qlc));

  // 3.57
  return median(qe, qmin, qmax);
}

static double
hyper_app_calc_rhs(hyper_app *app, const struct gkyl_array *qin, struct gkyl_array *rhs)
{
  int ndim = app->app0->ndim;
  int meqn = app->hyper_eqn.meqn;
  const struct gkyl_range *range = &app->app0->local;

  enum { I3M, I2M, IM, IP, I2P, I3P }; // interface is between IM and IP

  double *qproj_on_left[6], *qrec_with_right[6];
  do {
    double *ptr = (double *) gkyl_mem_buff_data(app->q_on_left_ev_buff);
    for (int i=0; i<6; ++i)
      qproj_on_left[i] = &ptr[meqn*i];

    ptr = (double *) gkyl_mem_buff_data(app->q_with_right_ev_buff);
    for (int i=0; i<6; ++i)
      qrec_with_right[i] = &ptr[meqn*i];
  } while(0);

  double amax = 0.0; // maximum absolute eigenvalue
  
  gkyl_array_clear_range(rhs, 0.0, range);
  for (int dir=0; dir<ndim; ++dir) {
    double dx = app->app0->grid.dx[dir];
    const double *qavg[6];

    long offsets[6];
    offsets[IP]  = get_offset(dir, 0, range);
    offsets[I2P] = get_offset(dir, 1, range);
    offsets[I3P] = get_offset(dir, 2, range);
    offsets[IM]  = get_offset(dir, -1, range);
    offsets[I2M] = get_offset(dir, -2, range);
    offsets[I3M] = get_offset(dir, -3, range);

    int upper[GKYL_MAX_DIM] = { 0 };
    for (int d=0; d<range->ndim; ++d) upper[d] = range->upper[d];
    upper[dir] += 1; // one extra edge than cells
    struct gkyl_range range_ext;
    gkyl_range_init(&range_ext, range->ndim, range->lower, upper);

    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &range_ext);
    while (gkyl_range_iter_next(&iter)) {
      long loc = gkyl_range_idx(range, iter.idx);

      // use arithmetic averages for use in characterisitic limiting
      double qifc[meqn];
      for (int m=0; m<meqn; ++m) qifc[m] = 0.5*(qavg[IM][m]+qavg[IP][m]);
      
      for (int i=0; i<6; ++i) {
        qavg[i] = gkyl_array_cfetch(qin, loc+offsets[i]);
        proj_on_left_ev(app, dir, qifc, qavg[i], qproj_on_left[i]);
      }

      double qul[meqn], qur[meqn];
      app->recon_fn(meqn, qproj_on_left[I3M], qproj_on_left[I2M], qproj_on_left[IM],
        qproj_on_left[IP], qproj_on_left[I2P], qproj_on_left[I3P],
        qul, qur);

      for (int m=0; m<meqn; ++m) {
        qur[m] = mp_limiter(qur[m],
          qavg[I3P][m], qavg[I2P][m], qavg[IP][m],
          qavg[IM][m], qavg[I2M][m]);
        
        qul[m] = mp_limiter(qul[m],
          qavg[I3M][m], qavg[I2M][m], qavg[IM][m],
          qavg[IP][m], qavg[I2P][m]);
      }

      // qle is left of edge (right edge of left cell), qre right of
      // edge (left edge of right cell)
      double *qle = gkyl_array_fetch(app->qr, loc+offsets[IM]);
      double *qre = gkyl_array_fetch(app->ql, loc+offsets[IP]);

      recwith_right_ev(app, dir, qifc, qul, qle);
      recwith_right_ev(app, dir, qifc, qur, qre);

      double *amdq_p = gkyl_array_fetch(app->amdq, loc+offsets[IM]);
      double *apdq_p = gkyl_array_fetch(app->apdq, loc+offsets[IP]);

      double a = app->hyper_eqn.fluct[dir](app->hyper_eqn.eqn_ctx, qle, qre, apdq_p, amdq_p);
      amax = fmax(amax, a);
    }

    double fle[meqn], fre[meqn];
    // Update RHS with contribution from flux jumps. Note this loop is
    // over interior cells
    gkyl_range_iter_init(&iter, range);
    while (gkyl_range_iter_next(&iter)) {
      long loc = gkyl_range_idx(range, iter.idx);

      const double *qle = gkyl_array_cfetch(app->ql, loc);
      const double *qre = gkyl_array_cfetch(app->qr, loc);
      
      app->hyper_eqn.flux[dir](app->hyper_eqn.eqn_ctx, qle, fle);
      app->hyper_eqn.flux[dir](app->hyper_eqn.eqn_ctx, qre, fre);

      const double *amdq_p = gkyl_array_cfetch(app->amdq, loc);
      const double *apdq_p = gkyl_array_cfetch(app->apdq, loc);

      double *rhs_p = gkyl_array_fetch(rhs, loc);
      for (int m=0; m<meqn; ++m)
        rhs_p[m] += -(fre[m]-fle[m])/dx - (apdq_p[m]+amdq_p[m])/dx;
      
    }
  }

  return amax;
}

struct gkyl_update_status
hyper_app_update(hyper_app *app, double dt)
{

  return (struct gkyl_update_status) {
    .dt_actual = dt,
    .dt_suggested = dt,
    .success = true
  };
}


void
hyper_app_release(struct hyper_app *app)
{
  app_0_release(app->app0);
  
  gkyl_mem_buff_release(app->q_on_left_ev_buff);
  gkyl_mem_buff_release(app->q_with_right_ev_buff);
  gkyl_array_release(app->ql);
  gkyl_array_release(app->qr);
  gkyl_array_release(app->apdq);
  gkyl_array_release(app->amdq);
  gkyl_free(app);
}
