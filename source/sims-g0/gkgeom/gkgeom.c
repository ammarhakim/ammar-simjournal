#include <gkgeom.h>

#define SQ(x) ((x) * (x))

// Adaptive quadrature
double wp34s(double (*func)(double, void *), void *ctx, double a, double b, int n, double eps);
// Root finder
double ridders(double (*func)(double,void*), void *ctx, double x1, double x2, double f1, double f2, double eps);

struct RdRdZ_sol {
  int nsol;
  double R[2], dRdZ[2];
};

static inline struct RdRdZ_sol
calc_RdR_p1(const double *psi, double psi0, double Z, double xc[2], double dx[2])
{
  struct RdRdZ_sol sol = { .nsol = 0 };

  double y = (Z-xc[1])/(dx[1]*0.5);
  
  double rnorm = (-(1.732050807568877*psi[2]*y)/(3.0*psi[3]*y+1.732050807568877*psi[1]))+(2.0*psi0)/(3.0*psi[3]*y+1.732050807568877*psi[1])-(1.0*psi[0])/(3.0*psi[3]*y+1.732050807568877*psi[1]) ;

  if ((-1<=rnorm) && (rnorm < 1)) {
    double drdznorm = -(3.0*(2.0*psi[3]*psi0-1.0*psi[0]*psi[3]+psi[1]*psi[2]))/SQ(3.0*psi[3]*y+1.732050807568877*psi[1]) ;
    
    sol.nsol = 1;
    sol.R[0] = rnorm*dx[0]*0.5 + xc[0];
    sol.dRdZ[0] = drdznorm*dx[0]/dx[1];
  }
  return sol;
}

static inline struct RdRdZ_sol
calc_RdR_p2(const double *psi, double psi0, double Z, double xc[2], double dx[2])
{
  struct RdRdZ_sol sol = { .nsol = 0 };
  double y = (Z-xc[1])/(dx[1]*0.5);

  double aq = 2.904737509655563*psi[6]*y+1.677050983124842*psi[4]; 
  double bq = 2.904737509655563*psi[7]*SQ(y)+1.5*psi[3]*y-0.9682458365518543*psi[7]+0.8660254037844386*psi[1]; 
  double cq = 1.677050983124842*psi[5]*SQ(y)-0.9682458365518543*psi[6]*y+0.8660254037844386*psi[2]*y-1.0*psi0-0.5590169943749475*psi[5]-0.5590169943749475*psi[4]+0.5*psi[0]; 
  double delta2 = bq*bq - 4*aq*cq;

  if (delta2 > 0) {
    double r1, r2;
    double delta = sqrt(delta2);
    // compute both roots
    if (bq>=0) {
      r1 = (-bq-delta)/(2*aq);
      r2 = 2*cq/(-bq-delta);
    }
    else {
      r1 = 2*cq/(-bq+delta);
      r2 = (-bq+delta)/(2*aq);
    }

    int sidx = 0;
    if ((-1<=r1) && (r1 < 1)) {
      sol.nsol += 1;
      sol.R[sidx] = r1*dx[0]*0.5 + xc[0];

      double x = r1;
      double C = 5.809475019311126*psi[7]*x*y+3.354101966249685*psi[5]*y+2.904737509655563*psi[6]*SQ(x)+1.5*psi[3]*x-0.9682458365518543*psi[6]+0.8660254037844386*psi[2]; 
      double A = 2.904737509655563*psi[7]*SQ(y)+5.809475019311126*psi[6]*x*y+1.5*psi[3]*y+3.354101966249685*psi[4]*x-0.9682458365518543*psi[7]+0.8660254037844386*psi[1];
      sol.dRdZ[sidx] = -C/A*dx[0]/dx[1];
      
      sidx += 1;
    }
    if ((-1<=r2) && (r2 < 1)) {
      sol.nsol += 1;
      sol.R[sidx] = r2*dx[0]*0.5 + xc[0];

      double x = r2;
      double C = 5.809475019311126*psi[7]*x*y+3.354101966249685*psi[5]*y+2.904737509655563*psi[6]*SQ(x)+1.5*psi[3]*x-0.9682458365518543*psi[6]+0.8660254037844386*psi[2]; 
      double A = 2.904737509655563*psi[7]*SQ(y)+5.809475019311126*psi[6]*x*y+1.5*psi[3]*y+3.354101966249685*psi[4]*x-0.9682458365518543*psi[7]+0.8660254037844386*psi[1];
      sol.dRdZ[sidx] = -C/A*dx[0]/dx[1];
      
      sidx += 1;      
    }
  }
  return sol;
}


struct gkgeom_app {
  char name[128];

  struct gkyl_rect_grid grid;
  struct gkyl_basis basis;
  struct gkyl_range local, local_ext;
  struct gkyl_array *psidg;

  double psi_min, psi_max;
  
  struct RdRdZ_sol (*calc_roots)(const double *psi, double psi0, double Z, double xc[2], double dx[2]);
  
  evalf_t psi;
  void *ctx;

  // stats
  long nfunc_calls;
};

static struct gkyl_array*
mkarr(size_t nc, size_t sz)
{
  return gkyl_array_new(GKYL_DOUBLE, nc, sz);
}

gkgeom_app*
gkgeom_app_new(const struct gkgeom_inp *inp)
{
  struct gkgeom_app *app = gkyl_malloc(sizeof(*app));
  strcpy(app->name, inp->name);

  gkyl_rect_grid_init(&app->grid, 2, inp->lower, inp->upper, inp->cells);
  gkyl_cart_modal_serendip(&app->basis, 2, inp->polyOrder);

  int nghost[GKYL_MAX_CDIM] = { 0, 0 };
  gkyl_create_grid_ranges(&app->grid, nghost, &app->local_ext, &app->local);

  // initilize psi(R,Z)
  app->psidg = mkarr(app->basis.num_basis, app->local_ext.volume);

  if (inp->use_proj_on_basis) {
    gkyl_proj_on_basis *pob = gkyl_proj_on_basis_new(&app->grid,
      &app->basis, inp->polyOrder+2, 1, inp->psi, inp->ctx);
    gkyl_proj_on_basis_advance(pob, 0.0, &app->local, app->psidg);
    gkyl_proj_on_basis_release(pob);
  }
  else {
    gkyl_eval_on_nodes *eon = gkyl_eval_on_nodes_new(&app->grid,
      &app->basis, 1, inp->psi, inp->ctx);
    gkyl_eval_on_nodes_advance(eon, 0.0, &app->local, app->psidg);
    gkyl_eval_on_nodes_release(eon);
  }

  char fileNm[1024]; // buffer for file name
  do {
    const char *fmt = "%s_psi.gkyl";
    snprintf(fileNm, sizeof fileNm, fmt, app->name);  
    gkyl_grid_sub_array_write(&app->grid, &app->local, app->psidg, fileNm);
  } while (0);

  app->calc_roots = calc_RdR_p1;
  if (inp->polyOrder == 2)  
    app->calc_roots = calc_RdR_p2;

  app->nfunc_calls = 0;
  
  return app;
}

int
gkgeom_app_R_psiz(const gkgeom_app *app, double psi, double Z, int nmaxroots,
  double *R, double *dR)
{
  double Zlower = app->grid.lower[1], dZ = app->grid.dx[1];

  int zcell = app->local.lower[1] + (int) floor((Z-Zlower)/dZ);
  zcell = zcell > app->local.upper[1] ? app->local.upper[1] : zcell;

  struct gkyl_range rangeR;
  gkyl_range_deflate(&rangeR, &app->local, (int[]) { 0, 1 }, (int[]) { 0, zcell });

  struct gkyl_range_iter riter;
  gkyl_range_iter_init(&riter, &rangeR);

  int idx[2] = { 0, zcell };
  double dx[2] = { app->grid.dx[0], app->grid.dx[1] };

  int sidx = 0;
  
  // loop over all R cells to find psi crossing
  while (gkyl_range_iter_next(&riter)) {
    long loc = gkyl_range_idx(&rangeR, riter.idx);
    const double *psih = gkyl_array_cfetch(app->psidg, loc);

    // compute R(psi,Z) and dR(psi,z)/dZ in this cell
    double xc[2];
    idx[0] = riter.idx[0];
    gkyl_rect_grid_cell_center(&app->grid, idx, xc);

    struct RdRdZ_sol sol = app->calc_roots(psih, psi, Z, xc, dx);
    if (sol.nsol > 0) {
      for (int s=0; s<sol.nsol; ++s) {
        R[sidx] = sol.R[s];
        dR[sidx] = sol.dRdZ[s];
        sidx += 1;
      }
    }
    if (sidx > nmaxroots) break;
  }
  return sidx;
}

struct countour_ctx {
  const gkgeom_app *app;
  double psi, last_R;
  int ncall;
};

static inline double
countour_func(double Z, void *ctx)
{
  struct countour_ctx *c = ctx;
  c->ncall += 1;
  double R[2] = { 0 }, dR[2] = { 0 };
  int nr = gkgeom_app_R_psiz(c->app, c->psi, Z, 2, R, dR);

  double dRdZ = 0;
  
  if (nr > 0) {
    dRdZ = dR[0];
    c->last_R = R[0];      
    if (fabs(R[1]-c->last_R)  < fabs(R[0]-c->last_R)) {
      dRdZ = dR[1];
      c->last_R = R[1];
    }
  }
  
  return nr>0 ? sqrt(1+dRdZ*dRdZ) : 0.0;
}

double
gkgeom_app_integrate_psi_contour(const gkgeom_app *app, double rclose,
  double zmin, double zmax, double psi)
{
  struct countour_ctx c = { .app = app, .psi = psi, .ncall = 0, .last_R = rclose };
  double res = wp34s(countour_func, &c, zmin, zmax, 10, 1e-10);
  ((gkgeom_app *)app)->nfunc_calls += c.ncall;
  return res;
}

struct arc_length_ctx {
  const gkgeom_app *app;
  double psi, rclose, zmin, arcL;
};

static inline double
arc_length_func(double Z, void *ctx)
{
  struct arc_length_ctx *actx = ctx;
  double psi = actx->psi, rclose = actx->rclose, zmin = actx->zmin, arcL = actx->arcL;
  return gkgeom_app_integrate_psi_contour(actx->app, rclose, zmin, Z, psi) - arcL;
}

static inline double
choose_closest(int ref, double R[2], double out[2])
{
  return fabs(R[0]-ref) < fabs(R[1]-ref) ? out[0] : out[1];
}

static void
write_nodal_coordinates(const gkgeom_app *app, struct gkyl_range *nrange,
  struct gkyl_array *nodes)
{
  double lower[3] = { 0.0, 0.0, 0.0 };
  double upper[3] = { 1.0, 1.0, 1.0 };
  int cells[3];
  for (int i=0; i<nrange->ndim; ++i)
    cells[i] = gkyl_range_shape(nrange, i);
  
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2, lower, upper, cells);

  char fileNm[1024]; // buffer for file name
  const char *fmt = "%s_node_coords.gkyl";
  snprintf(fileNm, sizeof fileNm, fmt, app->name);  
  gkyl_grid_sub_array_write(&grid, nrange, nodes, fileNm);
}

void
gkgeom_app_calcgeom(gkgeom_app *app, const struct gkyl_rect_grid *cgrid,
  int poly_order, struct gkyl_array *mapc2p)
{
  int nodes[3] = { 1, 1, 1 };
  if (poly_order == 1)
    for (int d=0; d<cgrid->ndim; ++d)
      nodes[d] = cgrid->cells[d]+1;
  if (poly_order == 2)
    for (int d=0; d<cgrid->ndim; ++d)
      nodes[d] = 2*cgrid->cells[d]+1;

  struct gkyl_range nrange;
  gkyl_range_init_from_shape(&nrange, cgrid->ndim, nodes);
  struct gkyl_array *mc2p = gkyl_array_new(GKYL_DOUBLE, cgrid->ndim, nrange.volume);

  double dtheta = cgrid->dx[0], dphi = cgrid->dx[1], dalpha = cgrid->dx[2];
  double theta_lo = cgrid->lower[0], phi_lo = cgrid->lower[1], alpha_lo = cgrid->lower[2];

  double dx_fact = poly_order == 1 ? 1 : 0.5; // USE THE ACTUAL NODE LOCATIONS!
  dtheta *= dx_fact; dphi *= dx_fact; dalpha *= dx_fact;

  enum { TH_IDX, PH_IDX, AL_IDX }; // arrangement of computational coordinates
  enum { R_IDX, Z_IDX }; // arrangement of physical coordinates

  double R[2] = { 0 }, dR[2] = { 0 };
  double Rmax = app->grid.upper[0];

  struct arc_length_ctx arc_ctx = { .app = app };

  int cidx[2] = { 0 };
  for (int ip=nrange.lower[PH_IDX]; ip<=nrange.upper[PH_IDX]; ++ip) {

    double zmin = app->grid.lower[1];
    double zmax = app->grid.upper[1];
    
    double psi_curr = phi_lo + ip*dphi;
    double arcL = gkgeom_app_integrate_psi_contour(app, Rmax, zmin, zmax, psi_curr);

    double delta_arcL = arcL/(poly_order*cgrid->cells[TH_IDX]);

    cidx[PH_IDX] = ip;

    // set node coordinates of first node
    cidx[TH_IDX] = nrange.lower[TH_IDX];
    double *mc2p_n = gkyl_array_fetch(mc2p, gkyl_range_idx(&nrange, cidx));
    mc2p_n[Z_IDX] = zmin;
    int nr = gkgeom_app_R_psiz(app, psi_curr, zmin, 2, R, dR);
    mc2p_n[R_IDX] = choose_closest(Rmax, R, R); // CHOOSE CLOSEST TO LAST R!!

    // set node coordinates of rest of nodes
    double arcL_curr = 0.0;
    for (int it=nrange.lower[TH_IDX]+1; it<nrange.upper[TH_IDX]; ++it) {
      arcL_curr += delta_arcL;

      arc_ctx.psi = psi_curr;
      arc_ctx.rclose = Rmax;
      arc_ctx.zmin = zmin;
      arc_ctx.arcL = arcL_curr;

      double z_curr = ridders(arc_length_func, &arc_ctx, zmin, zmax, 0, arcL, 1e-10);
      int nr = gkgeom_app_R_psiz(app, psi_curr, z_curr, 2, R, dR);
      double r_curr = choose_closest(Rmax, R, R);

      cidx[TH_IDX] = it;
      double *mc2p_n = gkyl_array_fetch(mc2p, gkyl_range_idx(&nrange, cidx));
      mc2p_n[Z_IDX] = z_curr;
      mc2p_n[R_IDX] = r_curr;
    }

    // set node coordinates of last node
    cidx[TH_IDX] = nrange.upper[TH_IDX];
    mc2p_n = gkyl_array_fetch(mc2p, gkyl_range_idx(&nrange, cidx));
    mc2p_n[Z_IDX] = zmax;
    nr = gkgeom_app_R_psiz(app, psi_curr, zmax, 2, R, dR);
    mc2p_n[R_IDX] = choose_closest(Rmax, R, R);
  }

  // write nodal coordinates for debugging
  write_nodal_coordinates(app, &nrange, mc2p);
  
  gkyl_array_release(mc2p);
  
}

long
gkgeom_app_nroots(const gkgeom_app *app)
{
  return app->nfunc_calls;
}

void
gkgeom_app_release(gkgeom_app *app)
{
  gkyl_array_release(app->psidg);
  gkyl_free(app);
}
