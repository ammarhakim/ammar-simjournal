#include <gkgeom.h>

#define SQ(x) ((x) * (x))

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
      sol.R[sidx++] = r1*dx[0]*0.5 + xc[0];
    }
    if ((-1<=r2) && (r2 < 1)) {
      sol.nsol += 1;
      sol.R[sidx++] = r2*dx[0]*0.5 + xc[0];
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
  
  return app;
}

int
gkgeom_app_R_psiz(const gkgeom_app *app, double psi, double Z, int nmaxroots,
  double *R, double *dR)
{
  double Zlower = app->grid.lower[1], dZ = app->grid.dx[1];

  int zcell = app->local.lower[1] + (int) floor((Z-Zlower)/dZ);

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

void
gkgeom_app_release(gkgeom_app *app)
{
  gkyl_array_release(app->psidg);
  gkyl_free(app);
}
