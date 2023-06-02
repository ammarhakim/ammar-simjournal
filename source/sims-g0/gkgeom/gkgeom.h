#pragma once

#include <gkylzero.h>

typedef struct gkgeom_app gkgeom_app;

struct gkgeom_inp {
  char name[128];
  
  int polyOrder;
  double lower[2], upper[2];
  int cells[2];

  bool use_proj_on_basis;

  evalf_t psi; // psi(R,Z)
  void *ctx; // eval context
};

// Create new app object
gkgeom_app *gkgeom_app_new(const struct gkgeom_inp *inp);

// Compute the geometry on specified computational grid. The output
// mapc2p array is a poly_oder DG representation of the mapping. The
// array must be pre-allocated.
void gkgeom_app_calcgeom(gkgeom_app *app, const struct gkyl_rect_grid *cgrid,
  int poly_order, struct gkyl_array *mapc2p);

// These are methods for debugging the geometry calculator

// Compute R given psi and Z (not to be called directly; just for
// testing). Output is in the R array. Number of roots is
// returned. The output array must be sufficiently large to hold the
// roots found.
int gkgeom_app_R_psiz(const gkgeom_app *app, double psi, double Z, int nmaxroots,
  double *R, double *dR);

// Integrate the psi contour from zmin to zmax. rclose is the r
// coordinate of a point close to the initial zmin
double gkgeom_app_integrate_psi_contour(const gkgeom_app *app, double rclose,
  double zmin, double zmax, double psi);

// Number of calls to roots computed
long gkgeom_app_nroots(const gkgeom_app *app);

// Release app
void gkgeom_app_release(gkgeom_app *app);
