#pragma once

#include <gkylzero.h>

typedef struct gkgeom_app gkgeom_app;

struct gkgeom_inp {
  char name[128]; // name of app: used as output prefix
  
  int polyOrder; // polynomial order
  double lower[2], upper[2]; // lower, upper bounds
  int cells[2]; // config-space cells

  evalf_t psi; // psi(R,Z)
  void *ctx; // eval context
};

// Create new app object
gkgeom_app *gkgeom_app_new(const struct gkgeom_inp *inp);

// Compute the geometry
void gkgeom_app_calcgeom(gkgeom_app *app);

// Compute R given psi and Z (not to be called directly; just for
// testing). Output is in the R array
void gkgeom_app_R_psiz(const gkgeom_app *app, double psi, double Z, double R[2], double dR[2]);

// Release app
void gkgeom_app_release(gkgeom_app *app);
