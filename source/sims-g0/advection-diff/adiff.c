#include <adiff.h>

struct adiff_app {
  int nframe; // number of output frames to write
  double tend; // time to run simulation
  double cfl; // CFL condition  
  
  // grid and range for cell-centered quantities
  struct gkyl_rect_grid grid;
  struct gkyl_range range;

  // grid and range for nodal quantities
  struct gkyl_rect_grid nc_grid;
  struct gkyl_range nc_range;

  // range for nodes to update
  struct gkyl_range nc_up_range;

  double alpha; // diffusion coefficient
  struct gkyl_array *xvel; // x-component of velocity
  struct gkyl_array *yvel; // y-component of velocity
};

adiff_app *
adiff_app_new(const struct adiff_app_inp *inp)
{
  adiff_app *app = gkyl_malloc(sizeof *app);

  app->nframe = inp->nframe;
  app->tend = inp->tend;
  app->cfl = inp->cfl;

  // create grid and range to on cell-centered grid  
  gkyl_rect_grid_init(&app->grid, 2, inp->lower, inp->upper, inp->cells);
  gkyl_range_init_from_shape(&app->range, app->grid.ndim, app->grid.cells);

  return app;
}

void
adiff_app_release(adiff_app *app)
{
  
  gkyl_free(app);
}
