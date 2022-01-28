#include "gkyl_alloc.h"
#include "gkyl_range.h"
#include "gkyl_rect_grid.h"
#include "gkyl_util.h"
#include <math.h>
#include <stdio.h>
#include <gkylzero.h>

#define SQ(x) ((x) * (x))

// App input
struct app_inp {
  int ndim; // space dimensions

  double lower[3], upper[3]; // lower/upper bounds of domain
  int cells[3]; // number of cells in each direction

  int nframe; // number of data-frames to write
  double tend; // end time for simulation

};

// Incompressible Navier-Stokes app
struct ns_app {
  struct gkyl_rect_grid grid; // nodal grid
  struct gkyl_range range; // nodal range
  double dx[3]; // cell-size

  struct gkyl_rect_grid xe_grid;
  struct gkyl_range xe_range;
};

// App API

/**
 * Construct new app to solve incompressible N-S equations.
 *
 * @param inp Input to app.
 * @return new N-S app
 */
struct ns_app* ns_app_new(struct app_inp inp);

/**
 * Release resources allocated by app.
 *
 * @param ns App to release.
 */
void ns_app_release(struct ns_app *ns);

// App definition
struct ns_app*
ns_app_new(struct app_inp inp)
{
  struct ns_app *ns = gkyl_malloc(sizeof(*ns));

  double lower[3], upper[3];
  int ndim = inp.ndim, cells[3];

  double dx[3];
  for (int d=0; d<ndim; ++d)
    ns->dx[d] = dx[d] = (inp.upper[d]-inp.lower[d])/inp.cells[d];

  for (int d=0; d<ndim; ++d) {
    cells[d] = inp.cells[d]+1; // 1 more node than cells
    lower[d] = inp.lower[d]-0.5*dx[d];
    upper[d] = inp.upper[d]+0.5*dx[d];
  }

  // construct grid and range over nodes
  gkyl_rect_grid_init(&ns->grid, ndim, lower, upper, cells);
  gkyl_range_init_from_shape(&ns->range, ndim, cells);

  return ns;
}

void
ns_app_release(struct ns_app *ns)
{
  gkyl_free(ns);
}

// rens app input (eventually from input file)
struct app_inp
get_app_inp(int argc, char *argv[])
{
  return (struct app_inp) {
    .cells = { 64, 64 },
    .lower = { 0, 0 },
    .upper = { 1.0, 1.0 },

    .nframe = 20,
    .tend = 2.0,
  };
}

int
main(int argc, char *argv[])
{
  struct app_inp inp = get_app_inp(argc, argv);
  struct ns_app *ns = ns_app_new(inp);


  ns_app_release(ns);

  return 0;
}
