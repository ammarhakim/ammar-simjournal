#include "gkyl_array.h"
#include "gkyl_array_ops.h"
#include "gkyl_array_rio.h"
#include "gkyl_elem_type.h"
#include "gkyl_fv_proj.h"
#include "gkyl_range.h"

#include <math.h>
#include <stdio.h>
#include <gkylzero.h>

#define SQ(x) ((x)*(x))

void
src_func(double t, const double * restrict xn, double* restrict fout, void *ctx)
{
  // source is a gaussian centered at (0.25, 0.5)
  
  double xc = 0.25, yc = 0.5;
  double sig = 0.05;
  double x = xn[0], y = xn[1];
  double r2 = SQ(x-xc) + SQ(y-yc);
  fout[0] = exp(-r2/(2.0*SQ(sig)));
}

void
init_source(struct gkyl_rect_grid grid, struct gkyl_range range,
  struct gkyl_array *S)
{
  gkyl_fv_proj *fv_proj = gkyl_fv_proj_new(&grid, 2, 1, src_func, 0);
  gkyl_fv_proj_advance(fv_proj, 0.0, &range, S);
  gkyl_fv_proj_release(fv_proj);
}

int
main(void)
{
  int NX = 64, NY = 64;
  int LX = 1.0, LY  = 1.0;

  /// create grid
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 2,
    (double[]) { 0.0, 0.0 }, // coordinate of lower-left corner
    (double[]) { LX, LY }, // coorindate of upper-right corner
    (int[]) { NX, NY } // number of cells in each direction
  );

  double dx = grid.dx[0], dy = grid.dx[0];

  // grid on which the nodes of 'grid' are cell-centers: this is
  // needed as G0 and pgkyl assume that we only work with
  // cell-centered fields
  struct gkyl_rect_grid nc_grid;
  gkyl_rect_grid_init(&nc_grid, 2,
    (double[]) { grid.lower[0]-0.5*dx, grid.lower[1]-0.5*dy },
    (double[]) { grid.upper[0]+0.5*dx, grid.upper[1]+0.5*dy },
    (int[]) { NX+1, NY+1 }
  );

  // create range to on cell-centered grid
  struct gkyl_range range;
  gkyl_range_init_from_shape(&range, grid.ndim, grid.cells);

  // create range to represent nodes
  struct gkyl_range nc_range;
  gkyl_range_init_from_shape(&nc_range, nc_grid.ndim, nc_grid.cells);

  // allocate temperature field on nodes
  struct gkyl_array *T = gkyl_array_new(GKYL_DOUBLE, 1, nc_range.volume);
  struct gkyl_array *Tnew = gkyl_array_new(GKYL_DOUBLE, 1, nc_range.volume);
  
  gkyl_array_clear(T, 0.0); gkyl_array_clear(Tnew, 0.0);
  
  // allocate source on nodes
  struct gkyl_array *S = gkyl_array_new(GKYL_DOUBLE, 1, nc_range.volume);

  // initialize source
  init_source(nc_grid, nc_range, S);

  // write source to file
  gkyl_grid_sub_array_write(&nc_grid, &nc_range, S, "source.gkyl");

  // allocate heat-flux (vector) field on cells
  struct gkyl_array *Q = gkyl_array_new(GKYL_DOUBLE, grid.ndim, range.volume);
  gkyl_array_clear(Q, 0.0);


  // release resources
  gkyl_array_release(T);
  gkyl_array_release(Tnew);
  gkyl_array_release(S);
  gkyl_array_release(Q);

  return 0;
}
