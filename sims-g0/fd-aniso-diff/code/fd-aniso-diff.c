#include "gkyl_array.h"
#include "gkyl_array_ops.h"
#include "gkyl_elem_type.h"
#include "gkyl_range.h"
#include <stdio.h>
#include <gkylzero.h>

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

  // create range to represent nodes (note: there is one extra node in
  // each direction than cells)
  struct gkyl_range node_range;
  gkyl_range_init(&node_range, grid.ndim,
    (int[]) { 0, 0 }, // index of lower-left node
    (int[]) { NX, NY } // index of upper-right node
  );

  // create sub-range of nodes to update
  struct gkyl_range up_node_range;
  gkyl_sub_range_init(&up_node_range, &node_range,
    (int[]) { 1, 1 },
    (int[]) { NX-1, NY-1 }
  );
  
  // allocate temperature field on nodes
  struct gkyl_array *T = gkyl_array_new(GKYL_DOUBLE, 1, node_range.volume);
  struct gkyl_array *Tnew = gkyl_array_new(GKYL_DOUBLE, 1, node_range.volume);

  // set initial value of temperature to 0
  gkyl_array_clear(T, 0.0);
  gkyl_array_clear(Tnew, 0.0);


  // release resources
  gkyl_array_release(T);
  gkyl_array_release(Tnew);

  return 0;
}
