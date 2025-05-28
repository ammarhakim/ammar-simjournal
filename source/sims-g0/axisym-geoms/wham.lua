-- Input file for mirrorgridgen Tool.

psiRZ_fname = "wham_hires.geqdsk_psi.gkyl"
include_axis = false -- should we include the r=0 axis?
write_psi_cubic = true -- write the bicubic interpolation to psi

-- field-line coordinate to use. one of:
-- sqrt_psi_cart_z
-- psi_cart_z
field_line_coordinate = psi_cart_z

-- lower and upper extents of computational space grid
-- (psi, phi, z)
lower = { 2.0e-6, 0.0, -2.0 }
upper = { 1.0e-3, 2*math.pi, 2.0 }
cells = { 10, 16, 64 }

