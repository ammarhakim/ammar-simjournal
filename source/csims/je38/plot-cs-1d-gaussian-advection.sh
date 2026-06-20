#!/usr/bin/sh
pgkyl cs-1d-gaussian-advection-nx-128_f_0.gkyl -l "Exact" cs-1d-gaussian-advection-nx-16_f_1.gkyl -l "16" cs-1d-gaussian-advection-nx-32_f_1.gkyl -l "32" cs-1d-gaussian-advection-nx-64_f_1.gkyl -l "64" cs-1d-gaussian-advection-nx-128_f_1.gkyl -l "128" pl -f0 --title "Advection of Gaussian in 1D: Upwind Scheme" --xlabel "X" --ylabel "\$f\$" --saveas "cs-1d-gaussian-advection-nx-scan.png" &

pgkyl cs-1d-gaussian-advection-nx-c128_f_0.gkyl -l "Exact" cs-1d-gaussian-advection-nx-c32_f_1.gkyl -l "32" cs-1d-gaussian-advection-nx-c64_f_1.gkyl -l "64" cs-1d-gaussian-advection-nx-c128_f_1.gkyl -l "128" pl -f0 --title "Advection of Gaussian in 1D: Central Scheme" --xlabel "X" --ylabel "\$f\$" --saveas "cs-1d-gaussian-advection-nx-cscan.png" &

pgkyl cs-1d-gaussian-advection-nx-l2-c64-a_diag.gkyl -l "\$\\Delta t\$" cs-1d-gaussian-advection-nx-l2-c64-b_diag.gkyl -l "\$\\Delta t/2\$" cs-1d-gaussian-advection-nx-l2-c64-c_diag.gkyl -l "\$\\Delta t/4\$" sel -c1 pl -f0 --xlabel "Time" --ylabel "\$L_2\$ Error" --title "Error in \$L_2\$ Norm for Central Scheme" --saveas "cs-1d-gaussian-advection-l2.png"

#pgkyl cs-1d-gaussian-advection-nx-16_f_0.gkyl cs-1d-gaussian-advection-nx-16_f_1.gkyl ev "f[0] f[1] - sq" integ 0 pr &
#pgkyl cs-1d-gaussian-advection-nx-32_f_0.gkyl cs-1d-gaussian-advection-nx-32_f_1.gkyl ev "f[0] f[1] - sq" integ 0 pr &
#pgkyl cs-1d-gaussian-advection-nx-64_f_0.gkyl cs-1d-gaussian-advection-nx-64_f_1.gkyl ev "f[0] f[1] - sq" integ 0 pr &
#pgkyl cs-1d-gaussian-advection-nx-128_f_0.gkyl cs-1d-gaussian-advection-nx-128_f_1.gkyl ev "f[0] f[1] - sq" integ 0 pr &
