#!/usr/bin/sh
pgkyl cs-2d-rotflow-c_f_1.gkyl -l "Central" cs-2d-rotflow-u_f_1.gkyl -l "Upwind" pl -ad -b --xlabel "X" --ylabel "Y" --title "Advection in Rotating Flow. T=\$2\\pi\$" --saveas="csd-2d-rotflow-2d.png" &

pgkyl cs-2d-rotflow-c_f_0.gkyl -l "Exact" cs-2d-rotflow-c_f_1.gkyl -l "Central" cs-2d-rotflow-u_f_1.gkyl -l "Upwind" sel --z0 0.25 pl -f0 --xlabel "Y" --ylabel "f" --title "Advection in Rotating Flow. T=\$2\\pi\$" --saveas="csd-2d-rotflow-2d-lineout.png" &
