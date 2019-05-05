pgkyl -f t2-lbo-shock_elc_vthSq_20.bp -f t2-lbo-shock_ion_vthSq_20.bp -l "\$T_i\$" interp -p 2 -b ms dataset -i0 ev -l "\$T_e\$" "f0 1836.2 /"  dataset -i1,2 sel --c0 15.0:30.0 pl -f0
