pgkyl -f t4-lbo-shock_elc_M0_10.bp -l "\$e^-\$"  -f t4-lbo-shock_ion_M0_10.bp -l "\$H^-\$" interp -p 2 -b ms pl -f0 --xscale 10.0 -x "\$x/\lambda_D\$" -y "Number Density" &

pgkyl -f t4-lbo-shock_field_10.bp -l "\$e^-\$"  interp -p 2 -b ms sel -c0 pl  --xscale 10.0 -x "\$x/\lambda_D\$" -y "\$E_x\$" &

pgkyl -f t4-lbo-shock_elc_M0_10.bp -f t4-lbo-shock_ion_M0_10.bp interp -p 2 -b ms ev "f1 f0 -" pl --xscale 10.0 -x "\$x/\lambda_D\$" -y "\$\varrho_c\$"
