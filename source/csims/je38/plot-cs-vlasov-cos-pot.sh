#!/usr/bin/sh

pgkyl cs-vlasov-cos-pot_f_4.gkyl -l "4" cs-vlasov-cos-pot_f_8.gkyl -l "8" cs-vlasov-cos-pot_f_12.gkyl -l "12" cs-vlasov-cos-pot_f_16.gkyl -l "16" pl -b --xlabel "X" --ylabel "V" --title "Charged Particles in a \$\\cos(x)\$ Potential" --saveas="cs-vlasov-cos-pot-u3.png"
