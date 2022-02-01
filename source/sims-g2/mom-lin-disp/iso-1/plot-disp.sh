pgkyl -f iso-1-waves_frequencies.bp val2coord -x0 -y 3::2 pl -s -f0 --no-legend --xlabel "k" --ylabel "$\omega_r$"  --hashtag --saveas iso-1_real --markersize=2 -t "Cold plasma dispersion $\mathbf{B}=(1, 0, 0)$" &

pgkyl -f iso-1-waves_frequencies.bp val2coord -x0 -y 4::2 pl -s -f0 --no-legend --xlabel "k" --ylabel "$\omega_i$"  --hashtag --saveas iso-1_imag --markersize=2 &
