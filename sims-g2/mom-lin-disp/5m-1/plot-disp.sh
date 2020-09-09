pgkyl -f 5m-1-waves_frequencies.bp val2coord -x0 -y 3::2 pl -s -f0 --no-legend --xlabel "k" --ylabel "$\omega_r$"  --hashtag --saveas 5m-1_real --markersize=2 -t "Five-moment dispersion $\mathbf{B}=(1/2, 0, 1)$" &

pgkyl -f 5m-1-waves_frequencies.bp val2coord -x0 -y 4::2 pl -s -f0 --no-legend --xlabel "k" --ylabel "$\omega_i$"  --hashtag --saveas 5m-1_imag --markersize=2 &
