pgkyl -f ../10m-ecdi-3/10m-ecdi-3_frequencies.bp -l '$m_i/m_e = 400$' val2 -x0 -y 4::2 pl -s --markersize=3 -f0 --xscale 0.2 --yscale 1.0 --hashtag --xlabel='$kV_{E\times B}/\omega_{ce}$' --ylabel='$\gamma/\omega_{ce}$' --saveas="10m-ecdi-growth" &

pgkyl -f ../10m-ecdi-3/10m-ecdi-3_frequencies.bp -l '$m_i/m_e = 400$' val2 -x 3::2 -y 4::2 pl -s --xlabel='$\omega_r/\omega_{ce}$' --ylabel='$\gamma/\omega_{ce}$' --markersize=3 -f0 --hashtag --saveas="10m-ecdi-complex" &
