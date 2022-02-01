pgkyl -f ../iso-ecdi-1/iso-ecdi-1_frequencies.bp -l 'isothermal' -f ../10m-ecdi-1/10m-ecdi-1_frequencies.bp -l 'ten-moment' val2 -x0 -y 4::2 pl -s --markersize=3 -f0 --xscale 0.2 --yscale 1.0 --hashtag --xlabel='$kV_{E\times B}/\omega_{ce}$' --ylabel='$\gamma/\omega_{ce}$' --saveas="10m-ecdi-growth" &

pgkyl -f ../10m-ecdi-1/10m-ecdi-1_frequencies.bp -l '$m_i/m_e = 400$' val2 -x 3::2 -y 4::2 pl -s --xlabel='$\omega_r/\omega_{ce}$' --ylabel='$\gamma/\omega_{ce}$' --markersize=3 -f0 --hashtag --saveas="10m-ecdi-complex" &

pgkyl -f ../10m-ecdi-1/10m-ecdi-1_frequencies.bp -l 'ten-moment' val2 -x0 -y 3::2 pl -s --markersize=2 -f0 --xscale 0.2 --yscale 1.0 --hashtag --xlabel='$kV_{E\times B}/\omega_{ce}$' --ylabel='$\omega_r/\omega_{ce}$'  &
