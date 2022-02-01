# comparison plot
pgkyl -f ../iso-2/iso-2-buneman_frequencies.bp -l '$m_i/m_e=25$' -f ../iso-3/iso-3-buneman_frequencies.bp -l '$m_i/m_e=200$' -f ../iso-4/iso-4-buneman_frequencies.bp -l '$m_i/m_e=1836.2$' val2 -x0 -y 4::2 pl -s --markersize=3 -f0 --xlabel='$k V_0/\omega_{pe}$' --ylabel='$\gamma/\omega_{pe}$' -t "Buneman instability growth rate" --hashtag --saveas="iso-bune-cmp" &

# real frequency
pgkyl -f ../iso-4/iso-4-buneman_frequencies.bp -l '$m_i/m_e=1836.2$' val2 -x0 -y 3::2 pl -s --markersize=3 -f0 --xlabel='$k V_0/\omega_{pe}$' --ylabel='$\omega_r/\omega_{pe}$' -t "Buneman oscillation frequency" --hashtag &

# imag frequency
pgkyl -f ../iso-4/iso-4-buneman_frequencies.bp -l '$m_i/m_e=1836.2$' val2 -x0 -y 4::2 pl -s --markersize=3 -f0 --xlabel='$k V_0/\omega_{pe}$' --ylabel='$\omega_r/\omega_{pe}$' --hashtag &

# complex frequency
pgkyl -f ../iso-4/iso-4-buneman_frequencies.bp -l '$m_i/m_e=1836.2$' val2 -x 3::2 -y 4::2 pl -s --markersize=3 -f0 --xlabel='$\omega_r/\omega_{pe}$' --ylabel='$\gamma/\omega_{pe}$' --hashtag --saveas="iso-bune-complex" &
