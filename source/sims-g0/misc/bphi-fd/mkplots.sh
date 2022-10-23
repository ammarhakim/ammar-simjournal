pgkyl Bphi_1000.gkyl pl -a --title "\$B_\phi\$" --xlabel "R" --ylabel "Z" &
pgkyl J_1000.gkyl sel -c0 pl -a --title "\$J_r\$" --xlabel "R" --ylabel "Z" &
pgkyl J_1000.gkyl sel -c1 pl -a --title "\$J_z\$" --xlabel "R" --ylabel "Z" &

pgkyl J_1000.gkyl Bphi_1000.gkyl ev "f[0][1] f[1] * -1 *" pl -a --title "\$F_r\$" --xlabel "R" --ylabel "Z" &
pgkyl J_1000.gkyl Bphi_1000.gkyl ev "f[0][0] f[1] *" pl -a --title "\$F_z\$" --xlabel "R" --ylabel "Z" &

pgkyl J_1000.gkyl pl -l -a --sdensity 0.5 --xlim 0,0.05 --ylim 0,0.01 --title "\$J\$" --xlabel "R" --ylabel "Z" &
