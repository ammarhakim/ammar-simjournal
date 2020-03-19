pgkyl -f w5-recovery-dg-diff_f_100.bp -f w5-recovery-dg-diff_exactSol.bp ev "f0[0] f1[0] - sq f0[1] f1[1] - sq f0[2] f1[2] - sq f0[3] f1[3] - sq + + +" integrate 0 integrate 1 ev "f0 sqrt" info
