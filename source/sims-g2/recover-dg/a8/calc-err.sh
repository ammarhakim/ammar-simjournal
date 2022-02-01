pgkyl -f a8-advection-recover-dg_f_0.bp -f a8-advection-recover-dg_f_1.bp ev "f0[0] f1[0] - sq f0[1] f1[1] - sq f0[2] f1[2] - sq + +" integrate 0 ev "f0 sqrt" info
