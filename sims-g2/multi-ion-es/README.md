# Multi-ion ES shock problems

- f1: Fluid simulation of a stationary shock

- n1: LBO Kinetic simulation corresponding to f1 (MFP = 0.5 on LX=10)
- n2: LBO Kinetic simulation corresponding to f1 (MFP = 0.25 on LX=10)
- n3: LBO Kinetic simulation corresponding to f1 (MFP = 0.1 on LX=10)

- b1: Same as n1, except with BGK operator
- b2: Same as n2, except with BGK operator
- b3: Same as n3, except with BGK operator

The following input file is not right. Ions and electrons have same
thermal velocity (which is probably ok) but then the mfp are
different. Should instead setup a stationary shock for ions and then
electrons will stick to them due to ES forces.

- t1: Two species (electron and ions) corresponding to n1