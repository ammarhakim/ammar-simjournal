# Multi-ion ES shock problems

- f1: Fluid simulation of a stationary shock

- n1: LBO Kinetic simulation corresponding to f1 (MFP = 0.5 on LX=10)
- n2: LBO Kinetic simulation corresponding to f1 (MFP = 0.25 on LX=10)
- n3: LBO Kinetic simulation corresponding to f1 (MFP = 0.1 on LX=10)

- b1: Same as n1, except with BGK operator
- b2: Same as n2, except with BGK operator
- b3: Same as n3, except with BGK operator

- t1: Two species (electron and ions) mfp = 0.1 (LX = 10)
- t2: Same as t1, except larger domain mfp = 0.1 (LX = 50)
- t3: Same as t2 except mfp = 0.25 (LX = 50)
- t3: Same as t2 except temperature is higher (mfp/lambda < 1) (LX = 50)