Linear dispersion solver test

Isothermal tests
----------------

- iso-1: Cold plasma waves; transverse propagation
- iso-2: Buneman instability, mi/me = 25
- iso-3: Buneman instability, mi/me = 200
- iso-4: Buneman instability, mi/me = 1836.2
- iso-ecdi-1: ECDI, mi/me = 400
- iso-ecdi-2: ECDI, mi/me = 1836.2
- iso-ecdi-3: ECDI, mi/me = 400, but both species are cold
- iso-weibel-1: Weibel. ud = 0.1. Cold fluids
- iso-weibel-2: Weibel. ud = 0.1. vt = ud.

5m tests
--------

- 5m-1: Same as iso-1, except with 5M model. Finite temperature (T = 0.1)
- 5m-ecdi-1: same as iso-ecdi-1 (ECDI, mi/me = 400)

10m tests
---------

- 10m-1: Same as iso-1, except with 10M model. Finite temperature (T = 0.1)
- 10m-ecdi-1: same as iso-ecdi-1 (ECDI, mi/me = 400)
- 10m-ecdi-3: using 10M electrons but isothermal ions instead. (shows
  the same second harmonic driven growing mode also).
- 10m-weibel-1: Weibel. ud = 0.1. elcTemp = 1e-6
- 10m-weibel-2: Weibel. ud = 0.1. ud = vt.
- 10m-weibel-3: Weibel. ud = 0.1. ud = 1.5*vt
