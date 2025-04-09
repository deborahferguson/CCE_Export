Cactus Code Thorn CCE_Export
Author(s)    :  Deborah Ferguson
                Noora Ghadiri
Maintainer(s):  Deborah Ferguson
Licence      : XXX
--------------------------------------------------------------------------

1. Purpose

This thorn exports the metric, lapse, and shift data (as well as their radial and time derivatives) on spherical shells, decomposed using scalar spherical harmonics (up through $\ell=8$).

Specifically, they are exported in the format compatible to be fed into the SpECTRE CCE module. A tutorial for running the SpECTRE CCE module can be found here: https://spectre-code.org/tutorial_cce.html

A separate h5 file is created for each desired extraction radius. The data is then stored in the following datasets:

```
gxx.dat, gxy.dat, gxz.dat, gyy.dat, gyz.dat, gzz.dat
Drgxx.dat, Drgxy.dat, Drgxz.dat, Drgyy.dat, Drgyz.dat, Drgzz.dat
Dtgxx.dat, Dtgxy.dat, Dtgxz.dat, Dtgyy.dat, Dtgyz.dat, Dtgzz.dat
Shiftx.dat, Shifty.dat, Shiftz.dat
DrShiftx.dat, DrShifty.dat, DrShiftz.dat
DtShiftx.dat, DtShifty.dat, DtShiftz.dat
Lapse.dat
DrLapse.dat
DtLapse.dat
```
The first column of each of these datasets is the coordinate time. Each pair of columns following is the real and imaginary components of each of the spherical harmonic modes in m-changes-fastest order:
```
"time", "Re(0,0)", "Im(0,0)", "Re(1,-1)", "Im(1,-1)", "Re(1,0)", "Im(1,0)", "Re(1,1)", "Im(1,1)", "Re(2,-2)", "Im(2,-2)", "Re(2,-1)", "Im(2,-1)", "Re(2,0)", "Im(2,0)", "Re(2,1)", "Im(2,1)", "Re(2,2)", "Im(2,2)", ...
```
The above is also included as a legend for each dataset.