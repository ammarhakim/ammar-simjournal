#include <math.h>
#include <float.h>

static inline double dsign(double x) { return x >= 0 ? 1 : -1; }

// Assume f1 and f2 have different signs: no checks are made
double
ridders(double (*func)(double,void*), void *ctx,
  double x1, double x2, double f1, double f2, double eps)
{
  double xl = x1, xr = x2, fl = f1, fr = f2;
  double res = DBL_MAX;

  int nev = 0, nitr = 0, iterating = 1;
  while(iterating) {
    double xm = 0.5*(xl+xr);
    double fm = func(xm, ctx); nev += 1;
    double W = sqrt(fm*fm - fl*fr);
    if (W == 0) { res = xm;  break; }
    
    double xnew = xm + (xm-xl)*dsign(fl-fr)*fm/W;
    if (fabs(xnew-res) < eps) break;
    res = xnew;
    double fnew = func(xnew, ctx); nev += 1;
    if (fnew == 0.0) break;

    if (fm*fnew < 0) {
      xl = xm; fl = fm;
      xr = res; fr = fnew;
    }
    else if (fl*fnew < 0) {
      xr = res; fr = fnew;      
    }
    else if (fr*fnew < fr) {
      xl = res; fl = fnew;
    }

    if (fabs(xr-xl) < eps) iterating = 0;
    nitr += 1;
  }
  return res;
}
