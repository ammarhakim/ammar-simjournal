#include <math.h>
#include <float.h>

static inline double sgn(double x) { return x >= 0 ? 1 : -1; }
static inline double babs(double a, double b) { return b>=0 ? fabs(a) : -fabs(a); }

// Assume f1 and f2 have different signs: no checks are made
double
ridders(double (*func)(double,void*), void *ctx,
  double x1, double x2, double f1, double f2, double eps)
{
  double xl = x1, xr = x2, fl = f1, fr = f2;
  double res = DBL_MAX;
  
  int iterating = 1;
  while(iterating) {
    double xm = 0.5*(xl+xr);
    double fm = func(xm, ctx);
    double sr = sqrt(fm*fm - fl*fr);

    if (sr == 0) return res;
    
    double xnew = xm + (xm-xl)*sgn(fl-fr)*fm/sr;
    if (fabs(xnew-res) < eps) return res;
    res = xnew;
    double fnew = func(xnew, ctx);
    if (fnew == 0.0) return res;

    if (babs(fm, fnew) != fm) {
      xl = xm; fl = fm;
      xr = res; fr = fnew;
    }
    else if (babs(fl,fnew) != fl) {
      xr = res; fr = fnew;      
    }
    else if (babs(fr,fnew) != fr) {
      xl = res; fl = fnew;
    }
    if (fabs(xr-xl) < eps) iterating = 0;
  }
  
  return res;
}
