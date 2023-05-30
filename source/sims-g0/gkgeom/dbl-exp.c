#include <math.h>
#include <stdio.h>

double
wp34s(double (*func)(double,void*), void *ctx, double a, double b, int n, double eps)
{
  double thr = 10*sqrt(eps); // too generous for larger eps, e.g. eps=1e-9
  double c = (a+b)/2; // center (mean)
  double d = (b-a)/2; // half distance
  double s = func(c, ctx);
  double fp = 0, fm = 0;
  double p, e, v, h = 2;
  double tmax = log(2/M_PI * log((d < 1 ? 2*d : 2) / eps));
  int k = 0; // level
  do {
    double q, t;
    int j = 1;
    v = s*d*M_PI/2*h; // last sum
    p = 0;
    h /= 2;
    t = h;
    do {
      double ch = cosh(t);
      double ecs = cosh(M_PI/2 * sqrt(ch*ch - 1)); // = cosh(pi/2*sinh(t))
      double w = 1/(ecs*ecs);
      double r = sqrt(ecs*ecs - 1)/ecs;
      double x = d*r;
      if (c+x > a) {
        double y = func(c+x, ctx);
        if (isfinite(y))
          fp = y;
      }
      if (c-x < b) {
        double y = func(c-x, ctx);
        if (isfinite(y))
          fm = y;
      }
      q = ch*w*(fp+fm);
      p += q;
      j += 1+(k>0);
      t = j*h;
    } while (t <= tmax && fabs(q) > eps*fabs(p));
    s += p;
    ++k;
  } while (s && fabs(2*fabs(p) - fabs(s)) >= fabs(thr*s) && k <= n);
  s *= d*M_PI/2*h;
  e = fabs(v-s);
  if (10*e >= fabs(s)) {
    e += fabs(s);
    s = 0;
  }
  return s; // result with estimated absolute error e
}
