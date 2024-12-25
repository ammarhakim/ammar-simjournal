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

// integrate function f, range a..b, max levels n, error tolerance eps
double qthsh(double (*func)(double,void*), void *ctx, double a, double b, int n, double eps)
{
  const double tol = 10*eps;
  double c = (a+b)/2; // center (mean)
  double d = (b-a)/2; // half distance
  double s = func(c, ctx);
  double e, v, h = 2;
  int k = 0;
  if (n <= 0) // use default levels n=6
    n = 6; // 6 is “optimal”, 7 just as good taking longer
  if (eps <= 0) // use default eps=1E-9
    eps = 1E-9;
  do {
    double p = 0, q, fp = 0, fm = 0, t, eh;
    h /= 2;
    t = eh = exp(h);
    if (k > 0)
      eh *= eh;
    do {
      double u = exp(1/t-t);
      // = exp(-2*sinh(j*h)) = 1/exp(sinh(j*h))^2
      double r = 2*u/(1+u);
      // = 1 - tanh(sinh(j*h))
      double w = (t+1/t)*r/(1+u); // = cosh(j*h)/cosh(sinh(j*h))^2
      double x = d*r;
      if (a+x > a) {
        // if too close to a then reuse previous fp
        double y = func(a+x,ctx);
        if (isfinite(y))
          fp = y;
        // if f(x) is finite, add to the local sum
      }
      if (b-x < b) {
        // if too close to b then reuse previous fm
        double y = func(b-x,ctx);
        if (isfinite(y))
          fm = y;
        // if f(x) is finite, add to the local sum
      }
      q = w*(fp+fm);
      p += q;
      t *= eh;
    } while (fabs(q) > eps*fabs(p));
    v = s-p;
    s += p;
    ++k;
  } while (fabs(v) > tol*fabs(s) && k <= n);
  e = fabs(v)/(fabs(s)+eps);
  return d*s*h; // result with estimated relative error e
}
