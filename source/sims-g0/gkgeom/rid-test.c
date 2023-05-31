#include <stdio.h>
#include <math.h>

double ridders(double (*func)(double,void*), void *ctx,
  double x1, double x2, double f1, double f2, double eps);

static inline double sq(double x) { return x*x; }

struct f_ctx { int ncall; };

static inline double f(double x, void *ctx)
{
  struct f_ctx *fc = ctx;
  fc->ncall++;
  return sin(x)-0.5;
}

int
main(void)
{
  double xl = -1.0, xr = 1.0;
  struct f_ctx fc = { };
  
  double fl = f(xl, &fc), fr = f(xr, &fc);
  double res = ridders(f, &fc, xl, xr, fl, fr, 1e-12);
  printf("res = %lg. ncall = %d\n", res, fc.ncall);
  
  return 0;
}
