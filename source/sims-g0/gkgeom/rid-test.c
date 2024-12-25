#include <stdio.h>
#include <math.h>

double ridders(double (*func)(double,void*), void *ctx,
  double x1, double x2, double f1, double f2, double eps);

double qthsh(double (*func)(double,void*), void *ctx, double a, double b, int n, double eps);
double
wp34s(double (*func)(double,void*), void *ctx, double a, double b, int n, double eps);

static inline double sq(double x) { return x*x; }

struct f_ctx { int ncall; };


static inline double f(double x, void *ctx)
{
  struct f_ctx *fc = ctx;
  fc->ncall++;
  return sin(x)-0.5;
}

double
func_12(double x, void *ctx)
{
  struct f_ctx *fc = ctx;
  fc->ncall++;  
  return 1/sqrt(sin(M_PI*x));
}

double
func_13(double x, void *ctx)
{
  struct f_ctx *fc = ctx;
  fc->ncall++;
  return pow(sin(M_PI*x), -0.8);
}

double
func_cir(double x, void *ctx)
{
  struct f_ctx *fc = ctx;
  fc->ncall++;
  return sqrt(1+x*x/(1-x*x));
}

int
main(void)
{
  double xl = -1.0, xr = 1.0;
  struct f_ctx fc = { };
  
  double fl = f(xl, &fc), fr = f(xr, &fc);
  double res = ridders(f, &fc, xl, xr, fl, fr, 1e-12);
  printf("res = %lg. ncall = %d\n", res, fc.ncall);

  {
    fc.ncall = 0;
    double res_q = qthsh(func_12, &fc, 0.0, 1.0, 10, 1e-9);
    printf("res_q = %g. Exact = %g. Rel error = %g. Ncall = %d\n",
      res_q, 1.669253683348149, (res_q-1.669253683348149)/1.669253683348149, fc.ncall);
  }

  {
    fc.ncall = 0;
    double res_q = wp34s(func_12, &fc, 0.0, 1.0, 10, 1e-15);
    printf("res_q = %g. Exact = %g. Rel error = %g. Ncall = %d\n",
      res_q, 1.669253683348149, (res_q-1.669253683348149)/1.669253683348149, fc.ncall);
  }

  {
    fc.ncall = 0;
    double res_q = qthsh(func_13, &fc, 0.0, 1.0, 10, 1e-15);
    printf("res_q = %g. Exact = %g. Rel error = %g. Ncall = %d\n",
      res_q, 3.604250526330095, (res_q-3.604250526330095)/3.604250526330095, fc.ncall);
  }

  {
    fc.ncall = 0;
    double res_q = wp34s(func_cir, &fc, -1.0, 1.0, 10, 1e-15);
    printf("res_q = %g. Exact = %g. Rel error = %g. Ncall = %d\n",
      res_q, M_PI, (res_q-M_PI)/M_PI, fc.ncall);
  }  
  
  return 0;
}
