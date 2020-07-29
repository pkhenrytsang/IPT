#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

struct data
{
  double *t;
  double *y;
  size_t n;
};

/* model function: A/(B+w) a * exp( -1/2 * [ (t - b) / c ]^2 ) */
double tailfunction(const double A, const double B, const double w)
{
  return A/(w+B);
}

int func_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  struct data *d = (struct data *) params;
  double A = gsl_vector_get(x, 0);
  double B = gsl_vector_get(x, 1);

  for (int i = 0; i < d->n; ++i)
  {
    double wi = d->w[i];
    double yi = d->y[i];
    double y = tailfunction(A,B,wi);

    gsl_vector_set(f, i, yi - y);
  }

  return GSL_SUCCESS;
}

int func_df (const gsl_vector * x, void *params, gsl_matrix * J)
{
  struct data *d = (struct data *) params;
  double A = gsl_vector_get(x, 0);
  double B = gsl_vector_get(x, 1);

  for (int i = 0; i < d->n; ++i)
  {
    double wi = d->w[i];
		double deno = B+wi;
    gsl_matrix_set(J, i, 0, 1/deno );  //derivative respect to A
    gsl_matrix_set(J, i, 1, -A/deno/deno );  //derivative respect to B
  }

  return GSL_SUCCESS;
}


int func_fvv (const gsl_vector * x, const gsl_vector * v,void *params, gsl_vector * fvv)
{
  struct data *d = (struct data *) params;
  double A = gsl_vector_get(x, 0);
  double B = gsl_vector_get(x, 1);
  double vA = gsl_vector_get(v, 0);
  double vB = gsl_vector_get(v, 1);

  for (int i = 0; i < d->n; ++i)
  {
    double wi = d->w[i];
    double deno = B+wi;
    
    //second derivatives
    double DAB = -1/deno/deno;
    double DBB = 2*A/deno/deno/deno;
    
    double sum;

    sum = 2.0 * vA * vB * DAB +
                vB * vB * DBB;

    gsl_vector_set(fvv, i, sum);
  }

  return GSL_SUCCESS;
}

void
callback(const size_t iter, void *params,
         const gsl_multifit_nlinear_workspace *w)
{
  gsl_vector *f = gsl_multifit_nlinear_residual(w);
  gsl_vector *x = gsl_multifit_nlinear_position(w);
  double avratio = gsl_multifit_nlinear_avratio(w);
  double rcond;

  (void) params; /* not used */

  /* compute reciprocal condition number of J(x) */
  gsl_multifit_nlinear_rcond(&rcond, w);

  fprintf(stderr, "iter %2zu: A = %.4f, B = %.4f, |A|/|v| = %.4f cond(J) = %8.4f, |f(x)| = %.4f\n",
          iter,
          gsl_vector_get(x, 0),
          gsl_vector_get(x, 1),
          avratio,
          1.0 / rcond,
          gsl_blas_dnrm2(f));
}

void
solve_system(gsl_vector *x, gsl_multifit_nlinear_fdf *fdf,
             gsl_multifit_nlinear_parameters *params)
{
  const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
  const size_t max_iter = 200;
  const double xtol = 1.0e-8;
  const double gtol = 1.0e-8;
  const double ftol = 1.0e-8;
  const size_t n = fdf->n;
  const size_t p = fdf->p;
  gsl_multifit_nlinear_workspace *work =
    gsl_multifit_nlinear_alloc(T, params, n, p);
  gsl_vector * f = gsl_multifit_nlinear_residual(work);
  gsl_vector * y = gsl_multifit_nlinear_position(work);
  int info;
  double chisq0, chisq, rcond;

  /* initialize solver */
  gsl_multifit_nlinear_init(x, fdf, work);

  /* store initial cost */
  gsl_blas_ddot(f, f, &chisq0);

  /* iterate until convergence */
  gsl_multifit_nlinear_driver(max_iter, xtol, gtol, ftol,
                              callback, NULL, &info, work);

  /* store final cost */
  gsl_blas_ddot(f, f, &chisq);

  /* store cond(J(x)) */
  gsl_multifit_nlinear_rcond(&rcond, work);

  gsl_vector_memcpy(x, y);

  /* print summary */

  fprintf(stderr, "NITER         = %zu\n", gsl_multifit_nlinear_niter(work));
  fprintf(stderr, "NFEV          = %zu\n", fdf->nevalf);
  fprintf(stderr, "NJEV          = %zu\n", fdf->nevaldf);
  fprintf(stderr, "NAEV          = %zu\n", fdf->nevalfvv);
  fprintf(stderr, "initial cost  = %.12e\n", chisq0);
  fprintf(stderr, "final cost    = %.12e\n", chisq);
  fprintf(stderr, "final x       = (%.12e, %.12e, %12e)\n",
          gsl_vector_get(x, 0), gsl_vector_get(x, 1), gsl_vector_get(x, 2));
  fprintf(stderr, "final cond(J) = %.12e\n", 1.0 / rcond);

  gsl_multifit_nlinear_free(work);
}
