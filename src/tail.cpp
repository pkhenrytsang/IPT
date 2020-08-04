#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "routines.h"

struct data
{
  double* w;
  double* y;
  size_t n;
};

/* model function: A/(B+w) */
double reDeltatailfunc(const double A, const double B, const double w)
{
  return A/(w+B);
}


/* model function: C+A/(B+w) */
double reSigmatailfunc(const double A, const double B, const double C, const double w)
{
  return C+A/(w+B);
}

int func_Sigmatail (const gsl_vector * x, void *params, gsl_vector * f)
{
  struct data *d = (struct data *) params;
  double A = gsl_vector_get(x, 0);
  double B = gsl_vector_get(x, 1);
  double C = gsl_vector_get(x, 2);

  for (int i = 0; i < d->n; ++i)
  {
    double wi = d->w[i];
    double yi = d->y[i];
    double y = reSigmatailfunc(A,B,C,wi);

    gsl_vector_set(f, i, yi - y);
  }

  return GSL_SUCCESS;
}

int func_Deltatail (const gsl_vector * x, void *params, gsl_vector * f)
{
  struct data *d = (struct data *) params;
  double A = gsl_vector_get(x, 0);
  double B = gsl_vector_get(x, 1);
  double C = gsl_vector_get(x, 2);

  for (int i = 0; i < d->n; ++i)
  {
    double wi = d->w[i];
    double yi = d->y[i];
    double y = reDeltatailfunc(A,B,wi);

    gsl_vector_set(f, i, yi - y);
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

  fprintf(stderr, "iter %2zu: A = %.4f, B = %.4f, C = %.4f, |A|/|v| = %.4f cond(J) = %8.4f, |f(x)| = %.4f\n",
          iter,
          gsl_vector_get(x, 0),
          gsl_vector_get(x, 1),
          gsl_vector_get(x, 2),
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
  const double xtol = 1.0e-14;
  const double gtol = 1.0e-14;
  const double ftol = 1.0e-14;
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
  fprintf(stderr, "final x       = (%.12e, %.12e)\n",
          gsl_vector_get(x, 0), gsl_vector_get(x, 1));
  fprintf(stderr, "final cond(J) = %.12e\n", 1.0 / rcond);

  gsl_multifit_nlinear_free(work);
}


int
main (void)
{
  int N;
  complex<double>* Delta;
  double* omega;
  
  {
  	int n,m;
  	double** output;
		ReadFunc("Sigma.out", n, m, output);  
		N=n;
		omega = new double[N]; //Grid
		Delta = new complex<double>[N];	//Bath
		
		for (int i = 0; i<n; i++) {
			omega[i] = output[0][i];
			Delta[i] = complex<double>(output[1][i],output[2][i]);
		}
		for (int i=0; i<m; i++)
			delete [] output[i];
		delete [] output;
  }


  const size_t n = 400;  /* number of data points to fit */
  const size_t p = 3;    /* number of model parameters */
  
  gsl_vector *f = gsl_vector_alloc(n);
  gsl_vector *x = gsl_vector_alloc(p);
  
  gsl_multifit_nlinear_fdf fdf;
  gsl_multifit_nlinear_parameters fdf_params =
    gsl_multifit_nlinear_default_parameters();
    
  struct data fit_data;
  size_t i;

	double* w = new double[n];
	double* y = new double[n];
  //fit_data.w = malloc(n * sizeof(double));
  //fit_data.y = malloc(n * sizeof(double));
  
  for (i = 0; i < n/2; ++i)
  {
      w[i] = omega[i];
      y[i] = real(Delta[i]);
  }
  
  for (i = 0; i < n/2; ++i)
  {
      w[i+n/2] = omega[N-1-i];
      y[i+n/2] = real(Delta[N-1-i]);
  }
  
  fit_data.w = w;
  fit_data.y = y;
  fit_data.n = n;
    
  /* define function to be minimized */
  fdf.f = func_Sigmatail;
  fdf.df = NULL;
  fdf.fvv = NULL;
  fdf.n = n;
  fdf.p = p;
  fdf.params = &fit_data;

  /* starting point */
  gsl_vector_set(x, 0, 1.0);
  gsl_vector_set(x, 1, 1.0);
  gsl_vector_set(x, 2, 1.0);

  fdf_params.trs = gsl_multifit_nlinear_trs_lmaccel;
  solve_system(x, &fdf, &fdf_params);
  
  /* print data and model */
  {
    double A = gsl_vector_get(x, 0);
    double B = gsl_vector_get(x, 1);
    double C = gsl_vector_get(x, 2);
    
    for (i = 0; i < n; ++i)
      {
        double wi = fit_data.w[i];
        double yi = fit_data.y[i];
        double fi = reSigmatailfunc(A, B,C, wi);

        printf("%f %f %f\n", wi, yi, fi);
      }
  }

  gsl_vector_free(f);
  gsl_vector_free(x);
  
  delete [] omega;
  delete [] Delta;
  delete [] w;
  delete [] y;

  return 0;
}
