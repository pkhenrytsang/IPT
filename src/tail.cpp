#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <complex>

struct fit_params
{
  size_t max_iter; //Cap on number of iterations for curve fit
  size_t ntail; // number of points on the tail to fit
  double xtol; // xtol for lm
  double gtol; // gtol for lm
  double ftol; // ftol for lm
  bool verbose; // verbose output
  bool quiet; // do not produce any output
};

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

int func_Deltatail (const gsl_vector * x, void *params, gsl_vector * f)
{
  struct data *d = (struct data *) params;
  double A = gsl_vector_get(x, 0);
  double B = gsl_vector_get(x, 1);

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

	struct fit_params *fp = (struct fit_params *) params;
	bool verbose = fp->verbose;
	bool quiet = fp->quiet;

  /* compute reciprocal condition number of J(x) */
  gsl_multifit_nlinear_rcond(&rcond, w);

	if (verbose && !quiet){
  fprintf(stderr, "iter %2zu: A = %.4f, B = %.4f, |A|/|v| = %.4f cond(J) = %8.4f, |f(x)| = %.4f\n",
          iter,
          gsl_vector_get(x, 0),
          gsl_vector_get(x, 1),
          avratio,
          1.0 / rcond,
          gsl_blas_dnrm2(f));
  }
}

void
solve_system(gsl_vector *x, gsl_multifit_nlinear_fdf *fdf,
             gsl_multifit_nlinear_parameters *params,struct fit_params * fparams)
{
  const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
  const size_t max_iter = fparams->max_iter;
  const double xtol = fparams->xtol;
  const double gtol = fparams->gtol;
  const double ftol = fparams->ftol;
  bool verbose = fparams->verbose;
  bool quiet = fparams->quiet;
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
                              callback, fparams, &info, work);

  /* store final cost */
  gsl_blas_ddot(f, f, &chisq);

  /* store cond(J(x)) */
  gsl_multifit_nlinear_rcond(&rcond, work);

  gsl_vector_memcpy(x, y);

  /* print summary */

	if (verbose && !quiet){
  fprintf(stderr, "NITER         = %zu\n", gsl_multifit_nlinear_niter(work));
  fprintf(stderr, "NFEV          = %zu\n", fdf->nevalf);
  fprintf(stderr, "NJEV          = %zu\n", fdf->nevaldf);
  fprintf(stderr, "NAEV          = %zu\n", fdf->nevalfvv);
  fprintf(stderr, "initial cost  = %.12e\n", chisq0);
  fprintf(stderr, "final cost    = %.12e\n", chisq);
  fprintf(stderr, "final x       = (%.12e, %.12e)\n",
          gsl_vector_get(x, 0), gsl_vector_get(x, 1));
  fprintf(stderr, "final cond(J) = %.12e\n", 1.0 / rcond);
  }

  gsl_multifit_nlinear_free(work);
}

int fit_tail(double* omega,std::complex<double>* Delta, size_t N,double &A1,double &B1,double &A2,double &B2, struct fit_params * params)
{

  const size_t n = params->ntail;  /* number of data points to fit */
  const size_t p = 2;    /* number of model parameters */
  
  gsl_vector *f1 = gsl_vector_alloc(n);
  gsl_vector *x1 = gsl_vector_alloc(p);
  gsl_vector *f2 = gsl_vector_alloc(n);
  gsl_vector *x2 = gsl_vector_alloc(p);
  
  gsl_multifit_nlinear_fdf fdf1;
  gsl_multifit_nlinear_fdf fdf2;
  gsl_multifit_nlinear_parameters fdf_params1 =
    gsl_multifit_nlinear_default_parameters();
  gsl_multifit_nlinear_parameters fdf_params2 =
    gsl_multifit_nlinear_default_parameters();
    
  struct data fit_data1;
  struct data fit_data2;

	double* w1 = new double[n];
	double* y1 = new double[n];
	double* w2 = new double[n];
	double* y2 = new double[n];
  
  for (int i = 0; i < n; ++i)
  {
      w1[i] = omega[i];
      y1[i] = std::real(Delta[i]);
  }
  
  for (int i = 0; i < n; ++i)
  {
      w2[i] = omega[N-1-i];
      y2[i] = std::real(Delta[N-1-i]);
  }
  
  fit_data1.w = w1;
  fit_data1.y = y1;
  fit_data1.n = n;
  
  fit_data2.w = w2;
  fit_data2.y = y2;
  fit_data2.n = n;
    
  /* define function to be minimized */
  fdf1.f = func_Deltatail;
  fdf1.df = NULL;
  fdf1.fvv = NULL;
  fdf1.n = n;
  fdf1.p = p;
  fdf1.params = &fit_data1;
  
  fdf2.f = func_Deltatail;
  fdf2.df = NULL;
  fdf2.fvv = NULL;
  fdf2.n = n;
  fdf2.p = p;
  fdf2.params = &fit_data2;

  /* starting point */
  gsl_vector_set(x1, 0, 1.0);
  gsl_vector_set(x1, 1, 1.0);
  gsl_vector_set(x2, 0, 1.0);
  gsl_vector_set(x2, 1, 1.0);

  fdf_params1.trs = gsl_multifit_nlinear_trs_lmaccel;
  fdf_params2.trs = gsl_multifit_nlinear_trs_lmaccel;
  solve_system(x1, &fdf1, &fdf_params1,params);
  solve_system(x2, &fdf2, &fdf_params2,params);
  
  /* print data and model */
  A1 = gsl_vector_get(x1, 0);
  B1 = gsl_vector_get(x1, 1);
  
  A2 = gsl_vector_get(x2, 0);
  B2 = gsl_vector_get(x2, 1);
  
  
    
  gsl_vector_free(f1);
  gsl_vector_free(x1);
  gsl_vector_free(f2);
  gsl_vector_free(x2);
  
  delete [] w1;
  delete [] y1;
  delete [] w2;
  delete [] y2;

  return 0;
}
