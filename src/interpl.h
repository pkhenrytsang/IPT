#include <cstddef> 

class dinterpl{

  public:

    dinterpl(const double x_array[], const double y_array[], size_t size);
    ~dinterpl();

  private:
    size_t  cache;        /* cache of index   */
    size_t  miss_count;   /* keep statistics  */
    size_t  hit_count;

    double* cspline_a;
    double* cspline_b;
    double* cspline_c;
    double* cspline_d;

    double* x_array;
    double* y_array;

    size_t size;


};


/* general interpolation object */
typedef struct {
  double* a;
  double* b;
  double* c;
  double* d;
  //double  xmin;
  //double  xmax;
  //size_t  size;
} cubic_spline;

typedef struct {
  size_t  cache;        /* cache of index   */
  size_t  miss_count;   /* keep statistics  */
  size_t  hit_count;
} interp_cache;

/*
  Caching for binary search
*/
interp_cache * interp_cache_alloc (void);
void interp_cache_free (interp_cache * a);

/*
  Linear Interpolation
*/
double linear_interpl_eval (double x , const double x_array[], const double y_array[], size_t size, interp_cache * cache);

/*
  Cubic Spline
*/
cubic_spline * cspline_init(const double x_array[], const double y_array[], size_t size);
void cspline_free(cubic_spline * spline);

double cspline_binary_eval(double x, const double x_array[], size_t size, cubic_spline * spline, interp_cache * cache);
