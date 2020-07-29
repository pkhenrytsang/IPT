#include <complex>
using namespace std;

class Grid
{
  public:
    Grid(const double omega_in[],size_t N);
    ~Grid();
    
    size_t N;

    double n;
    double n0;
    double mu;
    double mu0;

    double* omega;		//omega grid
    double* fermi;		//fermi function
    double* Ap;			//spectral functions
    double* Am;
    double* P1;			//polarizations
    double* P2;
    complex<double>* SOCSigma;	//Second order contribution in sigma
    complex<double>* Sigma;	//Sigma interpolating between exact limiting cases
    complex<double>* G;		//Greens function on real axis
    complex<double>* Delta;	//Bath
    complex<double>* G0;	//auxillary Green's function
    double* DOS;		//quasi-particle density of states (typical DOS in TMT)
    double* NIDOS;		//non-interacting density of states
    double* DOSmed;		//medium DOS in TMT, can be used as an auxiallry DOS in other cases

    void PrintResult(const char* ResultFN); //Debug version - basically prints everything
    void PrintResult(const char* Gffile,const char* Sigfile); //Non-debug version to save disk space
    bool ReadFromFile(const char* ResultFN);
    

  private:
    void Initialize(const double omega_in[],size_t N);
    void ReleaseMemory();
};
