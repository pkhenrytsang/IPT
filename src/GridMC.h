/*
	The class grid stores the grid that IPT read and writes, called by SIAM class.
	original code by Jaksha Vuchichevicc https://github.com/JaksaVucicevic/DMFT
	in the original code this class is named "Result"
*/


#include <complex>
using namespace std;

class GridMC
{
  public:
  
  	//Constructor-destructor
    GridMC(const double omega_in[],size_t N);
    ~GridMC();
    
    //Grid
    int N; //length of grid
    double* omega;		//omega grid
    
    //Input
    complex<double>* Delta;	//Bath
    
    //Output
    complex<double>* Sigma;	//Sigma interpolating between exact limiting cases
    complex<double>* G;		//Greens function on real axis
    
    //Occupation
    double n,n0;
    
    void PrintFullResult(const char* ResultFN); //Debug version - basically prints everything
    void PrintResult(const char* Gffile,const char* Sigfile); //Non-debug version to save disk space
    
  	//Intermediate steps
  	double* dos0;
  	double* dos;
    complex<double>* SOCSigma;	//Second order contribution in sigma
    
  private:
    void ReleaseMemory();
};
