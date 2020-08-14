//**********************************************************//
//              SIAM at Arbitrary Filling                   //
//                by Pak Ki Henry Tsang                     //
//        Adapted from code by Jaksha Vuchichevicc          //
//**********************************************************//

#include <iostream>
#include <complex>


/*

*/

//Buffer size for printouts  (It should be sufficient unless AmoebaMaxIts is increased to absurdly large values)
#define BUFFERSIZE 65536

//Maximum subdivision for adaptive quadrature
#define QUADLIMIT 200

//epsrel for adaptive quadrature
#define QUADACCR 1e-9

//
#define CLIPVAL 1e-5


class GridMC;
class dinterpl;

using namespace std;

//Siam params

struct siammcparams{
	double mu,mu0,U,T,epsilon,eta;
	double KKAccr,Accr;
	
	//Monte Carlo
	int M;
	double W;
	
	//Tail
	double * L;
	double * R;
	int fitorder;
	bool SymmetricCase,Fixmu0,usecubicspline,tailcorrection;
	
	bool CheckSpectralWeight;
	bool verbose;
	int AmoebaMaxIts;
	double AmoebaScanStart,AmoebaScanEnd,AmoebaScanStep;
	bool AmoebaForceScanAndPrintOut;
};

//======================= SIAM Class ==========================================//

class SIAMMC
{
  public:
  
    //--Constructors/destructors--//
    SIAMMC(const double omega[],size_t N, void * params);
    ~SIAMMC();
    
    //--------RUN SIAM--------//
    int Run(const complex<double> Delta[]); 
    
    //--Result IO--//
    void PrintResult(); //print to "Gf.out" and "Sig.out"
    void PrintResult(const char* Gffile,const char* Sigfile);
    void PrintFullResult(const char* ResultFN);
    
    //--Buffer IO--//
    const char* buffer();
    void PrintBuffer(const char* FN);
    void PrintBuffer(const char* FN,bool quiet);

    //No meaning so far
    int status;
  private:
  
  	//MPI
  	int rank,size;
  	
  	//Integration
  	int M;
		double W;
  
  	dinterpl * dos0spline;
  	dinterpl * dosspline;
  	dinterpl * imSOCSigmaspline;
  
  	//verbose output (currently not useful as we print all output)
  	bool verbose;
  	
  	//print buffer
  	char * ibuffer;
  	
  	//parameters
  	struct siammcparams *p;
  	
    //--Grid--//
    GridMC* g;
    int N;
    
    //--impurity parameters--//
    double U;			//on-site repulsion
    double T;			//temperature
    double epsilon;		//impurity energy level

    //--bath parameters--// 
    double mu;			//global chemical potential
    double mu0;			//fictious chemical potential
    
    //---BROADENING---//
    double eta;
    complex<double> ieta;

    //----tail correction----//
    bool tailcorrection;
    int fitorder;
    double *L;
    double *R;
    
		double getn0corr();
		double getwG0corr();

    //----kramers kronig----//
    double KKAccr;
    bool usecubicspline;
    
     //--get functions--//
    double get_n(double dos[]);

		inline double Ap(double om);
		inline double Am(double om);
		
		//Cubic-spline version
		inline double Apc(double om);
		inline double Amc(double om);

    //--get procedures--//
    void get_dos0();
    void get_SOCSigma();
    double get_b();
    void get_Sigma();
    void get_G();

		//--Clipping of ImSigma and ImG--//
    bool ClipOff(complex<double> &X);
    bool Clipped;
    double ClippingValue = CLIPVAL;

    //--- SIAM solver ---//
    void SolveSiam(double* V);
    void Amoeba(double accr, double* V);	//amoeba method for mu0 search.
  
    //-- mu0 search --//
    bool SymmetricCase;
    bool HalfFilling;
    bool Fixmu0;
    
    int MAX_ITS; 
    
    double Accr;
    
    double AmoebaScanStart;	//before amoeba starts, the equation is solved roughly (with accuracy AmobeScanStep) by scanning from AmoebaScanStart to AmoebaScanEnd.
    double AmoebaScanEnd; 	//make sure AmoebaScanStart and AmoebaScanEnd are far enough apart (when U or W is large).
    double AmoebaScanStep;
    int AmoebaMaxIts;		//maximum number of Amoeba iterations
    
    bool AmoebaForceScanAndPrintOut;	//output n, n0, n-n0 and result when scanning for mu0 candidate
    
    
    //------ OPTIONS -------//
    bool CheckSpectralWeight;   //if true program prints out spectral weights of G and G0 after each iteration
    
    
   
};
