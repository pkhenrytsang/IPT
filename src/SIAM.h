//**********************************************************//
//              SIAM at Arbitrary Filling                   //
//           original by Jaksha Vuchichevicc                //
//           rewritten by Pak Ki Henry Tsang                //
//**********************************************************//

#include <iostream>
#include <complex>

//Buffer size for printouts (not sure if enough!)
#define BUFFERSIZE 8192

//Maximum subdivision for adaptive quadrature
#define QUADLIMIT 200

//epsrel for adaptive quadrature
#define QUADACCR 1e-9

//
#define CLIPVAL 1e-5


class Grid;

using namespace std;

//Siam params

struct siamparams{
	double mu,mu0,U,T,epsilon,eta;
	double KKAccr,Accr;
	double A1,A2,B1,B2;
	bool SymmetricCase,Fixmu0,usecubicspline,tailcorrection;
	
	bool CheckSpectralWeight;
	bool verbose;
	int AmoebaMaxIts;
	double AmoebaScanStart,AmoebaScanEnd,AmoebaScanStep;
	bool AmoebaForceScanAndPrintOut;
};

//======================= SIAM Class ==========================================//

class SIAM
{
  public:
  
    //--Constructors/destructors--//
    SIAM(const double omega[],size_t N, void * params);
    ~SIAM();
    
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
  
  	//verbose output (currently not useful as we print all output)
  	bool verbose;
  	
  	//print buffer
  	char* ibuffer;
  
  	//parameters
  	struct siamparams *p;
  	
    //--Grid--//
    Grid* g;
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
    double A1,B1,A2,B2;
    
		double getn0corr(double A1,double B1);
		double getwG0corr(double A1,double B1,double A2,double B2);

    //----kramers kronig----//
    double KKAccr;
    bool usecubicspline;
    
     //--get functions--//
    double get_fermi(int i);
    double get_n(complex<double> X[]);

    //--get procedures--//
    void get_fermi();
    void get_G0();
    void get_As();
    void get_Ps();  
    void get_SOCSigma();
    double get_b();
    void get_Sigma();
    void get_G();

		//--Clipping of ImSigma and ImG--//
    bool ClipOff(complex<double> &X);
    bool Clipped;
    double ClippingValue = CLIPVAL;

    //--- SIAM solver ---//
    void SolveSiam(complex<double>* V);
    void Amoeba(double accr, complex<double>* V);	//amoeba method for mu0 search.
  
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
