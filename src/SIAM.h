//**********************************************************//
//              SIAM at Arbitrary Filling                   //
//                    by Jaksha Vuchichevicc                //
//**********************************************************//

#include <iostream>
#include <complex>

class Grid;
class Logging;
class Params;

using namespace std;

//======================= SIAM Class ==========================================//

class SIAM
{
  private:

    void Defaults();
    string ParamsFN;

    Grid* g;
    Logging* logs;
    
    bool usecubicspline;

    //--impurity parameters--//
    double U;			//on-site repulsion
    double T;			//temperature
    double epsilon;		//impurity energy level
    
		double getn0corr(double A1,double B1);
		double getwG0corr(double A1,double B1,double A2,double B2);
    //---BROADENING---//
    double eta;
    complex<double> ieta;
    
    //--bath parameters--// 
    double mu;			//global chemical potential
    double mu0;			//fictious chemical potential
    
    //----lattice---------//
    double KKAccr;
    
    //--don't touch this---//
    bool SymmetricCase;
    bool HalfFilling;
    bool Fixmu0;
    
    //-- Broyden solver options--//
    double Accr;
    int MAX_ITS; 

    //-- mu0 search --//
    
    //--storage arrays--//
    size_t N;

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

    bool ClipOff(complex<double> &X);
    bool Clipped;

    //--- SIAM solver ---//
    void SolveSiam(complex<double>* V);
    void Amoeba(double accr, complex<double>* V);	//amoeba method for mu0 search. not applicable when MPT corrections are used (TODO: generalize this method)
  
    double AmoebaScanStart;	//before amoeba starts, the equation is solved roughly (with accuracy AmobeScanStep) by scanning from AmoebaScanStart to AmoebaScanEnd.
    double AmoebaScanEnd; 	//make sure AmoebaScanStart and AmoebaScanEnd are far enough apart (when U or W is large).
    double AmoebaScanStep;
    int AmoebaMaxIts;		//maximum number of Amoeba iterations
    bool AmoebaForceScanAndPrintOut;	//output n, n0, n-n0 and result when scanning for mu0 candidate
    
  public:
  
    //--Constructors/destructors--//
    SIAM();  
    SIAM(Params* params);
    ~SIAM();
  
    //------ OPTIONS -------//
    bool CheckSpectralWeight;   //if true program prints out spectral weights of G and G0 after each iteration
    
    double get_T(){return T;};
    double get_mu(){return mu;};
    double get_epsilon(){return epsilon;};
    
    //Tail
    double A1,B1,A2,B2;
    bool tailcorrection;
    
    
    //--------RUN SIAM--------//
    bool Run(Grid* r,Logging* logs); 
   
};
