#include "SIAM.h"
#include "routines.h"
#include "Grid.h"
#include <cstring>

//GSL Libraries for Adaptive Cauchy Integration
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

using namespace std;

//---Tail correction---//

double tailfunction(double om, void *params);

struct tailparams
{
	double A,B;
	double mu0,eta;
};

double tailfunction(double om, void *params)
{
  struct tailparams *p= (struct tailparams *)params;
  double A = p->A;
  double B = p->B;
  double mu0 = p->mu0;
  double eta = p->eta;
  return -eta/(pow(om+mu0-A/(om+B),2)+pow(eta,2));
}


//================== Constructors/Destructors ====================//

SIAM::SIAM(const double omega[],size_t N, void * params)
{
	//This loads all necessary parameters into the class object
  p = (struct siamparams *)params;
  
  this -> verbose = p->verbose;
  //Initialize Grid
  this->N = N;
  g = new Grid(omega,N);
  
  //Initialize class internal buffer
  ibuffer = new char[BUFFERSIZE];
  sprintf(ibuffer,"----- Initializing SIAM solver -----\n");
  
	//impurity parameters
  this->U = p->U;
  this->T = p->T;
  this->epsilon = p->epsilon;
  this->mu = p->mu;
	this->mu0 = p->mu0;
	  
  //PH symmetry
  this->SymmetricCase =  p->SymmetricCase;
  this->Fixmu0 =  p->Fixmu0;
  
  this->Accr = p->Accr;
  this->AmoebaScanStart = p->AmoebaScanStart;
  this->AmoebaScanEnd = p->AmoebaScanEnd;
  this->AmoebaScanStep = p->AmoebaScanStep;
  this->AmoebaMaxIts = p->AmoebaMaxIts;
  this->AmoebaForceScanAndPrintOut = p->AmoebaForceScanAndPrintOut;
  
  //Kramers Kronig
  this->KKAccr = p->KKAccr;
  this->usecubicspline =  p->usecubicspline;
  
  //broadening of G0
  this->eta = p->eta;
  ieta = complex<double>(0.0,eta);
  
  //G0 integral tail correction
  this->tailcorrection =  p->tailcorrection;
  this->A1 =  p->A1;
  this->A2 =  p->A2;
  this->B1 =  p->B1;
  this->B2 =  p->B2;

  //options
  this->CheckSpectralWeight = p->CheckSpectralWeight; //default false
}

SIAM::~SIAM()
{
	delete g; //Grid
	delete [] ibuffer; //Buffer
}

//========================= RUN SIAM WITH FIXED Mu ==========================//

int SIAM::Run(const complex<double> Delta[]) //output
{ 
	//Read Delta into grid
	for (int i=0;i<N;i++){
		g->Delta[i] = Delta[i];
	}

	//Initialize some variables
  Clipped = false;
  g->n = 0.5;
  if (SymmetricCase) mu0=0;//-epsilon-U*g->n; //If forcing PH symmetry we have mu0=0
  
	//Obtain the fermi function
  get_fermi();
  
  //print something to class internal buffer
  sprintf(ibuffer + strlen(ibuffer),"SIAM::run::start SIAM solver\n");
  
  sprintf(ibuffer + strlen(ibuffer),"SIAM::run::%s mu=%f U=%f T=%f epsilon=%f eta=%f\n"
  , (SymmetricCase) ? "Symmetric" : "Asymmetric", mu, U, T, epsilon,eta);
  
  //----- initial guess for mu0------// 

  complex<double>* V = new complex<double>[1];
  V[0] = mu0;

  //------ SOLVE SIAM ---------//
  if (SymmetricCase or Fixmu0) {//mu0 and n are known => there's no solving of system of equations
    SolveSiam(V);
  }
  else //Search for mu0
  {
		V[0] = 0.0;
		Amoeba(Accr, V); 
  }
  
  delete [] V;
  //----------------------------//

  //output spectral weight if opted
  if (CheckSpectralWeight)
  {
  	double wG = -imag(TrapezIntegral(N, g->G, g->omega))/M_PI;
  	double wG0 = -imag(TrapezIntegral(N, g->G0, g->omega))/M_PI;
  	
  	if (tailcorrection) { //Correct the spectral weight using tail-correction
  		wG0+=getwG0corr(A1,B1,A2,B2);
  	}
  	
		sprintf(ibuffer + strlen(ibuffer),"SIAM::run::Spectral weight G: %f\n",wG);
		sprintf(ibuffer + strlen(ibuffer),"SIAM::run::Spectral weight G0: %f\n",wG0);
  }
  
  //print occupation
	sprintf(ibuffer + strlen(ibuffer), "SIAM::run::mu=%f\nSIAM::n=%f\n",mu,g->n);

  return 0;
}

//=================================== FUNCTIONS ===================================//

inline double fermi_func(double omega,double T)
{
	return 1.0 / ( 1.0 + exp(omega/T ) );
}

void SIAM::get_fermi()
{  
  #pragma omp parallel for
  for (int i=0; i<N; i++) g->fermi[i] = fermi_func(g->omega[i],T);
}

void SIAM::get_G0()
{
  #pragma omp parallel for
  for (int i=0; i<N; i++) 
    g->G0[i] = complex<double>(1.0)
               / ( complex<double>(g->omega[i] + mu0, eta) - g->Delta[i] ); 
}

double SIAM::get_n(complex<double> X[])
{
  double* dos = new double[N];
  #pragma omp parallel for
  for (int i=0; i<N; i++) 
    dos[i]=-(1/M_PI)*imag(X[i])*g->fermi[i];
  
  double n = TrapezIntegral(N, dos, g->omega);
  delete [] dos;
  return n; 
}

void SIAM::get_As() 
{
  #pragma omp parallel for
  for (int i=0; i<N; i++)
  {
    g->Ap[i] = -imag(g->G0[i]) * fermi_func(g->omega[i],T) / M_PI;
    g->Am[i] = -imag(g->G0[i]) * (1.0 - fermi_func(g->omega[i],T)) / M_PI;
  }
}

double SIAM::get_b()
{ //we used mu0 as (mu0 - epsilon - U*n) in G0, now we're correcting that
  if (!SymmetricCase)
    return ( (1.0 - 2.0 * g->n) * U - mu + (mu0 + epsilon + U * g->n) ) 
           /             ( g->n * (1.0 - g->n) * U*U );
  else return 0;
}

void SIAM::get_Sigma()
{
 
  if (!SymmetricCase)
  { 
    double b = get_b();    
    #pragma omp parallel for
    for (int i=0; i<N; i++) 
      g->Sigma[i] =  U*g->n + g->SOCSigma[i] 
                              / ( 1.0 - b * g->SOCSigma[i] );
    
  }
  else
  { 
    #pragma omp parallel for
    for (int i=0; i<N; i++) 
      g->Sigma[i] =  U * g->n + g->SOCSigma[i];
  }

}

//---------------- Get G -------------------------------//

void SIAM::get_G()
{
  #pragma omp parallel for
  for (int i=0; i<N; i++) 
  {    
    g->G[i] =  1.0
               / (g->omega[i] + mu - epsilon - g->Delta[i] - g->Sigma[i]) ;
  }
  
  //We can't parallelize this due to Clipped being shared amongst the cores - it'll slow things down.
  for (int i=0; i<N; i++) 
  {
    if (ClipOff(g->G[i])) Clipped = true;
  }
  if (Clipped) sprintf(ibuffer + strlen(ibuffer),"SIAM::run::(Warning) !!Clipping G!!\n");
}

//------------------------------------------------------//


double SIAM::getn0corr(double A1,double B1)
{

		sprintf(ibuffer + strlen(ibuffer),"SIAM::run::correcting n0 integral tail using A1=%f B1=%f\n",A1,B1);
		
  	gsl_set_error_handler_off();
    struct tailparams params;
    gsl_function F;
    params.A = A1;
    params.B = B1;
    params.eta = eta;
	  params.mu0 = mu0;
  	F.function = &tailfunction;
  	F.params = &params;
  	
  	
    const double epsabs = 0, epsrel = QUADACCR; // requested errors
    double result; // the integral value
    double error; // the error estimate
    
    size_t limit = QUADLIMIT;// work area size
    gsl_integration_workspace *ws = gsl_integration_workspace_alloc (limit);
    
		int S = gsl_integration_qagil(&F, g->omega[0], epsabs, epsrel, limit, ws, &result, &error);
		
    gsl_integration_workspace_free (ws);
    
    double corr = -result/M_PI;
    
		sprintf(ibuffer + strlen(ibuffer),"SIAM::run::corr = %f\n",corr);
		
    if (corr > 1e-1){
    	corr=0.0;
    	sprintf(ibuffer + strlen(ibuffer),"SIAM::run::(Warning) n0 correction is too large, setting correction to 0\n");
    }
    
    return corr;
}

double SIAM::getwG0corr(double A1,double B1,double A2,double B2)
{
		double corr1,corr2;
		gsl_set_error_handler_off();
		
		sprintf(ibuffer + strlen(ibuffer),"SIAM::run::correcting G0dos integral tail using A1=%f B1=%f A2=%f B2=%f\n",A1,B1,A2,B2);
		
		{
		  struct tailparams params;
		  gsl_function F;
		  params.A = A1;
		  params.B = B1;
		  params.eta = eta;
			params.mu0 = mu0;
			F.function = &tailfunction;
			F.params = &params;
			
			
		  const double epsabs = 0, epsrel = QUADACCR; // requested errors
		  double result; // the integral value
		  double error; // the error estimate
		  
		  size_t limit = QUADLIMIT;// work area size
		  gsl_integration_workspace *ws = gsl_integration_workspace_alloc (limit);
		  
			int S = gsl_integration_qagil(&F, g->omega[0], epsabs, epsrel, limit, ws, &result, &error);
		  gsl_integration_workspace_free (ws);
		  corr1 = -result/M_PI;
		}
		{
		  struct tailparams params2;
		  gsl_function F2;
		  params2.A = A2;
		  params2.B = B2;
		  params2.eta = eta;
			params2.mu0 = mu0;
			F2.function = &tailfunction;
			F2.params = &params2;
			
			
		  const double epsabs = 0, epsrel = QUADACCR; // requested errors
		  double result2; // the integral value
		  double error2; // the error estimate
		  
		  size_t limit = QUADLIMIT;// work area size
		  gsl_integration_workspace *ws2 = gsl_integration_workspace_alloc (limit);
			int S = gsl_integration_qagiu(&F2, g->omega[N-1], epsabs, epsrel, limit, ws2, &result2, &error2);
		  gsl_integration_workspace_free (ws2);
		  corr2 = -result2/M_PI;
		}
				
		sprintf(ibuffer + strlen(ibuffer),"SIAM::run::corr1 = %f corr2 = %f\n",corr1,corr2);
    
    return corr1+corr2;
}


void SIAM::SolveSiam(complex<double>* V)
{
  mu0 = real(V[0]);

  //--------------------//
  get_G0();
  
	
	g->n0 = get_n(g->G0);
	
	if (tailcorrection){
		double n0corr = getn0corr(A1,B1);
  	g->n0 += n0corr;
	}
	
  get_As();
  get_Ps();
  get_SOCSigma();
  
  g->n = g->n0; 
	
  get_Sigma();   
  get_G();
  
  g->n = get_n(g->G);
  //--------------------//

  V[0] = mu0 + (g->n - g->n0); //we need to satisfy (get_n(G) == n)
}

void SIAM::Amoeba(double accr, complex<double>* V)
{
  //x here stands for mu0
  
  double x_start = AmoebaScanStart;
  double x_end   = AmoebaScanEnd;
  double x_step  = AmoebaScanStep;

  int sign_old=0;
  double x;
  bool found = false;
  int try_count = 0;
  double x_candidate;
  double x_best=0, diff_best=10e+100;
  
  sprintf(ibuffer + strlen(ibuffer),"SIAM::run::Amoeba::Start search for mu0\n");
  
  while( (not found) and (try_count<1) ) 
  {
     for(x=x_start; x<x_end; x+=x_step)
     {
       V[0] = x;
       SolveSiam(V);
       
       double x_res=real(V[0]);
       
			 sprintf(ibuffer + strlen(ibuffer),"SIAM::run::Amoeba::mu0: %.15f n(G0)-n(G): %.2le step: %.2le true diff: %.3f n(G): %.3f g->n: %.3f\n",
                              x, x-x_res, x_step, get_n(g->G0)- get_n(g->G),get_n(g->G),g->n);

       if (sign_old==0) 
       { sign_old = ((x-x_res)>=0) ? 1 : -1;//int_sign(x-x_res);
         continue;
       }

       int sign =  ((x-x_res)>=0) ? 1 : -1;//int_sign(x - x_res);
       if (abs(x-x_res) < diff_best) { x_best=x; diff_best = abs(x-x_res); };
       if ((sign_old!=sign) and (not found))
       {  x_candidate = x-x_step;
          found = true; 
          break; 
       }
    }
    
    try_count++;
    
    if (not found) { 
    	x_start *=2.0; x_end *= 2.0; x_step *= 2.0; 
    	sprintf(ibuffer + strlen(ibuffer),"SIAM::run::Amoeba::mu0 candidate NOT found! now scanning a wider range...\n"); }
  } 
 
  
  if (not found)
  {
  
     sprintf(ibuffer + strlen(ibuffer),"SIAM::run::Amoeba::mu0 candidate NOT found! setting mu0 to to best choice: mu0_best: %f diff: %.2le\n",x_best,diff_best);
     
     V[0] = x_best;
     SolveSiam(V);
  }
  else
  {
    sprintf(ibuffer + strlen(ibuffer),"SIAM::run::Amoeba::mu0 candidate found! proceeding with aomeba...\n");  
    
    x = x_candidate;
    x_step *= 0.5;
    x += x_step;

    bool converged = false;
    int it = 0;
    while( (not converged) and (it<=AmoebaMaxIts) )
    {
			it ++;
      V[0] = x;
      SolveSiam(V);
      double x_res=real(V[0]);
      
		  sprintf(ibuffer + strlen(ibuffer),"SIAM::run::Amoeba::it: %d mu0: %.15f n(G0)-n(G): %.2le step: %le n0: %.3f n: %.3f\n", it, x, x-x_res, x_step, get_n(g->G0), get_n(g->G));
		  
      converged = ( abs(x-x_res) < accr );
      int sign =  ((x-x_res)>=0) ? 1 : -1;//int_sign(x - x_res);
      if (sign_old==sign)
         x_step = abs(x_step);
      else
         x_step = -abs(x_step); 
      x_step *= 0.5;
      x += x_step;
    }
    if (converged) sprintf(ibuffer + strlen(ibuffer),"SIAM::run::Amoeba::desired accuracy reached!\n");
  }
  sprintf(ibuffer + strlen(ibuffer),"SIAM::run::Amoeba::--- Amoeba DONE ---\n");
}

//-----------------------Miscellaneous---------------------------------//

bool SIAM::ClipOff(complex<double> &X)
{
  if (imag(X)>0) 
  {
    X = complex<double>(real(X),-ClippingValue);
    return true;
  }
  else
    return false;
}

//--IO--//
void SIAM::PrintFullResult(const char* ResultFN)
{
	g->PrintFullResult(ResultFN);
}

void SIAM::PrintResult()
{
	g->PrintResult("Gf.out","Sig.out");
}

void SIAM::PrintResult(const char* Gffile,const char* Sigfile)
{
	g->PrintResult(Gffile,Sigfile);
}

const char* SIAM::buffer() { return ibuffer; }

void SIAM::PrintBuffer(const char* FN) {

	FILE * flog;
  flog = fopen(FN, "a");
  if (flog == NULL) { 
  	char buff[100]; 
  	snprintf(buff, sizeof(buff), "-- ERROR -- Failed to open file %s : ",FN); 
  	perror(buff); };
  fprintf(flog,ibuffer);
  fclose(flog);
}

void SIAM::PrintBuffer(const char* FN,bool quiet) {

	if (!quiet) printf(ibuffer);
	
	FILE * flog;
  flog = fopen(FN, "a");
  if (flog == NULL) { 
  	char buff[100]; 
  	snprintf(buff, sizeof(buff), "-- ERROR -- Failed to open file %s : ",FN); 
  	perror(buff); };
  fprintf(flog,ibuffer);
  fclose(flog);
}
