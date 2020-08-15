/*
	IPT Anderson Impurity Solver CPU+GPU code
	by Pak Ki Henry Tsang
	Adapted from code by Jaksha Vuchichevicc https://github.com/JaksaVucicevic/DMFT
	
	In CPU mode this code is linked against SIAM.cpu.cpp
	In GPU mode this code is linked against SIAM.cu
*/

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
	double * LR;
	double mu0,eta;
	int fitorder;
};

double tailfunction(double om, void *params)
{
  struct tailparams *p= (struct tailparams *)params;
  double mu0 = p->mu0;
  double eta = p->eta;
  double reDelta = 0.0;
  int fitorder = p->fitorder;
  for (int i=0;i<fitorder;i++) reDelta += p->LR[i]/pow(om,2*i+1);
  return -eta/(pow(om+mu0-reDelta,2) +pow(eta,2));
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
  this->L = p->L;
  this->R = p->R;
  this->fitorder = p->fitorder;

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

  double* V = new double[1];
  V[0] = mu0;

  //------ SOLVE SIAM ---------//
  if (SymmetricCase or Fixmu0) {//mu0 and n are known => there's no solving of system of equations
    SolveSiam(V);
  }
  else //Search for mu0
  {
		V[0] = 0.0;
		//Amoeba(Accr, V); 
		Hybrid_Bisection(Accr,V);
  }
  
  delete [] V;
  //----------------------------//

  //output spectral weight if opted
  if (CheckSpectralWeight)
  {
  	double wG = -imag(TrapezIntegral(N, g->G, g->omega))/M_PI;
  	double wG0 = -imag(TrapezIntegral(N, g->G0, g->omega))/M_PI;
  	
  	if (tailcorrection) { //Correct the spectral weight using tail-correction
  		wG0+=getwG0corr();
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


double SIAM::getn0corr()
{

		if (verbose) sprintf(ibuffer + strlen(ibuffer),"SIAM::run::correcting n0 integral tail\n");
		//for (int i=0;i<fitorder;i++) sprintf(ibuffer + strlen(ibuffer),"SIAM::run::L[%d] = %f\n",i,L[i]);
		
  	gsl_set_error_handler_off();
    struct tailparams params;
    gsl_function F;
    params.LR = L;
    params.eta = eta;
	  params.mu0 = mu0;
	  params.fitorder = fitorder;
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
    
		if (verbose) sprintf(ibuffer + strlen(ibuffer),"SIAM::run::corr = %f\n",corr);
		
    if (corr > 1e-1){
    	corr=0.0;
    	sprintf(ibuffer + strlen(ibuffer),"SIAM::run::(Warning) n0 correction is too large, setting correction to 0\n");
    }
    
    return corr;
}

double SIAM::getwG0corr()
{
		double corr1,corr2;
		gsl_set_error_handler_off();
		
		if (verbose) sprintf(ibuffer + strlen(ibuffer),"SIAM::run::correcting G0dos integral tail\n");
		for (int i=0;i<fitorder;i++) sprintf(ibuffer + strlen(ibuffer),"SIAM::run::L[%d] = %f\n",i,L[i]);
		for (int i=0;i<fitorder;i++) sprintf(ibuffer + strlen(ibuffer),"SIAM::run::R[%d] = %f\n",i,R[i]);
		
		{
		  struct tailparams params;
		  gsl_function F;
		  params.LR = L;
		  params.eta = eta;
			params.mu0 = mu0;
	  	params.fitorder = fitorder;
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
		  params2.LR = R;
		  params2.eta = eta;
			params2.mu0 = mu0;
	  	params2.fitorder = fitorder;
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
				
		if (verbose) sprintf(ibuffer + strlen(ibuffer),"SIAM::run::corr1 = %f corr2 = %f\n",corr1,corr2);
    
    return corr1+corr2;
}


void SIAM::SolveSiam(double* V)
{
  mu0 = V[0];

  //--------------------//
  get_G0();
  
	
	g->n0 = get_n(g->G0);
	
	if (tailcorrection){
		double n0corr = getn0corr();
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

void SIAM::Hybrid_Bisection(double accr, double* V)
{
  double mu0L = AmoebaScanStart;
  double mu0R   = AmoebaScanEnd;
  int max_tries = AmoebaMaxIts;
  int it=0;
  
  double dmu0;
  double mu0M;
  
	double x,x_prev;
	double x_res;
  double diff(0.0),diff_prev(0.0);
  
  bool bisection = true;
  bool converged = false;
  
  sprintf(ibuffer + strlen(ibuffer),"SIAM::run::Hybrid_Bisection::Start search for mu0\n");
  
  mu0M = mu0R-mu0L;
  
  while (it<max_tries and !converged){
  
  	dmu0 = mu0R-mu0L;
  	mu0M = (mu0R+mu0L)*0.5;
  	
		if (bisection){ mu0 = mu0M; }
		else{ 
		mu0 = mu0-diff*(x_prev-x)/(diff_prev-diff);
		if (mu0>mu0R or mu0<mu0L) {mu0 = mu0M; bisection=true;}
		}

		V[0] = mu0 ;
		
		SolveSiam(V);
		
		x_prev = x;
		diff_prev= diff;
		
		x=mu0;
		x_res=V[0];
    diff = x-x_res;
    
	 	sprintf(ibuffer + strlen(ibuffer),"SIAM::run::Hybrid_Bisection::%s::it: %d mu0: %.15f n(G0)-n(G): %.3le dmu0: %.3le n: %.3le n0: %.3le\n",
                              (bisection) ? "bisection" : "interpolate",it,x, diff, dmu0,g->n,g->n0);
                              
    if (it!=0) bisection = !bisection; //Interchange method
    
  	if (abs(diff)<accr){ converged=true; }
  	else{ //shrink region
  		mu0L = (diff<0) ? mu0 : mu0L;
  		mu0R = (diff<0) ? mu0R : mu0;
  	}
  	it++;
	}
	
  if (converged) sprintf(ibuffer + strlen(ibuffer),"SIAM::run::Hybrid_Bisection::desired accuracy reached!\n");
  sprintf(ibuffer + strlen(ibuffer),"SIAM::run::Hybrid_Bisection::--- Hybrid_Bisection DONE ---\n");
}

void SIAM::Amoeba(double accr, double* V)
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
       
       double x_res=V[0];
       
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
      double x_res=V[0];
      
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
  fprintf(flog,"%s",ibuffer);
  fclose(flog);
}

void SIAM::PrintBuffer(const char* FN,bool quiet) {

	if (!quiet) printf("%s",ibuffer);
	
	FILE * flog;
  flog = fopen(FN, "a");
  if (flog == NULL) { 
  	char buff[100]; 
  	snprintf(buff, sizeof(buff), "-- ERROR -- Failed to open file %s : ",FN); 
  	perror(buff); };
  fprintf(flog,"%s",ibuffer);
  fclose(flog);
}
