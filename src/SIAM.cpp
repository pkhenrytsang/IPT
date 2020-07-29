#include "SIAM.h"
#include "routines.h"
#include "Grid.h"
#include "Params.h"
#include "log.h"
#include <omp.h>
//GSL Libraries for Adaptive Cauchy Integration
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

//================== Constructors/DEstructors ====================//

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

void SIAM::Defaults()
{
  U = 2.0;
  
  T = 0.05;
  
  epsilon = 0;
  
  SymmetricCase = true;
  Fixmu0 = false;
  
  usecubicspline = false;
  
  tailcorrection = false;
  
  A1=0.0;
  A2=0.0;
  B1=1.0;
  B2=1.0;
  
  KKAccr = 1e-9;

  //broadening
  eta = 1e-3;
   
  //options
  CheckSpectralWeight = false; //default false
  
  Accr = 1e-9;
  AmoebaScanStart = -1;
  AmoebaScanEnd = 1;
  AmoebaScanStep = 0.2;
  AmoebaMaxIts = 20;
  AmoebaForceScanAndPrintOut = false;
}

SIAM::SIAM()
{
  Defaults();
}

SIAM::SIAM(Params* params)
{
	
  Defaults();
  //this->ParamsFN.assign(ParamsFN);
  //cout << "-- INFO -- SIAM: Params File Name set to:" << this->ParamsFN << endl;

  //Params params(ParamsFN);
	//this->params = 
  params->ReadParam(U,"'U'");
  params->ReadParam(T,"'T'");
  params->ReadParam(epsilon,"'epsilon'");
  params->ReadParam(eta,"'eta'");
  params->ReadParam(CheckSpectralWeight, "'CheckSpectralWeight'");
  
  //Away from half filling
  params->ReadParam(SymmetricCase,"'SymmetricCase'"); //Force half-filling?
  params->ReadParam(Fixmu0,"'Fixmu0'"); //Force half-filling?
  params->ReadParam(usecubicspline,"'UseCubicSpline'");
  params->ReadParam(Accr,"'Accr'");
  params->ReadParam(KKAccr,"'KKAccr'");
  params->ReadParam(AmoebaScanStart,"'AmoebaScanStart'");
  params->ReadParam(AmoebaScanEnd,"'AmoebaScanEnd'");
  params->ReadParam(AmoebaScanStep,"'AmoebaScanStep'");
  params->ReadParam(AmoebaMaxIts,"'AmoebaMaxIts'");
  params->ReadParam(AmoebaForceScanAndPrintOut,"'AmoebaForceScanAndPrintOut'");
  
  ieta = complex<double>(0.0,eta);
}

SIAM::~SIAM()
{}

//========================= RUN SIAM WITH FIXED Mu ==========================//

bool SIAM::Run(Grid* g,Logging* logs) //output
{  
  this->g = g;
  this->logs = logs;
  N = g->N;
  get_fermi();

  Clipped = false;
  
  mu = g->mu;
  mu0 = g->mu0;
  g->n = 0.5;//default
  if (SymmetricCase) mu0=0;
  
  logs->print("SIAM::start SIAM solver\n");
  {
		char buff[200];
		snprintf(buff, sizeof(buff), "SIAM::%s\nSIAM::mu=%f\nSIAM::U=%f\nSIAM::T=%f\nSIAM::epsilon=%f\nSIAM::eta=%f\n", (SymmetricCase) ? "Symmetric" : "Asymmetric", g->mu, U, T, epsilon,eta);
		logs->print(buff);
  }
  
  //----- initial guess ------// 

  complex<double>* V = new complex<double>[1];
  V[0] = mu0;

  //------ SOLVE SIAM ---------//
  if (SymmetricCase or Fixmu0) {//mu0 and n are known => there's no solving of system of equations
    SolveSiam(V);
  }
  else 
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
  	
  	if (tailcorrection) wG0+=getwG0corr(A1,B1,A2,B2);
  	//printf("A1=%f B1=%f A2=%f B2=%f\n",A1,B1,A2,B2);
  	
		char buff1[100];
		char buff2[100];
		snprintf(buff1, sizeof(buff1),"SIAM::Spectral weight G: %f\n",wG);
		snprintf(buff2, sizeof(buff2),"SIAM::Spectral weight G0: %f\n",wG0);
		logs->print(buff1);
		logs->print(buff2);
  }

  g->mu0 = mu0;
  {
		char buff[100];
		snprintf(buff, sizeof(buff), "SIAM::mu=%f\nSIAM::n=%f\n",g->mu,g->n);
		logs->print(buff);
  }

	//Calculates the DOS
  #pragma omp parallel for
  for (int i=0; i<N; i++)
    g->DOS[i] = - imag(g->G[i]) / M_PI;

  return false;
}

//=================================== FUNCTIONS ===================================//

void SIAM::get_fermi()
{  
  #pragma omp parallel for
  for (int i=0; i<N; i++) g->fermi[i] = 1.0 / ( 1.0 + exp( g->omega[i]/T ) );
}

inline double fermi_func(double omega,double T)
{
	return 1.0 / ( 1.0 + exp( omega/T ) );
}


void SIAM::get_G0()
{
  #pragma omp parallel for
  for (int i=0; i<N; i++) 
    g->G0[i] = complex<double>(1.0)
               / ( complex<double>(g->omega[i] + mu0, eta)
                   - g->Delta[i] ); 
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
  //printf("         SIAM::get_b : MPT_B = %.3f, MPT_B0 = %.3f\n");  
  if (!SymmetricCase)
    return ( (1.0 - 2.0 * g->n) * U - g->mu + (mu0 + epsilon + U * g->n) ) 
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
  //#pragma omp parallel for
  for (int i=0; i<N; i++) 
  {    
    g->G[i] =  1.0
               / (g->omega[i] + g->mu - epsilon - g->Delta[i] - g->Sigma[i]) ;
    
    if (ClipOff(g->G[i])) Clipped = true;
  }
  
  if (Clipped) logs->print("SIAM::!!!!Clipping G!!!!\n");
}

//------------------------------------------------------//



double SIAM::getn0corr(double A1,double B1)
{
  	gsl_set_error_handler_off();
    struct tailparams params;
    gsl_function F;
    params.A = A1;
    params.B = B1;
    params.eta = eta;
	  params.mu0 = mu0;
  	F.function = &tailfunction;
  	F.params = &params;
  	
  	
    const double epsabs = 0, epsrel = Accr; // requested errors
    double result; // the integral value
    double error; // the error estimate
    
    size_t limit = 200;// work area size
    gsl_integration_workspace *ws = gsl_integration_workspace_alloc (limit);
    
		int S = gsl_integration_qagil(&F, g->omega[0], epsabs, epsrel, limit, ws, &result, &error);
		//logs->print("SIAM::tail integration for n0\n");
    gsl_integration_workspace_free (ws);
    
    double corr = -result/M_PI;
    
    if (corr > 1e-2) corr=0.0;
    
    return corr;
}

double SIAM::getwG0corr(double A1,double B1,double A2,double B2)
{
		double corr1,corr2;
		gsl_set_error_handler_off();
		{
		  struct tailparams params;
		  gsl_function F;
		  params.A = A1;
		  params.B = B1;
		  params.eta = eta;
			params.mu0 = mu0;
			F.function = &tailfunction;
			F.params = &params;
			
			
		  const double epsabs = 0, epsrel = Accr; // requested errors
		  double result; // the integral value
		  double error; // the error estimate
		  
		  size_t limit = 400;// work area size
		  gsl_integration_workspace *ws = gsl_integration_workspace_alloc (limit);
		  
			int S = gsl_integration_qagil(&F, g->omega[0], epsabs, epsrel, limit, ws, &result, &error);
			//logs->print("SIAM::tail integration for n0\n");
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
			
			
		  const double epsabs = 0, epsrel = Accr; // requested errors
		  double result2; // the integral value
		  double error2; // the error estimate
		  
		  size_t limit = 400;// work area size
		  gsl_integration_workspace *ws2 = gsl_integration_workspace_alloc (limit);
			int S = gsl_integration_qagiu(&F2, g->omega[N-1], epsabs, epsrel, limit, ws2, &result2, &error2);
			//logs->print("SIAM::tail integration for n0\n");
		  gsl_integration_workspace_free (ws2);
		  corr2 = -result2/M_PI;
		}
		
		//printf("corr1 = %f corr2 = %f\n",corr1,corr2);
    
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
	
	char buff[100];
	snprintf(buff, sizeof(buff),"n0corr = %f\n",n0corr);
	logs->print(buff);
	
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
  logs->print("SIAM::Start search for mu0\n");
  while( (not found) and (try_count<1) ) 
  {
     FILE* ScanFile;  
     if (AmoebaForceScanAndPrintOut)
     {
       char ScanFN[50];
       sprintf(ScanFN, "scan.eps%.3f",epsilon);
       ScanFile = fopen(ScanFN,"w");
     }
     for(x=x_start; x<x_end; x+=x_step)
     {
       V[0] = x;
       SolveSiam(V);
       
       if (AmoebaForceScanAndPrintOut)
       {  char FN[50];
          sprintf(FN,"siam.eps%.3f.mu0_%.3f",epsilon, mu0);
          g->PrintResult(FN);
       }
      
       double x_res=real(V[0]);
       {
				 char buff[200];
				 snprintf(buff, sizeof(buff), "SIAM::Amoeba::mu0: %.15f n(G0)-n(G): %.2le step: %.2le true diff: %.3f n(G): %.3f g->n: %.3f\n",
                              x, x-x_res, x_step, get_n(g->G0)- get_n(g->G),get_n(g->G),g->n);
				 logs->print(buff);
			 }
       //printf("         mu0: %.15f n(G0)-n(G): %.2le step: %.2le true diff: %.3f n(G): %.3f g->n: %.3f\n",
       //                       x, x-x_res, x_step, get_n(g->G0)- get_n(g->G),get_n(g->G),g->n);

       if (AmoebaForceScanAndPrintOut)
         fprintf(ScanFile,"%.15le %.15le %.15le %.15le\n", x, x-x_res, g->n, g->n0);

       if (sign_old==0) 
       { sign_old = ((x-x_res)>=0) ? 1 : -1;//int_sign(x-x_res);
         continue;
       }

       int sign =  ((x-x_res)>=0) ? 1 : -1;//int_sign(x - x_res);
       if (abs(x-x_res) < diff_best) { x_best=x; diff_best = abs(x-x_res); };
       if ((sign_old!=sign) and (not found))
       {  x_candidate = x-x_step;
          found = true; 
          if (not AmoebaForceScanAndPrintOut) break; 
       }
    }
    try_count++;
    if (not found) { x_start *=2.0; x_end *= 2.0; x_step *= 2.0; logs->print("SIAM::Amoeba::mu0 candidate NOT found! now scanning a wider range...\n"); }
    if (AmoebaForceScanAndPrintOut) fclose(ScanFile);
  } 
 
  
  if (not found)
  {
		 char buff[100];
     snprintf(buff, sizeof(buff), "SIAM::Amoeba::mu0 candidate NOT found! setting mu0 to to best choice: mu0_best: %f diff: %.2le\n",x_best,diff_best);
		 logs->print(buff);
     V[0] = x_best;
     SolveSiam(V);
  }
  else
  {
    logs->print("SIAM::Amoeba::mu0 candidate found! proceeding with aomeba...\n");  
    x = x_candidate;
    x_step *= 0.5;
    x += x_step;

    bool converged = false;
    int it = 0;
    while( (not converged) and (it<=AmoebaMaxIts) )
    { it ++;
      V[0] = x;
      SolveSiam(V);
      double x_res=real(V[0]);
      {
		  char buff[200];
		  snprintf(buff, sizeof(buff), "SIAM::Amoeba::it: %d mu0: %.15f n(G0)-n(G): %.2le step: %le n0: %.3f n: %.3f\n", it, x, x-x_res, x_step, get_n(g->G0), get_n(g->G));
		  logs->print(buff);
		  }
      //printf("         it: %d mu0: %.15f n(G0)-n(G): %.2le step: %le n0: %.3f n: %.3f\n", it, x, x-x_res, x_step, get_n(g->G0), get_n(g->G));
      converged = ( abs(x-x_res) < accr );
      int sign =  ((x-x_res)>=0) ? 1 : -1;//int_sign(x - x_res);
      if (sign_old==sign)
         x_step = abs(x_step);
      else
         x_step = -abs(x_step); 
      x_step *= 0.5;
      x += x_step;
    }
    if (converged) logs->print("SIAM::Amoeba::Amoeba: desired accuracy reached!\n");
  }
  logs->print("SIAM::Amoeba::--- Amoeba DONE ---\n");
}

//-----------------------Miscellaneous---------------------------------//

bool SIAM::ClipOff(complex<double> &X)
{
  if (imag(X)>0) 
  {
    X = complex<double>(real(X),-1e-5);
    return true;
  }
  else
    return false;
}

