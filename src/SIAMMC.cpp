/*
	IPT Anderson Impurity Solver Monte-Carlo code
	by Pak Ki Henry Tsang
	adapted from code by Jaksha Vuchichevicc https://github.com/JaksaVucicevic/DMFT
	
	This code is designed to compute 
*/

#include "SIAMMC.h"
#include "GridMC.h"
#include "dinterpl.h"
#include <ctime>
#include <mpi.h>
#include <stdlib.h>
#include "mkl_vsl.h"

//GSL Libraries for Adaptive Cauchy Integration
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#define BRNG VSL_BRNG_MCG31
//#define BRNG VSL_BRNG_SFMT19937

// Struct


struct tailparams
{
	double * LR;
	double mu0,eta;
	int fitorder;
};



struct qparams
    {
      double* omega;
      double* Y;
      double T;
      int N;
      dinterpl* spline;
    };
    
// inline function

inline double fermi_func(double omega,double T){
	return 1/ ( 1.0 + exp(omega/T ) );
}

// Integrands

double imSOCSigmafc(double om, void *params);
double imSOCSigmafl(double om, void *params);

double nfc(double om, void *params);
double nfl(double om, void *params);

double imSOCSigmafc(double om, void *params)
{
  struct qparams *p= (struct qparams *)params;
  return p->spline->cspline_eval(om);
}

double imSOCSigmafl(double om, void *params)
{
  struct qparams *p= (struct qparams *)params;
  //printf("om: %f , N : %d\n",om,p->N);
  //double val = dinterpl::linear_eval(om,p->omega, p->Y , p->N);
  //printf("%f\n",val);
  //return val;
  return dinterpl::linear_eval(om,p->omega, p->Y , p->N);
}

double nfl(double om, void *params)
{
  struct qparams *p= (struct qparams *)params;
  return dinterpl::linear_eval(om,p->omega, p->Y , p->N)*fermi_func(om,p->T);
}

double nfc(double om, void *params)
{
  struct qparams *p= (struct qparams *)params;
  return p->spline->cspline_eval(om)*fermi_func(om,p->T);
}

double tailfunction(double om, void *params);

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

SIAMMC::SIAMMC(const double omega[],size_t N, void * params)
{
	//This loads all necessary parameters into the class object
  p = (struct siammcparams *)params;
  
  this -> verbose = p->verbose;
  //Initialize Grid
  this->N = N;
  g = new GridMC(omega,N);
  
  //Integration
  this->W = p->W;
  this->M = p->M;
  
  //Initialize class internal buffer
  ibuffer = new char[BUFFERSIZE];
  sprintf(ibuffer,"----- Initializing SIAMMC solver -----\n");
  
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

SIAMMC::~SIAMMC()
{
	delete g; //Grid
	delete [] ibuffer; //Buffer
}


//========================= RUN SIAM WITH FIXED Mu ==========================//

int SIAMMC::Run(const complex<double> Delta[]) //output
{ 
	//MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	//printf("SIAM from rank %d\n",rank);

	//Make a copy of Delta
	if (!SymmetricCase){
		for (int i=0;i<N;i++) g->Delta[i] = Delta[i];
	}
	else{ //Enforce symmetry
		for (int i=0;i<N/2;i++) g->Delta[i] = complex<double>(-real(Delta[N-1-i]),imag(Delta[N-1-i]));
		for (int i=0;i<N/2;i++) g->Delta[i+N/2] = Delta[i+N/2];
	}
	
	//printf("Finish copying Delta from rank %d\n",rank);
	//this->Delta = Delta;

	//Initialize some variables
  Clipped = false;
  g->n = 0.5;
  if (SymmetricCase) mu0=0;//-epsilon-U*g->n; //If forcing PH symmetry we have mu0=0
  
  //print something to class internal buffer
  sprintf(ibuffer + strlen(ibuffer),"SIAMMC::run::start SIAMMC solver from rank %d\n",rank);
  
  sprintf(ibuffer + strlen(ibuffer),"SIAMMC::run::%s mu=%f U=%f T=%f epsilon=%f eta=%f\n"
  , (SymmetricCase) ? "Symmetric" : "Asymmetric", mu, U, T, epsilon,eta);
  
  //----- initial guess for mu0------// 

  double * V = new double[1];
  V[0] = mu0;

  //------ SOLVE SIAM ---------//
  if (SymmetricCase or Fixmu0) {//mu0 and n are known => there's no solving of system of equations
  	//printf("Solving SIAM from rank %d\n",rank);
    SolveSiam(V);
  	//printf("Finished Solving SIAM from rank %d\n",rank);
  }
  else //Search for mu0
  {
		V[0] = 0.0;
		Amoeba(Accr, V);
  }
  
  delete [] V;
  //----------------------------//

  //output spectral weight if opted
  /*
	if (CheckSpectralWeight)
	{
		double wG = -imag(TrapezIntegral(N, g->G, g->omega))/M_PI;
		double wG0 = -imag(TrapezIntegral(N, g->G0, g->omega))/M_PI;
		
		if (tailcorrection) { //Correct the spectral weight using tail-correction
			wG0+=getwG0corr();
		}
		
		sprintf(ibuffer + strlen(ibuffer),"SIAMMC::run::Spectral weight G: %f\n",wG);
		sprintf(ibuffer + strlen(ibuffer),"SIAMMC::run::Spectral weight G0: %f\n",wG0);
	}
	*/
	//print occupation
	sprintf(ibuffer + strlen(ibuffer), "SIAMMC::run::mu=%f\nSIAMMC::n=%f\n",mu,g->n);
		 
  return 0;
}

//=================================== Self energy ===================================//



void SIAMMC::get_dos0()
{	//Obtain eta broadened imG0 for interpolation
  for (int i=0; i<N; i++) {
    complex<double> G0 = 1/( complex<double>(g->omega[i] + mu0, eta) - g->Delta[i] ); 
    g->dos0[i] = -imag(G0)/M_PI;
  }
}

inline double SIAMMC::Ap(double om){ //inline?
	return dinterpl::linear_eval(om,g->omega,g->dos0,N)*fermi_func(om,T);}

inline double SIAMMC::Am(double om){ //inline?
	return dinterpl::linear_eval(om,g->omega,g->dos0,N)*(1-fermi_func(om,T));}
	
	
inline double SIAMMC::Apc(double om){ //inline?
	return dos0spline->cspline_eval(om)*fermi_func(om,T);}

inline double SIAMMC::Amc(double om){ //inline?
	return dos0spline->cspline_eval(om)*(1-fermi_func(om,T));}

void SIAMMC::get_SOCSigma(){
  double *imSOCSigmaS = new double[N]; //Slaves
  double *reSOCSigmaS = new double[N]; //Slaves
  
  double *imSOCSigma = new double[N]; //Result
  double *reSOCSigma = new double[N]; //Result
  
  double Vol = 4*W*W; //Integration Volume
	
	//Begin Monte Carlo integration
	//Right now we do just uniform distribution
  for (int i=0;i<N;i++) {
		imSOCSigmaS[i]=0.0;
		reSOCSigmaS[i]=0.0;
		imSOCSigma[i]=0.0;
		reSOCSigma[i]=0.0;
  }//Set zero
  
	
	double * w1s = new double[M];
	double * w2s = new double[M];
	double * sigma_local = new double[N];
	double * sigma_global = new double[N];
	double ** fw = new double*[N];
	for (int i=0;i<N;i++) fw[i] = new double[M];
	
	VSLStreamStatePtr stream;
	
	vslNewStream( &stream, BRNG, time(NULL) + rank );
	vdRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream, M, w1s, -W, W );
	vslDeleteStream( &stream );
	
	vslNewStream( &stream, BRNG, time(NULL) * rank );
	vdRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream, M, w2s, -W, W );  
	vslDeleteStream( &stream );
  
  
  /*
  	This particular loop order j->i is the most optimum as we have
  	(2M+2M)+(2NM) = 4M + 2NM functional evaluations of Ap/Am
  	
  	The other way round i->j one evaluate 
  	
  	N(2M+2M+2M) = 8NM functional evaluations, which is ~4 times as much evaluations for large N!! 
  	
  	We use the same set of random p1 and p2 for all N omega, otherwise computational cost will be
  	6MN which is ~3 times as long for large N!
  */
  
  if (!usecubicspline){
		for (int j=0;j<M;j++){ //M random points
			double p1 = Ap(w2s[j])*Am(w2s[j]-w1s[j]); //2M functional evaluations
			double p2 = Am(w2s[j])*Ap(w2s[j]-w1s[j]); //2M functional evaluations
			#pragma ivdep
			for (int i=0;i<N;i++){
				double fint = (Ap(g->omega[i]-w1s[j])*p2 + Am(g->omega[i]-w1s[j])*p1); //2N functional evaluations
				fw[i][j] = fint;
				imSOCSigmaS[i]+=fint;
			}
		}
  }
  else{
		for (int j=0;j<M;j++){ //M random points
			double p1 = Ap(w2s[j])*Am(w2s[j]-w1s[j]);
			double p2 = Am(w2s[j])*Ap(w2s[j]-w1s[j]);
			#pragma ivdep
			for (int i=0;i<N;i++){
				double fint = (Apc(g->omega[i]-w1s[j])*p2 + Amc(g->omega[i]-w1s[j])*p1);
				fw[i][j] = fint;
				imSOCSigmaS[i]+=fint;
			}
		}
  }

  //Collect results through MPI interface
  MPI_Allreduce(imSOCSigmaS, imSOCSigma, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	
  for (int i=0;i<N;i++){
  	#pragma ivdep
		for (int j=0;j<M;j++){
			sigma_local[i] += (fw[i][j]-imSOCSigma[i])*(fw[i][j]-imSOCSigma[i]);	
		}
		sigma_local[i] = sigma_local[i]/((double) M * size - 1);
  }	
	
  MPI_Allreduce(sigma_local, sigma_global, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
	for (int i=0;i<N;i++) imSOCSigma[i] = -U*U*imSOCSigma[i]*M_PI*Vol/( (double) M * size );
	for (int i=0;i<N;i++) sigma_global[i] = pow(-U*U*M_PI*Vol/( (double) M * size ),2)*sigma_global[i];
	
	
  //if (rank==0) printf("imSOCSigma[N/2] = %f , sigma_global[N/2] = %f\n",imSOCSigma[N/2],sigma_global[N/2]);
  
	delete [] w1s;
  delete [] w2s;
  delete [] sigma_local;
  delete [] sigma_global;
  for (int i=0;i<N;i++) delete [] fw[i];
  delete [] fw;
	
	if (SymmetricCase){	
		for (int i=0;i<N/2;i++) {
			imSOCSigma[i] = (imSOCSigma[N-1-i]+imSOCSigma[i])*0.5; 
		}
		for (int i=0;i<N/2;i++) {
			imSOCSigma[N-1-i] = imSOCSigma[i]; 
		}
	}
	
	gsl_set_error_handler_off();
	
	if (usecubicspline) imSOCSigmaspline = new dinterpl(g->omega, imSOCSigma , N);
	
	//Run Kramers Kronig using MPI by splitting integration region
	for (int i=1; i<N-1; i++)
	{ 
		if (i%size != rank) continue;
		
	  const double a = g->omega[0], b = g->omega[N-1]; // limits of integration
	  const double epsabs = 0, epsrel = KKAccr; // requested errors
	  double result; // the integral value
	  double error; // the error estimate

	  double c = g->omega[i];

	  struct qparams params;
	  gsl_function F;
	  
    if (!usecubicspline){
		  params.omega = g->omega;
		  params.Y = imSOCSigma;
		  params.N = N;
    	F.function = &imSOCSigmafl;
    }
    else{	
    	params.spline = imSOCSigmaspline;
    	F.function = &imSOCSigmafc;
    }
    F.params = &params;
	  
	  
	  size_t limit = QUADLIMIT;// work area size
	  gsl_integration_workspace *ws = gsl_integration_workspace_alloc (limit);

	  gsl_integration_qawc (&F, a, b , c , epsabs, epsrel, limit, ws, &result, &error);
	  
		//printf("Finished integration from rank %d\n",rank);

	  gsl_integration_workspace_free (ws);

		reSOCSigmaS[i] = result/M_PI;
	}
	//Collect Results from MPI workers
	MPI_Allreduce(reSOCSigmaS, reSOCSigma, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	
	if (SymmetricCase){	
		for (int i=0;i<N/2;i++) {
			reSOCSigma[i] = (-reSOCSigma[N-1-i]+reSOCSigma[i])*0.5; 
		}
		for (int i=0;i<N/2;i++) {
			reSOCSigma[N-1-i] = -reSOCSigma[i]; 
		}
	}
	
	//Store as complex number array
	for (int i=0;i<N;i++) g->SOCSigma[i] = complex<double>(reSOCSigma[i],imSOCSigma[i]);
	//End points
	g->SOCSigma[0] = g->SOCSigma[1];
	g->SOCSigma[N-1] = g->SOCSigma[N-2];
	delete [] imSOCSigma;
  delete [] reSOCSigma;
  delete [] imSOCSigmaS;
  delete [] reSOCSigmaS;
  delete imSOCSigmaspline;
}

double SIAMMC::get_b()
{ //we used mu0 as (mu0 - epsilon - U*n) in G0, now we're correcting that
  if (!SymmetricCase)
    return ( (1.0 - 2.0 * g->n) * U - mu + (mu0 + epsilon + U * g->n) ) 
           /             ( g->n * (1.0 - g->n) * U*U );
  else return 0;
}

void SIAMMC::get_Sigma()
{
 
  if (!SymmetricCase)
  { 
    double b = get_b();    
    for (int i=0; i<N; i++) 
      g->Sigma[i] =  U*g->n + g->SOCSigma[i] 
                              / ( 1.0 - b * g->SOCSigma[i] );
    
  }
  else
  {
    for (int i=0; i<N; i++) 
      g->Sigma[i] =  U * g->n + g->SOCSigma[i];
  }

}

double SIAMMC::get_n(double dos[])
{

	gsl_set_error_handler_off();
	
	const double a = g->omega[0], b = g->omega[N-1]; // limits of integration
  const double epsabs = 0, epsrel = KKAccr; // requested errors
  double result; // the integral value
  double error; // the error estimate

  struct qparams params;
  gsl_function F;
  
  if (usecubicspline) dosspline = new dinterpl(g->omega, dos, N);
  if (!usecubicspline){
	  params.omega = g->omega;
	  params.Y = dos;
	  params.N = N;
	  params.T = T;
  	F.function = &nfl;
  }
  else{	
  	params.spline = dosspline;
  	F.function = &nfc;
  }
  F.params = &params;
  
  size_t limit = QUADLIMIT;// work area size
  gsl_integration_workspace *ws = gsl_integration_workspace_alloc (limit);

	gsl_integration_qag(&F, a, b, epsabs, epsrel, limit, 6, ws, &result, &error);
	
  gsl_integration_workspace_free (ws);

	delete dosspline;
  return result; 
}

//---------------- Get G -------------------------------//

void SIAMMC::get_G()
{
  for (int i=0; i<N; i++) g->G[i] =  1.0/(g->omega[i] + mu - epsilon - g->Delta[i] - g->Sigma[i]) ;
  
  //We can't parallelize this due to Clipped being shared amongst the cores - it'll slow things down.
  for (int i=0; i<N; i++) {if (ClipOff(g->G[i])) Clipped = true;}
  if (Clipped) sprintf(ibuffer + strlen(ibuffer),"SIAMMC::run::(Warning) !!Clipping G!!\n");

	for (int i=0; i<N; i++) g->dos[i] = -imag(g->G[i])/M_PI;
}

//------------------------------------------------------//


double SIAMMC::getn0corr()
{

		sprintf(ibuffer + strlen(ibuffer),"SIAMMC::run::correcting n0 integral tail\n");
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
    
		sprintf(ibuffer + strlen(ibuffer),"SIAMMC::run::corr = %f\n",corr);
		
    if (corr > 1e-1){
    	corr=0.0;
    	sprintf(ibuffer + strlen(ibuffer),"SIAMMC::run::(Warning) n0 correction is too large, setting correction to 0\n");
    }
    
    return corr;
}

double SIAMMC::getwG0corr()
{
		double corr1,corr2;
		gsl_set_error_handler_off();
		
		sprintf(ibuffer + strlen(ibuffer),"SIAMMC::run::correcting G0dos integral tail\n");
		for (int i=0;i<fitorder;i++) sprintf(ibuffer + strlen(ibuffer),"SIAMMC::run::L[%d] = %f\n",i,L[i]);
		for (int i=0;i<fitorder;i++) sprintf(ibuffer + strlen(ibuffer),"SIAMMC::run::R[%d] = %f\n",i,R[i]);
		
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
				
		sprintf(ibuffer + strlen(ibuffer),"SIAMMC::run::corr1 = %f corr2 = %f\n",corr1,corr2);
    
    return corr1+corr2;
}


void SIAMMC::SolveSiam(double* V)
{
	MPI_Bcast(V,1,MPI_DOUBLE,0,MPI_COMM_WORLD); //Start MPI : Sync V
  mu0 = V[0];

  //Obtain dos0=-imag(G0)/pi
	get_dos0();
	if (usecubicspline) dos0spline = new dinterpl(g->omega, g->dos0, N);
	//printf("Getting DOS0 from rank %d\n",rank);
	
	//Get SOCSigma (MPI) - this require only g->dos0 and g->omega
  get_SOCSigma();
  
	//Obtain n0
	g->n0 = get_n(g->dos0);
	
	if (tailcorrection){
		double n0corr = getn0corr();
		g->n0 += n0corr;
	}
	
	//Set n = n0 to compute Sigma
	g->n = g->n0; 
	//printf("Getting Sigma from rank %d\n",rank);
	
	//Compute self energy
	get_Sigma();   
	get_G();
	
	g->n = get_n(g->dos);
	//--------------------//

	V[0] = mu0 + (g->n - g->n0); //we need to satisfy (get_n(G) == n)
	if (usecubicspline) delete dos0spline;
	MPI_Bcast(V,1,MPI_DOUBLE,0,MPI_COMM_WORLD); //Start MPI : Sync V
}

void SIAMMC::Amoeba(double accr, double* V)
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
  double n,n0;
  
  sprintf(ibuffer + strlen(ibuffer),"SIAMMC::run::Amoeba::Start search for mu0\n");
  
  while( (not found) and (try_count<1) ) 
  {
     for(x=x_start; x<x_end; x+=x_step)
     {
       V[0] = x;
       
       SolveSiam(V); //Run impurity Solver
       
       double x_res=V[0];
       
			 //sprintf(ibuffer + strlen(ibuffer),"SIAMMC::run::Amoeba::mu0: %.15f n(G0)-n(G): %.2le step: %.2le true diff: %.3f n(G): %.3f g->n: %.3f\n",x, x-x_res, x_step, get_n(g->G0)- get_n(g->G),get_n(g->G),g->n);
			 sprintf(ibuffer + strlen(ibuffer),"SIAMMC::run::Amoeba::mu0: %.15f n(G0)-n(G): %.2le step: %.2le true diff: %.3f n(G): %.3f\n",x, x-x_res, x_step, n0-n,n);

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
    	sprintf(ibuffer + strlen(ibuffer),"SIAMMC::run::Amoeba::mu0 candidate NOT found! now scanning a wider range...\n"); }
  } 
 
  
  if (not found)
  {
  
     sprintf(ibuffer + strlen(ibuffer),"SIAMMC::run::Amoeba::mu0 candidate NOT found! setting mu0 to to best choice: mu0_best: %f diff: %.2le\n",x_best,diff_best);
     
     V[0] = x_best;
     SolveSiam(V); //Run impurity Solver
  }
  else
  {
    sprintf(ibuffer + strlen(ibuffer),"SIAMMC::run::Amoeba::mu0 candidate found! proceeding with amoeba...\n");  
    
    x = x_candidate;
    x_step *= 0.5;
    x += x_step;

    bool converged = false;
    int it = 0;
    while( (not converged) and (it<=AmoebaMaxIts) )
    {
			it ++;
      V[0] = x;
      SolveSiam(V); //Run impurity Solver
      double x_res=V[0];
      
		  //sprintf(ibuffer + strlen(ibuffer),"SIAMMC::run::Amoeba::it: %d mu0: %.15f n(G0)-n(G): %.2le step: %le n0: %.3f n: %.3f\n", it, x, x-x_res, x_step, get_n(g->G0), get_n(g->G));
		  sprintf(ibuffer + strlen(ibuffer),"SIAMMC::run::Amoeba::it: %d mu0: %.15f n(G0)-n(G): %.2le step: %le n0: %.3f n: %.3f\n", it, x, x-x_res, x_step, n0, n);
		  
      converged = ( abs(x-x_res) < accr );
      int sign =  ((x-x_res)>=0) ? 1 : -1;//int_sign(x - x_res);
      if (sign_old==sign)
         x_step = abs(x_step);
      else
         x_step = -abs(x_step); 
      x_step *= 0.5;
      x += x_step;
    }
    if (converged) sprintf(ibuffer + strlen(ibuffer),"SIAMMC::run::Amoeba::desired accuracy reached!\n");
  }
  sprintf(ibuffer + strlen(ibuffer),"SIAMMC::run::Amoeba::--- Amoeba DONE ---\n");
}

//-----------------------Miscellaneous---------------------------------//

bool SIAMMC::ClipOff(complex<double> &X)
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
void SIAMMC::PrintFullResult(const char* ResultFN)
{
	if (rank==0) g->PrintFullResult(ResultFN);
}

void SIAMMC::PrintResult()
{
	if (rank==0) g->PrintResult("Gf.out","Sig.out");
}

void SIAMMC::PrintResult(const char* Gffile,const char* Sigfile)
{
	if (rank==0) g->PrintResult(Gffile,Sigfile);
}

const char* SIAMMC::buffer() { return ibuffer; }

void SIAMMC::PrintBuffer(const char* FN) {
	if (rank==0){
		FILE * flog;
		flog = fopen(FN, "a");
		if (flog == NULL) { 
			char buff[100]; 
			snprintf(buff, sizeof(buff), "-- ERROR -- Failed to open file %s : ",FN); 
			perror(buff); };
		fprintf(flog,"%s",ibuffer);
		fclose(flog);
	}
}

void SIAMMC::PrintBuffer(const char* FN,bool quiet) {
	if (rank==0){
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
}
