/*
	IPT impurity Solver Interface
	by Pak Ki Henry Tsang
	Affiliation: National High Magnetic Field Laboratory, Tallahassee, FL, USA
	August 2020
	
	Email: henrytsang222@gmail.com
	
	DISCLAIMER: original code by Jaksha Vuchichevicc https://github.com/JaksaVucicevic/DMFT
	
	Notes: 
	1. This is purely an impurity solver meant to be used by an external interface
	2. The design is aimed to simplify use and all one need to provide is the hybridization function for the Anderson Impurity Model to solve, as well as a PARAMS file. This is inspired by K. Haule's CTQMC code
	3. Use trapezoidal rule for most integrals at the moment, but there isn't a lot of error estimation.
	4. Both CPU and GPU versions are controlled by this file.
	5. by default, the output is stored in "Gf.out" and "Sig.out" unless debug mode is enabled, which output and input will be saved to output.res and input.res respectively
	6. PARAMS file must be supplied completely and accurately - else program will abort. For some fields "default" is acceptable.
	7. Internally, SIAM solver class is "standalone" in the sense that it is the only class needed to run the impurity solver code.
*/


#include <stdlib.h>
#include <argp.h>
#include "routines.h"
#include "SIAM.h"
#include "Params.h"
#include "dinterpl.h"
#include <cmath>
#include <complex>
#include <fstream>

using namespace std;


const char *argp_program_version =
  "IPT-cuda Aug-2020";
const char *argp_program_bug_address =
  "<henrytsang222@gmail.com>";
  
/* TODO Program documentation. */
static char doc[] =
  "Iterative Pertubation Theory (IPT) Solver for Anderson Impurity Model on the real axis";
  
/* INPUT DESCRIPTION */
static char args_doc[] = "[INPUT]";

/* Options */
static struct argp_option options[] = {
  {"params",  'p',  "FILE",      0,  "Read params from FILE instead of standard input" },
  {"log",  'l',  "FILE",      0,  "Output log to FILE instead of standard output" },
  {"grid",  'w', "FILE",      0,  "Read grid from FILE" },
  {"Gf",  'g', "FILE",      0,  "Output Gf to FILE" },
  {"debug",  'd', 0,      0,  "Output intermediate steps" },
  {"verbose",  'v', 0,      0,  "Produce verbose output" },
  {"quiet",    'q', 0,      0,  "Don't produce any output" },
  {"silent",   's', 0,      OPTION_ALIAS },
  {"output",   'o', "FILE", 0,  "Output to FILE instead of standard output" },
  { 0 }
};


/* Used by main to communicate with parse_opt. */
struct arguments
{
  char *args[1];                /* arg1 */
  int silent, verbose,debug,Gf;
  char *Sigma_file,*Gf_file,*grid_file,*log_file,*params_file;
};

/* Parse a single option. */
static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
  /* Get the input argument from argp_parse, which we
     know is a pointer to our arguments structure. */
  struct arguments *arguments = (struct arguments*)(state)->input;

  switch (key)
    {
    case 'q': case 's':
      arguments->silent = 1;
      break;
    case 'v':
      arguments->verbose = 1;
      break;
    case 'd':
      arguments->debug = 1;
      break;
    case 'o':
      arguments->Sigma_file = arg;
      break;
    case 'g':
      arguments->Gf_file = arg;
      break;
    case 'w':
      arguments->grid_file = arg;
      break;

    case ARGP_KEY_ARG:
      if (state->arg_num >= 1)
        /* Too many arguments. */
        argp_usage (state);

      arguments->args[state->arg_num] = arg;

      break;

    case ARGP_KEY_END:
      if (state->arg_num < 1)
        /* Not enough arguments. */
        argp_usage (state);
      break;

    default:
      return ARGP_ERR_UNKNOWN;
    }
  return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc };




int main(int argc, char* argv[])
{

  struct arguments arguments;

  /* Default values. */
  arguments.debug = 0;
  arguments.silent = 0;
  arguments.verbose = 0;
  arguments.Sigma_file = "\0";
  arguments.Gf_file = "\0";
  arguments.grid_file = "\0";
  arguments.log_file = "IPT.log";
  arguments.params_file = "PARAMS";

  argp_parse (&argp, argc, argv, 0, 0, &arguments);
  
  //IPT parameters
  bool quiet = arguments.silent;
  bool verbose = arguments.verbose;
  bool defaultgrid = true;
  bool debug = arguments.debug;
  //string gridfile,deltafile,Gffile,Sigfile;
  string logfile = arguments.log_file;
  string paramsfile = arguments.params_file;
  char* gridfile = new char[BUFFERSIZE];
  char* deltafile = arguments.args[0];//new char[BUFFERSIZE];
  char* Gffile = new char[BUFFERSIZE];
  char* Sigfile = new char[BUFFERSIZE];
  bool tailcorrection = true;

	//Init Log file
	{
		FILE* LOGFILE = fopen(logfile.c_str(),"w");
		if ( LOGFILE==NULL )
		{
		  fprintf(LOGFILE,"-- ERROR -- Cannot open log file!\n");
		  if (!quiet) printf("-- ERROR -- Cannot open log file!\n");
		  return -1;
		}
		fclose(LOGFILE);
		if (!quiet && verbose) printf("-- INFO -- Program log will be output to %s\n",logfile.c_str());
	}
	
	//Init Params readin
  Params params(paramsfile.c_str(),logfile.c_str(),!verbose || quiet);
  
  //Output-Gf
  bool delete_Gf = true;
  if (arguments.Gf_file=="\0")
  { 
    char* pstring = "'Gf-file'";
		int S = params.ReadParam(Gffile,pstring);
		
		string str(Gffile);
		if (str=="default") Gffile = "Gf.out";
		
		
		
  	FILE* flog = fopen(logfile.c_str(), "a");
		fprintf(flog,"-- INFO -- Gf will be output to %s\n",Gffile);
		fclose(flog);
  }
  else {Gffile = arguments.Gf_file; delete_Gf = false;}
  if (!quiet && verbose) printf("-- INFO -- Gf will be output to %s\n",Gffile);
  
  //Output-Sigma
  bool delete_Sigma = true;
  if (arguments.Sigma_file=="\0")
  { 
    char* pstring = "'Sig-file'";
		int S = params.ReadParam(Sigfile,pstring); 
		
		string str(Sigfile);
		if (str=="default") Gffile = "Sig.out";
		
  	FILE* flog = fopen(logfile.c_str(), "a");
		fprintf(flog,"-- INFO -- Sigma will be output to %s\n",Sigfile);
		fclose(flog);
  }
  else {Sigfile = arguments.Sigma_file; delete_Sigma = false;}
	if (!quiet && verbose) printf("-- INFO -- Sigma will be output to %s\n",Sigfile);
  
  //Grid
  int n,m,N;
  double** output;
  double* omega;
  double* omega_in;
  complex<double>* Delta;
  double* reDelta_in;
  double* imDelta_in;
  bool delete_grid = false;
  {
    char* pstring = "'Grid'";
    if (arguments.grid_file=="\0"){
			int S = params.ReadParam(gridfile,pstring);
			delete_grid = true;
    }
    else if (arguments.grid_file=="default"){
    	gridfile = "default";
		}
		else{
			gridfile = arguments.grid_file;
		}
		
		string str(gridfile);
		if (str=="default"){
			defaultgrid = true; 
			
			if (!quiet && verbose) printf("-- INFO -- Use grid supplied by Delta (default)\n");
			
			FILE* flog = fopen(logfile.c_str(), "a");
			fprintf(flog,"-- INFO -- Use grid supplied by Delta (default)\n");
			fclose(flog);
		}
		else{
			//check if omega file is readable
			FILE * fgrid = fopen(gridfile, "r");
			if (fgrid == NULL) { 
				char buff[BUFFERSIZE]; 
				snprintf(buff, sizeof(buff), "-- ERROR -- Failed to open file %s :",gridfile); 
				perror(buff);
				
				FILE* flog = fopen(logfile.c_str(), "a");
				fprintf(flog,"-- ERROR -- Failed to open file : %s\n",gridfile); 
				fclose(flog);
				
				return -1;
			}
			fclose(fgrid);
			
			{
				if (!quiet && verbose) printf("-- INFO -- Reading in grid from file : %s\n",gridfile);
				FILE* flog = fopen(logfile.c_str(), "a");
				fprintf(flog,"-- INFO -- Reading in grid from file : %s\n",gridfile);
				fclose(flog);
			}
			
			//ReadFunc(gridfile, n, m, output);
			ReadFunc(gridfile, N, omega);
			//N = n;
			//omega = new double[N]; //Grid
			{
				if (!quiet && verbose) printf("-- INFO -- %d omega points read in from %s\n",N,gridfile); 
				FILE* flog = fopen(logfile.c_str(), "a");
				fprintf(flog,"-- INFO -- %d omega points read in from %s\n",N,gridfile); 
				fclose(flog);
			}
			defaultgrid = false;
			
			//read in omega
			//for (int i = 0; i<n; i++)
			//	omega[i] = output[0][i];
				
			//Clear memory	
			//for (int i=0; i<m; i++)
			//	delete [] output[i];
			//delete [] output;
		}
	}	
	
	//Delta
	{
		if (!quiet && verbose) printf("-- INFO -- deltafile set to %s\n",deltafile);
		FILE* flog = fopen(logfile.c_str(), "a");
		fprintf(flog,"-- INFO -- deltafile set to %s\n",deltafile);
		fclose(flog);
	}
	/*
  {
    char* pstring = "'Delta'";
		int S = params.ReadParam(deltafile,pstring); 
		
		string str(deltafile);
		if (str=="default") deltafile = "Delta.inp"; 
  }
  */
  
	
  //Check if Delta file is readable
  FILE * fDelta = fopen(deltafile, "r");
  if (fDelta == NULL) { 
  	char buff[BUFFERSIZE]; 
  	snprintf(buff, sizeof(buff), "-- ERROR -- Failed to open file %s : ",deltafile); 
  	perror(buff); 
		
		FILE* flog = fopen(logfile.c_str(), "a");
		fprintf(flog,"-- ERROR -- Failed to open file %s\n",deltafile); 
		fclose(flog);
  	return -1;
  }
  fclose(fDelta);
  
  
  //Read-in Delta
  //ReadFunc(deltafile, n, m, output);  
	//N=n;
	if (defaultgrid == true){ //Use grid provided by Delta file
		//memory for omega is not yet assigned in this case
		ReadFunc(deltafile, N,Delta,omega); 
		//omega = new double[N]; //omega-grid
  	//Delta = new complex<double>[N];	//Bath
		//for (int i = 0; i<n; i++) {
		//	omega[i] = output[0][i];
		//	Delta[i] = complex<double>(output[1][i],output[2][i]);
		//}
	}
	else{ //Read-in Delta and interpolate to the omega-grid
		int N_in;
		ReadFunc(deltafile, N_in,reDelta_in,imDelta_in,omega_in); 
		Delta = new complex<double>[N];	//Bath
		for (int i = 0; i<N; i++)
		  Delta[i] = complex<double>( 
				dinterpl::linear_eval (omega[i], omega_in, reDelta_in, N_in) ,
				dinterpl::linear_eval (omega[i], omega_in, imDelta_in, N_in) );
  	delete [] omega_in;
  	delete [] reDelta_in;
	}
  //for (int i=0; i<m; i++)
  //	delete [] output[i];
  //delete [] output;


	if (!quiet && verbose) printf("-- INFO -- %d points of Delta read in\n",N);
	{
	FILE* flog = fopen(logfile.c_str(), "a");
	fprintf(flog,"-- INFO -- %d points of Delta read in\n",N);
	fclose(flog);
	}
	
  //Initialize SIAM solver
  struct siamparams solverparams;
  
  solverparams.verbose = verbose;
  
	//impurity parameters
  params.ReadParam(solverparams.U,"'U'");
  params.ReadParam(solverparams.T,"'T'");
  params.ReadParam(solverparams.epsilon,"'epsilon'");
  params.ReadParam(solverparams.mu,"'mu'");
	params.ReadParam(solverparams.mu0,"'mu0'");
	
  //broadening of G0
  params.ReadParam(solverparams.eta,"'eta'");
	  
  //PH symmetry
  params.ReadParam(solverparams.SymmetricCase,"'SymmetricCase'");
  params.ReadParam(solverparams.Fixmu0,"'Fixmu0'");
  
  params.ReadParam(solverparams.Accr,"'Accr'");
  params.ReadParam(solverparams.AmoebaScanStart,"'AmoebaScanStart'");
  params.ReadParam(solverparams.AmoebaScanEnd,"'AmoebaScanEnd'");
  params.ReadParam(solverparams.AmoebaScanStep,"'AmoebaScanStep'");
  params.ReadParam(solverparams.AmoebaMaxIts,"'AmoebaMaxIts'");
  params.ReadParam(solverparams.AmoebaForceScanAndPrintOut,"'AmoebaForceScanAndPrintOut'");
  
  //Kramers Kronig
  params.ReadParam(solverparams.KKAccr,"'KKAccr'");
  params.ReadParam(solverparams.usecubicspline,"'UseCubicSpline'");
  
  
  //G0 integral tail correction
  params.ReadParam(solverparams.A1,"'A1'");
  params.ReadParam(solverparams.A2,"'A2'");
  params.ReadParam(solverparams.B1,"'B1'");
  params.ReadParam(solverparams.B2,"'B2'");
  params.ReadParam(solverparams.tailcorrection,"'TailCorrection'");

  //options
  params.ReadParam(solverparams.CheckSpectralWeight,"'CheckSpectralWeight'"); //default false
  
	if (!quiet && verbose) printf("-- INFO -- Finish initializing paramters\n");
	{
	FILE* flog = fopen(logfile.c_str(), "a");
	fprintf(flog,"-- INFO -- Finish initializing paramters\n");
	fclose(flog);
	}
  
  //Initialize the solver using supplied arguments
  SIAM Solver(omega,N,&solverparams);
	
	if (debug==true) Solver.PrintFullResult("input.res"); //for debug
	if (!quiet && verbose) printf("-- INFO -- Finish printing full input to input.res\n");
	
  //run the impurity solver
  Solver.Run(Delta); 
  
  Solver.PrintBuffer(logfile.c_str(),quiet);
  
	Solver.PrintResult(Gffile,Sigfile);
	
  if (debug==true) Solver.PrintFullResult("output.res"); //for debug
	if (!quiet && verbose) printf("-- INFO -- Finish printing full output to output.res\n");
  
  //Free heap
  
  delete [] omega;
  delete [] Delta;
  //delete [] omega_in;
  //delete [] reDelta_in;
  
  delete [] imDelta_in;
  
  if (delete_Gf) delete [] Gffile;
  
  if (delete_Sigma) delete [] Sigfile;
  
  //delete [] deltafile;
  
  
  if (delete_grid) delete [] gridfile;
  
  return 0;
}
