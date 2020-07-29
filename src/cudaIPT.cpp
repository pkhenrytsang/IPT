/*
	IPT impurity Sovler
	by Pak Ki Henry Tsang
	
	original code by Jaksa
	
	Design notes: 
	1. This is purely an impurity solver meant to be used by an external interface
	2. The design is aimed to simplify use and all one need to provide is the hybridization function for the Anderson Impurity Model to solve, as well as a PARAMS file. This is inspired by K. Haule's CTQMC code
	3. No arguments parsing here, yet.
	4. Use trapezoidal rule for most integrals at the moment.
	5. Both CPU and GPU versions should be included. (Only GPU at this moment...)
	6. Broyden should be implemented on the external interface
	7. by default, the output is stored in "Gf.out" and "Sigma.out" unless debug mode is enabled, which output and input will be saved to output.res and input.res respectively
	8. TODO : allow code to choose where to output results
*/


#include "Grid.h"
#include "routines.h"
#include "SIAM.h"
#include "Params.h"
#include "dinterpl.h"
#include "log.h"
#include <cmath>
#include <complex>
#include <string>
#include <fstream>

#define CHARSIZE 50
using namespace std;



int main(int argc, char* argv[])
{
  
  //IPT parameters
  bool verbose = true;
  bool defaultgrid = true;
  bool debug = false;
  //string gridfile,deltafile,Gffile,Sigfile;
  char* gridfile = new char[CHARSIZE];
  char* deltafile = new char[CHARSIZE];
  char* Gffile = new char[CHARSIZE];
  char* Sigfile = new char[CHARSIZE];
  bool tailcorrection = true;
  
  //Logging
  Logging logs("IPT.log",verbose);

  Params params("PARAMS",&logs);
    
  //verbose
  { 
    char* pstring = "'Verbose'";
		int S = params.ReadParam(logs.verbose,pstring); 
		if (S==-1) { 
			char buff[100]; snprintf(buff, sizeof(buff), "(%d) params %s not found in Params file\n",S,pstring); 
			logs.verbose = true; //default
			logs.print(buff);
			logs.print("verbose set to true\n");
		}
  }
  
  //debug
  { 
    char* pstring = "'Debug'";
		int S = params.ReadParam(debug,pstring); 
		if (S==-1) { 
			char buff[100]; snprintf(buff, sizeof(buff), "(%d) params %s not found in Params file\n",S,pstring); 
			logs.print(buff);
			debug = false; //default
			logs.print("debug set to false\n");
		}
  }
  
  //Output-Gf
  { 
    char* pstring = "'Gf-file'";
		int S = params.ReadParam(Gffile,pstring); 
		if (S==-1) { 
			char buff[100]; snprintf(buff, sizeof(buff), "(%d) params %s not found in Params file\n",S,pstring); 
			logs.print(buff,false);
			Gffile = "Gf.out"; //default
			//delGf = false;
			logs.print("Gf will be output to Gf.out\n");
		}
  }
  
  //Output-Sigma
  { 
    char* pstring = "'Sig-file'";
		int S = params.ReadParam(Sigfile,pstring); 
		if (S==-1) { 
			char buff[100]; snprintf(buff, sizeof(buff), "(%d) params %s not found in Params file\n",S,pstring); 
			logs.print(buff,false);
			Sigfile = "Sig.out"; //default
			logs.print("Sigma will be output to Sig.out\n");
		}
  }
  
  //Grid
  int n,m,N;
  double** output;
  double* omega;
  complex<double>* Delta;
  {
    char* pstring = "'Grid'";
		int S = params.ReadParam(gridfile,pstring); 
		if (S==-1) { 
			char buff[100]; snprintf(buff, sizeof(buff), "(%d) params %s not found in Params file\n",S,pstring); 
			logs.print(buff,false);
			defaultgrid = true;
			gridfile = "DEBUG";
		}
		else {
			string str(gridfile);
			if (str=="default"){
				defaultgrid = true; 
				logs.print("Use grid supplied by Delta (default)\n");
			}
			else{
				//check if omega file is readable
				FILE * fgrid = fopen(gridfile, "r");
				if (fgrid == NULL) { 
					char buff[CHARSIZE]; 
					snprintf(buff, sizeof(buff), "Failed to open file %s : ",gridfile); 
					perror(buff); 
					char buff2[CHARSIZE]; 
					snprintf(buff2, sizeof(buff2), "Failed to open file %s\n",gridfile); 
					logs.print(buff2,true);
					return -1;
				}
				fclose(fgrid);
				
				ReadFunc(gridfile, n, m, output);
				N = n;
				omega = new double[N]; //Grid
				Delta = new complex<double>[N];	//Bath
				char buff[100]; snprintf(buff, sizeof(buff), "%d omega points read in from %s\n",N,gridfile); 
				logs.print(buff);
				defaultgrid = false;
				
				//read in omega
				for (int i = 0; i<n; i++)
					omega[i] = output[0][i];
					
				//Clear memory	
				for (int i=0; i<m; i++)
					delete [] output[i];
				delete [] output;
			}
		}
	}	
	
	//Delta
  {
    char* pstring = "'Delta'";
		int S = params.ReadParam(deltafile,pstring); 
		if (S==-1) { 
			char buff[100]; snprintf(buff, sizeof(buff), "(%d) params %s not found in Params file\n",S,pstring); 
			logs.print(buff);
			deltafile = "Delta.inp";
			logs.print("deltafile set to Delta.inp\n");
		}
  }
  
  //TailCorrection
  { 
    char* pstring = "'TailCorrection'";
		int S = params.ReadParam(tailcorrection,pstring); 
		if (S==-1) { 
			char buff[100]; snprintf(buff, sizeof(buff), "(%d) params %s not found in Params file\n",S,pstring); 
			logs.print(buff);
			tailcorrection = false; //default
			logs.print("tailcorrection set to false (default)\n");
		}
  }
  
  double A1,B1,A2,B2; //tail info
  if (tailcorrection)
  { 
    char* pstring = "'A1'";
		int S = params.ReadParam(A1,pstring); 
		if (S==-1) { 
			char buff[100]; snprintf(buff, sizeof(buff), "(%d) params %s not found in Params file\n",S,pstring); 
			logs.print(buff);
			tailcorrection = false; //default
			logs.print("tailcorrection set to false (default)\n");
		}
  }
  if (tailcorrection)
  { 
    char* pstring = "'B1'";
		int S = params.ReadParam(B1,pstring); 
		if (S==-1) { 
			char buff[100]; snprintf(buff, sizeof(buff), "(%d) params %s not found in Params file\n",S,pstring); 
			logs.print(buff);
			tailcorrection = false; //default
			logs.print("tailcorrection set to false (default)\n");
		}
  }
  if (tailcorrection)
  { 
    char* pstring = "'A2'";
		int S = params.ReadParam(A2,pstring); 
		if (S==-1) { 
			char buff[100]; snprintf(buff, sizeof(buff), "(%d) params %s not found in Params file\n",S,pstring); 
			logs.print(buff);
			tailcorrection = false; //default
			logs.print("tailcorrection set to false (default)\n");
		}
  }
  if (tailcorrection)
  { 
    char* pstring = "'B2'";
		int S = params.ReadParam(B2,pstring); 
		if (S==-1) { 
			char buff[100]; snprintf(buff, sizeof(buff), "(%d) params %s not found in Params file\n",S,pstring); 
			logs.print(buff);
			tailcorrection = false; //default
			logs.print("tailcorrection set to false (default)\n");
		}
  }
  
  //Check if Delta file is readable
  FILE * fDelta = fopen(deltafile, "r");
  if (fDelta == NULL) { 
  	char buff[CHARSIZE]; 
  	snprintf(buff, sizeof(buff), "Failed to open file %s : ",deltafile); 
  	perror(buff); 
  	char buff2[CHARSIZE]; 
  	snprintf(buff2, sizeof(buff2), "Failed to open file %s\n",deltafile); 
  	logs.print(buff2,true);
  	return -1;
  }
  fclose(fDelta);
  
  
  
  //Read-in Delta
  ReadFunc(deltafile, n, m, output);  
	N=n;
	if (defaultgrid == true){ //Use grid provided by Delta file
		omega = new double[N]; //omega-grid
  	Delta = new complex<double>[N];	//Bath
		for (int i = 0; i<n; i++) {
			omega[i] = output[0][i];
			Delta[i] = complex<double>(output[1][i],output[2][i]);
		}
	}
	else{ //Read-in Delta and interpolate to the omega-grid
		for (int i = 0; i<N; i++)
		  Delta[i] = complex<double>( 
		  dinterpl::linear_eval (omega[i], output[0], output[1], n) ,
		  dinterpl::linear_eval (omega[i], output[0], output[2], n) );
	}
  for (int i=0; i<m; i++)
  	delete [] output[i];
  delete [] output;
  
  //Initialize grid class (it holds memory for everything)
  Grid grid(omega,N);
  
  for (int i=0;i<N;i++){
  	grid.Delta[i] = Delta[i];
  }
  
  //mu
  { 
    char* pstring = "'mu'";
		int S = params.ReadParam(grid.mu,pstring); 
		if (S==-1) { 
			char buff[100]; snprintf(buff, sizeof(buff), "(%d) params %s not found in Params file\nmu set to 0\n",S,pstring); 
			logs.print(buff);
			grid.mu=0;
		}
		else{
			char buff[100]; snprintf(buff, sizeof(buff), "PARAMS::mu set to %f\n",grid.mu); 
			logs.print(buff);
		}
  }
  //mu0
  { 
    char* pstring = "'mu0'";
		int S = params.ReadParam(grid.mu0,pstring); 
		if (S==-1) { 
			char buff[100]; snprintf(buff, sizeof(buff), "(%d) params %s not found in Params file\nmu0 set to 0\n\n",S,pstring); 
			logs.print(buff);
			grid.mu0=0;
		}
		else{
			char buff[100]; snprintf(buff, sizeof(buff), "PARAMS::mu0 set to %f\n",grid.mu0); 
			logs.print(buff);
		}
  }

  //Initialize SIAM solver
  SIAM siam(&params);
  //Feed in tail correction parameters
  siam.tailcorrection = tailcorrection;
  siam.A1 = A1;
  siam.A2 = A2;
  siam.B1 = B1;
  siam.B2 = B2;

  //Initial run
  //printf("--- Initial DMFT run ---\n");
	
	if (debug==true) grid.PrintResult("input.res"); //for debug
  
  siam.Run(&grid,&logs); //run the impurity solver
		
  
  grid.PrintResult(Gffile,Sigfile);
  if (debug==true) grid.PrintResult("output.res"); //for debug
  
  //Free heap
  
  delete [] omega;
  delete [] Delta;
  
  
  delete [] Gffile;
  
  delete [] Sigfile;
  
  delete [] deltafile;
  
  delete [] gridfile;
  
  return 0;
}
