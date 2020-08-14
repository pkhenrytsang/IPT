/*
	The class grid stores the grid that IPT read and writes, called by SIAM class.
*/

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "GridMC.h"

// Constructor and Destructors

GridMC::GridMC(const double omega_in[],size_t N)
{
	//This constructor assign all memory needed in this class
  
  //Grid
  this->N = N;
  
  omega = new double[N];
  
  //Input
  Delta = new complex<double>[N];
  //Output
  Sigma = new complex<double>[N];
  G = new complex<double>[N];
  
  //copy omega-grid
	for (int i=0;i<N;i++)	omega[i] = omega_in[i];

	//Intermediate steps
	dos = new double[N];
	dos0 = new double[N];
  SOCSigma = new complex<double>[N];
}

GridMC::~GridMC()
{
  ReleaseMemory();
}

void GridMC::ReleaseMemory()
{
  delete [] Sigma;
  delete [] G;
  
  delete [] SOCSigma;
  
  delete [] omega;
  delete [] Delta;
  
  delete [] dos;
  delete [] dos0;
}

void GridMC::PrintFullResult(const char* ResultFN) //this prints everything, including intermediate step
{ 
  FILE *f;
  f = fopen(ResultFN, "w");

  for (int i=0; i<N; i++)
  { 
     // loop through and store the numbers into the file
    fprintf(f, "%.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le\n", 
                   omega[i],					//1 
                   real(Delta[i]), imag(Delta[i]),			//2 3
                   dos0[i], //4
                   real(SOCSigma[i]), imag(SOCSigma[i]), 		//5 6
                   real(Sigma[i]), imag(Sigma[i]),			//7 8
                   real(G[i]), imag(G[i]), // 9 10
                   dos[i] //11
                   );	
                   
                
  }
  fclose(f);
}

void GridMC::PrintResult(const char* Gffile,const char* Sigfile) //Non-debug version to save disk space
{ 
  FILE *fGf;
  fGf = fopen(Gffile, "w");
  for (int i=0; i<N; i++) fprintf(fGf, "%.15le %.15le %.15le\n", omega[i], real(G[i]), imag(G[i]));  
  fclose(fGf);
  
  FILE *fSigma;
  fSigma = fopen(Sigfile, "w");
  for (int i=0; i<N; i++) fprintf(fSigma, "%.15le %.15le %.15le\n", omega[i], real(Sigma[i]), imag(Sigma[i]));  
  fclose(fSigma);
}

