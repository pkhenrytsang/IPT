#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "Grid.h"
#include "routines.h"

// Constructor and Destructors

Grid::Grid(const double omega_in[],size_t N)
{
  Initialize(omega_in,N);
}

Grid::~Grid()
{
  ReleaseMemory();
}

void Grid::Initialize(const double omega_in[],size_t N)
{
  this->N = N;
  
  omega = new double[N];

	for (int i=0;i<N;i++){
		omega[i] = omega_in[i]; //copy omega-grid
	}

  Delta = new complex<double>[N];
  fermi = new double[N];
  G0 = new complex<double>[N];
  Ap = new double[N];
  Am = new double[N];
  P1 = new double[N];
  P2 = new double[N];
  SOCSigma = new complex<double>[N];
  Sigma = new complex<double>[N];
  G = new complex<double>[N];
  DOS = new double[N];
  NIDOS = new double[N];
  DOSmed = new double[N];

  n=0.0;
  n0=0.0;
  mu=0.0;
  mu0=0.0;
  
}

void Grid::ReleaseMemory()
{
  delete [] omega;
  delete [] Delta;
  delete [] fermi;          
  delete [] G0;  
  delete [] Ap;          
  delete [] Am;
  delete [] P1;         
  delete [] P2;
  delete [] SOCSigma;
  delete [] Sigma;
  delete [] G;
  delete [] DOS;
  delete [] NIDOS;
  delete [] DOSmed;
}

void Grid::PrintResult(const char* ResultFN)
{ 
  FILE *f;
  f = fopen(ResultFN, "w+");
  
  fprintf(f,"# n = %le n0 = %le mu = %le mu0=%le\n",n,n0,mu,mu0);   

  for (int i=0; i<N; i++)
  { 
     // loop through and store the numbers into the file
    fprintf(f, "%.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le\n", 
                   omega[i], fermi[i],					//1 2
                   real(Delta[i]), imag(Delta[i]),			//3 4
                   real(G0[i]), imag(G0[i]), 				//5 6
                   Ap[i], Am[i], P1[i], P2[i],				//7 8 9 10 
                   real(SOCSigma[i]), imag(SOCSigma[i]), 		//11 12
                   real(Sigma[i]), imag(Sigma[i]),			//13 14
                   real(G[i]), imag(G[i]),				//15 16
                   DOS[i], NIDOS[i], DOSmed[i]);			//17 18 19 
                   
                
  }
  fclose(f);
  //PrintResult(const char* Gffile,const char* Sigfile);
}

void Grid::PrintResult(const char* Gffile,const char* Sigfile) //Non-debug version to save disk space
{ 
  FILE *fGf;
  fGf = fopen("Gf.out", "w");
  fprintf(fGf,"# n = %le n0 = %le mu = %le mu0=%le\n",n,n0,mu,mu0);
  for (int i=0; i<N; i++) fprintf(fGf, "%.15le %.15le %.15le\n", omega[i], real(G[i]), imag(G[i]));  
  fclose(fGf);
  
  FILE *fSigma;
  fSigma = fopen("Sigma.out", "w");
  fprintf(fSigma,"# n = %le n0 = %le mu = %le mu0=%le\n",n,n0,mu,mu0);
  for (int i=0; i<N; i++) fprintf(fSigma, "%.15le %.15le %.15le\n", omega[i], real(Sigma[i]), imag(Sigma[i]));  
  fclose(fSigma);
}

bool Grid::ReadFromFile(const char* ResultFN)
{ 
  FILE *f;
  f = fopen(ResultFN, "r");
  if (f==NULL) { return false; }

  char rstLine[1000];
  fgets ( rstLine, 1000, f );

  char * pch;
  printf ("rstline: %s\n",rstLine);
  pch = strtok (rstLine,"=");
  int counter=1;
  while (pch != NULL)
  { 
    switch (counter)
    {
      case 2: n=atof(pch); printf ("%d: %s, atof: %le \n",counter,pch,n); break;
      case 3: n0=atof(pch); printf ("%d: %s, atof: %le \n",counter,pch,n0); break;
      case 4: mu=atof(pch); printf ("%d: %s, atof: %le \n",counter,pch,mu); break;
      case 5: mu0=atof(pch); printf ("%d: %s, atof: %le \n",counter,pch,mu0); break;
    }
          
    pch = strtok (NULL, "=");
    counter++; 
  }
  
  for (int i=0; i<N; i++)
  { double o, fer, rd, id, rg0, ig0, ap, am, p1, p2, rsocs, isocs, rs, is, rg, ig, dos, nidos, dosmed;
     // loop through and store the numbers into the file
    fscanf(f, "%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le\n", 
                   &o, &fer,			//1 2
                   &rd, &id,			//3 4
                   &rg0, &ig0, 			//5 6
                   &ap, &am, &p1, &p2,		//7 8 9 10 
                   &rsocs, &isocs, 		//11 12
                   &rs, &is,			//13 14
                   &rg, &ig,			//15 16
                   &dos, &nidos, &dosmed);	//17 18 19 

    omega[i] = o;
    fermi[i] = fer;					
    Delta[i] = complex<double>(rd,id);			
    G0[i] = complex<double>(rg0,ig0); 				
    Ap[i] = ap; Am[i]=am; P1[i]=p1; P2[i]=p2;		
    SOCSigma[i] = complex<double>(rsocs,isocs); 	
    Sigma[i] = complex<double>(rs,is);			
    G[i] = complex<double>(rg,ig);			
    DOS[i] = dos; NIDOS[i] = nidos; DOSmed[i]=dosmed;	                                  
  }
  fclose(f);
  
  return true;
}

