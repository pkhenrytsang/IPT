/*
	Some routines with regards to file IO and trapezoidal integration
*/


#include <cstdio>
#include "routines.h"
#include <iostream>
#include <cmath>
#include <omp.h>

using namespace std;

/*

	Mathematical Routines

*/

void PrintFunc3D(const char* FileName, int N, complex<double>** Y, double* X)
{
  FILE *f;
  f = fopen(FileName, "w");
  for (int i=0; i<N; i+=20)
    for (int j=0; j<N; j+=20)
      fprintf(f,"%.15le %.15le %.15le %.15le\n", X[i], X[j], real(Y[i][j]), imag(Y[i][j]));
  fclose(f);
}

void PrintFunc(const char* FileName, int N, int M, double** Y, double* X)
{
  FILE *f;
  f = fopen(FileName, "w");
  for (int i=0; i<N; i++)
  {
    fprintf(f, "%.15le", X[i]);
    for (int j=0; j<M; j++)
    {
       // loop through and store the numbers into the file
       fprintf(f, "%.15le", Y[i][j] );
    }
   fprintf(f, "\n");
  }
  fclose(f);
}

void PrintFunc(const char* FileName, int N, complex<double>* Y, double* X)
{
  FILE *f;
  f = fopen(FileName, "w");
  for (int i=0; i<N; i++)
    fprintf(f,"%.15le %.15le %.15le\n", X[i], real(Y[i]), imag(Y[i]));
  fclose(f);
}

void PrintFunc(const char* FileName, int N, complex<double>* Y)
{
  FILE *f;
  f = fopen(FileName, "w");
  for (int i=0; i<N; i++)
    fprintf(f,"%d %.15le %.15le\n", i, real(Y[i]), imag(Y[i]));
  fclose(f);
}

void PrintFunc(const char* FileName, int N, double* Y)
{
  FILE *f;
  f = fopen(FileName, "w");
  for (int i=0; i<N; i++)
    fprintf(f,"%d %.15le \n", i, Y[i]);
  fclose(f);
}


void PrintFunc(const char* FileName, int N, double* Y, double* X)
{
  FILE *f;
  f = fopen(FileName, "w");
  for (int i=0; i<N; i++)
    fprintf(f,"%.15le %.15le\n", X[i], Y[i]);
  fclose(f);
}

void PrintFunc(const char* FileName, std::vector< complex<double> > Y, std::vector<double> X)
{
  FILE *f;
  f = fopen(FileName, "w");
  for (int i=0; i<X.size(); i++)
    fprintf(f,"%.15le %.15le %.15le\n", X[i], real(Y[i]), imag(Y[i]));
  fclose(f);
}

void GetDimensions(const char* FileName, int &N, int &M)
{
 //----open file---//
  FILE *f;
  f = fopen(FileName, "r");
  if (f==NULL) perror ("Error opening file");

  //----count rows----//
  int i=0;
  char str[100000];
  while (!feof(f))
  {
    fgets ( str, 100000, f );
    i++;
  }
  N=i-1;

  //----count columns---//
  i=1;
  int j=0;
  while (str[i] != '\0')
  {  if ((str[i]!=' ')and(str[i+1]!='\0')and(str[i-1]==' ')) j++;
     i++;
  }
  M=j+1;

  //---close file-----//
  fclose(f);
}

void ReadFunc(const char* FileName, int &N, int &M, double** &X)
{
 GetDimensions(FileName, N, M);
 //printf("N: %d M: %d\n", N, M);

 X = new double*[M];
 for (int i=0; i<M; i++)
   X[i] = new double[N];

 FILE *f;
 f = fopen(FileName, "r");

 for (int i=0; i<N; i++)
   for (int j=0; j<M; j++)
   { double Dummy;
     fscanf(f, "%le", &Dummy);
     X[j][i]=Dummy;
   }
 fclose(f);
}

void ReadFunc(const char* FileName, int &N, complex<double>* &Y, double* &X, bool PurelyReal)
{ // reads file formatted like 
  //  X  ReY(X) ImY(X)
  //----open file---//
  FILE *f;
  f = fopen(FileName, "r");
  if (f==NULL) perror ("Error opening file");
  
  //----count rows----//
  int i=0;
  char str[1000];  
  int prelines = 0;
  while (!feof(f))
  {  
    fgets ( str, 1000, f ); 
    if (str[0]=='\n') continue;
  
    string sline(str);   
    if ( ( sline.find("nan") != string::npos )
          or
         (str[0]=='#')
       )
       { prelines++; continue; };
    i++;
  }
  N=i-1;
  printf("N: %d, prelines: %d \n", N, prelines);
  fclose(f);
 
  X = new double[N];
  Y = new complex<double>[N];

  f = fopen(FileName, "r");

  for (int i=0; i<prelines; i++)
    fgets ( str, 1000, f ); 
   
  for (int i=0; i<N; i++)
  { double Dummy1,Dummy2,Dummy3;
    fscanf(f, "%le", &Dummy1);
    fscanf(f, "%le", &Dummy2);
    if (not PurelyReal) fscanf(f, "%le", &Dummy3);
    X[i]=Dummy1;
    Y[i]=complex<double>(Dummy2, (PurelyReal) ? 0.0 : Dummy3);
  }
  fclose(f);
}

void ReadFunc(const char* FileName, int &N, double* &Y, double* &X)
{ // reads file formatted like 
  //  X  ReY(X) ImY(X)
  //----open file---//
  FILE *f;
  f = fopen(FileName, "r");
  if (f==NULL) perror ("Error opening file");
  
  //----count rows----//
  int i=0;
  char str[1000];  
  int prelines = 0;
  while (!feof(f))
  {  
    fgets ( str, 1000, f ); 
    if (str[0]=='\n') continue;
  
    string sline(str);   
    if ( ( sline.find("nan") != string::npos )
          or
         (str[0]=='#')
       )
       { prelines++; continue; };
    i++;
  }
  N=i-1;
  printf("N: %d, prelines: %d \n", N, prelines);
  fclose(f);
 
  X = new double[N];
  Y = new double[N];

  f = fopen(FileName, "r");

  for (int i=0; i<prelines; i++)
    fgets ( str, 1000, f ); 
   
  for (int i=0; i<N; i++)
  { double Dummy1,Dummy2;
    fscanf(f, "%le", &Dummy1);
    fscanf(f, "%le", &Dummy2);

    X[i]=Dummy1;
    Y[i]=Dummy2;
  }
  fclose(f);
}

void ReadFunc(const char* FileName, int N, double* Y, double* X)
{ // reads file formatted like 
  //  X  ReY(X) ImY(X)
  //----open file---//
  FILE *f;
  f = fopen(FileName, "r");
  if (f==NULL) { printf("Error opening file"); return; }

  for (int i=0; i<N; i++)
  { double Dummy1,Dummy2;
    fscanf(f, "%le", &Dummy1);
    fscanf(f, "%le", &Dummy2);

    X[i]=Dummy1;
    Y[i]=Dummy2;
  }
  fclose(f);
}

bool FileExists(const char* FN){ FILE* f = fopen(FN,"r"); if (f==NULL) return false; else { fclose(f); return true; } };

//---------------- vectors and matrices--------------------//


//------------------ integral routine ---------------------//


double TrapezIntegral(int N, double Y[], double X[])
{
  
  double sum = Y[0]*(X[1]-X[0]) + Y[N-1]*(X[N-1]-X[N-2]);
  for (int i=1; i<N-1; i++)
    sum+=Y[i]*(X[i+1]-X[i-1]);
  return sum*0.5;
}

complex<double> TrapezIntegral(int N, complex<double> Y[], double X[])
{
  complex<double> sum = Y[0]*complex<double>(X[1]-X[0]) +
                         Y[N-1]*complex<double>(X[N-1]-X[N-2]);
  for (int i=1; i<N-1; i++)
    sum+=Y[i]*complex<double>(X[i+1]-X[i-1]);

  return sum*0.5;
}
