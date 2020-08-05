/*
	Some routines with regards to file IO and trapezoidal integration
*/


#include <cstdio>
#include "routines.h"
#include <iostream>
#include <cmath>

using namespace std;

/*

	Mathematical Routines

*/

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


void ReadFunc(const char* FileName, int &N, double* &X)
{ // reads file formatted like 
  //  X
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

  f = fopen(FileName, "r");

  for (int i=0; i<prelines; i++)
    fgets ( str, 1000, f ); 
   
  for (int i=0; i<N; i++)
  { double Dummy1;
    fscanf(f, "%le", &Dummy1);

    X[i]=Dummy1;
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

void ReadFunc(const char* FileName, int &N, double* &reY, double* &imY, double* &X)
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
  reY = new double[N];
  imY = new double[N];

  f = fopen(FileName, "r");

  for (int i=0; i<prelines; i++)
    fgets ( str, 1000, f ); 
   
  for (int i=0; i<N; i++)
  { double Dummy1,Dummy2,Dummy3;
    fscanf(f, "%le", &Dummy1);
    fscanf(f, "%le", &Dummy2);
    fscanf(f, "%le", &Dummy3);

    X[i]=Dummy1;
    reY[i]=Dummy2;
    imY[i]=Dummy3;
  }
  fclose(f);
}

void ReadFunc(const char* FileName, int &N, complex<double>* &Y, double* &X)
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
  //printf("N: %d, prelines: %d \n", N, prelines);
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
    fscanf(f, "%le", &Dummy3);


    X[i]=Dummy1;
    Y[i]=complex<double>(Dummy2,Dummy3);
    //printf("%f %f %f\n",Dummy1,Dummy2,Dummy3);
  }
  fclose(f);
}


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
