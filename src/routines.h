#include <complex>
#include <vector>

using namespace std;

//--- integral ---//

double TrapezIntegral(int N, double Y[], double X[]);
complex<double> TrapezIntegral(int N, complex<double> Y[], double X[]);

//======================== IO =======================================//

void ReadFunc(const char* FileName, int &N, int &M, double** &X);

void ReadFunc(const char* FileName, int& N, double* &Y, double* &X);
void ReadFunc(const char* FileName, int &N, double* &Y, double* &X);
void ReadFunc(const char* FileName, int &N, complex<double>* &Y, double* &X);


                                                                                                
