#include <complex>
#include <vector>

using namespace std;

//--- integral ---//

double TrapezIntegral(int N, double Y[], double X[]);
complex<double> TrapezIntegral(int N, complex<double> Y[], double X[]);

//======================== IO =======================================//

void PrintFunc(const char* FileName, int N, int M, double** Y, double* X);
void PrintFunc(const char* FileName, int N, complex<double>* Y, double* X);
void PrintFunc(const char* FileName, int N, complex<double>* Y);
void PrintFunc(const char* FileName, int N, double* Y);
void PrintFunc(const char* FileName, int N, double* Y, double* X);
void PrintFunc(const char* FileName, std::vector< complex<double> > Y, std::vector<double> X);
void PrintFunc3D(const char* FileName, int N, complex<double>** Y, double* X);
void ReadFunc(const char* FileName, int &N, int &M, double** &X);
void ReadFunc(const char* FileName, int N, double* Y, double* X);
bool FileExists(const char* FN);
                                                                                                
