#include <cstdio>

using namespace std;

class Logging{
   public:
     Logging(char* filename,bool verbose);
     ~Logging();
     FILE* flog;
     bool verbose;
     void print(const char* buffer);
     void print(const char* buffer,bool verbose);
   private:
     void lprintf(FILE ** file,const char* buffer,bool verbose);
     const int BUFFSIZE = 100;
};
