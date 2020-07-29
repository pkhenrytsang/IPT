#include <cstdio>
#include "log.h"

using namespace std;

Logging::Logging(char* filename,bool verbose)
{
	this->verbose = verbose;
  flog = fopen(filename, "w");
  if (flog == NULL) { 
  	char buff[BUFFSIZE]; 
  	snprintf(buff, sizeof(buff), "Failed to open file %s : ",filename); 
  	perror(buff); };
  fprintf(flog,"Begin IPT logging\n");
}

Logging::~Logging() {fclose(flog);}

void Logging::print(const char* buffer){
	lprintf(&flog,buffer,verbose);
}

void Logging::print(const char* buffer,bool verbose){
	lprintf(&flog,buffer,verbose);
}

void Logging::lprintf(FILE ** file,const char* buffer,bool verbose){
	fprintf(*file,buffer);
	if (verbose==true) printf(buffer);
}
