#include <cstdio>
#include <string>
#include <iostream>
#include <sstream>
#include "Params.h"
#include "log.h"
using namespace std;


Params::Params(const char* InputFN,Logging* logs)
{
  SetInputFN(InputFN);
  this->logs=logs;
}

Params::~Params()
{
  
}

void Params::SetInputFN(const char* InputFN)
{ 
  this->InputFN.assign(InputFN);
  FILE* InputFile = fopen(this->InputFN.c_str(),"r");
  if ( InputFile==NULL )
    printf("-- WARNING -- Params: Params File does not exist!\n");
  else
  {  //cout << "-- INFO -- Params: Params File name set to: " << this->InputFN << endl;
     fclose(InputFile);
  }
}

template <typename T> int Params::ReadArray(int N, T* Param,const char* ParamName)
{
  char* line = ReadParam(ParamName);
  if (line==NULL) return -1;
  stringstream ss;
  ss << line;
  for(int i = 0; i < N; i++)
    ss >> Param[i];
  return 0;
}

int Params::ReadArray(int N, double* Param, const char* ParamName)
{
  int Err;
  Err = ReadArray<double>(N, Param, ParamName);
  return Err;
}

int Params::ReadArray(int N, int* Param, const char* ParamName)
{
  int Err;
  Err = ReadArray<int>(N, Param, ParamName);
  return Err;
}

int Params::ReadParam(int& Param, const char* ParamName)
{
  char* line = ReadParam(ParamName);
  if (line==NULL) return -1;

  if ( sscanf(line,"%d", &Param) == EOF ) 
  {  
     printf("-- ERROR -- Params: Param %s can not be read\n",ParamName);
     return -1;    
  }
  return 0;
}

int Params::ReadParam(double& Param, const char* ParamName)
{
  char* line = ReadParam(ParamName);
  if (line==NULL) return -1;

  if ( sscanf(line,"%le", &Param) == EOF ) 
  {  
     printf("-- ERROR -- Params: Param %s can not be read\n",ParamName);
     return -1;    
  }
  return 0;
}

int Params::ReadParam(char& Param, const char* ParamName)
{
  char* line = ReadParam(ParamName);
  if (line==NULL) return -1;

  Param = line[0]; 
  return 0;
}

int Params::ReadParam(char* Param, const char* ParamName)
{
  char* line = ReadParam(ParamName);
  if (line==NULL) return -1;

  sscanf(line,"%s",Param);
  return 0;

}

int Params::ReadParam(bool& Param, const char* ParamName)
{
  char* line = ReadParam(ParamName);
  if (line==NULL) return -1;

  if (line[0]=='T')
    Param = true;
  else 
    if (line[0]=='F')
      Param = false;
    else
    {
			char buff[100];
			snprintf(buff, sizeof(buff),"-- ERROR -- Params: Param %s can not be read\n",ParamName);
			logs->print(buff);
      return -1; 
    }
  return 0;
}


char* Params::ReadParam(const char* ParamName)
{ 
  FILE* InputFile = fopen(InputFN.c_str(),"r");
  if ( InputFile==NULL ) return NULL;
    
  char* line = new char[128];
  while ( fgets(line, 128, InputFile ) != NULL ) 
  { if (line[0]=='#') continue;
    string sline(line);
    if ( sline.find(ParamName) != string::npos )
      return line;
  }
  delete [] line;
  fclose (InputFile);
  {
	char buff[100];
	snprintf(buff, sizeof(buff),"-- INFO -- Params: Param %s not found in Params File\n",ParamName);
	logs->print(buff);
	}
  return NULL;
}
