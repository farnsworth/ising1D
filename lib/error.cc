
#include "error.hh"
#include <iostream>
#include <string>
using namespace std;

int ierr=0;
string before="  ";

void error(const char* file, const char* function, const int line,const char* message)
{
  cerr << "ERROR: " << message << endl;
  cerr << "Tracking: " << endl;
  cerr << "function " << function << ", file " << file << ", line " << line << endl;
  // cerr << "from function " << function;
  // cerr << " (contained in " << file << ", line " << line << ")" << endl;
  ierr = 1;
  before="  ";
}

void warning(const char* file, const char* function, const int line,const char* message)
{
  cerr << "WARNING: " << message << endl;
  cerr << "function " << function << ", file " << file << ", line " << line << endl;
  ierr = -1;
}

void message(const char* file, const char* function, const int line,const char* message)
{
  cout << "MESSAGE: " << message << endl;
  cout << "function " << function << ", file " << file << ", line " << line << endl;
}


void error_tracking(const char* file, const char* function, const int line)
{
  cerr << before << "function " << function << ", file " << file << ", line " << line << endl;
  before = before+"  ";
}
