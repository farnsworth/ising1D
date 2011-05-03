
#include "ising1D.hh"
#include "error.hh"
#include "io.hh"
#include "matrices.hh"

#include <iostream>
#include <cmath>
#include <complex>
using namespace std;


int main(int argc, char *argv[])
{

  if (argc > 1){
    double eps=1.0e-12;
    double err=0.0,dtemp,dtemp2;
    string filename = argv[1];

    in_file ifile(filename);

    quench gigi( &ifile );
    state s(gigi.get_size());

    dtemp = 0.0;
    while (1){
      dtemp += gigi.get_calpha2(&s);
      if (s.islast()) break;
      ++s;
    }
    err = abs(dtemp-1.0);
    if (err>eps) exit(1);
    cout << "max error " << err << endl;
  }
  else{
    _ERROR_("no file name given");
    return -1;
  }
  
  exit(0);
}
