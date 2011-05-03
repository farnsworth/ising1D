
#include "ising1D.hh"
#include "error.hh"
#include "io.hh"
#include "matrices.hh"

#include <iostream>
#include <cmath>
#include <complex>
using namespace std;

double occupation_time_evolution(int mu,quench * qin)
{
  int size = qin->get_size();
  matrix< complex<double> > temp(size,size);
  double result = 0.0;
  
  temp = qin->system->VV->transpose() * *qin->UUt + qin->system->UU->transpose() * *qin->VVt;

  for (int i=0;i<size;++i){
    result += norm( temp(mu,i) );
  }

  return result;

}


int main(int argc, char *argv[])
{

  if (argc > 1){
    double eps=1.0e-12;
    double err=0.0,dtemp;
    string filename = argv[1];

    in_file ifile(filename);

    quench gigi( &ifile );
    gigi.set_gge_occupations();
    gigi.set_time_evolution(100.0);
  
    for (int i =0;i<gigi.get_size();++i){
      dtemp = abs(gigi.gge[i]-occupation_time_evolution(i, &gigi) );
      if (dtemp > eps)
	exit(1);
      if (dtemp > err)
	err = dtemp;
    }
    cout << "max error " << err << endl; 
  }
  else{
    _ERROR_("no file name given");
    return -1;
  }
  
  exit(0);
}
