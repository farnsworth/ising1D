
#include "ising1D.hh"
#include "error.hh"
#include "io.hh"
#include "matrices.hh"

#include <iostream>
#include <cmath>
#include <complex>
using namespace std;

double delta(int i, int j)
{
  if (i==j) return 1.0;
  return 0.0;
}


int main(int argc, char *argv[])
{

  if (argc > 1){
    double eps=1.0e-12;
    string filename = argv[1];

    in_file ifile(filename);

    ising1D gigi( &ifile );
    int size = gigi.get_size();
    matrix<double> temp(size,size);


    temp = gigi.UU->daga() * *gigi.UU + \
      gigi.VV->daga() * *gigi.VV;
    for (int i=0;i<size;++i)
      for (int j=0;j<size;++j)
	if ( abs(temp(i,j)-delta(i,j)) > eps){
	  temp.print();
	  exit(1);
	}
    
    temp = gigi.VV->transpose() * *gigi.UU + \
      gigi.UU->transpose() * *gigi.VV;
    for (int i=0;i<size;++i)
      for (int j=0;j<size;++j)
	if ( abs(temp(i,j)) > eps){
	  temp.print();
	  exit(2);
	}
    
    temp = gigi.UU->daga() * gigi.VV->conjugate() + \
      gigi.VV->daga() * gigi.UU->conjugate();
    for (int i=0;i<size;++i)
      for (int j=0;j<size;++j)
	if ( abs(temp(i,j)) > eps){
	  temp.print();
	  exit(3);
	}


    /*gigi.set_time_evolution(t);*/
  }
  else{
    _ERROR_("no file name given");
    return -1;
  }

  exit(0);
  
  _ERROR_TRACKING_(-1);
}
