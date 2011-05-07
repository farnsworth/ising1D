
#include "ising1D.hh"
#include "error.hh"
#include "io.hh"
#include "matrices.hh"

#include <iostream>
#include <cmath>
#include <complex>
using namespace std;

complex<double> delta(int i, int j)
{
  if (i==j) return complex<double>(1.0,0.0);
  return complex<double>(0.0,0.0);
}


int main(int argc, char *argv[])
{

  if (argc > 1){
    double eps=1.0e-8;
    double err=0.0,dtemp;
    string filename = argv[1];

    in_file ifile(filename);

    quench gigi( &ifile );
    int size = gigi.get_size();
    matrix< complex<double> > temp(size,size);
    
    gigi.set_time_evolution(100.0);
    
    
    temp = gigi.UUt->daga() * *gigi.UUt +	\
      gigi.VVt->daga() * *gigi.VVt;
    for (int i=0;i<size;++i)
      for (int j=0;j<size;++j){
	dtemp = abs(temp(i,j)-delta(i,j));
	if ( dtemp > eps){
	  temp.print();
	  exit(1);
	}
        if ( dtemp > err)
          err = dtemp;
      }
    temp = gigi.VVt->transpose() * *gigi.UUt +	\
      gigi.UUt->transpose() * *gigi.VVt;
    for (int i=0;i<size;++i)
      for (int j=0;j<size;++j){
	dtemp = abs(temp(i,j));
	if ( dtemp > eps){
	  temp.print();
	  exit(2);
	}
        if ( dtemp > err)
          err = dtemp;
      }
    
    temp = gigi.UUt->daga() * gigi.VVt->conjugate() +	\
      gigi.VVt->daga() * gigi.UUt->conjugate();
    for (int i=0;i<size;++i)
      for (int j=0;j<size;++j){
	dtemp = abs(temp(i,j));
	if (dtemp > eps){
	  temp.print();
	  exit(3);
	}
        if ( dtemp > err)
          err = dtemp;
      }
    cout << "max error " << err << endl;
  }
  else
    _ERROR_("no file name given",-1);
  
  exit(0);
}
