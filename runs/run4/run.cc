
#include "ising1D.hh"
#include "error.hh"
#include "matrices.hh"
#include "io.hh"
#include "common.hh"
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
using namespace std;

int main(int argc, char *argv[])
{
  if (argc > 1){

    double t;
    string name,data,fileout="file.out";

    string filename = argv[1];

    in_file ifile(filename);
    ifile.find_tag("parameters");
    while ( ifile.read_data(name,data)>0){
      if (name=="fileout")
	fileout = data;
    }

    out_file ofile(fileout);

    ofile.copyfile( &ifile );

    loop<double> time_loop( &ifile, "time_loop");

    quench gigi( &ifile );
    gigi.set_gge_occupations();
    int size = gigi.get_size();

    tmag tmag1( gigi.system );
    tmag1.set_spvs();
    double mgge = tmag1.get_ensemble_average( gigi.gge )/double(gigi.get_size());

    localtmag ltmag(size/2,gigi.system);
    ltmag.set_spvs();
    double lgge = ltmag.get_ensemble_average( gigi.gge );

    cidagacj cicj1(size/2,size/2+5,gigi.system);
    cicj1.set_spvs();
    double cicj1gge = cicj1.get_ensemble_average( gigi.gge );

    cidagacj cicj2(size/2,size/2+10,gigi.system);
    cicj2.set_spvs();
    double cicj2gge = cicj2.get_ensemble_average( gigi.gge );
    
    matrix<double> outdata(time_loop.get_steps(),11);

#pragma omp parallel firstprivate(time_loop)
    {
      time_loop.splitOMP();
      complex<double> temp;
      matrix< complex<double> > UU(size,size),VV(size,size);

      while (time_loop.again()){
      	t = time_loop.get_val();
	int ii = time_loop.get_index();
	
	outdata(ii,0) = t;

	gigi.set_time_evolution(t,&UU,&VV);

	outdata(ii,1) = tmag1.get_time_evolution( &UU, &VV)/double(size);
	outdata(ii,2) = mgge;

	outdata(ii,3) = ltmag.get_time_evolution( &UU, &VV);
	outdata(ii,4) = lgge;

	temp = cicj1.get_time_evolution( &UU, &VV);
	outdata(ii,5) = real(temp);
	outdata(ii,6) = imag(temp);
	outdata(ii,7) = cicj1gge;
      
	temp = cicj2.get_time_evolution( &UU, &VV );
	outdata(ii,8) = real(temp);
	outdata(ii,9) = imag(temp);
	outdata(ii,10) = cicj2gge;
	
      	time_loop.next();
      }
    }
    ofile << setprecision(15) << setw(25);
    ofile << outdata;
  }
  else
    _ERROR_("no file name given",-1);
  
  _ERROR_TRACKING_(-1);
}
