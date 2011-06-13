
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

    int ndata = 5;
    int delta = 1;
    rho rhov(size/2, ndata*delta,gigi.system);
    rhov.set_ensemble_average(gigi.gge);
    double rhogge[ndata];
    

    for (int j=0;j<ndata;++j){
      int dist = (j+1)*delta;
      rhogge[j] = rhov.get_ensemble_average( dist );
    }

    // no calculation for l=0 of \rho^xx
    matrix<double> outdata(time_loop.get_steps(),2*ndata+1);

#pragma omp parallel firstprivate(time_loop,rhov)
    {
      time_loop.splitOMP();
      
      matrix< complex<double> > UU(size,size),VV(size,size);
      
      while (time_loop.again()){
      	t = time_loop.get_val();
	int ii = time_loop.get_index();
	
	outdata(ii,0) = t;

	gigi.set_time_evolution(t,&UU,&VV);
	rhov.set_time_evolution( &UU, &VV );


	for (int j=0;j<ndata;++j){
	  int dist = (j+1)*delta; 
	  outdata(ii,j*2+1) = rhov.get_time_evolution(dist);
	  outdata(ii,j*2+2) = rhogge[j];
	}

      	time_loop.next();
      }
    }
    ofile << setprecision(15) << setw(25);
    ofile << outdata;

    // check
    int n = outdata.get_nrow();
    for (int i=0;i<n-1;++i)
      if (outdata(i,0)>outdata(i+1,0))
	cout << "errore 1" << i << endl;

  }
  else
    _ERROR_("no file name given",-1);
  
  _ERROR_TRACKING_(-1);
}
