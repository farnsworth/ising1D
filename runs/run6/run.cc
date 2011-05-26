
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
    string name,data,fileout1="file1.out",fileout2="file2.out";

    string filename = argv[1];

    in_file ifile(filename);
    ifile.find_tag("parameters");
    while ( ifile.read_data(name,data)>0){
      if (name=="fileout1")
	fileout1 = data;
      if (name=="fileout2")
	fileout2 = data;
    }

    out_file ofile1(fileout1);
    out_file ofile2(fileout2);

    ofile1.copyfile( &ifile );
    ofile2.copyfile( &ifile );

    loop<double> time_loop( &ifile, "time_loop");

    quench gigi( &ifile );
    gigi.set_gge_occupations();
    int size = gigi.get_size();

    int ndata = 10;
    int delta = 2;
    rho rhov(size/2, ndata*delta,gigi.system);
    rhov.set_ensemble_average(gigi.gge);
    rhozz * rhozzv[ndata];
    double rhozzgge[ndata];
    double rhogge[ndata];
    

    for (int j=0;j<ndata;++j){
      int dist = (j+1)*delta;
      rhogge[j] = rhov.get_ensemble_average( dist );

      rhozzv[j] = new rhozz( size/2, dist, gigi.system );
      rhozzv[j]->set_ensemble_average();
      rhozzgge[j] = rhozzv[j]->get_ensemble_average( gigi.gge );
    }

    // no calculation for l=0 of \rho^xx
    matrix<double> outdata1(time_loop.get_steps(),2*ndata+1);
    matrix<double> outdata2(time_loop.get_steps(),2*ndata+1);

#pragma omp parallel firstprivate(time_loop)
    {
      time_loop.splitOMP();
      
      matrix< complex<double> > UU(size,size),VV(size,size);
      
      while (time_loop.again()){
      	t = time_loop.get_val();
	int ii = time_loop.get_index();
	
	outdata1(ii,0) = t;
	outdata2(ii,0) = t;

	gigi.set_time_evolution(t,&UU,&VV);

	rhov.set_time_evolution( &UU, &VV );


	for (int j=0;j<ndata;++j){
	  int dist = (j+1)*delta; 
	  outdata1(ii,j*2+1) = rhov.get_time_evolution(dist);
	  outdata1(ii,j*2+2) = rhogge[j];

	  outdata2(ii,j*2+1) = rhozzv[j]->get_time_evolution( &UU, &VV);
	  outdata2(ii,j*2+2) = rhozzgge[j];
	}

      	time_loop.next();
      }
    }
    ofile1 << setprecision(15) << setw(25);
    ofile1 << outdata1;
    ofile2 << setprecision(15) << setw(25);
    ofile2 << outdata2;

    // check
    int n = outdata2.get_nrow();
    for (int i=0;i<n-1;++i){
      if (outdata1(i,0)>outdata1(i+1,0))
	cout << "errore 1" << i << endl;
      if (outdata2(i,0)>outdata2(i+1,0))
	cout << "errore 2" << i << endl;
      }

    for (int j=0;j<ndata;++j)
      delete rhozzv[j];
  }
  else
    _ERROR_("no file name given",-1);
  
  _ERROR_TRACKING_(-1);
}
