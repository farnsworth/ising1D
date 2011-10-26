
#include "ising1D.hh"
#include "error.hh"
#include "io.hh"
#include "common.hh"
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
using namespace std;

int delta;
#include "generators.cc"

int main(int argc, char *argv[])
{
  if (argc > 1){

    FPType t;
    string name,data,fileout="file.out";
    string filename = argv[1];

    in_file ifile(filename);
    ifile.find_tag("parameters");
    while ( ifile.read_data(name,data)>0){
      if (name=="fileout")
	fileout = data;
      else if (name=="delta")
	istringstream(data) >> delta;
    }

    out_file ofile1(fileout);
    ofile1.copyfile( &ifile );

    loop<FPType> time_loop( &ifile, "time_loop");

    quench gigi( &ifile,&genJJ,&genHHpeaks,&genJJ,&genHH);

    gigi.set_gge_occupations();

    int size = gigi.get_size();

    int ndata = 5;
    int every = 5;
    localtmag **ltmagvect;
    ltmagvect = new localtmag*[ndata];
    FPType *ltmaggge;
    ltmaggge = new FPType[ndata];

    for (int i=0;i<ndata;++i){
      int temp = size/2+i*every;
      ltmagvect[i] = new localtmag( temp, gigi.system );
      ltmagvect[i]->set_spvs();
      ltmaggge[i] = ltmagvect[i]->get_ensemble_average( gigi.gge );
    }

    ofile1 << setprecision(15) << setw(25);


#pragma omp parallel firstprivate(time_loop)
    {
      time_loop.splitOMP();
      
      matrix< complex<double> > UU(size,size),VV(size,size);
      
      while (time_loop.again()){

      	t = time_loop.get_val();
	int ii = time_loop.get_index();
	ofile1 << t << "\t";
	gigi.set_time_evolution(t,&UU,&VV);

	for (int i = 0;i<ndata;++i){
	  ofile1 << ltmagvect[i]->get_time_evolution( &UU, &VV) << "\t";
	  ofile1 << ltmaggge[i] << "\t";
	}

	ofile1 << endl;
	time_loop.next();
      }
    }
    delete [] ltmaggge;
    for (int i=0;i<ndata;++i){
      delete ltmagvect[i];
    }
    delete [] ltmagvect;
  }
  else
    _ERROR_("no file name given",-1);
  
  _ERROR_TRACKING_(-1);
}
