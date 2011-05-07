
#include "ising1D.hh"
#include "error.hh"
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
    int nsites=3;
    int sites[nsites];
    double ggeltmag[nsites];
    localtmag *ltmag[nsites];

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

    tmag tmag1( gigi.system);
    tmag1.set_spvs();
    double mgge = tmag1.get_ensemble_average( gigi.gge )/double(gigi.get_size());

    sites[0] = 0;
    sites[1] = 50;
    sites[2] = gigi.get_size()/2;


    for (int isite=0;isite<nsites;++isite){
      ltmag[isite] = new localtmag(sites[isite],gigi.system);
      ltmag[isite]->set_spvs();
      ggeltmag[isite] = ltmag[isite]->get_ensemble_average( gigi.gge );
      }

    if (ierr > 0){
      cout << "error" << endl;
      exit(1);
    }
    
    ofile << setprecision(15) << setw(25);

    while (time_loop.next()){
      t = time_loop.get_val();
      ofile << t << "\t";

      gigi.set_time_evolution(t);
      
      for (int isite=0;isite<nsites;++isite){
	ofile << ltmag[isite]->get_time_evolution( gigi.UUt,gigi.VVt ) << "\t";
	ofile << ggeltmag[isite] << "\t";
      }
      
      ofile << tmag1.get_time_evolution(gigi.UUt, gigi.VVt)/double(gigi.get_size()) << "\t";
      ofile << mgge << endl;
    }
    
    for (int isite=0;isite<nsites;++isite){
      delete ltmag[isite];
      ltmag[isite] = NULL;
    }
  }
  else
    _ERROR_("no file name given",-1);
  
  _ERROR_TRACKING_(-1);
}
