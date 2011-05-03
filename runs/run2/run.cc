
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
    int maxdist;
    int fSite = -1;
    double t;
    string name,data,fileout="file.out";

    string filename = argv[1];

    in_file ifile(filename);
    ifile.find_tag("parameters");
    while ( ifile.read_data(name,data)>0){
      if (name=="time")
	istringstream(data) >> t;
      else if (name=="istart")
	istringstream(data) >> fSite;
      else if (name=="fileout")
	fileout = data;
    }
    out_file ofile(fileout);

    ofile.copyfile( &ifile );

    loop<double> time_loop( &ifile, "time_loop");
    loop<int> dist_loop( &ifile, "dist_loop");

    maxdist = dist_loop.get_max_val();

    quench gigi( &ifile );

    gigi.set_gge_occupations();

    tmag tmag1(gigi.system);
    tmag1.set_spvs();
    double mgge = tmag1.get_ensemble_average( gigi.gge )/double(gigi.get_size());


    if (fSite<0)
      fSite= gigi.get_size()/2;
    else
      cout << "you have had manual istart: " << fSite << endl;
    rho rho1( fSite, maxdist, gigi.system );
    
    rho1.set_ensemble_average( gigi.gge );

    ofile << setprecision(10) << setw(20);

    while (time_loop.next()){
      t = time_loop.get_val();
      ofile << t << "\t";

      gigi.set_time_evolution(t);
      ofile << tmag1.get_time_evolution(gigi.UUt, gigi.VVt)/double(gigi.get_size()) << "\t";
      ofile << mgge;
      
      rho1.set_time_evolution( gigi.UUt, gigi.VVt );
      while (dist_loop.next()){
	int idist = dist_loop.get_val();
	ofile << "\t" << rho1.get_time_evolution(idist) << "\t";
	ofile << rho1.get_ensemble_average(idist);
      }
      ofile << endl;
      dist_loop.restart();
    }
  }
  else{
    _ERROR_("no file name given");
    return -1;
  }
  
  _ERROR_TRACKING_(-1);
}
