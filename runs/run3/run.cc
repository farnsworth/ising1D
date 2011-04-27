
#include "ising1D.hh"
#include "random.hh"
#include "error.hh"
#include "io.hh"
#include "common.hh"

#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
#include <complex>
using namespace std;

int main(int argc, char *argv[])
{
  //string name,data;

  int seed = time(NULL);
  rand_init( &seed );


  if (argc > 1){
    double time;
    string name,data,fileout="file.out";

    string filename = argv[1];

    in_file ifile(filename);
    ifile.find_tag("parameters");
    while ( ifile.read_data(name,data)>0){
      if (name=="time")
	istringstream(data) >> time;
      else if (name=="fileout")
	fileout = data;
    }
    
    out_file ofile(fileout);

    ofile.copyfile( &ifile );

    loop dloop( &ifile, "loop");

    quench gigi( &ifile );
    gigi.set_gge_occupations();

    
    int fSite= gigi.get_size()/2;
    rho rho1(fSite,dloop.steps,gigi.system);
    
    rho1.set_ensemble_average(gigi.gge);

    gigi.set_time_evolution(time);
    rho1.set_time_evolution( gigi.UUt, gigi.VVt );

    ofile << setprecision(10) << setw(20);
    while (dloop.next()){
      int idist = dloop.get_index();
      ofile << idist << "\t";
      ofile << rho1.get_time_evolution(idist) << "\t";
      ofile << rho1.get_ensemble_average(idist) << endl;
    }

  }
  else{
    _ERROR_("no file name given");
    return -1;
  }
  
  _ERROR_TRACKING_;
}
