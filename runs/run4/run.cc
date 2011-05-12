
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
    
    ofile << setprecision(15) << setw(25);

    complex<double> temp;

    while (time_loop.next()){
      t = time_loop.get_val();
      ofile << t << "\t";

      gigi.set_time_evolution(t);

      ofile << tmag1.get_time_evolution(gigi.UUt, gigi.VVt)/double(gigi.get_size()) << "\t";
      ofile << mgge << "\t";

      ofile << ltmag.get_time_evolution(gigi.UUt, gigi.VVt) << "\t" << lgge << "\t";

      temp = cicj1.get_time_evolution(gigi.UUt, gigi.VVt);
      ofile << real(temp) << "\t" << imag(temp) << "\t" << cicj1gge << "\t";
      
      temp = cicj2.get_time_evolution(gigi.UUt, gigi.VVt);
      ofile << real(temp) << "\t" << imag(temp) << "\t" << cicj2gge << endl;

    }
  }
  else
    _ERROR_("no file name given",-1);
  
  _ERROR_TRACKING_(-1);
}
