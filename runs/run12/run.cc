
#include "ising1D.hh"
#include "error.hh"
#include "io.hh"
#include "common.hh"
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <cstdlib>
using namespace std;

int delta;

FPType genHH0(int i,ising1D* system)
{
  return system->get_h();
}


FPType genJJ0(int i,ising1D* system)
{
  return system->get_J();
}


FPType genHH(int i, ising1D* system)
{
  int size = system->get_size();
  //  return system->get_h() + system->get_epsilon()*exp(double( pow(i - system->get_size()/2, 2 ) )/(2*sigma*sigma) )
  if ( abs(i-size/2) <= delta )
    return system->get_h()+system->get_epsilon();
  else
    return system->get_h();
}


FPType genJJ(int i,ising1D* system)
{
  return system->get_J();
}



int main(int argc, char *argv[])
{
  if (argc > 1){

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

    out_file ofile(fileout);
    ofile.copyfile( &ifile );

    quench gigi( &ifile,&genJJ,&genHH,&genJJ,&genHH);

    gigi.set_gge_occupations();

    int size = gigi.get_size();

    localtmag *obstmp;

    ofile << setprecision(15) << setw(25);

    for (int isite=0;isite<size;++isite){
      obstmp = new localtmag( isite, gigi.system );
      obstmp->set_spvs( );
      ofile << isite-size/2 << "\t" <<  obstmp->get_ensemble_average( gigi.gge ) << endl;
      delete obstmp;
    }    
  }
  else
    _ERROR_("no file name given",-1);
  
  _ERROR_TRACKING_(-1);
}
