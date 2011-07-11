
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

    out_file ofile(fileout);
    ofile.copyfile( &ifile );

    ising1D gigi( &ifile, "system" , &genJJ, &genHH);

    int size = gigi.get_size();

    for (int isite=0; isite<size; ++isite){
      localtmag ltmag( isite, &gigi);
      ltmag.set_spvs();
      ofile << isite << "\t" << ltmag.get_gsv() << endl;
    }
  }
  else
    _ERROR_("no file name given",-1);
  
  _ERROR_TRACKING_(-1);
}
