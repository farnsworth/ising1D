
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
    string name,data,fileout1="file_timeevolution.out",fileout2="file_snapshots.out";
    string filename = argv[1];

    in_file ifile(filename);
    ifile.find_tag("parameters");
    while ( ifile.read_data(name,data)>0){
      if (name=="fileout1")
	fileout1 = data;
      else if (name=="fileout2")
	fileout2 = data;
      else if (name=="delta")
	istringstream(data) >> delta;
    }

    out_file ofile1(fileout1);
    ofile1.copyfile( &ifile );
    out_file ofile2(fileout2);
    ofile2.copyfile( &ifile );

    loop<FPType> time_loop( &ifile, "time_loop");

    quench gigi( &ifile,&genJJ,&genHH,&genJJ,&genHH);

    gigi.set_gge_occupations();

    int size = gigi.get_size();

    localtmag **ltmagvect;
    ltmagvect = new localtmag*[size];
    FPType *ltmaggge;
    ltmaggge = new FPType[size];

    for (int isite=0;isite<size;++isite){
      ltmagvect[isite] = new localtmag( isite, gigi.system );
      ltmagvect[isite]->set_spvs();
      ltmaggge[isite] = ltmagvect[isite]->get_ensemble_average( gigi.gge );
    }

   
    ofile1 << setprecision(15) << setw(25);
    ofile2 << setprecision(15) << setw(25);

    int ndata = 5;
    int every = 5;
    matrix<FPType> outdata1(time_loop.get_steps(),2*ndata+1);
    matrix<FPType> outdata2(size,2);


#pragma omp parallel firstprivate(time_loop)
    {
      time_loop.splitOMP();
      
      matrix< complex<double> > UU(size,size),VV(size,size);
      
      while (time_loop.again()){

      	t = time_loop.get_val();
	int ii = time_loop.get_index();
	outdata1(ii,0) = t;

	gigi.set_time_evolution(t,&UU,&VV);

	for (int isite = 0;isite<size;++isite){
	  outdata2(isite,0) = ltmagvect[isite]->get_time_evolution( &UU, &VV);
	  outdata2(isite,1) = ltmaggge[isite];
	  if ( isite >= size/2){
	    int temp = isite-size/2;
	    if (( temp % every == 0 ) && ( temp/every < ndata )){
	      outdata1(ii,2*(temp/every)+1) = outdata2(isite,0);
	      outdata1(ii,2*(temp/every)+2) = ltmaggge[isite];
	    }
	  }
	}
	ofile2 << outdata2;
	ofile2 << endl;
	ofile2 << endl;

	time_loop.next();
      }
      ofile1 << outdata1;
    }
    delete [] ltmaggge;
    for (int isite=0;isite<size;++isite){
      delete ltmagvect[isite];
    }
    delete [] ltmagvect;
  }
  else
    _ERROR_("no file name given",-1);
  
  _ERROR_TRACKING_(-1);
}
