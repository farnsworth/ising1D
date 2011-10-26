
#include "io.hh"
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
using namespace std;

double pi = 4.0 * atan(1.0);


double epsilon(double k,double h){
  return 2.0*abs(cos(k)-h);
}


double kFromL(int l, int n, int Ntot){
  return l*(pi/(n*Ntot));
}


double corr(int l, int n, int Ntot){
  if ( (l%Ntot) == 0){
    if ( (l%(2*n)) == 0)
      return 0.5;
    else{
      double temp = kFromL(l,n,Ntot/2.0);
      return sin(temp*n)/(2.0*n*sin(temp));
    }
  }
  else
    return 0.0;
}


double f(double t,int j1,int j2, int l1, int l2, double h, int n, int Ntot){
  double alpha;
  alpha = (epsilon(kFromL(l1,n,Ntot),h) - epsilon(kFromL(l2,n,Ntot),h))*t;
  alpha -= kFromL(l1,n,Ntot)*(j1-0.5-n/2.0-n*Ntot);
  alpha += kFromL(l2,n,Ntot)*(j2-0.5-n/2.0-n*Ntot);
  return cos(alpha)*corr(l1-l2,n,Ntot)/(n*Ntot);
}


double te(double t, int j1, int j2, double h, int n, int Ntot){
  double res = 0.0;
  for (int l1=-n*Ntot;l1<n*Ntot;++l1){
    for (int l2=-n*Ntot;l2<n*Ntot;++l2){
      res += f(t,j1,j2,l1,l2,h,n,Ntot);
    }
  }
  return res;
}


int main(int argc, char *argv[]){
  int n,Ntot,j1,j2,ntime;
  double t,h,delta;

  string filename = argv[1],fileout="file.out",name,data;

  in_file ifile(filename);

  ifile.find_tag("parameters");

  while ( ifile.read_data(name,data)>0){
    if (name=="fileout")
      fileout = data;
    else if (name=="h")
      istringstream(data) >> h;
    else if (name=="j1")
      istringstream(data) >> j1;
    else if (name=="j2")
      istringstream(data) >> j2;
    else if (name=="n")
      istringstream(data) >> n;
    else if (name=="Ntot")
      istringstream(data) >> Ntot;
    else if (name=="delta")
      istringstream(data) >> delta;
    else if (name=="ntime")
      istringstream(data) >> ntime;
  }

  out_file ofile(fileout);
  ofile.copyfile( &ifile );

  for (int itime=0;itime<ntime;++itime){
    ofile << itime*delta << "\t" << te( itime*delta, j1, j2, h, n, Ntot) << endl;
  }

}
