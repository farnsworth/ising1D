
#include "ising1D_o.hh"
#include "error.hh"
#include "common.hh"
#include <iostream>
#include <cmath>
using namespace std;

#define ERROR(a) error(__FILE__,__FUNCTION__,__LINE__,a)
#define WARNING(a) warning(__FILE__,__FUNCTION__,__LINE__,a)


ising1D_o::ising1D_o(int size, bool symmetry, double h, double J)
{
  cout << "Initialization"<< endl;
  if (size%2 == 1)
    {
    ERROR("odd number of particle");
    return;
    }
  this->size = size;

  if (symmetry)
    {
      this->nquantum = size/2;
    }
  else
    {
      this->nquantum = size;
    }

  this->symmetry = symmetry;
  this->J = J;
  this->h = h;
  this->kpoints = this->get_kpoints();
  //  this->energy.set( this->get_spenergy, this );
}

double* ising1D_o::get_kpoints()
{
  double delta = (double) 2 * M_PI / (double) this->size;
  double* result = new double[this->nquantum];
  
  result[0] = (symmetry) ? delta/ 2.0 : -M_PI + delta/ 2.0;

  for (int i=1; i< this->nquantum; i+=1)
    result[i] = result[i-1] + delta;

  return result;
}


ising1D_o::~ising1D_o()
{
  delete [] this->kpoints;
}


void ising1D_o::print()
{
  cout << "size = " << size << "\nJ = " << J << "\nh = " << h << "\nsymmetry = "<< symmetry << endl;
}


local_obs_o::local_obs_o( ising1D_o* system) : local_obs<double>(system->nquantum)
{
  this->occupancy = (system->symmetry) ? 2.0 : 1.0;
  this->system = system;
}



energy::energy( ising1D_o* system_in) : local_obs_o(system_in)
{
  _gsv = 0.0;
  for (int i=0; i<_size; i+=1){
    _gsv -= (occupancy/(double) 2)*get_value(system->kpoints[i]);
    _spvs[i] = get_value(system->kpoints[i]);
  }
}



double energy::get_value( double k)
{
  double J = this->system->J;
  double h = this->system->h;
  return 2.0 * sqrt( J*J + h*h - 2.0*h*J*cos(k) );
}



sigmaz::sigmaz( ising1D_o* system_in ) : local_obs_o(system_in)
{
  _gsv = 0.0;
  for (int i=0; i<_size; i+=1){
    _gsv -= (occupancy/(double) 2)*get_value(system->kpoints[i]);
    _spvs[i] = get_value(system->kpoints[i]);
  }
  _gsv = _gsv/(double) system->size;
}


double sigmaz::get_value( double k)
{
  double J = this->system->J;
  double h = this->system->h;
  double e = 2.0 * sqrt( J*J + h*h - 2.0*h*J*cos(k) );
  return 4.0*(J*cos(k) - h)/e;
}


imAiAj::imAiAj( ising1D_o* system_in, int d_in ) : local_obs_o(system_in)
{
  d = d_in;
  _gsv = 0.0;
  for (int i=0; i<_size; i+=1)
    _spvs[i] = get_value(system->kpoints[i]);
}

double imAiAj::get_value( double k )
{
  double J = this->system->J;
  double h = this->system->h;
  if ((d==0)||( this->system->symmetry))
    return 0.0;
  else  
    return 2.0*sin(abs(k)*d);
}


imBiBj::imBiBj( ising1D_o* system_in, int d_in ) : local_obs_o(system_in)
{
  d = d_in;
  _gsv = 0.0;
  for (int i=0; i<_size; i+=1)
    _spvs[i] = get_value(system->kpoints[i]);
}

double imBiBj::get_value( double k )
{
  double J = this->system->J;
  double h = this->system->h;
  if ((d==0)||( this->system->symmetry))
    return 0.0;
  else  
    return -2.0*sin(abs(k)*d);
}


BiAj::BiAj( ising1D_o* system_in, int d_in ) : local_obs_o(system_in)
{
  d = d_in;
  _gsv = 0.0;
  for (int i=0; i<_size; i+=1){
    _gsv -= (occupancy/(double) 2)*get_value(system->kpoints[i]);
    _spvs[i] = get_value(system->kpoints[i]);
  }
  _gsv = _gsv/(double) system->size;
}

double BiAj::get_value( double k )
{
  double J = this->system->J;
  double h = this->system->h;
  double e = 2.0 * sqrt( J*J + h*h - 2.0*h*J*cos(k) );

  return ( - cos( k*(double)(d-1) ) + h * cos(k* (double)d ) )/e;
}
