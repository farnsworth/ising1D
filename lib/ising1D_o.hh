
#ifndef __ISING1Do_hh__
#define __ISING1Do_hh__

#include "common.hh"

// ordered 1D ising model object
class ising1D_o{
public:
  ising1D_o(int size,bool symmetry=false, double h=0.0, double J=1.0);
  ~ising1D_o();

  // length of the chain
  int size;

  // number of single particle quantum numbers
  // equal to size if there isn't symmetry k,-k
  // and equal to size/2 in presence of this symmetry
  int nquantum;

  // parameters of the hamiltonian:
  // H = - J \sum_i \sigma_i^x \sigma_{i+1}^x - h \sum_i \sigma^z_i
  double h,J;

  // true if we impose the symmetry k,-k
  bool symmetry;

  // momentums of the quasiparticles
  double* kpoints;

  // print the parameters that define the system
  void print();


private:
  // get the kpoints (after initialization they are in the array kpoints)
  double* get_kpoints();
};

class local_obs_o : public local_obs<double>{
public:
  local_obs_o( ising1D_o* );
protected:
  ising1D_o* system;
  double occupancy;
};

class energy : public local_obs_o{
protected:
  double get_value( double );
public:
  energy( ising1D_o* );
};


class sigmaz : public local_obs_o{
protected:
  double get_value( double );
public:
  sigmaz( ising1D_o* );
};


class imAiAj : public local_obs_o{
protected:
  double get_value( double );
  int d;
public:
  imAiAj( ising1D_o*, int );
};


class imBiBj : public local_obs_o{
protected:
  double get_value( double );
  int d;
public:
  imBiBj( ising1D_o*, int );
};


class BiAj : public local_obs_o{
protected:
  double get_value( double );
  int d;
public:
  BiAj( ising1D_o*, int );
  };

#endif
