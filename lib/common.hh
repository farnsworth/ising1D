#ifndef __COMMON_hh__
#define __COMMON_hh__

#include "matrices.hh"
#include "io.hh"
#include "error.hh"
#include <string>

double fermi_dist(double, double);

class state {
public:
  state( int );
  ~state();
  bool* conf;
  void clear();
  void next();
  void random();
  void print();
  void from_int( unsigned int );
  void flip(int);
  bool islast();
  friend ostream& operator<<( ostream& , const state &s );
private:
  int size;
};

ostream& operator<<( ostream& , const state &s );
void operator++(state&);

template<class T>
class obs{
public:
  obs(int size);
  //  virtual double get_from_state(bool* state) = 0;
  // virtual double get_ensemble_average(double* nk) = 0;
  virtual T get_time_evolution( matrix< complex<double> > *UUt, matrix< complex<double> > * VVt ) = 0;
protected:
  int _size;
  double _gsv;
};

template<class T>
class local_obs : public obs<T> {
public:
  local_obs( int );
  local_obs( double*, int );
  ~local_obs();
  double get_from_state(state*);
  double get_ensemble_average(double* nk);
  void set_gsv();
  virtual T get_time_evolution( matrix< complex<double> > *UUt, matrix< complex<double> > * VVt )=0;
  double get_spv( int );
protected:
  double* _spvs;
};


template<class T>
class loop {
public:
  loop( in_file *file, const string name="loop");
  loop( int steps, T delta, T initval  );
  bool next();
  T get_val();
  T get_max_val();
  int get_index();
  void read( in_file* file, const string name="loop");
  void restart();
  loop<T>& operator=( const loop<T> &source);
private:
  int _steps;
  int _index;
  T _delta;
  T _initval;
};


#include "common.hpp"
#endif
