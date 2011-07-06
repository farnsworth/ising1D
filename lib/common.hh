#ifndef __COMMON_hh__
#define __COMMON_hh__

#ifdef _OPENMP
#include <omp.h>
#endif
#include "matrices.hh"
#include "io.hh"
#include "error.hh"
#include <string>

typedef double FPType;

FPType fermi_dist(FPType, FPType);

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
  obs(const obs<T>& source);
  //  virtual double get_from_state(bool* state) = 0;
  // virtual double get_ensemble_average(double* nk) = 0;
  virtual T get_time_evolution( matrix< complex<FPType> > *UUt, matrix< complex<FPType> > * VVt ) = 0;
  FPType get_gsv();
protected:
  int _size;
  FPType _gsv;
};

template<class T>
class local_obs : public obs<T> {
public:
  local_obs( int );
  local_obs( FPType*, int );
  local_obs(const local_obs<T>& source);
  ~local_obs();
  FPType get_from_state(state*);
  FPType get_ensemble_average(FPType* nk);
  void set_gsv();
  virtual T get_time_evolution( matrix< complex<FPType> > *UUt, matrix< complex<FPType> > * VVt )=0;
  FPType get_spv( int );
  FPType * getSpvsPointer();
protected:
  FPType* _spvs;
};


template<class T>
class loop {
public:
  loop( in_file *file, const string name="loop");
  loop( int steps, T delta, T initval  );
  void next();
  bool again();
  T get_val( int index=-1);
  T get_max_val();
  int get_index();
  int get_steps();
  void read( in_file* file, const string name="loop");
  void restart();
  void splitOMP();
  loop<T>& operator=( const loop<T> &source);
private:
  int _sindex;
  int _steps;
  int _index;
  T _delta;
  T _initval;
};


template<class T>
class histogramm {
public:
  histogramm( T xmin , T xmax, int nbin );
  ~histogramm();
  void putData( T x);
  void putData( T *x, int ndata);
  template <class T2> friend ostream& operator<<( ostream& , const histogramm<T2> & );
private:
  T _xmin;
  T _xmax;
  T _xdelta;
  int _nbin;
  int *_hist;
  int _ndata;
};


template<class T>
class wHistogramm {
public:
  wHistogramm( T xmin , T xmax, int nbin );
  ~wHistogramm();
  void putData( T x, T w);
  void putData( T *x, T *w,int ndata);
  template <class T2> friend ostream& operator<<( ostream& , const wHistogramm<T2> & );
private:
  T _xMin;
  T _xMax;
  T _xDelta;
  int _nBin;
  T *_wHist;
  T _totWeight;
};


#include "common.hpp"
#endif
