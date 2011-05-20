#ifndef __GEN_MATRICES_hh__

#define __GEN_MATRICES_hh__

#include "error.hh"
#ifdef BLAS
#include <vecLib/cblas.h>
#endif
#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
using namespace std;

extern "C" void dsyev_(char *, char *, int *, double *, int *, double *, double *, int *, int * );
extern "C" void dsyevd_(char *, char *, int *, double *, int *, double *, double *, int *, int*, int*, int*);
extern "C" void dgetrf_(int *,int *, double *, int *, int *, int * );
extern "C" void zgetrf_(int *,int *, complex<double> *, int *, int *, int * );


template <class T>
class matrix {
public:
  matrix( int,int );
  matrix(const matrix<T>& source); /*< copy constructor */
  ~matrix();
  void print();

  T& operator() ( const int, const int );
  matrix<T>& operator=( const matrix<T> &source);
  matrix<T>& operator=(const T);

  template <class T2> friend ostream& operator<<( ostream& , const matrix<T2> & );
  template <class T2> friend matrix<T2> operator+(const matrix<T2> &, const matrix<T2> &);
  template <class T2> friend matrix<T2> operator-(const matrix<T2> &, const matrix<T2> &);
  template <class T2> friend matrix<T2> operator*(const matrix<T2> &, const matrix<T2> &);
  template <class T2> friend matrix<T2> operator*(const matrix<T2> &, const T2 );
  template <class T2> friend matrix<T2> operator*(const T2 , const matrix<T2> &);

  template <class T2> friend matrix< complex<T2> > operator*(const matrix< complex<T2> > &c1, const matrix<T2> &c2);

  template <class T2> friend matrix< complex<T2> > operator*(const matrix<T2> &c1, const matrix< complex<T2> > &c2);

  T* diagonalize( bool evect=false );
  T det();
  matrix<T> transpose();
  matrix<T> conjugate();
  matrix<T> daga();
  int get_nrow() const;
  int get_ncol() const;

private:
  T* pointer;
  int nrow;
  int ncol;
  int index( const int irow, const int icol ) const;
};

#include "matrices.hpp"

#endif
