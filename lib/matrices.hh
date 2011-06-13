#ifndef __GEN_MATRICES_hh__

#define __GEN_MATRICES_hh__

#include "error.hh"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#ifdef BLAS
#ifdef __APPLE__
#include <vecLib/cblas.h>
#else
#include <gsl/gsl_cblas.h>
#endif
#endif
using namespace std;

extern "C" void ssyev_(char *, char *, int *, float *, int *, float *, float *, int *, int * );
extern "C" void dsyev_(char *, char *, int *, double *, int *, double *, double *, int *, int * );
extern "C" void dsyevd_(char *, char *, int *, double *, int *, double *, double *, int *, int*, int*, int*);
extern "C" void sgetrf_(int *,int *, float *, int *, int *, int * );
extern "C" void dgetrf_(int *,int *, double *, int *, int *, int * );
extern "C" void zgetrf_(int *,int *, complex<double> *, int *, int *, int * );


template <class T>
class matrix {
public:
  matrix( int,int );
  matrix(const matrix<T>& source); /*< copy constructor */
  ~matrix();
  void print() const;

  T& operator() ( const int, const int ) const;
  matrix<T>& operator=( const matrix<T> &source);
  matrix<T>& operator=(const T);

  T* diagonalize( bool evect=false );
  T det();
  matrix<T> transpose() const;
  matrix<T> conjugate() const;
  matrix<T> daga() const;
  int get_nrow() const;
  int get_ncol() const;

  template <class T2> friend ostream& operator<<( ostream& , const matrix<T2> & );
  template <class T2> friend matrix<T2> operator+(const matrix<T2> &, const matrix<T2> &);
  template <class T2> friend matrix<T2> operator-(const matrix<T2> &, const matrix<T2> &);
  template <class T2> friend matrix<T2> operator*(const matrix<T2> &, const matrix<T2> &);
  template <class T2> friend matrix<T2> operator*(const matrix<T2> &, const T2 );
  template <class T2> friend matrix<T2> operator*(const T2 , const matrix<T2> &);


  template <class T2> friend matrix<T2> gemm( const matrix<T2> &a, const char trana, const matrix<T2> &b, const char tranb);
  template <class T2> friend matrix< complex<T2> > gemm(const matrix< complex<T2> > &a, const char trana, const matrix<T2> &b, const char tranb);
  template <class T2> friend matrix< complex<T2> > gemm(const matrix<T2> &a, const char trana, const matrix< complex<T2> > &b, const char tranb);


private:
  T* pointer;
  int nrow;
  int ncol;
  int index( const int irow, const int icol ) const;
};

template<class T>
matrix<T> gemm( const matrix<T> &a, const char trana, const matrix<T> &b, const char tranb );

template<class T>
matrix< complex<T> > gemm( const matrix< complex<T> > &a, const char trana, const matrix<T> &b, const char tranb );

template<class T>
matrix< complex<T> > gemm( const matrix<T> &a, const char trana, const matrix< complex<T> > &b, const char tranb );


#include "matrices.hpp"

#endif
