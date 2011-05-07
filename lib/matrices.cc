
#include "matrices.hh"
#include <complex>
#include <iostream>

template <>
double matrix<double>::det( )
{
  if (nrow != ncol){
    _ERROR_("you are trying to compute the determiant of a non square matrix",0.0);
    return 0.0;
  }
  int lda = ncol;
  int info;
  int sign = 0;
  int* pivot = new int[ncol];
  double result = 1.0;

  dgetrf_( &ncol, &ncol , pointer, &lda, pivot, &info );
  
  for (int i=0;i<ncol;++i){
    result = result*(*this)(i,i);
    /* the index start from 1 and not from 0 */
    if ( pivot[i] != i+1) ++sign;
  }

  result = (sign%2 == 0) ? result : -result;

  delete [] pivot;
  return result;
}


template <>
complex<double> matrix< complex<double> >::det( )
{
  if (nrow != ncol)
    _ERROR_("you are trying to compute the determiant of a non square matrix",complex<double>(0.0,0.0));

  int lda = ncol;
  int info;
  int sign = 0;
  int* pivot = new int[ncol];
  complex<double> result = 1.0;

  zgetrf_( &ncol, &ncol , pointer, &lda, pivot, &info );
  
  for (int i=0;i<ncol;++i){
    result = result*(*this)(i,i);
    if ( pivot[i] != i+1) ++sign;
  }

  result = (sign%2 == 0) ? result : -result;

  delete [] pivot;
  return result;
}

#ifdef USUAL_ALG
template <>
double* matrix<double>::diagonalize( bool evect )
{
  if (nrow != ncol)
    _ERROR_("you are trying to diagonalize a non square matrix",NULL);

  double* eigenvalues = new double[nrow];

  char jobz = (evect) ? 'V': 'N';
  char uplo = 'U';
  int lda = ncol;
  int info;
  //  int lwork = 3*ncol-1;
  //int lwork = 3*ncol-1;
  //double* work = new double[lwork];
  // change pointer because it destroy the initial matrix

  int temp = -1;
  double twork;
  dsyev_( &jobz, &uplo, &ncol, pointer, &lda, eigenvalues , &twork, &temp, &info );

  int lwork = int(twork);
  double* work = new double[lwork];

  dsyev_( &jobz, &uplo, &ncol, pointer, &lda, eigenvalues , work, &lwork, &info );

  delete [] work;
  return eigenvalues;
}
#endif


#ifdef DK_ALG
template <>
double* matrix<double>::diagonalize( bool evect )
{
  if (nrow != ncol){
    _ERROR_("you are trying to diagonalize a non square matrix",NULL);
    return NULL;
  }

  _WARNING_("You are using devide and konqueror algorithm to diagonalize the matrix");

  double* eigenvalues = new double[nrow];

  char jobz = (evect) ? 'V': 'N';
  char uplo = 'U';
  int lda = ncol;
  int info;

  //int lwork = (jobz=='N') ? 2*nrow+1 : 1+6*nrow + 2*nrow*nrow;
  //int liwork = (jobz=='N') ? 1 : 3+5*nrow;

  int lwork = -1;
  int liwork = -1;
  int  cliwork;
  double  clwork;
  
  dsyevd_( &jobz, &uplo, &ncol, pointer, &lda, eigenvalues, &clwork, &lwork, &cliwork, &liwork, &info );

  liwork = cliwork;
  lwork = int(clwork);

  double *work = new double[lwork];
  int * iwork = new int[liwork];

  dsyevd_( &jobz, &uplo, &ncol, pointer, &lda, eigenvalues, work, &lwork, iwork, &liwork, &info );
  
  delete [] work;
  delete [] iwork;
  
  return eigenvalues;
}
#endif


template <>
matrix< complex<double> > matrix< complex<double> >::conjugate()
{
  matrix< complex<double> > result = matrix(nrow,ncol);

  for (int i=0;i<nrow;++i)
    for (int j=0;j<ncol;++j){
      result(i,j) = conj((*this)(i,j));
    }
  return result;
}
