
#include "matrices.hh"
#include <complex>

template <>
double matrix<double>::det( )
{
  if (nrow != ncol){
    _ERROR_("you are trying to compute the determiant of a non square matrix");
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
  if (nrow != ncol){
    _ERROR_("you are trying to compute the determiant of a non square matrix");
    return complex<double>(0.0,0.0);
  }

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


template <>
double* matrix<double>::diagonalize( bool evect )
{
  if (nrow != ncol){
    _ERROR_("you are trying to diagonalize a non square matrix");
    return NULL;
  }

  double* eigenvalues = new double[nrow];

  char jobz = (evect) ? 'V': 'N';
  char uplo = 'U';
  int lda = ncol;
  int info;
  int lwork = 3*ncol-1;
  double* work = new double[lwork];

  // change pointer because it destroy the initial matrix
  dsyev_(&jobz, &uplo, &ncol, pointer, &lda, eigenvalues , work, &lwork, &info );
  
  delete [] work;
  
  return eigenvalues;
}

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
