
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




#ifdef BLAS
/*SUBROUTINE SGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*  Purpose
*  =======
*
*  SGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X**T,
*
*
*  Arguments
*  ==========
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n',  op( A ) = A.
*              TRANSA = 'T' or 't',  op( A ) = A**T.
*              TRANSA = 'C' or 'c',  op( A ) = A**T.
*
*  TRANSB - CHARACTER*1.
*           On entry, TRANSB specifies the form of op( B ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSB = 'N' or 'n',  op( B ) = B.
*              TRANSB = 'T' or 't',  op( B ) = B**T.
*              TRANSB = 'C' or 'c',  op( B ) = B**T.
*
*  M      - INTEGER.
*           On entry,  M  specifies  the number  of rows  of the  matrix
*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N  specifies the number  of columns of the matrix
*           op( B ) and the number of columns of the matrix C. N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry,  K  specifies  the number of columns of the matrix
*           op( A ) and the number of rows of the matrix op( B ). K must
*           be at least  zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by m  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - REAL             array of DIMENSION ( LDB, kb ), where kb is
*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  n by k  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  BETA   - REAL            .
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - REAL             array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n  matrix
*           ( alpha*op( A )*op( B ) + beta*C ).
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit. */

template <>
matrix<float> operator*(const matrix<float> &a,const matrix<float> &b)
{
  matrix<float> result(a.nrow,b.ncol);

  cout << "you are using float cblas routine";
  
  if (a.ncol != b.nrow)
    _ERROR_("you are multipling  two incompatible matrices",result);
  
  enum CBLAS_TRANSPOSE transa = CblasNoTrans,transb = CblasNoTrans;
  enum CBLAS_ORDER order = CblasRowMajor;
  int M = a.nrow, N = b.ncol, K = a.ncol;
  int LDA = M, LDB = K, LDC = M;
  float ALPHA = 1.0, BETA = 0.0;

  cblas_sgemm( order, transa, transb, M, N, K, ALPHA, a.pointer, LDA, b.pointer, LDB,BETA, result.pointer, LDC );
  
  return result;
}


template <>
matrix<double> operator*(const matrix<double> &a,const matrix<double> &b)
{
  cout << "you are using float cblas routine";
}


template <>
matrix< complex<float> > operator*(const matrix< complex<float> > &c1,const matrix< complex<float> > &c2)
{
  cout << "you are using float cblas routine";
}


template <>
matrix< complex<double> > operator*(const matrix< complex<double> > &c1,const matrix< complex<double> > &c2)
{
  cout << "you are using float cblas routine";
}
#endif
