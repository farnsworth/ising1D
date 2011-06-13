
#include "matrices.hh"
#include <complex>
#include <iostream>

template <>
float matrix<float>::det( )
{
  if (nrow != ncol){
    _ERROR_("you are trying to compute the determiant of a non square matrix",0.0);
    return 0.0;
  }
  int lda = ncol;
  int info;
  int sign = 0;
  int* pivot = new int[ncol];
  float result = 1.0;

  sgetrf_( &ncol, &ncol , pointer, &lda, pivot, &info );
  
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
float* matrix<float>::diagonalize( bool evect )
{
  if (nrow != ncol)
    _ERROR_("you are trying to diagonalize a non square matrix",NULL);

  float* eigenvalues = new float[nrow];

  char jobz = (evect) ? 'V': 'N';
  char uplo = 'U';
  int lda = ncol;
  int info;
  //  int lwork = 3*ncol-1;
  //int lwork = 3*ncol-1;
  //double* work = new double[lwork];
  // change pointer because it destroy the initial matrix

  int temp = -1;
  float twork;
  ssyev_( &jobz, &uplo, &ncol, pointer, &lda, eigenvalues , &twork, &temp, &info );

  int lwork = int(twork);
  float* work = new float[lwork];

  ssyev_( &jobz, &uplo, &ncol, pointer, &lda, eigenvalues , work, &lwork, &info );

  delete [] work;
  return eigenvalues;
}



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
matrix< complex<double> > matrix< complex<double> >::conjugate() const
{
  matrix< complex<double> > result = matrix(nrow,ncol);

  for (int i=0;i<nrow;++i)
    for (int j=0;j<ncol;++j){
      result(i,j) = conj((*this)(i,j));
    }
  return result;
}




#ifdef BLAS
template <>
matrix<float> gemm( const matrix<float> &a, const char transa, const matrix<float> &b, const char transb)
{

  // cout << "you are using float gemm routine" << endl;

  enum CBLAS_ORDER order = CblasRowMajor;  
  enum CBLAS_TRANSPOSE btransa,btransb;

  btransa = ( (transa=='N')||(transa=='n') ) ? CblasNoTrans : CblasTrans;
  btransb = ( (transb=='N')||(transb=='n') ) ? CblasNoTrans : CblasTrans;

  int M = ( (transa=='N')||(transa=='n') ) ? a.nrow : a.ncol;
  int N = ( (transb=='N')||(transb=='n') ) ? b.ncol : b.nrow;
  int K = ( (transa=='N')||(transa=='n') ) ? a.ncol : a.nrow;
  int Kcheck = ( (transb=='N')||(transb=='n') ) ? b.nrow : b.ncol;

  matrix<float> result(M,N);

  if ( (K != Kcheck) )
    _ERROR_("incompatible matrices",result);

  int LDA = ( (transa=='N')||(transa=='n') ) ? K : M;
  int LDB = ( (transb=='N')||(transb=='n') ) ? N : K;
  int LDC = N;

  float ALPHA = 1.0, BETA = 0.0;


  cblas_sgemm( order, btransa, btransb, M, N, K, ALPHA, a.pointer, LDA, b.pointer, LDB,BETA, result.pointer, LDC );
  
  return result;
}


template <>
matrix<double> gemm( const matrix<double> &a, const char transa, const matrix<double> &b, const char transb)
{

  // cout << "you are using double gemm routine" << endl;

  enum CBLAS_ORDER order = CblasRowMajor;  
  enum CBLAS_TRANSPOSE btransa,btransb;

  btransa = ( (transa=='N')||(transa=='n') ) ? CblasNoTrans : CblasTrans;
  btransb = ( (transb=='N')||(transb=='n') ) ? CblasNoTrans : CblasTrans;

  int M = ( (transa=='N')||(transa=='n') ) ? a.nrow : a.ncol;
  int N = ( (transb=='N')||(transb=='n') ) ? b.ncol : b.nrow;
  int K = ( (transa=='N')||(transa=='n') ) ? a.ncol : a.nrow;
  int Kcheck = ( (transb=='N')||(transb=='n') ) ? b.nrow : b.ncol;

  matrix<double> result(M,N);

  if ( (K != Kcheck) )
    _ERROR_("incompatible matrices",result);

  int LDA = ( (transa=='N')||(transa=='n') ) ? K : M;
  int LDB = ( (transb=='N')||(transb=='n') ) ? N : K;
  int LDC = N;

  double ALPHA = 1.0, BETA = 0.0;

  cblas_dgemm( order, btransa, btransb, M, N, K, ALPHA, a.pointer, LDA, b.pointer, LDB,BETA, result.pointer, LDC );
  
  return result;
}


template <>
matrix< complex<float> > gemm( const matrix< complex<float> > &a, const char transa, const matrix< complex<float> > &b, const char transb )
{
  // cout << "you are using complex float gemm routine" << endl;

  enum CBLAS_ORDER order = CblasRowMajor;  
  enum CBLAS_TRANSPOSE btransa,btransb;

  if ( ( (transa=='N')||(transa=='n') ) )
    btransa = CblasNoTrans;
  else if ( ( (transa=='T')||(transa=='t') ) )
    btransa = CblasTrans;
  else
    btransa = CblasConjTrans;

  if ( ( (transb=='N')||(transb=='n') ) )
    btransb = CblasNoTrans;
  else if ( ( (transb=='T')||(transb=='t') ) )
    btransb = CblasTrans;
  else
    btransb = CblasConjTrans;

  int M = ( (transa=='N')||(transa=='n') ) ? a.nrow : a.ncol;
  int N = ( (transb=='N')||(transb=='n') ) ? b.ncol : b.nrow;
  int K = ( (transa=='N')||(transa=='n') ) ? a.ncol : a.nrow;
  int Kcheck = ( (transb=='N')||(transb=='n') ) ? b.nrow : b.ncol;

  matrix< complex<float> > result(M,N);
  if ( (K != Kcheck) )
    _ERROR_("incompatible matrices",result);

  int LDA = ( (transa=='N')||(transa=='n') ) ? K : M;
  int LDB = ( (transb=='N')||(transb=='n') ) ? N : K;
  int LDC = N;

  complex<float> ALPHA = 1.0, BETA = 0.0;

  cblas_cgemm( order, btransa, btransb, M, N, K, &ALPHA, a.pointer, LDA, b.pointer, LDB,&BETA, result.pointer, LDC );
  
  return result;
}


template <>
matrix< complex<double> > gemm( const matrix< complex<double> > &a, const char transa, const matrix< complex<double> > &b, const char transb )
{
  //  cout << "you are using complex double gemm routine" << endl;

  enum CBLAS_ORDER order = CblasRowMajor;  
  enum CBLAS_TRANSPOSE btransa,btransb;

  if ( ( (transa=='N')||(transa=='n') ) )
    btransa = CblasNoTrans;
  else if ( ( (transa=='T')||(transa=='t') ) )
    btransa = CblasTrans;
  else
    btransa = CblasConjTrans;

  if ( ( (transb=='N')||(transb=='n') ) )
    btransb = CblasNoTrans;
  else if ( ( (transb=='T')||(transb=='t') ) )
    btransb = CblasTrans;
  else
    btransb = CblasConjTrans;

  int M = ( (transa=='N')||(transa=='n') ) ? a.nrow : a.ncol;
  int N = ( (transb=='N')||(transb=='n') ) ? b.ncol : b.nrow;
  int K = ( (transa=='N')||(transa=='n') ) ? a.ncol : a.nrow;
  int Kcheck = ( (transb=='N')||(transb=='n') ) ? b.nrow : b.ncol;

  matrix< complex<double> > result(M,N);
  if ( (K != Kcheck) )
    _ERROR_("incompatible matrices", result);

  int LDA = ( (transa=='N')||(transa=='n') ) ? K : M;
  int LDB = ( (transb=='N')||(transb=='n') ) ? N : K;
  int LDC = N;

  complex<double> ALPHA = 1.0, BETA = 0.0;

  cblas_zgemm( order, btransa, btransb, M, N, K, &ALPHA, a.pointer, LDA, b.pointer, LDB,&BETA, result.pointer, LDC );

  return result;
}

#endif
