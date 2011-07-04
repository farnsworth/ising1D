
template <class T>
matrix<T>::matrix( int nrow_in, int ncol_in )
{
  nrow = nrow_in;
  ncol = ncol_in;
  pointer = new T[nrow*ncol];
}

template <class T>
matrix<T>::matrix(const matrix<T>& source)
{
  nrow = source.nrow;
  ncol = source.ncol;
  pointer = new T[nrow*ncol];
  *this = source;
}


template <class T>
int matrix<T>::index( const int irow, const int icol ) const
{
  if ((irow>=nrow)||(icol>=ncol))
    _ERROR_("matrices index are outside of range",0);
  return irow*ncol+icol;
}

template <class T>
int matrix<T>::get_nrow() const
{
  return nrow;
}

template <class T>
int matrix<T>::get_ncol() const
{
  return ncol;
}

template <class T>
T& matrix<T>::operator() ( const int irow, const int icol ) const
{
  return pointer[index(irow,icol)];
}


template <class T>
matrix<T>::~matrix()
{
  delete [] pointer;
}


template <class T>
matrix<T> matrix<T>::transpose() const
{
  matrix<T> result = matrix(ncol,nrow);
  
  for (int i=0;i<nrow;++i)
    for (int j=0;j<ncol;++j){
      result(j,i) = (*this)(i,j);
    }
  return result;
}

/* in general it is an identity, for complex numbers see matrices.cc*/
template <class T>
matrix<T> matrix<T>::conjugate() const
{
  matrix<T> result = matrix(nrow,ncol);
  
  for (int i=0;i<nrow;++i)
    for (int j=0;j<ncol;++j){
      result(i,j) = (*this)(i,j);
    }
  return result;
}


template <class T>
matrix<T> matrix<T>::daga() const
{
  return (this->conjugate()).transpose();
}


template <class T>
inline matrix<T>& matrix<T>::operator=( const T source)
{
  for (int i=0;i<nrow*ncol;++i)
    this->pointer[i] = source;

  return *this;
}


template <class T>
matrix<T>& matrix<T>::operator=( const matrix<T> &source)
{
  if (this == &source){
    _WARNING_("you are assigning a matrix to itself");
    return *this;
  }
  if ( (nrow != source.nrow)||(ncol != source.ncol))
    _ERROR_("you are assigning two matrix with different size",*this);

  for (int i=0;i<nrow*ncol;++i)
    this->pointer[i] = source.pointer[i];

  return *this;
}

template <class T>
ostream& operator<<( ostream& out, const matrix<T> &m )
{
  for (int irow=0;irow<m.nrow;++irow){
    for (int icol=0;icol<m.ncol;++icol)
      out << m.pointer[m.index(irow,icol)] << "\t";
    out << endl;
  }
  return out;
}



template <class T>
matrix<T> operator+(const matrix<T> &c1, const matrix<T> &c2)
{
  matrix<T> result(c1.nrow,c1.ncol);

  if ( (c1.nrow != c2.nrow)||(c1.ncol != c2.ncol)){
    _ERROR_("you are summing two matrices with different size",result);
  }

  for (int i=0;i<c1.nrow*c1.ncol;++i)
    result.pointer[i] = c1.pointer[i] + c2.pointer[i];

  return result;
}


template <class T>
matrix<T> operator-(const matrix<T> &c1, const matrix<T> &c2)
{
  matrix<T> result(c1.nrow,c1.ncol);
  
  if ( (c1.nrow != c2.nrow)||(c1.ncol != c2.ncol))
    _ERROR_("you are summing two matrices with different size",result);

  for (int i=0;i<c1.nrow*c1.ncol;++i)
    result.pointer[i] = c1.pointer[i] - c2.pointer[i];

  return result;
}


template <class T>
matrix<T> operator*(const matrix<T> &c1,const matrix<T> &c2)
{
  return gemm(c1,'N',c2,'N');
}


template <class T>
matrix<T> gemm( const matrix<T> &a, const char transa, const matrix<T> &\
	   b, const char transb)
{
  int M = ( (transa=='N')||(transa=='n') ) ? a.nrow : a.ncol;
  int N = ( (transb=='N')||(transb=='n') ) ? b.ncol : b.nrow;
  int K = ( (transa=='N')||(transa=='n') ) ? a.ncol : a.nrow;
  int Kcheck = ( (transb=='N')||(transb=='n') ) ? b.nrow : b.ncol;

  matrix<T> result(M,N);
  /* check for dimension errors */
  if ( (K != Kcheck) )
    _ERROR_("incompatible matrices",result);

  matrix<T> c1(M,K),c2(K,N);
  
  if ( ( (transa=='N')||(transa=='n') ) )
    c1 = a;
  else if ( ( (transa=='T')||(transa=='t') ) )
    c1 = a.transpose();
  else
    c1 = a.daga();


  if ( ( (transb=='N')||(transb=='n') ) )
    c2 = b;
  else if ( ( (transb=='T')||(transb=='t') ) )
    c2 = b.transpose();
  else
    c2 = b.daga();

  result = 0.0;
  
  for (int irow=0;irow<M;++irow)
    for (int icol=0;icol<N;++icol)
      for (int ik=0;ik<K;++ik)
	result(irow,icol) += c1.pointer[c1.index(irow,ik)]*c2.pointer[c2.index(ik,icol)];

  return result;
}




template <class T>
matrix<T> operator*(const T alpha, const matrix<T> &source )
{
  matrix<T> result(source.nrow,source.ncol);

  for (int i=0;i<source.nrow*source.ncol;++i)
    result.pointer[i] = alpha*source.pointer[i];

  return result;
}


template <class T>
matrix<T> operator*(const matrix<T> &source, const T alpha )
{
  matrix<T> result(source.nrow,source.ncol);

  result = alpha*source;
  return result;
}


/* ------------------------------------------------------------------- */
/*      complex matrix time a matrix gives always a complex matrix     */

template <class T>
matrix< complex<T> > gemm(const matrix< complex<T> > &a, const char transa, const matrix<T> &b, const char transb)
{
  matrix< complex<T> > temp(b.nrow,b.ncol);
  for (int irow=0;irow<b.nrow;++irow)
    for (int icol=0;icol<b.ncol;++icol)
      temp(irow,icol) = b(irow,icol);
  
  return gemm(a,transa,temp,transb);
}

template <class T>
matrix< complex<T> > gemm( const matrix<T> &a, const char transa, const matrix< complex<T> > &b, const char transb )
{
  matrix< complex<T> > temp(a.nrow,a.ncol);
  for (int irow=0;irow<a.nrow;++irow)
    for (int icol=0;icol<a.ncol;++icol)
      temp(irow,icol) = a(irow,icol);
  return gemm(temp,transa,b,transb);
}
/* -------------------------------------------------------------------- */

template <class T>
void matrix<T>::print() const
{
  for (int irow=0;irow<nrow;++irow){
    for (int icol=0;icol<ncol;++icol)
      cout << setw(10) << setiosflags(ios::fixed) << setprecision(5) << pointer[irow*ncol+icol];
    cout << endl;
  }
}

template <class T>
T* matrix<T>::diagonalize( bool evect)
{
  _ERROR_("diagonalization is not yet implemented for this kind of matrices",NULL);
}

template <class T>
T matrix<T>::det()
{
  _ERROR_("determinant is not yet implemented for this kind of matrices",0);
}

template<>
float* matrix<float>::diagonalize( bool evect );

template<>
double* matrix<double>::diagonalize( bool evect );

template <>
matrix< complex<double> > matrix< complex<double> >::conjugate() const;

template<>
float matrix<float>::det();

template<>
double matrix<double>::det();

template <>
complex<double> matrix< complex<double> >::det();


#if defined(BLAS) || defined(MIXED) || defined(CUDA)
template<>
matrix<float> gemm( const matrix<float> &a, const char trana, const matrix<float> &b, const char tranb);

template<>
matrix<double> gemm( const matrix<double> &a, const char trana, const matrix<double> &b, const char tranb);

template<>
matrix<complex<float> > gemm( const matrix<complex<float> > &a, const char trana, const matrix<complex<float> > &b, const char tranb);

template<>
matrix<complex<double> > gemm( const matrix<complex<double> > &a, const char trana, const matrix<complex<double> > &b, const char tranb );
#endif
