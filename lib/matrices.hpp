
template <class T>
matrix<T>::matrix( int nrow_in, int ncol_in )
{
  nrow = nrow_in;
  ncol = ncol_in;
  pointer = new T[nrow*ncol];
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
T& matrix<T>::operator() ( const int irow, const int icol )
{
  return pointer[index(irow,icol)];
}


template <class T>
matrix<T>::~matrix()
{
  delete [] pointer;
}


template <class T>
matrix<T> matrix<T>::transpose()
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
matrix<T> matrix<T>::conjugate()
{
  matrix<T> result = matrix(nrow,ncol);
  
  for (int i=0;i<nrow;++i)
    for (int j=0;j<ncol;++j){
      result(i,j) = (*this)(i,j);
    }
  return result;
}


template <class T>
matrix<T> matrix<T>::daga()
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
  matrix<T> result(c1.nrow,c2.ncol);
  
  if (c1.ncol != c2.nrow)
    _ERROR_("you are multipling  two incompatible matrices",result);
  
  result = 0.0;
  
  for (int irow=0;irow<c1.nrow;++irow)
    for (int icol=0;icol<c2.ncol;++icol)
      for (int k=0;k<c1.ncol;++k)
	result(irow,icol) += c1.pointer[c1.index(irow,k)]*c2.pointer[c2.index(k,icol)];

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
matrix< complex<T> > operator*(const matrix< complex<T> > &c1, const matrix<T> &c2)
{
  matrix< complex<T> > result(c1.nrow,c2.ncol);

  if (c1.ncol != c2.nrow)
    _ERROR_("you are multipling  two incompatible matrices",result);
  
  result = complex<T>(0.0,0.0);
  
  for (int irow=0;irow<c1.nrow;++irow)
    for (int icol=0;icol<c2.ncol;++icol)
      for (int k=0;k<c1.ncol;++k)
	result(irow,icol) += c1.pointer[c1.index(irow,k)]*c2.pointer[c2.index(k,icol)];

  return result;
}

template <class T>
matrix< complex<T> > operator*( const matrix<T> &c1, const matrix< complex<T> > &c2 )
{
  matrix< complex<T> > result(c1.nrow,c2.ncol);

  if (c1.ncol != c2.nrow)
    _ERROR_("you are multipling  two incompatible matrices",result);
  
  result = complex<T>(0.0,0.0);
  
  for (int irow=0;irow<c1.nrow;++irow)
    for (int icol=0;icol<c2.ncol;++icol)
      for (int k=0;k<c1.ncol;++k)
	result(irow,icol) += c1.pointer[c1.index(irow,k)]*c2.pointer[c2.index(k,icol)];

  return result;
}

/* -------------------------------------------------------------------- */

template <class T>
void matrix<T>::print()
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
  _ERROR_("determinant is not yet implemented for this kind of matrices",);
}

template<>
extern double* matrix<double>::diagonalize( bool evect );

template <>
extern matrix< complex<double> > matrix< complex<double> >::conjugate();

template<>
extern double matrix<double>::det();

template <>
extern complex<double> matrix< complex<double> >::det();
