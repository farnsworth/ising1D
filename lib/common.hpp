template<class T>
obs<T>::obs( int size)
{
  this->_size = size;
}

template<class T>
obs<T>::obs( const obs<T>& source)
{
  _size = source._size;
  _gsv  = source._gsv;
}


template<class T>
local_obs<T>::local_obs( int size_in) : obs<T>(size_in)
{
  this->_spvs = NULL;
}


template<class T>
local_obs<T>::local_obs( FPType* val, int size_in) : obs<T>(size_in)
{
  this->_spvs = new T[this->_size];
  this->_gsv = 0.0;
  for (int i=0;i<size_in;++i){
    _spvs[i] = val[i];
    this->_gsv -= val[i];
  }
}

template<class T>
local_obs<T>::local_obs(const local_obs<T>& source) : obs<T>(source)
{
  if (source._spvs == NULL)
    _spvs = NULL;
  else{
    _spvs = new T[this->_size];
    for (int i=0;i<this->_size;++i)
      _spvs[i] = source._spvs[i];
  }
}

template<class T>
local_obs<T>::~local_obs()
{
  delete [] _spvs;
}

template<class T>
void local_obs<T>::set_gsv()
{
  this->_gsv = 0.0;
  for (int i=0;i<this->_size;++i)
    this->_gsv -= _spvs[i];
}

template<class T>
FPType local_obs<T>::get_spv(int num)
{
  return _spvs[num];
}


template<class T>
FPType local_obs<T>::get_from_state(state* s)
{
  FPType result=0.0;
  for (int i=0; i<this->_size; ++i)
    if (s->conf[i])
      result += 2.0*_spvs[i];
  return this->_gsv + result;
}

template<class T>
FPType local_obs<T>::get_ensemble_average(FPType* nk)
{
  FPType result=0.0;
  for (int i=0; i < this->_size; ++i)
    result += 2.0*nk[i]*_spvs[i];
  return this->_gsv + result;
}

template<class T>
FPType* local_obs<T>::getSpvsPointer(){
  return _spvs;
}

template<class T>
loop<T>::loop( int steps, T delta, T initval  )
{
  _index = 0;
  _sindex = 0;
  _steps = steps;
  _delta = delta;
  _initval = initval;
#ifdef DEBUG
  _ERROR_TRACKING_();
#endif
}


template<class T>
loop<T>::loop(in_file *file, const string name)
{
  _index = 0;
  _sindex = 0;
  _steps = 0;
  _delta = 1;
  _initval = 0;
  read( file, name );
#ifdef DEBUG
  _ERROR_TRACKING_();
#endif
}

template<class T>
void loop<T>::read(in_file *file, const string tagname)
{
  string name,data;
  
  if (file->find_tag(tagname)<0){
    _WARNING_("not able to find loop tag");
    return;
  }
  while (  file->read_data(name,data) > 0 ){
    if (name=="steps")
      istringstream(data) >> _steps;
    else if (name=="delta")
      istringstream(data) >> _delta;
    else if (name=="initval")
      istringstream(data) >> _initval;
  }
#ifdef DEBUG
  _ERROR_TRACKING_();
#endif
}


template<class T>
inline loop<T>& loop<T>::operator=( const loop<T> &source)
{
  this->_steps = source._steps;
  this->_index = source._index;
  this->_delta = source._delta;
  this->_initval = source._initval;

  return *this;
}

template<class T>
void loop<T>::splitOMP()
{
#ifdef _OPENMP
  int ithread = omp_get_thread_num();
  int nthread = omp_get_num_threads();
  int window = this->_steps/nthread;
  _sindex = ithread*window;
  _index = _sindex;
  _steps = (ithread==nthread-1) ? _steps - (nthread-1)*window : window;
#else
  _WARNING_("library not compiled with OMP futures");
#endif
}

template<class T>
bool loop<T>::again()
{
  if (_index < _sindex + _steps ) return true;
  return false;
}


template<class T>
void loop<T>::next()
{
  ++_index;
}


template<class T>
int loop<T>::get_index()
{
  return _index;
}

template<class T>
int loop<T>::get_steps()
{
  return _steps;
}

template<class T>
void loop<T>::restart()
{
  _index = _sindex;
}

template<>
double loop<double>::get_val(int index);

template<>
int loop<int>::get_val(int index);

template<>
float loop<float>::get_val(int index);

template<>
double loop<double>::get_max_val();

template<>
int loop<int>::get_max_val();

template<>
float loop<float>::get_max_val();


template<class T>
histogramm<T>::histogramm( T xmin , T xmax, int nbin ){
  _xmin = xmin;
  _xmax = xmax;
  _nbin = nbin;
  _xdelta = (_xmax-_xmin)/T(_nbin);
  _hist = new int[nbin];
  for (int ibin=0;ibin<_nbin;++ibin){
    _hist[ibin] = 0;
  }
  _ndata = 0;
}


template<class T>
histogramm<T>::~histogramm(){
  delete [] _hist;
}


template<class T>
void histogramm<T>::putData( T x){
  putData(&x,1);
}


template<class T>
void histogramm<T>::putData( T *x, int ndata){
  int ibin;
  _ndata += ndata;
  for (int idata=0;idata<ndata;++idata){
    ibin = int( ( x[idata] - _xmin ) / _xdelta );
    if ( (ibin>=0) and (ibin<_nbin) )
      _hist[ibin] += 1;
  }
}


template <class T>
ostream& operator<<( ostream& out, const histogramm<T> &h )
{
  T xbin = h._xmin+0.5*h._xdelta;
  for (int ibin=0;ibin<h._nbin;++ibin){
    out << xbin << "\t" << double(h._hist[ibin])/(double(h._ndata)*h._xdelta);
    out << endl;
    xbin += h._xdelta;
  }
  cout << "ndata" << h._ndata << endl;
  return out;
}


template<class T>
wHistogramm<T>::wHistogramm( T xmin , T xmax, int nbin ){
  _xMin = xmin;
  _xMax = xmax;
  _nBin = nbin;
  _xDelta = (_xMax-_xMin)/T(_nBin);
  _wHist = new T[nbin];
  for (int ibin=0;ibin<_nBin;++ibin){
    _wHist[ibin] = T(0);
  }
  _totWeight = T(0);
}


template<class T>
wHistogramm<T>::~wHistogramm(){
  delete [] _wHist;
}


template<class T>
void wHistogramm<T>::putData( T x, T w){
  putData(&x,&w,1);
}


template<class T>
void wHistogramm<T>::putData( T *x, T *w, int ndata){
  int ibin;
  for (int idata=0;idata<ndata;++idata){
    ibin = int( ( x[idata] - _xMin ) / _xDelta );
    if ( (ibin>=0) and (ibin<_nBin) ){
      _wHist[ibin] += w[idata];
      _totWeight += w[idata];
    }
  }
}


template <class T>
ostream& operator<<( ostream& out, const wHistogramm<T> &h )
{
  T xbin = h._xMin+0.5*h._xDelta;
  for (int ibin=0;ibin<h._nBin;++ibin){
    out << xbin << "\t" << h._wHist[ibin]/(h._totWeight * h._xDelta);
    out << endl;
    xbin += h._xDelta;
  }
  //cout << "ndata" << h._ndata << endl;
  return out;
}
