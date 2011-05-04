template<class T>
obs<T>::obs( int size)
{
  this->_size = size;
}


template<class T>
local_obs<T>::local_obs( int size_in) : obs<T>(size_in)
{
  this->_spvs = NULL;
}

template<class T>
local_obs<T>::local_obs( double* val, int size_in) : obs<T>(size_in)
{
  this->_spvs = new T[this->_size];
  this->_gsv = 0.0;
  for (int i=0;i<size_in;++i){
    _spvs[i] = val[i];
    this->_gsv -= val[i];
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
double local_obs<T>::get_spv(int num)
{
  return _spvs[num];
}


template<class T>
double local_obs<T>::get_from_state(state* s)
{
  double result=0.0;
  for (int i=0; i<this->_size; ++i)
    if (s->conf[i])
      result += 2.0*_spvs[i];
  return this->_gsv + result;
}

template<class T>
double local_obs<T>::get_ensemble_average(double* nk)
{
  double result=0.0;
  for (int i=0; i < this->_size; ++i)
    result += 2.0*nk[i]*_spvs[i];
  return this->_gsv + result;
}


template<class T>
loop<T>::loop(in_file *file, const string name)
{
  index = 0;
  steps = 0;
  delta = 1;
  initval = 0;
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
      istringstream(data) >> steps;
    else if (name=="delta")
      istringstream(data) >> delta;
    else if (name=="initval")
      istringstream(data) >> initval;
  }
#ifdef DEBUG
  _ERROR_TRACKING_();
#endif
}

template<class T>
bool loop<T>::next()
{
  if (index>=steps)
    return false;
  ++index;
  return true;
}


template<class T>
int loop<T>::get_index()
{
  return index;
}

template<class T>
void loop<T>::restart()
{
  index = 0;
}
