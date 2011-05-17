
#include "ising1D.hh"
#include "common.hh"
#include "matrices.hh"
#include "error.hh"
#include "io.hh"
#include "random.hh"

#include <iostream>
#include <cmath>
#include <string>
#include <complex>
#include <time.h>
using namespace std;


energy::energy(double * data,int size_in) : local_obs<double>( data,size_in ){
}


double energy::get_time_evolution( matrix< complex<double> > *UUt, matrix< complex<double> > *VVt )
{
  return 0.0;
}

localtmag::localtmag(int site, ising1D * system_in) : local_obs<double>( system_in->get_size() ){
  _system = system_in;
  _site = site;
}


double localtmag::_get_time_evolution(int site, matrix< complex<double> > *UUt, matrix< complex<double> > *VVt )
{
  int size = UUt->get_ncol();
  double temp = 0.0;
  
  for (int icol=0;icol<size;++icol)
    temp += norm( (*VVt)(site,icol) );
  
  return 2.0*temp-1.0;
}


double localtmag::get_time_evolution( matrix< complex<double> > *UUt, matrix< complex<double> > *VVt )
{
  return _get_time_evolution( _site, UUt, VVt );
}


void localtmag::set_spvs()
{
  _spvs = new double[_size];
  
  _gsv = -1.0;
  for (int mu=0;mu<_size;++mu){
    _spvs[mu] = (*_system->UU)(_site,mu)*(*_system->UU)(_site,mu)-(*_system->VV)(_site,mu)*(*_system->VV)(_site,mu);
    _gsv += 2.0*(*_system->VV)(_site,mu)*(*_system->VV)(_site,mu);
  }
}


/********* to modify ************/
cidagacj::cidagacj(int i, int j, ising1D * system_in) : local_obs< complex<double> >( system_in->get_size() ){
  _system = system_in;
  _i = i;
  _j = j;
}


complex<double> cidagacj::_get_time_evolution(int i,int j, matrix< complex<double> > *UUt, matrix< complex<double> > *VVt )
{
  int size = UUt->get_ncol();
  complex<double> temp = 0.0;
  
  for (int icol=0;icol<size;++icol)
    temp += (*VVt)(i,icol) * conj( (*VVt)(j,icol) );
  
  return temp;
}


complex<double> cidagacj::get_time_evolution( matrix< complex<double> > *UUt, matrix< complex<double> > *VVt )
{
  return _get_time_evolution( _i,_j, UUt, VVt );
}


void cidagacj::set_spvs()
{
  _spvs = new double[_size];
  
  _gsv = 0.0;
  for (int mu=0;mu<_size;++mu){
    _spvs[mu] = ((*_system->UU)(_i,mu)*(*_system->UU)(_j,mu)-(*_system->VV)(_i,mu)*(*_system->VV)(_j,mu))/2.0;
    _gsv += (*_system->VV)(_i,mu)*(*_system->VV)(_j,mu);
  }

}
/********************************/









tmag::tmag(ising1D * system_in) : local_obs<double>( system_in->get_size() ){
  _system = system_in;
}


double tmag::_get_time_evolution( matrix< complex<double> > *UUt, matrix< complex<double> > *VVt )
{
  int size = UUt->get_ncol();
  double temp = 0.0;
  
  for (int irow=0;irow<size;++irow)
    for (int icol=0;icol<size;++icol)
      //      temp += pow( abs( (*VVt)(irow,icol) ),2);
      temp += norm( (*VVt)(irow,icol) );
  
  return 2.0*temp-double(size);
}


double tmag::get_time_evolution( matrix< complex<double> > *UUt, matrix< complex<double> > *VVt )
{
  return _get_time_evolution( UUt, VVt );
}


void tmag::set_spvs()
{
  _spvs = new double[_size];
  
  for (int mu=0;mu<_size;++mu){
    _spvs[mu] = 0.0;
    for (int irow=0;irow<_size;++irow)
      _spvs[mu] += (*_system->VV)(irow,mu)*(*_system->VV)(irow,mu);
    _spvs[mu] = 1.0-2.0*_spvs[mu];
  }
  set_gsv();
}


BiAj::BiAj(int i_in, int j_in, ising1D * system_in) : local_obs< double >( system_in->get_size() ){
  if ((i_in < 0) || (j_in < 0) || (i_in > system_in->get_size()) || (j_in > system_in->get_size()) ){
    _ERROR_("indexes outside the chain",);
  }
  _system = system_in;
  _i = i_in;
  _j = j_in;
}


double BiAj::get_time_evolution( matrix< complex<double> > *UUt, matrix< complex<double> > *VVt )
{
  return _get_time_evolution(_i,_j,UUt,VVt );
}


double BiAj::_get_time_evolution( int i, int j, matrix< complex<double> > *UUt, matrix< complex<double> > *VVt )
{
  int size = UUt->get_ncol();
  complex<double> temp1,temp2;
  double temp;
  
  if (i==j)
    temp = -1.0;
  else
    temp = 0.0;

  for (int imu=0;imu<size;++imu){
    temp1 = (*VVt)(i,imu) - (*UUt)(i,imu);
    temp2 = conj((*VVt)(j,imu));
    temp += 2.0*real(temp1*temp2); 
  }

  return temp;
}


void BiAj::set_spvs()
{
  _spvs = new double[_size];
  double temp1,temp2;
  
  for (int imu=0;imu<_size;++imu){
    temp1 = (*_system->UU)(_i,imu) - (*_system->VV)(_i,imu);
    temp2 = (*_system->UU)(_j,imu) + (*_system->VV)(_j,imu);
    _spvs[imu] = temp1*temp2;
  }
  set_gsv();
}


AiAj::AiAj(int i_in, int j_in, ising1D * system_in) : local_obs< complex<double> >( system_in->get_size() ){
  if ((i_in < 0) || (j_in < 0) || (i_in > system_in->get_size()) || (j_in > system_in->get_size()) ){
    _ERROR_("indexes outside the chain",);
  }
  _system = system_in;
  _i = i_in;
  _j = j_in;
}


complex<double> AiAj::get_time_evolution( matrix< complex<double> > *UUt, matrix< complex<double> > *VVt )
{
  return _get_time_evolution(_i,_j,UUt,VVt );
}


complex<double> AiAj::_get_time_evolution( int i, int j, matrix< complex<double> > *UUt, matrix< complex<double> > *VVt )
{
  int size = UUt->get_ncol();
  complex<double> temp,temp1,temp2;
  
  temp = complex<double>(0.0,0.0);

  for (int imu=0;imu<size;++imu){
    temp1 = (*VVt)(i,imu) + (*UUt)(i,imu);
    temp2 = conj((*VVt)(j,imu) + (*UUt)(j,imu));
    temp += temp1*temp2;
  }

  return temp;
}


void AiAj::set_spvs()
{
  _spvs = new double[_size];
  for (int imu=0;imu<_size;++imu)
    _spvs[imu] = 0.0;
  set_gsv();
}


void AiAj::set_gsv()
{
  if (_i==_j)
    _gsv=1.0*double(_size);
  else
    _gsv=0.0;
}


BiBj::BiBj(int i_in, int j_in, ising1D * system_in) : local_obs< complex<double> >( system_in->get_size() ){
  if ((i_in < 0) || (j_in < 0) || (i_in > system_in->get_size()) || (j_in > system_in->get_size()) ){
    _ERROR_("indexes outside the chain",);
  }
  _system = system_in;
  _i = i_in;
  _j = j_in;
}


complex<double> BiBj::get_time_evolution( matrix< complex<double> > *UUt, matrix< complex<double> > *VVt )
{
  return _get_time_evolution(_i,_j,UUt,VVt );
}


complex<double> BiBj::_get_time_evolution( int i, int j, matrix< complex<double> > *UUt, matrix< complex<double> > *VVt )
{
  int size = UUt->get_ncol();
  complex<double> temp,temp1,temp2;
  
  temp = complex<double>(0.0,0.0);

  for (int imu=0;imu<size;++imu){
    temp1 = (*VVt)(i,imu) - (*UUt)(i,imu);
    temp2 = conj(-(*VVt)(j,imu) + (*UUt)(j,imu));
    temp += temp1*temp2;
  }

  return temp;
}

void BiBj::set_spvs()
{
  _spvs = new double[_size];
  for (int imu=0;imu<_size;++imu)
    _spvs[imu] = 0.0;
  set_gsv();
}

void BiBj::set_gsv()
{
  if (_i==_j)
    _gsv=-1.0*double(_size);
  else
    _gsv=0.0;
}

rho::rho( int i,int r, ising1D * system) : obs<double>( system->get_size())
{
  if ((i+r > system->get_size()) || (i < 0)){
    _ERROR_("indexes outside the chain",);
  }
  _system = system;
  _r = r;
  _i = i;
  _BiAjmatrix = NULL;
  _reduced_matrix = NULL;
  _full_matrix = NULL;
  set_BiAjpointers();
}


void rho::set_BiAjpointers()
{
  _BiAjmatrix = new matrix< BiAj* >(_r,_r);

  for( int irow=0;irow<_r;++irow)
    for ( int icol=0;icol<_r;++icol){
      (*_BiAjmatrix)(irow,icol) = new BiAj(_i+irow,_i+1+icol,_system);
      (*_BiAjmatrix)(irow,icol)->set_spvs();
    }
}


void rho::set_ensemble_average(double* nk)
{
  _reduced_matrix = new matrix<double>(_r,_r);
  
  for( int irow=0;irow<_r;++irow)
    for ( int icol=0;icol<_r;++icol)
      (*_reduced_matrix)(irow,icol) = (*_BiAjmatrix)(irow,icol)->get_ensemble_average(nk);
}


double rho::get_ensemble_average(int l)
{
  if (l<0) l=_r;
  if (l>_r) _ERROR_("distance greater than the initialized one",0.0);
 
  matrix<double> temp(l,l);
  double result;

  for( int irow=0;irow<l;++irow)
    for ( int icol=0;icol<l;++icol)
      temp(irow,icol) = (*_reduced_matrix)(irow,icol);

  result = temp.det();
  
  return result;
}

double rho::get_time_evolution( matrix< complex<double> > *UUt, matrix< complex<double> > * VVt )
{
  return 0.0;
}

void rho::set_time_evolution( matrix< complex<double> > *UUt, matrix< complex<double> > * VVt ){
  int i,j;
  _full_matrix = new matrix<complex<double> >(2*_r,2*_r);
  
  *_full_matrix = complex<double>(0.0,0.0);
  for (int irow=0;irow<2*_r;++irow){

    if (irow%2 == 0)
      i = _i + irow/2;
    else
      i = _i + (irow+1)/2;

    for (int icol=irow+1;icol<2*_r;++icol){
      if (icol%2 == 0)
	j = _i + icol/2;
      else
	j = _i + (icol+1)/2;
      
      if ((icol%2 == 0)&&(irow%2 == 0)){
	(*_full_matrix)(irow,icol) = BiBj::_get_time_evolution(i,j,UUt,VVt);
      }
      if ((irow%2 == 0)&&(icol%2 == 1)){
	(*_full_matrix)(irow,icol) = complex<double>(1.0,0.0)*BiAj::_get_time_evolution(i,j,UUt,VVt);
      }
      if ((irow%2 == 1)&&(icol%2 == 0)){
	(*_full_matrix)(irow,icol) = complex<double>(-1.0,0.0)*BiAj::_get_time_evolution(j,i,UUt,VVt);
      }
      if ((irow%2 == 1)&&(icol%2 == 1)){
	(*_full_matrix)(irow,icol) = AiAj::_get_time_evolution(i,j,UUt,VVt);
      }
      
      (*_full_matrix)(icol,irow) = - (*_full_matrix)(irow,icol);
    }
  }
}



double rho::get_time_evolution(int l)
{
  if (l<0) l=_r;
  if (l>_r) _ERROR_("distance greater than the initialized one",0.0);
 
  matrix< complex<double> > temp(2*l,2*l);
  complex<double> tResult;

  for( int irow=0;irow<2*l;++irow)
    for ( int icol=0;icol<2*l;++icol)
      temp(irow,icol) = (*_full_matrix)(irow,icol);
  
  //cout << "matrix";
  //temp.print();

  tResult = temp.det();
  return sqrt(real(tResult));
}



rho::~rho()
{
  for( int irow=0;irow<_r;++irow)
    for ( int icol=0;icol<_r;++icol)
      delete (*_BiAjmatrix)(irow,icol);
  delete _reduced_matrix;
  delete _BiAjmatrix;
  delete _full_matrix;
}



ising1D::ising1D( in_file* file, const string name, double (*JJgen)(int,ising1D*),double (*HHgen)(int,ising1D*) )
{
  size = _SIZE_;
  _h = _H_;
  _J = _J_;
  _gamma = _GAMMA_;
  _epsilon = _EPSILON_;
  _pbc = _PBC_;
  _seed = _SEED_;
  _from_center = _FROM_CENTER_;
  read( file, name );
  init( (*JJgen), (*HHgen) );
#ifdef DEBUG
  _ERROR_TRACKING_;
#endif
}


ising1D::ising1D(int size_in, double h, double J, double epsilon, double gamma, bool pbc, int seed)
{
  size = size_in;
  _J = J;
  _h = h;
  _gamma = gamma;
  _pbc = pbc;
  _epsilon = epsilon;
  _seed = seed;
  _from_center = _FROM_CENTER_;
  init( (&randomJJ),(&randomHH));
#ifdef DEBUG
  _ERROR_TRACKING_;
#endif
}


void ising1D::init( double (*JJgen)(int,ising1D*),double (*HHgen)(int,ising1D*))
{
  if (size<1)
    _ERROR_("Non valid value of size",);

  if (size%2 == 1)
    _WARNING_("Odd number of sites");

  if (_pbc)
    _MESSAGE_("You are using PBC");

  _hamiltonian = new matrix<double>(2*size,2*size);

  if (_seed < 0)
    _seed = time(NULL);

  rand_init( &_seed );
  
  UU = new matrix<double>(size,size);
  VV = new matrix<double>(size,size);
  _JJ = new double[size];
  _hh = new double[size];

  if (_from_center){
    _MESSAGE_("Disorder is generated from the center of the chain");
    for (int i=0;i<size/2;++i){
      //_hh[size/2 - i - 1] = _h * (1.0 + 2.0*_epsilon*(drand1()-0.5));
      //_hh[size/2 + i] = _h * (1.0 + 2.0*_epsilon*(drand1()-0.5));
      //_JJ[size/2 - i - 1] = _J * (1.0 + 2.0*_epsilon*(drand1()-0.5));
      //_JJ[size/2 + i] = _J * (1.0 + 2.0*_epsilon*(drand1()-0.5));
      _hh[size/2 - i - 1] = (*HHgen)(size/2 - i - 1,this);
      _hh[size/2 + i] = (*HHgen)(size/2 + i,this);
      _JJ[size/2 - i - 1] = (*JJgen)(size/2 - i - 1,this);
      _JJ[size/2 + i] = (*JJgen)(size/2 + i,this);
    }
  }
  else
    for (int i=0;i<size;++i){
      //_hh[i] = _h * (1.0 + 2.0*_epsilon*(drand1()-0.5));
      //_JJ[i] = _J * (1.0 + 2.0*_epsilon*(drand1()-0.5));
      _hh[i] = (*HHgen)(i,this);
      _JJ[i] = (*JJgen)(i,this);
    }

  for (int i=0;i<10;++i){
    cout << _hh[i] << "\t" << _JJ[i] << endl; 
  }  

  (*_hamiltonian) = get_hamiltonian();
  solve_diagonalization();

#ifdef DEBUG
  _ERROR_TRACKING_;
#endif
}



double randomHH(int i, ising1D* system){
  return system->get_h() * (1.0 + 2.0*system->get_epsilon()*(drand1()-0.5));
}


double randomJJ(int i, ising1D* system){
  return system->get_J() * (1.0 + 2.0*system->get_epsilon()*(drand1()-0.5));
}


int ising1D::get_size()
{
  return size;
}

double ising1D::get_J()
{
  return _J;
}

double ising1D::get_h()
{
  return _h;
}

double ising1D::get_epsilon()
{
  return _epsilon;
}


void ising1D::solve_diagonalization()
{
  double * tempeigenval = new double[2*size];
  double * eigenval = new double[size];
  matrix<double> temph(2*size,2*size);

  temph = get_hamiltonian();

  tempeigenval = temph.diagonalize(true);
  if (not _pbc)
    check( tempeigenval, &temph );
#ifdef DEBUG
  _ERROR_TRACKING_;
#endif

  for (int i=0;i<size;++i){
    eigenval[i] = tempeigenval[size+i];
    // uncomment to print the eigenvalues
    //cout << eigenval[i] << endl;
  }

  e = new energy(eigenval, size);

  for (int irow=0;irow<size;++irow)
    for (int icol=0;icol<size;++icol){
      (*VV)(irow,icol) = temph(size+icol,size+irow);
      (*UU)(irow,icol) = temph(size+icol,irow);
    }

  delete [] eigenval;
  delete [] tempeigenval;

#ifdef DEBUG
  _ERROR_TRACKING_;
#endif
  return;
}

void ising1D::check( double* eigenval,  matrix<double> * eigvect )
{
  bool docor = false;
  int size = eigvect->get_ncol()/2;
  double temp;

  for (int irow=0;irow<size;++irow)
    for (int icol=0;icol<size;++icol){
      temp = abs(abs( (*eigvect)(irow,icol)) - abs( (*eigvect)(2*size-1-irow,icol+size)) );
      if ( temp > 1.0e-11){
	if (irow == size-1)
	  docor = true;
	else
	  {
	    cout << "discrepancy " << temp << " val" << abs( (*eigvect)(irow,icol)) << endl;
	    _ERROR_("some error in the diagonalization",);
	  }
      }
      temp = abs( abs( (*eigvect)(size +irow,icol) ) - abs( (*eigvect)(size-irow-1, icol+size ))); 
      if ( temp >1.0e-10){
	if (irow == 0)
	  docor = true;
	else
	  {
	    cout << "discrepancy " << temp << " val" << abs( (*eigvect)(irow,icol)) << endl;
	    _ERROR_("some error in the diagonalization",);
	  }
      }
    }
    
  if (docor){
    _WARNING_("manual correction of eigenvectors");
    double gamma = ( (*eigvect)(size-1,size) + (*eigvect)(size,0) ) / ( (*eigvect)(size,size) - (*eigvect)(size-1,0) );
    double beta = 1.0/sqrt(1.0+gamma*gamma);
    double alpha = gamma*beta;
    double temp1,temp2;
    for (int i=0;i<2*size;++i){
      temp1 = (*eigvect)(size-1,i);
      temp2 = (*eigvect)(size,i);
      (*eigvect)(size,i) = alpha*temp1 + beta*temp2;
      (*eigvect)(size-1,i) = -beta*temp1 + alpha*temp2;
    }

    // calculation of the eigenvalue associated to the eigenvectors size-1 and size
    temp1 = 0.0;
    temp2 = 0.0;
    for (int j=0;j<2*size;++j){
      double t1=0.0,t2=0.0;
      for (int i=0;i<2*size;++i){
	t1 += (*_hamiltonian)(j,i) * (*eigvect)(size-1,i);
	t2 += (*_hamiltonian)(j,i) * (*eigvect)(size,i);
      }
      temp1 += t1*(*eigvect)(size-1,j);
      temp2 += t2*(*eigvect)(size,j);
    }

    double temp;
    if ( temp1>0 ){
      // inversion of the eigenvectors
      eigenval[size-1] = temp2;
      eigenval[size] = temp1;
      for (int i=0;i<2*size;++i){
	temp = (*eigvect)(size,i);
	(*eigvect)(size,i) = (*eigvect)(size-1,i);
	(*eigvect)(size-1,i) = temp;
      }
    }
    else{
      eigenval[size-1] = temp1;
      eigenval[size] = temp2; 
    }
  }
}


matrix<double> ising1D::get_hamiltonian()
{
  matrix<double> AA(size,size);
  matrix<double> BB(size,size);
  matrix<double> result(2*size,2*size);

  AA = get_matrix_A();
  BB = get_matrix_B();

  for (int irow=0;irow<size;++irow){
    for (int icol=0;icol<size;++icol){
      result(irow,icol) = AA(irow,icol);
      result(irow+size,icol+size) = - AA(irow,icol);
      result(irow,size+icol) = BB(irow,icol);
      result(irow+size,icol) = - BB(irow,icol);
    }
  }

#ifdef DEBUG
  _ERROR_TRACKING_;
#endif
  return result;
}


matrix<double> ising1D::get_matrix_A()
{
  matrix<double>  AA(size,size);
  
  AA = 0.0;
  
  for (int irow=0;irow<size;++irow){
    AA(irow,irow) = - _hh[irow];
    if (irow+1<size){
      AA(irow,irow+1) = - _JJ[irow]*0.5;
      AA(irow+1,irow) = - _JJ[irow]*0.5;
    }
  }
  
  if (_pbc){
    AA(size-1,0) += _JJ[size-1]*0.5;
    AA(0,size-1) += _JJ[size-1]*0.5;
  }

#ifdef DEBUG
  _ERROR_TRACKING_;
#endif
  return AA;
}


matrix<double> ising1D::get_matrix_B()
{
  matrix<double> BB(size,size);
  //double temp = (1.0-gamma)/2.0;
  double temp = _gamma/2.0;
  
  BB = 0.0;
  
  for (int irow=0;irow<size-1;++irow){
    BB(irow,irow+1) = - _JJ[irow]*temp;
    BB(irow+1,irow) =  _JJ[irow]*temp;
    }

  // supposing even NF
  if (_pbc){
    BB(size-1,0) += _JJ[size-1]*temp;
    BB(0,size-1) += - _JJ[size-1]*temp;
  }

#ifdef DEBUG
  _ERROR_TRACKING_;
#endif
  return BB;
}


ising1D::~ising1D()
{
  delete [] _JJ;
  delete [] _hh;
  delete _hamiltonian;
  delete UU;
  delete VV;
  delete e;
}


/*************************************************/
void ising1D::print()
{
  cout << "size = " << size << endl;
  cout << "h = " << _h << endl;
  cout << "J = " << _J << endl;
  cout << "gamma = " << _gamma << endl;
  cout << "epsilon = " << _epsilon << endl;
  cout << "pbc = " << _pbc << endl;
}


void ising1D::write(out_file* file)
{
  file->write_tag("system");
  file->write_parameter("size", size );
  file->write_parameter("h", _h );
  file->write_parameter("J", _J );
  file->write_parameter("epsilon", _epsilon );
  file->write_parameter("gamma", _gamma );
  file->write_parameter("pbc", _pbc );
  file->write_parameter("seed",_seed );
#ifdef DEBUG
  _ERROR_TRACKING_;
#endif
}

/*************************************************/

void ising1D::read( in_file* file,const string systemname)
{
  string name, data;

  if (file->find_tag(systemname)<0){
    _ERROR_("not able to find tag",);
  }
#ifdef DEBUG
  _ERROR_TRACKING_;
#endif
  while (  file->read_data(name,data)>0 ){
    if (name=="size")
      istringstream(data) >> size;
    else if (name=="h")
      istringstream(data) >> _h;
    else if (name=="J")
      istringstream(data) >> _J;
    else if (name=="gamma")
      istringstream(data) >> _gamma;
    else if (name=="epsilon")
      istringstream(data) >> _epsilon;
    else if (name=="pbc")
      istringstream(data) >> _pbc;
    else if (name=="seed")
      istringstream(data) >> _seed;
    else if (name=="from_center")
      istringstream(data) >> _from_center;
  }
#ifdef DEBUG
  _ERROR_TRACKING_;
#endif
}


quench::quench( int size_in, double h0, double h, double J, double epsilon, double gamma, bool pbc)
{
  system0 = new ising1D(size_in, h0, J, epsilon, gamma, pbc);
  system = new ising1D(size_in, h,  J, epsilon, gamma, pbc);
  init();
#ifdef DEBUG
  _ERROR_TRACKING_;
#endif
}


quench::~quench()
{
  delete system0;
  delete system;
  delete UUt;
  delete VVt;
}


quench::quench( in_file* file,double (*JJgen0)(int,ising1D*),double (*HHgen0)(int,ising1D*),double (*JJgen)(int,ising1D*),double (*HHgen)(int,ising1D*))
{
  system0 = new ising1D(file,"system0",JJgen0,HHgen0);
  system = new ising1D(file,"system",JJgen,HHgen);
  init();
#ifdef DEBUG
  _ERROR_TRACKING_;
#endif
}


void quench::init()
{
  size = system->size;
  gge = NULL;
  UUt = new matrix< complex<double> >(size,size);
  VVt = new matrix< complex<double> >(size,size);
#ifdef DEBUG
  _ERROR_TRACKING_;
#endif
}


int quench::get_size()
{
  return size;
}


void quench::set_gge_occupations()
{
  matrix<double> temp(size,size);
  gge = new double[size];
  
  temp = system->UU->transpose()* *(system0->VV) + system->VV->transpose()* *(system0->UU);
  
  for (int mu=0;mu<size;++mu){
    gge[mu] = 0.0;
    for (int icol=0;icol<size;++icol)
      gge[mu] += pow(abs(temp(mu,icol)),2);
    }

  //  for (int mu=0;mu<size;++mu)
  //cout << gge[mu] << endl;
}


void quench::set_time_evolution( const double time )
{
  matrix< complex<double> > te(size,size);
  matrix< complex<double> > tec(size,size);
  matrix< complex<double> > AA(size,size);
  matrix< complex<double> > BB(size,size);

  te = get_evolution_matrix(time);
  tec = te.conjugate();
  
  AA = *(system->UU) * te * system->UU->daga();
  AA = AA + (system->VV->conjugate() * tec * system->VV->transpose());
  BB = *(system->VV) * te * system->UU->daga();
  BB = BB + (system->UU->conjugate() * tec * system->VV->transpose());
  
  (*UUt) = AA * *(system0->UU) + BB.conjugate() * *(system0->VV);
  (*VVt) = BB * *(system0->UU) + AA.conjugate() * *(system0->VV);

  /* uncomment to check the unitarity of time evolution
     cout << "check" << endl;
     (UUt->daga() * *(UUt) + VVt->daga() * *(VVt)).print();
     (UUt->daga() * *(UUt) + VVt->daga() * *(VVt)).print();
  */
  
#ifdef DEBUG
  _ERROR_TRACKING_;
#endif
}


matrix< complex<double> > quench::get_evolution_matrix( const double time )
{
  matrix< complex<double> > temp(size,size);
  complex<double> imath(0.0,1.0);
  
  temp = complex<double>(0.0,0.0);
  for (int i=0;i<size;++i)
      temp(i,i) = exp(- 2.0*system->e->get_spv(i) * time * imath );

#ifdef DEBUG
  _ERROR_TRACKING_;
#endif
  return temp;
}

double quench::get_calpha2( state* s)
{
  matrix<double> u1(size,size),v1(size,size),temp(size,size);
 
  for (int icol=0;icol<size;++icol){
    if (s->conf[icol]){
      for (int irow=0;irow<size;++irow){
	u1(irow,icol) = (*system->VV)(irow,icol);
	v1(irow,icol) = (*system->UU)(irow,icol);
      }}
    else
      for (int irow=0;irow<size;++irow){
	u1(irow,icol) = (*system->UU)(irow,icol);
	v1(irow,icol) = (*system->VV)(irow,icol);
      }
  }

  temp = system0->UU->daga() * u1 + system0->VV->daga() * v1;
      
  return abs(temp.det());
}
