
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


energy::energy(FPType * data,int size_in) : local_obs<FPType>( data,size_in ){
}


FPType energy::get_time_evolution( matrix< complex<FPType> > *UUt, matrix< complex<FPType> > *VVt )
{
  return 0.0;
}

localtmag::localtmag(int site, ising1D * system_in) : local_obs<FPType>( system_in->get_size() ){
  _system = system_in;
  _site = site;
}


FPType localtmag::_get_time_evolution(int site, matrix< complex<FPType> > *UUt, matrix< complex<FPType> > *VVt )
{
  int size = UUt->get_ncol();
  FPType temp = 0.0;
  
  for (int icol=0;icol<size;++icol)
    temp += norm( (*VVt)(site,icol) );
  
  return 2.0*temp-1.0;
}


FPType localtmag::get_time_evolution( matrix< complex<FPType> > *UUt, matrix< complex<FPType> > *VVt )
{
  return _get_time_evolution( _site, UUt, VVt );
}


void localtmag::set_spvs()
{
  _spvs = new FPType[_size];
  
  _gsv = -1.0;
  for (int mu=0;mu<_size;++mu){
    _spvs[mu] = (*_system->UU)(_site,mu)*(*_system->UU)(_site,mu)-(*_system->VV)(_site,mu)*(*_system->VV)(_site,mu);
    _gsv += 2.0*(*_system->VV)(_site,mu)*(*_system->VV)(_site,mu);
  }
}


/********* to modify ************/
cidagacj::cidagacj(int i, int j, ising1D * system_in) : local_obs< complex<FPType> >( system_in->get_size() ){
  _system = system_in;
  _i = i;
  _j = j;
}


complex<FPType> cidagacj::_get_time_evolution(int i,int j, matrix< complex<FPType> > *UUt, matrix< complex<FPType> > *VVt )
{
  int size = UUt->get_ncol();
  complex<FPType> temp = 0.0;
  
  for (int icol=0;icol<size;++icol)
    temp += (*VVt)(i,icol) * conj( (*VVt)(j,icol) );
  
  return temp;
}


complex<FPType> cidagacj::get_time_evolution( matrix< complex<FPType> > *UUt, matrix< complex<FPType> > *VVt )
{
  return _get_time_evolution( _i,_j, UUt, VVt );
}


void cidagacj::set_spvs()
{
  _spvs = new FPType[_size];
  
  _gsv = 0.0;
  for (int mu=0;mu<_size;++mu){
    _spvs[mu] = ((*_system->UU)(_i,mu)*(*_system->UU)(_j,mu)-(*_system->VV)(_i,mu)*(*_system->VV)(_j,mu))/2.0;
    _gsv += (*_system->VV)(_i,mu)*(*_system->VV)(_j,mu);
  }

}
/********************************/









tmag::tmag(ising1D * system_in) : local_obs<FPType>( system_in->get_size() ){
  _system = system_in;
}


FPType tmag::_get_time_evolution( matrix< complex<FPType> > *UUt, matrix< complex<FPType> > *VVt )
{
  int size = UUt->get_ncol();
  FPType temp = 0.0;
  
  for (int irow=0;irow<size;++irow)
    for (int icol=0;icol<size;++icol)
      //      temp += pow( abs( (*VVt)(irow,icol) ),2);
      temp += norm( (*VVt)(irow,icol) );
  
  return 2.0*temp-FPType(size);
}


FPType tmag::get_time_evolution( matrix< complex<FPType> > *UUt, matrix< complex<FPType> > *VVt )
{
  return _get_time_evolution( UUt, VVt );
}


void tmag::set_spvs()
{
  _spvs = new FPType[_size];
  
  for (int mu=0;mu<_size;++mu){
    _spvs[mu] = 0.0;
    for (int irow=0;irow<_size;++irow)
      _spvs[mu] += (*_system->VV)(irow,mu)*(*_system->VV)(irow,mu);
    _spvs[mu] = 1.0-2.0*_spvs[mu];
  }
  set_gsv();
}


BiAj::BiAj(int i_in, int j_in, ising1D * system_in) : local_obs< FPType >( system_in->get_size() ){
  if ((i_in < 0) || (j_in < 0) || (i_in > system_in->get_size()) || (j_in > system_in->get_size()) ){
    _ERROR_("indexes outside the chain",);
  }
  _system = system_in;
  _i = i_in;
  _j = j_in;
}


FPType BiAj::get_time_evolution( matrix< complex<FPType> > *UUt, matrix< complex<FPType> > *VVt )
{
  return _get_time_evolution(_i,_j,UUt,VVt );
}


FPType BiAj::_get_time_evolution( int i, int j, matrix< complex<FPType> > *UUt, matrix< complex<FPType> > *VVt )
{
  int size = UUt->get_ncol();
  complex<FPType> temp1,temp2;
  FPType temp;
  
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
  _spvs = new FPType[_size];
  FPType temp1,temp2;
  
  for (int imu=0;imu<_size;++imu){
    temp1 = (*_system->UU)(_i,imu) - (*_system->VV)(_i,imu);
    temp2 = (*_system->UU)(_j,imu) + (*_system->VV)(_j,imu);
    _spvs[imu] = temp1*temp2;
  }
  set_gsv();
}


AiAj::AiAj(int i_in, int j_in, ising1D * system_in) : local_obs< complex<FPType> >( system_in->get_size() ){
  if ((i_in < 0) || (j_in < 0) || (i_in > system_in->get_size()) || (j_in > system_in->get_size()) ){
    _ERROR_("indexes outside the chain",);
  }
  _system = system_in;
  _i = i_in;
  _j = j_in;
}


complex<FPType> AiAj::get_time_evolution( matrix< complex<FPType> > *UUt, matrix< complex<FPType> > *VVt )
{
  return _get_time_evolution(_i,_j,UUt,VVt );
}


complex<FPType> AiAj::_get_time_evolution( int i, int j, matrix< complex<FPType> > *UUt, matrix< complex<FPType> > *VVt )
{
  int size = UUt->get_ncol();
  complex<FPType> temp,temp1,temp2;
  
  temp = complex<FPType>(0.0,0.0);

  for (int imu=0;imu<size;++imu){
    temp1 = (*VVt)(i,imu) + (*UUt)(i,imu);
    temp2 = conj((*VVt)(j,imu) + (*UUt)(j,imu));
    temp += temp1*temp2;
  }

  return temp;
}


void AiAj::set_spvs()
{
  _spvs = new FPType[_size];
  for (int imu=0;imu<_size;++imu)
    _spvs[imu] = 0.0;
  set_gsv();
}


void AiAj::set_gsv()
{
  if (_i==_j)
    _gsv=1.0*FPType(_size);
  else
    _gsv=0.0;
}


BiBj::BiBj(int i_in, int j_in, ising1D * system_in) : local_obs< complex<FPType> >( system_in->get_size() ){
  if ((i_in < 0) || (j_in < 0) || (i_in > system_in->get_size()) || (j_in > system_in->get_size()) ){
    _ERROR_("indexes outside the chain",);
  }
  _system = system_in;
  _i = i_in;
  _j = j_in;
}


complex<FPType> BiBj::get_time_evolution( matrix< complex<FPType> > *UUt, matrix< complex<FPType> > *VVt )
{
  return _get_time_evolution(_i,_j,UUt,VVt );
}


complex<FPType> BiBj::_get_time_evolution( int i, int j, matrix< complex<FPType> > *UUt, matrix< complex<FPType> > *VVt )
{
  int size = UUt->get_ncol();
  complex<FPType> temp,temp1,temp2;
  
  temp = complex<FPType>(0.0,0.0);

  for (int imu=0;imu<size;++imu){
    temp1 = (*VVt)(i,imu) - (*UUt)(i,imu);
    temp2 = conj(-(*VVt)(j,imu) + (*UUt)(j,imu));
    temp += temp1*temp2;
  }

  return temp;
}

void BiBj::set_spvs()
{
  _spvs = new FPType[_size];
  for (int imu=0;imu<_size;++imu)
    _spvs[imu] = 0.0;
  set_gsv();
}

void BiBj::set_gsv()
{
  if (_i==_j)
    _gsv=-1.0*FPType(_size);
  else
    _gsv=0.0;
}


rho::rho( int i,int r, ising1D * system) : obs<FPType>( system->get_size())
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


rho::rho(const rho& source) : obs<FPType>( source )
{
  /* ATTENTION */
  /* they point to the same system and same _BiAjmatrix elements (not matrix) */
  _system = source._system;
  _r = source._r;
  _i = source._i;

  if (source._BiAjmatrix == NULL)
    _BiAjmatrix = NULL;
  else{
    _BiAjmatrix = new matrix< BiAj* >(_r,_r);
    for( int irow=0;irow<_r;++irow)
      for ( int icol=0;icol<_r;++icol)
	(*_BiAjmatrix)(irow,icol) = new BiAj(*( (*source._BiAjmatrix)(irow,icol) ));
  }

  if (source._reduced_matrix == NULL)
    _reduced_matrix = NULL;
  else
    _reduced_matrix = new matrix<FPType>( *(source._reduced_matrix) );

  if (source._full_matrix == NULL)
    _full_matrix = NULL;
  else
    _full_matrix = new matrix< complex<FPType> >( *(source._full_matrix));
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


void rho::set_ensemble_average(FPType* nk)
{
  _reduced_matrix = new matrix<FPType>(_r,_r);
  
  for( int irow=0;irow<_r;++irow)
    for ( int icol=0;icol<_r;++icol)
      (*_reduced_matrix)(irow,icol) = (*_BiAjmatrix)(irow,icol)->get_ensemble_average(nk);
}


FPType rho::get_ensemble_average(int l)
{
  if (l<0) l=_r;
  if (l>_r) _ERROR_("distance greater than the initialized one",0.0);
 
  matrix<FPType> temp(l,l);
  FPType result;

  for( int irow=0;irow<l;++irow)
    for ( int icol=0;icol<l;++icol)
      temp(irow,icol) = (*_reduced_matrix)(irow,icol);

  result = temp.det();
  
  return result;
}

FPType rho::get_time_evolution( matrix< complex<FPType> > *UUt, matrix< complex<FPType> > * VVt )
{
  _ERROR_("wrong use of rho object",0.0);
}

void rho::set_time_evolution( matrix< complex<FPType> > *UUt, matrix< complex<FPType> > * VVt ){
  int i,j;

  if (_full_matrix == NULL)
    _full_matrix = new matrix<complex<FPType> >(2*_r,2*_r);
  
  *_full_matrix = complex<FPType>(0.0,0.0);
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
	(*_full_matrix)(irow,icol) = complex<FPType>(1.0,0.0)*BiAj::_get_time_evolution(i,j,UUt,VVt);
      }
      if ((irow%2 == 1)&&(icol%2 == 0)){
	(*_full_matrix)(irow,icol) = complex<FPType>(-1.0,0.0)*BiAj::_get_time_evolution(j,i,UUt,VVt);
      }
      if ((irow%2 == 1)&&(icol%2 == 1)){
	(*_full_matrix)(irow,icol) = AiAj::_get_time_evolution(i,j,UUt,VVt);
      }
      
      (*_full_matrix)(icol,irow) = - (*_full_matrix)(irow,icol);
    }
  }
}



FPType rho::get_time_evolution(int l)
{
  if (l<0) l=_r;
  if (l>_r) _ERROR_("distance greater than the initialized one",0.0);
 
  matrix< complex<FPType> > temp(2*l,2*l);
  complex<FPType> tResult;

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




rhozz::rhozz( int i,int r, ising1D * system) : obs<FPType>( system->get_size())
{
  if ((i > system->get_size()) || (i < 0)){
    _ERROR_("indexes outside the chain",);
  }
  if ( i+r >= system->get_size() ){
    _MESSAGE_("second site is outside the chain");
  }
  _system = system;
  _r = r;
  _i = i;
  int secondSite = _i+_r;

  if (secondSite >= system->get_size() ){
    secondSite = secondSite % system->get_size();
    _MESSAGE_("second site has been put inside");
    //_r = secondSite - _i;
  }

  _BlAl = new BiAj(i,i,system);
  _BmAm = new BiAj(secondSite, secondSite , system);
  _BmAl = new BiAj(secondSite, i, system);
  _BlAm = new BiAj(i, secondSite, system);
  _AlAm = new AiAj(i, secondSite, system);
  _BlBm = new BiBj(i, secondSite, system);
}

void rhozz::set_ensemble_average()
{
  _BlAl->set_spvs();
  _BmAm->set_spvs();
  _BmAl->set_spvs();
  _BlAm->set_spvs();
  _AlAm->set_spvs();
  _BlBm->set_spvs();
}

FPType rhozz::get_time_evolution( matrix< complex<FPType> > *UUt, matrix< complex<FPType> > *VVt )
{
  FPType result;
  
  result = _BlAl->get_time_evolution(UUt,VVt)* _BmAm->get_time_evolution(UUt,VVt);
  result -= _BmAl->get_time_evolution(UUt,VVt)* _BlAm->get_time_evolution(UUt,VVt);
  result -= real( _AlAm->get_time_evolution(UUt,VVt)* _BlBm->get_time_evolution(UUt,VVt) );

  return result;
}



FPType rhozz::get_ensemble_average( FPType *nk )
{
  FPType result;
  
  result = _BlAl->get_ensemble_average(nk) * _BmAm->get_ensemble_average(nk);
  result -= _BmAl->get_ensemble_average(nk) * _BlAm->get_ensemble_average(nk);
  if (_r==0)
    result += 1.0;
  
  return result;
}



rhozz::~rhozz()
{
  delete _BlAl;
  delete _BmAm;
  delete _BmAl;
  delete _BlAm;
  delete _AlAm;
  delete _BlBm;
}



ising1D::ising1D( in_file* file, const string name, FPType (*JJgen)(int,ising1D*),FPType (*HHgen)(int,ising1D*) )
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

  init();
  gendis( (*JJgen), (*HHgen) );

#ifdef DEBUG
  _ERROR_TRACKING_;
#endif
}


ising1D::ising1D(int size_in, FPType h, FPType J, FPType epsilon, FPType gamma, bool pbc, int seed)
{
  size = size_in;
  _J = J;
  _h = h;
  _gamma = gamma;
  _pbc = pbc;
  _epsilon = epsilon;
  _seed = seed;
  _from_center = _FROM_CENTER_;

  init();
  gendis( (&randomJJ),(&randomHH));

#ifdef DEBUG
  _ERROR_TRACKING_;
#endif
}


void ising1D::init()
{
  if (size<1)
    _ERROR_("Non valid value of size",);

  if (size%2 == 1)
    _WARNING_("Odd number of sites");

  if (_pbc)
    _MESSAGE_("You are using PBC");

  _hamiltonian = new matrix<FPType>(2*size,2*size);

  if (_seed < 0)
    _seed = time(NULL);

  rand_init( &_seed );
  
  UU = new matrix<FPType>(size,size);
  VV = new matrix<FPType>(size,size);
  _JJ = new FPType[size];
  _hh = new FPType[size];
  e = NULL;

#ifdef DEBUG
  _ERROR_TRACKING_;
#endif
}


void ising1D::gendis( FPType (*JJgen)(int,ising1D*),FPType (*HHgen)(int,ising1D*))
{
  if (_from_center){
    _MESSAGE_("Disorder is generated from the center of the chain");
    for (int i=0;i<size/2;++i){
      _hh[size/2 - i - 1] = (*HHgen)(size/2 - i - 1,this);
      _hh[size/2 + i] = (*HHgen)(size/2 + i,this);
      _JJ[size/2 - i - 1] = (*JJgen)(size/2 - i - 1,this);
      _JJ[size/2 + i] = (*JJgen)(size/2 + i,this);
    }
  }
  else
    for (int i=0;i<size;++i){
      _hh[i] = (*HHgen)(i,this);
      _JJ[i] = (*JJgen)(i,this);
    }

  (*_hamiltonian) = get_hamiltonian();

  solve_diagonalization();
#ifdef DEBUG
  _ERROR_TRACKING_;
#endif
}




FPType randomHH(int i, ising1D* system){
  return system->get_h() * (1.0 + 2.0*system->get_epsilon()*(drand1()-0.5));
}


FPType randomJJ(int i, ising1D* system){
  return system->get_J() * (1.0 + 2.0*system->get_epsilon()*(drand1()-0.5));
}


int ising1D::get_size()
{
  return size;
}

FPType ising1D::get_J()
{
  return _J;
}

FPType ising1D::get_h()
{
  return _h;
}

FPType ising1D::get_epsilon()
{
  return _epsilon;
}


void ising1D::solve_diagonalization()
{
  FPType * tempeigenval = new FPType[2*size];
  FPType * eigenval = new FPType[size];
  matrix<FPType> temph(2*size,2*size);

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

  if (e != NULL){
    delete e;
    e = NULL;
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

void ising1D::check( FPType* eigenval,  matrix<FPType> * eigvect )
{
  bool docor = false;
  int size = eigvect->get_ncol()/2;
  FPType temp;

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
    FPType gamma = ( (*eigvect)(size-1,size) + (*eigvect)(size,0) ) / ( (*eigvect)(size,size) - (*eigvect)(size-1,0) );
    FPType beta = 1.0/sqrt(1.0+gamma*gamma);
    FPType alpha = gamma*beta;
    FPType temp1,temp2;
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
      FPType t1=0.0,t2=0.0;
      for (int i=0;i<2*size;++i){
	t1 += (*_hamiltonian)(j,i) * (*eigvect)(size-1,i);
	t2 += (*_hamiltonian)(j,i) * (*eigvect)(size,i);
      }
      temp1 += t1*(*eigvect)(size-1,j);
      temp2 += t2*(*eigvect)(size,j);
    }

    FPType temp;
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


matrix<FPType> ising1D::get_hamiltonian()
{
  matrix<FPType> AA(size,size);
  matrix<FPType> BB(size,size);
  matrix<FPType> result(2*size,2*size);

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


matrix<FPType> ising1D::get_matrix_A()
{
  matrix<FPType>  AA(size,size);
  
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


matrix<FPType> ising1D::get_matrix_B()
{
  matrix<FPType> BB(size,size);
  //FPType temp = (1.0-gamma)/2.0;
  FPType temp = _gamma/2.0;
  
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


matrix< complex<FPType> > ising1D::temp( const FPType time )
{
  matrix< complex<FPType> > te = getEvolutionMatrix(time);
  matrix< complex<FPType> > tec = te.conjugate();
  matrix< complex<FPType> > AA(size,size);
  //matrix< complex<FPType> > BB(size,size);
  
  //AA = gemm( *(system->UU),'N',gemm(te,'N', *system->UU,'D'),'N');
  //  AA = AA + gemm(system->VV->conjugate(),'N',gemm(tec,'N',*system->VV,'T'),'N');
  AA = gemm(VV->conjugate(),'N',gemm(tec,'N',*VV,'T'),'N');
  return AA;
  //BB = gemm(*(system->VV),'N', gemm(te,'N',*system->UU,'D'),'N');
  //BB = BB + gemm( system->UU->conjugate(),'N',gemm(tec,'N',*system->VV,'T'),'N');
  

#ifdef DEBUG
  _ERROR_TRACKING_;
#endif
}

matrix< complex<FPType> > ising1D::getEvolutionMatrix( const FPType time )
{
  matrix< complex<FPType> > temp(size,size);
  complex<FPType> imath(0.0,1.0);
  FPType two=2.0;
  
  temp = complex<FPType>(0.0,0.0);
  for (int i=0;i<size;++i)
    temp(i,i) = exp(- two*e->get_spv(i) * time * imath );

#ifdef DEBUG
  _ERROR_TRACKING_;
#endif
  return temp;
}



quench::quench( int size_in, FPType h0, FPType h, FPType J, FPType epsilon, FPType gamma, bool pbc)
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


quench::quench( in_file* file,FPType (*JJgen0)(int,ising1D*),FPType (*HHgen0)(int,ising1D*),FPType (*JJgen)(int,ising1D*),FPType (*HHgen)(int,ising1D*))
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
  UUt = new matrix< complex<FPType> >(size,size);
  VVt = new matrix< complex<FPType> >(size,size);
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
  matrix<FPType> temp(size,size);
  gge = new FPType[size];
  
  // modified
  //  temp = system->UU->transpose()* *(system0->VV) + system->VV->transpose()* *(system0->UU);

  temp = gemm(*system->UU,'T',*system0->VV,'N')+gemm(*system->VV,'T',*system0->UU,'N');
  
  for (int mu=0;mu<size;++mu){
    gge[mu] = 0.0;
    for (int icol=0;icol<size;++icol)
      gge[mu] += pow(abs(temp(mu,icol)),2);
    }

  //  for (int mu=0;mu<size;++mu)
  //cout << gge[mu] << endl;
}


void quench::set_time_evolution( const FPType time, matrix<complex<FPType> > * UU, matrix<complex<FPType> > * VV )
{
  matrix< complex<FPType> > te = system->getEvolutionMatrix(time);
  matrix< complex<FPType> > tec = te.conjugate();
  matrix< complex<FPType> > AA(size,size);
  matrix< complex<FPType> > BB(size,size);
  
  AA = gemm( *(system->UU),'N',gemm(te,'N', *system->UU,'D'),'N');
  AA = AA + gemm(system->VV->conjugate(),'N',gemm(tec,'N',*system->VV,'T'),'N');
  BB = gemm(*(system->VV),'N', gemm(te,'N',*system->UU,'D'),'N');
  BB = BB + gemm( system->UU->conjugate(),'N',gemm(tec,'N',*system->VV,'T'),'N');
  
  (*UU) = gemm(AA,'N',*(system0->UU),'N') + gemm(BB.conjugate(),'N', *(system0->VV), 'N');
  (*VV) = gemm(BB,'N',*(system0->UU),'N') + gemm(AA.conjugate(),'N', *(system0->VV), 'N');

#ifdef DEBUG
  _ERROR_TRACKING_;
#endif
}

void quench::set_time_evolution( const FPType time )
{
  set_time_evolution(time, UUt,VVt);
}


//matrix< complex<FPType> > quench::get_evolution_matrix( const FPType time )
//{
//  matrix< complex<FPType> > temp(size,size);
//  complex<FPType> imath(0.0,1.0);
//  FPType two=2.0;
//  
//  temp = complex<FPType>(0.0,0.0);
//  for (int i=0;i<size;++i)
//    temp(i,i) = exp(- two*system->e->get_spv(i) * time * imath );
//
//#ifdef DEBUG
//  _ERROR_TRACKING_;
//#endif
//  return temp;
//}

FPType quench::get_calpha2( state* s)
{
  matrix<FPType> u1(size,size),v1(size,size),temp(size,size);
 
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
