
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
using namespace std;


energy::energy(double * data,int size_in) : local_obs<double>( data,size_in ){
}


double energy::get_time_evolution( matrix< complex<double> > *UUt, matrix< complex<double> > *VVt )
{
  return 0.0;
}



tmag::tmag(ising1D * system_in) : local_obs<double>( system_in->size ){
  _system = system_in;
}



double tmag::_get_time_evolution( matrix< complex<double> > *UUt, matrix< complex<double> > *VVt )
{
  int size = UUt->get_ncol();
  double temp = 0.0;
  
  for (int irow=0;irow<size;++irow)
    for (int icol=0;icol<size;++icol)
      temp += pow(abs( (*VVt)(irow,icol) ),2);
  
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
      _spvs[mu] += pow( (*_system->VV)(irow,mu),2);
    _spvs[mu] = 1.0-2.0*_spvs[mu];
  }
  set_gsv();
}


BiAj::BiAj(int i_in, int j_in, ising1D * system_in) : local_obs< double >( system_in->size ){
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


AiAj::AiAj(int i_in, int j_in, ising1D * system_in) : local_obs< complex<double> >( system_in->size ){
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
  
  temp = (0.0,0.0);

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


BiBj::BiBj(int i_in, int j_in, ising1D * system_in) : local_obs< complex<double> >( system_in->size ){
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
  
  temp = (0.0,0.0);

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

rho::rho( int i,int r, ising1D * system) : obs<double>( system->size)
{
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
  if (l>_r) _ERROR_("distance greater than the initialized one");
 
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
  
  *_full_matrix = (0.0,0.0);
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
  if (l>_r) _ERROR_("distance greater than the initialized one");
 
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



ising1D::ising1D( in_file* file, const string name)
{
  size = _SIZE_;
  h = _H_;
  J = _J_;
  gamma = _GAMMA_;
  epsilon = _EPSILON_;
  pbc = _PBC_;
  read( file, name );
  init();
#ifdef DEBUG
  _ERROR_TRACKING_;
#endif
}


ising1D::ising1D(int size_in, double h_in, double J_in, double epsilon_in, double gamma_in, bool pbc_in)
{
  size = size_in;
  J = J_in;
  h = h_in;
  gamma = gamma_in;
  pbc = pbc_in;
  epsilon = epsilon_in;
  init();
#ifdef DEBUG
  _ERROR_TRACKING_;
#endif
}


void ising1D::init()
{
  if (size<1)
    _ERROR_("Non valid value of size");

  if (pbc)
    _WARNING_("You are using PBC");

  hamiltonian = new matrix<double>(2*size,2*size);
  
  UU = new matrix<double>(size,size);
  VV = new matrix<double>(size,size);
  JJ = new double[size];
  hh = new double[size];
  
  for (int i=0;i<size;++i){
    hh[i] = h + 2.0*epsilon*(drand1()-0.5);
    JJ[i] = J + 2.0*epsilon*(drand1()-0.5);
  }
  (*hamiltonian) = get_hamiltonian();
  solve_diagonalization();
#ifdef DEBUG
  _ERROR_TRACKING_;
#endif
}


void ising1D::solve_diagonalization()
{
  double * tempeigenval = new double[2*size];
  double * eigenval = new double[size];

  matrix<double> temph(2*size,2*size);

  temph = get_hamiltonian();

  tempeigenval = temph.diagonalize(true);

  for (int i=0;i<size;++i){
    eigenval[i] = tempeigenval[2*size-i-1];
    
    /* uncomment to print the eigenvalues
      cout << eigenval[i] << endl;
    */
  }

  e = new energy(eigenval, size);

  for (int irow=0;irow<size;++irow)
    for (int icol=0;icol<size;++icol){
      (*VV)(irow,icol) = temph(icol,irow);
      /* attention: the order of eigenvectors is inverted for e>0
	 respect to our notation */
      (*UU)(irow,icol) = temph(icol,irow+size);
    }

  /* uncomment to check the definition of U and V  
  for (int irow=0;irow<size;++irow)
    for (int icol=0;icol<size;++icol){
      temph(irow,icol) = (*UU)(irow,icol) ;
      temph(irow+size,icol+size) = (*UU)(irow,icol);
      temph(irow+size,icol) = (*VV)(irow,icol) ;
      temph(irow,icol+size) = (*VV)(irow,icol);
    }
  (temph.transpose() * *hamiltonian * temph).print();
  */

  delete [] eigenval;
  delete [] tempeigenval;

#ifdef DEBUG
  _ERROR_TRACKING_;
#endif
  return;
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
    AA(irow,irow) = -hh[irow];
    if (irow+1<size){
      AA(irow,irow+1) = -JJ[irow]*0.5;
      AA(irow+1,irow) = -JJ[irow]*0.5;
    }
  }
  
  if (pbc){
    AA(size-1,0) += JJ[size-1]*0.5;
    AA(0,size-1) += JJ[size-1]*0.5;
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
  double temp = gamma/2.0;
  
  BB = 0.0;
  
  for (int irow=0;irow<size-1;++irow){
    BB(irow,irow+1) = -JJ[irow]*temp;
    BB(irow+1,irow) =  JJ[irow]*temp;
    }

  // supposing even NF
  if (pbc){
    BB(size-1,0) += JJ[size-1]*temp;
    BB(0,size-1) += -JJ[size-1]*temp;
  }

#ifdef DEBUG
  _ERROR_TRACKING_;
#endif
  return BB;
}


ising1D::~ising1D()
{
  delete [] JJ;
  delete [] hh;
  delete hamiltonian;
  delete UU;
  delete VV;
  delete e;
}


void ising1D::print()
{
  cout << "size = " << size << endl;
  cout << "h = " << h << endl;
  cout << "J = " << J << endl;
  cout << "gamma = " << gamma << endl;
  cout << "epsilon = " << epsilon << endl;
  cout << "pbc = " << pbc << endl;
}


void ising1D::write(out_file* file)
{
  file->write_tag("system");
  file->write_parameter("size", size );
  file->write_parameter("h", h );
  file->write_parameter("J", J );
  file->write_parameter("epsilon", epsilon );
  file->write_parameter("gamma", gamma );
  file->write_parameter("pbc", pbc );
#ifdef DEBUG
  _ERROR_TRACKING_;
#endif
}


void ising1D::read( in_file* file,const string systemname)
{
  string name, data;

  if (file->find_tag(systemname)<0){
    _ERROR_("not able to find tag");
  }
#ifdef DEBUG
  _ERROR_TRACKING_;
#endif
  while (  file->read_data(name,data)>0 ){
    if (name=="size")
      istringstream(data) >> size;
    else if (name=="h")
      istringstream(data) >> h;
    else if (name=="J")
      istringstream(data) >> J;
    else if (name=="gamma")
      istringstream(data) >> gamma;
    else if (name=="epsilon")
      istringstream(data) >> epsilon;
    else if (name=="pbc")
      istringstream(data) >> pbc;
  }
#ifdef DEBUG
  _ERROR_TRACKING_;
#endif
}


quench::quench( int size_in, double h0_in, double h_in, double J_in, double epsilon_in, double gamma_in, bool pbc_in)
{
  system0 = new ising1D(size_in, h0_in, J_in, epsilon_in, gamma_in, pbc_in);
  system = new ising1D(size_in, h_in,  J_in, epsilon_in, gamma_in, pbc_in);
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


quench::quench( in_file* file)
{
  system0 = new ising1D(file,"system0");
  system = new ising1D(file,"system");
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
