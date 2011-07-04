//! ising1D objects definitions
/*!
 * ising1D.hh
 * Purpose: definitons of objects that describes the general ising model
 * with transverse field (general means without imposing order) and its
 * observables
 * @author Simone Ziraldo
 * @version 0.1 19/04/2011
*/
#ifndef __ISING1D_hh__
#define __ISING1D_hh__


#include "matrices.hh"
#include "common.hh"
#include "io.hh"
#include <complex>
#include <string>


#define _SIZE_ 0 /**< int value. Default size. */
#define _H_ 0.0  /**< double value. Default transverse magnetization amplitude. */
#define _J_ 1.0  /**< double value. Default coupling amplitude. */
#define _GAMMA_ 1.0 /**< double value. Default x-y coupling. */
#define _EPSILON_ 0.0 /**< double value. Default disorder. */
#define _PBC_ true /**< boolean value. Default periodic boundary conditions. */
#define _SEED_ -1 /**< default seed value. */
#define _FROM_CENTER_ false /**< boolean value. Starting point of generating disorder */

class ising1D;

//!  Energy observable.
/*!
  Local observable that contains the eigenstate of the Hamiltonian.
 */
class energy : public local_obs<FPType>{
public:
  energy( FPType *, int );
  FPType get_time_evolution( matrix< complex<FPType> > *UUt, matrix< complex<FPType> > *VVt );
};


class tmag : public local_obs<FPType>{
public:
  tmag( ising1D *);
  FPType get_time_evolution( matrix< complex<FPType> > *UUt, matrix< complex<FPType> > *VVt );
  static FPType _get_time_evolution( matrix< complex<FPType> > *UUt, matrix< complex<FPType> > *VVt );
  void set_spvs();
private:
  ising1D * _system;
};


class localtmag : public local_obs<FPType>{
public:
  localtmag( int, ising1D * );
  FPType get_time_evolution( matrix< complex<FPType> > *UUt, matrix< complex<FPType> > *VVt );
  static FPType _get_time_evolution( int,matrix< complex<FPType> > *UUt, matrix< complex<FPType> > *VVt );
  void set_spvs();
private:
  ising1D * _system;
  int _site;
};


class cidagacj : public local_obs< complex<FPType> >{
public:
  cidagacj( int i, int j, ising1D * );
  complex<FPType> get_time_evolution( matrix< complex<FPType> > *UUt, matrix< complex<FPType> > *VVt );
  static complex<FPType> _get_time_evolution( int i,int j,matrix< complex<FPType> > *UUt, matrix< complex<FPType> > *VVt );
  void set_spvs();
private:
  ising1D * _system;
  int _i;
  int _j;
};


class AiAj : public local_obs< complex<FPType> >{
public:
  AiAj(int i, int j, ising1D * );
  complex<FPType> get_time_evolution( matrix< complex<FPType> > *UUt, matrix< complex<FPType> > *VVt );
  static complex<FPType> _get_time_evolution(int i, int j, matrix< complex<FPType> > *UUt, matrix< complex<FPType> > *VVt );
  void set_spvs();
  void set_gsv();
private:
  ising1D * _system;
  int _i,_j;
};



class BiBj : public local_obs< complex<FPType> >{
public:
  BiBj(int i, int j, ising1D * );
  complex<FPType> get_time_evolution( matrix< complex<FPType> > *UUt, matrix< complex<FPType> > *VVt );
  static complex<FPType> _get_time_evolution(int i, int j, matrix< complex<FPType> > *UUt, matrix< complex<FPType> > *VVt );
  void set_spvs();
  void set_gsv();
private:
  ising1D * _system;
  int _i,_j;
  };



class BiAj : public local_obs< FPType >{
public:
  BiAj(int i, int j, ising1D * );
  FPType get_time_evolution( matrix< complex<FPType> > *UUt, matrix< complex<FPType> > *VVt );
  static FPType _get_time_evolution(int i, int j, matrix< complex<FPType> > *UUt, matrix< complex<FPType> > *VVt );
  void set_spvs();
private:
  ising1D * _system;
  int _i,_j;
};



class rho : public obs<FPType>{
public:
  rho(int i, int r, ising1D *);
  ~rho();
  rho( const rho& source ); /*< copy constructor */

  void set_ensemble_average(FPType* nk);
  FPType get_ensemble_average(int l=-1);

  void set_time_evolution( matrix< complex<FPType> > *UUt, matrix< complex<FPType> > * VVt );
  FPType get_time_evolution( matrix< complex<FPType> > *UUt, matrix< complex<FPType> > * VVt );
  FPType get_time_evolution(int l=-1);
private:
  ising1D *_system;
  //  quench  * _quench;
  int _r,_i;
  matrix< BiAj* > *_BiAjmatrix;
  void set_BiAjpointers();
  matrix<FPType> *_reduced_matrix;
  matrix< complex<FPType> > *_full_matrix;
};


class rhozz : public obs<FPType>{
public:
  rhozz(int i, int r, ising1D *);
  ~rhozz();

  void set_ensemble_average();
  FPType get_ensemble_average(FPType* nk);
  FPType get_time_evolution( matrix< complex<FPType> > *UUt, matrix< complex<FPType> > * VVt );
private:
  ising1D *_system;
  int _r,_i;
  BiAj *_BlAl,*_BmAm,*_BmAl,*_BlAm;
  AiAj *_AlAm;
  BiBj *_BlBm;
};





FPType randomHH(int i, ising1D* system);
FPType randomJJ(int i, ising1D* system);

// general 1D ising model object
class ising1D{
public:  
  ising1D( int size_in, FPType h_in=_H_, FPType J_in=_J_, FPType epsilon_in = _EPSILON_, FPType gamma_in=_GAMMA_, bool pbc_in=_PBC_,int seed=_SEED_);

  ising1D( in_file *file, const string name="system", FPType (*JJgen)( int,ising1D* )=&randomJJ, FPType (*HHgen)( int,ising1D* )=&randomHH);
  ~ising1D();

  matrix<FPType> *UU,*VV;

  void print();
  void write( out_file *);
  void read( in_file *, const string);
  int get_size();
  FPType get_epsilon();
  FPType get_h();
  FPType get_J();
  energy* e;
  void gendis( FPType (*JJgen)( int,ising1D* )=&randomJJ, FPType (*HHgen)( int,ising1D* )=&randomHH);

  matrix< complex<FPType> > getEvolutionMatrix(const FPType );

  matrix< complex<FPType> > temp( const FPType time );

  friend class quench;

private:
  void init();
  matrix<FPType>  get_hamiltonian();
  void solve_diagonalization();
  void check( FPType* eigenval, matrix<FPType> * eigvect );
  matrix<FPType>  get_matrix_A();
  matrix<FPType>  get_matrix_B();
  int _seed;
  bool _pbc;
  FPType _h,_J,_gamma;
  FPType _epsilon;
  FPType *_hh;
  FPType *_JJ;
  // matrix that represent the Hamiltonian
  matrix<FPType> *_hamiltonian;
  // length of the chain
  int size;
  bool _from_center; /**< generate disorder from the center of the chain or not */
  bool _ext;
};


class quench{
public:
  quench( int size_in, FPType h0_in, FPType h_in, FPType J_in=_J_, FPType epsilon_in = _EPSILON_, FPType gamma_in=_GAMMA_, bool pbc_in=_PBC_);
  quench( in_file*,FPType (*JJgen0)( int,ising1D* )=&randomJJ, FPType (*HHgen0)( int,ising1D* )=&randomHH,FPType (*JJgen)( int,ising1D* )=&randomJJ, FPType (*HHgen)( int,ising1D* )=&randomHH );
  ~quench();

  matrix< complex<FPType> > *UUt,*VVt;
  FPType* gge;
  void set_time_evolution( const FPType );
  void set_time_evolution( const FPType , matrix< complex<FPType> > *, matrix< complex<FPType> > * );
  void set_gge_occupations( );
  int get_size();
  FPType get_calpha2( state* );

  ising1D * system0;
  ising1D * system;
private:
  int size;
  void init();
  //matrix< complex<FPType> > get_evolution_matrix(const FPType );
};

#endif
