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
class energy : public local_obs<double>{
public:
  energy( double *, int );
  double get_time_evolution( matrix< complex<double> > *UUt, matrix< complex<double> > *VVt );
private:
};


class tmag : public local_obs<double>{
public:
  tmag( ising1D *);
  double get_time_evolution( matrix< complex<double> > *UUt, matrix< complex<double> > *VVt );
  static double _get_time_evolution( matrix< complex<double> > *UUt, matrix< complex<double> > *VVt );
  void set_spvs();
private:
  ising1D * _system;
};


class localtmag : public local_obs<double>{
public:
  localtmag( int, ising1D * );
  double get_time_evolution( matrix< complex<double> > *UUt, matrix< complex<double> > *VVt );
  static double _get_time_evolution( int,matrix< complex<double> > *UUt, matrix< complex<double> > *VVt );
  void set_spvs();
private:
  ising1D * _system;
  int _site;
};


class cidagacj : public local_obs< complex<double> >{
public:
  cidagacj( int i, int j, ising1D * );
  complex<double> get_time_evolution( matrix< complex<double> > *UUt, matrix< complex<double> > *VVt );
  static complex<double> _get_time_evolution( int i,int j,matrix< complex<double> > *UUt, matrix< complex<double> > *VVt );
  void set_spvs();
private:
  ising1D * _system;
  int _i;
  int _j;
};


class AiAj : public local_obs< complex<double> >{
public:
  AiAj(int i, int j, ising1D * );
  complex<double> get_time_evolution( matrix< complex<double> > *UUt, matrix< complex<double> > *VVt );
  static complex<double> _get_time_evolution(int i, int j, matrix< complex<double> > *UUt, matrix< complex<double> > *VVt );
  void set_spvs();
  void set_gsv();
private:
  ising1D * _system;
  int _i,_j;
};



class BiBj : public local_obs< complex<double> >{
public:
  BiBj(int i, int j, ising1D * );
  complex<double> get_time_evolution( matrix< complex<double> > *UUt, matrix< complex<double> > *VVt );
  static complex<double> _get_time_evolution(int i, int j, matrix< complex<double> > *UUt, matrix< complex<double> > *VVt );
  void set_spvs();
  void set_gsv();
private:
  ising1D * _system;
  int _i,_j;
  };



class BiAj : public local_obs< double >{
public:
  BiAj(int i, int j, ising1D * );
  double get_time_evolution( matrix< complex<double> > *UUt, matrix< complex<double> > *VVt );
  static double _get_time_evolution(int i, int j, matrix< complex<double> > *UUt, matrix< complex<double> > *VVt );
  void set_spvs();
private:
  ising1D * _system;
  int _i,_j;
};


class rho : public obs<double>{
public:
  rho(int i, int r, ising1D *);
  ~rho();

  void set_ensemble_average(double* nk);
  double get_ensemble_average(int l=-1);

  void set_time_evolution( matrix< complex<double> > *UUt, matrix< complex<double> > * VVt );
  double get_time_evolution( matrix< complex<double> > *UUt, matrix< complex<double> > * VVt );
  double get_time_evolution(int l=-1);
private:
  ising1D *_system;
  //  quench  * _quench;
  int _r,_i;
  matrix< BiAj* > *_BiAjmatrix;
  void set_BiAjpointers();
  matrix<double> *_reduced_matrix;
  matrix< complex<double> > *_full_matrix;
};

double randomHH(int i, ising1D* system);
double randomJJ(int i, ising1D* system);

// general 1D ising model object
class ising1D{
public:
  
  ising1D( int size_in, double h_in=_H_, double J_in=_J_, double epsilon_in = _EPSILON_, double gamma_in=_GAMMA_, bool pbc_in=_PBC_,int seed=_SEED_);

  ising1D( in_file *file, const string name="system", double (*JJgen)( int,ising1D* )=&randomJJ, double (*HHgen)( int,ising1D* )=&randomHH);
  ~ising1D();

  matrix<double> *UU,*VV;

  void print();
  void write( out_file *);
  void read( in_file *, const string);
  int get_size();
  double get_epsilon();
  double get_h();
  double get_J();
  energy* e;

  friend class quench;

private:
  void init( double (*JJgen)( int,ising1D* ), double (*HHgen)( int,ising1D* ));
  matrix<double>  get_hamiltonian();
  void solve_diagonalization();
  void check( double* eigenval, matrix<double> * eigvect );
  matrix<double>  get_matrix_A();
  matrix<double>  get_matrix_B();
  int _seed;
  bool _pbc;
  double _h,_J,_gamma;
  double _epsilon;
  double *_hh;
  double *_JJ;
  // matrix that represent the Hamiltonian
  matrix<double> *_hamiltonian;
  // length of the chain
  int size;
  bool _from_center; /**< generate disorder from the center of the chain or not */
  bool _ext;
};


class quench{
public:
  quench( int size_in, double h0_in, double h_in, double J_in=_J_, double epsilon_in = _EPSILON_, double gamma_in=_GAMMA_, bool pbc_in=_PBC_);
  quench( in_file*,double (*JJgen0)( int,ising1D* )=&randomJJ, double (*HHgen0)( int,ising1D* )=&randomHH,double (*JJgen)( int,ising1D* )=&randomJJ, double (*HHgen)( int,ising1D* )=&randomHH );
  ~quench();

  matrix< complex<double> > *UUt,*VVt;
  double* gge;
  void set_time_evolution( const double );
  void set_time_evolution( const double , matrix< complex<double> > *, matrix< complex<double> > * );
  void set_gge_occupations( );
  int get_size();
  double get_calpha2( state* );

  ising1D * system0;
  ising1D * system;
private:
  int size;
  void init();
  matrix< complex<double> > get_evolution_matrix(const double );
};

#endif
