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


class ising1D;

//!  Energy observable.
/*!
  Local observable that contains the eigenstate of the Hamiltonian.
 */
class energy : public local_obs<double>{
public:
  energy( double *, int );
  double get_time_evolution( matrix< complex<double> > *UUt, matrix< complex<double> > *VVt );
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


// general 1D ising model object
class ising1D{
public:
  // length of the chain
  int size;
  // parameters of the hamiltonian:
  // H = - \sum_i J_i \sigma_i^x \sigma_{i+1}^x - \sum_i h_i \sigma^z_i
  // hh[i] = h + epsilon*rand()
  double *hh;
  double *JJ;
  double h,J,gamma;
  double epsilon;
  bool pbc;
  
  ising1D( int size_in, double h_in=_H_, double J_in=_J_, double epsilon_in = _EPSILON_, double gamma_in=_GAMMA_, bool pbc_in=_PBC_,int seed=_SEED_);

  ising1D( in_file *file, const string name="system");
  ~ising1D();

  // matrix that represent the Hamiltonian
  matrix<double> *hamiltonian;
  matrix<double> *UU,*VV;
  energy* e;

  void print();
  void write( out_file *);
  void read( in_file *, const string);

private:
  void init();
  matrix<double>  get_hamiltonian();
  void solve_diagonalization();
  matrix<double>  get_matrix_A();
  matrix<double>  get_matrix_B();
  int _seed;
};


class quench{
public:
  quench( int size_in, double h0_in, double h_in, double J_in=_J_, double epsilon_in = _EPSILON_, double gamma_in=_GAMMA_, bool pbc_in=_PBC_);
  quench( in_file* );
  ~quench();

  matrix< complex<double> > *UUt,*VVt;
  double* gge;
  void set_time_evolution( const double );
  void set_gge_occupations( );
  int get_size();

  ising1D * system0;
  ising1D * system;
protected:
  int size;
  void init();
  matrix< complex<double> > get_evolution_matrix(const double );
};


#endif
