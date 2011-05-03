
#include "common.hh"
#include "error.hh"
#include "io.hh"
#include "random.hh"
#include <iostream>
#include <cmath>
#include <string>
using namespace std;

double fermi_dist(double energy, double temp)
{
  double beta = 1/temp;
  return 1.0/(exp(beta*energy)+1.0);
}


state::state( int size_in )
{
  size = size_in;
  conf = new bool[size];
  clear();
}


state::~state()
{
  delete [] conf;
}

void state::clear()
{
  for (int i=0;i<size;++i)
    conf[i] = false;
}

void state::print()
{
  cout << "state ";
  for (int i=0;i<this->size;++i)
    cout << this->conf[i];
  cout << endl;
}

void state::random()
{
  for (int i=0;i<size;++i)
    conf[i] = int (2.0*drand1());
}

void state::flip( int n)
{
  if ((n>=size) || (n<0))
    _ERROR_("you are flipping a site outside the chain");
  else
    conf[n] = !conf[n];
}

void state::next()
{
  int i=0;
  int rest=0;
  
  do {
    if (i==size){
      _ERROR_("You are producing too big state");
      break;
    }
    if (conf[i])
      rest = 1;
    else
      rest = 0;

    conf[i] = !conf[i];
    ++i;
  } while (rest != 0);
}


bool state::islast()
{ 
  for (int i=0;i<size;++i)
    if (not conf[i]) return false;
  return true;
}


void state::from_int(unsigned int num)
{
  clear();
  int i = 0;
  while (num != 0){
    if (i==size){
      _WARNING_("you are generating too many states");
      break;
    }
    conf[i] = num % 2;
    num = (num - num%2)/2;
    ++i;
  }
}

void operator++( state& s)
{
  s.next();
}

ostream& operator<<(ostream& out, const state &s)
{
  for (int i=0;i<s.size;++i)
    out << s.conf[i];
  return out;
}

template<>
double loop<double>::get_val()
{
  return double(index)*delta+initval;
}

template<>
int loop<int>::get_val()
{
  return int(index)*delta+initval;
}

template<>
float loop<float>::get_val()
{
  return float(index)*delta+initval;
}

template<>
double loop<double>::get_max_val()
{
  return double(steps)*delta+initval;
}

template<>
int loop<int>::get_max_val()
{
  return int(steps)*delta+initval;
}

template<>
float loop<float>::get_max_val()
{
  return float(steps)*delta+initval;
}
