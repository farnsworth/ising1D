
#ifndef __IO_hh__

#define __IO_hh__

#include "error.hh"
#include <fstream>
#include <string>
using namespace std;

#define SEP " "

class out_file;

/* input file */
class in_file{
public:
  in_file( const string filename );
  ~in_file();
  int find_tag( const string );
  int read_data(string& , string&);
  friend class out_file;
private:
  ifstream * file;
};


/* output file */
class out_file{
public:
  out_file( const string filename );
  ~out_file();
  void write_tag( const string );
  template <class T> void write_parameter(const string name, const T val );
  void copyfile( in_file* ifile, const string = "#" );
  template <class T> friend out_file& operator<<(out_file& ,const T &data);
  // necessary to work with endl, that is a function
  friend out_file& operator<<(out_file& , ostream& (*pf)(ostream&));
private:
  ofstream * file;
};

out_file& operator<<(out_file& , ostream& (*pf)(ostream&));

template <class T>
void out_file::write_parameter( const string name, const T val )
{
  *file << name << SEP << val << endl;
}

template<class T>
out_file& operator<<(out_file& obj, const T & data)
{
  *obj.file << data;
  return obj;
}

#endif
