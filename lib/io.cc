
#include "io.hh"
#include "error.hh"
#include <iostream>
#include <string>
#include <fstream>
using namespace std;

in_file::in_file( const string filename )
{
  file = new ifstream(filename.c_str());
}

in_file::~in_file()
{
  file->close();
}



int in_file::find_tag(const string tag)
{
  string line;
  size_t found;

  file->seekg (0, ios::beg);

  while ( file->good() ){
    getline (*file,line);
    found = line.find("["+tag+"]");
    if (int(found)>=0){
      return int(found);
    }
  }
  file->clear();
  return -1;
}


int in_file::read_data( string &name, string &val )
{
  size_t found;
  string line;

  while ( file->good() ){
    getline (*file,line);
    found = line.find("[");
    if (int(found)>=0){
      return -1;
    }
    found = line.find(TAB);
    if (int(found)<0)
      continue;
    name = line.substr(0,found);
    name = name.substr(name.find_first_not_of(' '),name.find_last_not_of(' ')+1);
    val = line.substr(found+1,line.length()-found);
    val = val.substr(val.find_first_not_of(' '),val.find_last_not_of(' ')+1);
    return 1;
  }

  file->clear();
  return -1;
}


out_file::out_file( const string filename )
{
  file = new ofstream(filename.c_str());
}

out_file::~out_file()
{
  file->close();
}



void out_file::write_tag( const string tag )
{
  string temp;
  temp = "["+tag+"]";
  *file << temp << endl;
}


void out_file::copyfile(in_file * ifile, const string comment)
{
  string line;
  ifile->file->seekg (0, ios::beg);
  while ( ifile->file->good() ){
    getline (*ifile->file,line);
    *file << comment << line << endl;
  }
  ifile->file->clear();
}

out_file& operator<<(out_file& obj, ostream& (*pf)(ostream&))
{
  *obj.file << pf;
  return obj;
}
