#ifndef WANNIER_PARSER
#define WANNIER_PARSER

#include <string>
#include <map>
#include <array>
#include <sstream>
#include <vector>
#include <tuple>
#include <iostream>
#include <fstream>
#include <limits>
#include <iomanip>
#include <cassert>

using namespace std;


inline int safe_stoi(const std::string& s ){
    int output;
    try{ output = stoi(s) ;  }
    catch(std::exception const & e)
    {
        cerr<<"error in conversion performed by: " << e.what() <<" in function read_wannier_file"<<endl;
        exit(-1);
    }
    return output;
};

vector< vector<string> >read_disorder_file(const string disorder_filename);

tuple<int, vector<string> > read_wannier_file(const string wannier_filename);

vector< tuple<string, array<double, 3> > > read_xyz_file(const string xyz_filename);

array< array<double,3> , 3 >  read_unit_cell_file(const string uc_filename);

#endif
