#ifndef GT_IO_H
#define GT_IO_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <set>
#include <map>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cassert>

// BOOST
#include <boost/rational.hpp>

using namespace std;
using namespace boost;

/************ MACRO DEFINITIONS ****************/

#define PRINT(A) cout << endl << #A << " : " << (A) << endl;

#define PRINTLN(A) cout << endl << #A << " ...\n" << (A) << endl;

/************ TYPE DEFINITIONS ****************/

typedef rational<int> Rational;

typedef vector<int> intVector;
typedef vector<intVector> intMatrix;

typedef vector<Rational> rationalVector;
typedef vector<rationalVector> rationalMatrix;

typedef vector<double> doubleVector;
typedef vector<doubleVector> doubleMatrix;

typedef complex<double> Complex;
typedef complex<double> complexdouble;
typedef vector<complexdouble> complexVector;
typedef vector<complexVector> complexMatrix;

typedef complex<rational<int> > complexrational;
typedef vector<complexrational> complexrationalVector;
typedef vector<complexrationalVector> complexrationalMatrix;


/************ ERROR HANDLING ****************/

class Error {
public:
  string text;
  Error(string text1) : text(text1) { };
};


/************ FUNCTION DECLARATIONS ****************/

bool to_bool(const string & str);

template <typename T>
vector<T> convert_string_to_vector(string vec);

template <typename T>
vector<vector<T> > read_matrix_from_disk(ifstream & file);

template <typename T>
vector<vector<T> > read_matrix_from_disk(const string & filename);

vector<string> return_next_lines(ifstream & file, unsigned offset);

vector<string> return_next_lines(ifstream & file, const string & regexp);

ostream & operator<<(ostream & os, const intVector & vec);

ostream & operator<<(ostream & os, const rationalVector & vec);

ostream & operator<<(ostream & os, const doubleVector & vec);

ostream & operator<<(ostream & os, const complexVector & vec); 

ostream & operator<<(ostream & os, const vector<complex<rational<int> > > & vec); 

ostream & operator<<(ostream & os, const intMatrix & A);

ostream & operator<<(ostream & os, const rationalMatrix & A);

ostream & operator<<(ostream & os, const doubleMatrix & A);

ostream & operator<<(ostream & os, const vector<doubleMatrix> & A);

ostream & operator<<(ostream & os, const complexMatrix & A);

ostream & operator<<(ostream & os, const vector<vector<complex<rational<int> > > > & A);

ostream & operator<<(ostream & os, const set<unsigned> & A);

ostream & operator<<(ostream & os, const set<int> & A);

ostream & operator<<(ostream & os, const set<rationalVector> & A);

ostream & operator<<(ostream & os, const set<intVector> & A);

ostream & operator<<(ostream & os, const set<doubleVector> & A);

ostream & operator<<(ostream & os, const multiset<doubleVector> & A);

ostream & operator<<(ostream & os, const multiset<intVector> & A);

ostream & operator<<(ostream & os, const multiset<rationalVector> & A);

ostream & operator<<(ostream & os, const map<intVector,int>  & A);

ostream & operator<<(ostream & os, const map<intVector,intVector>  & A);

ostream & operator<<(ostream & os, const map<doubleVector,doubleVector>  & A);

ostream & operator<<(ostream & os, const vector<string> & A);

ostream & operator<<(ostream & os, const vector<unsigned> & A);

ostream & operator<<(ostream & os, const set<string> & A);

ostream & operator<<(ostream & os, const vector<vector<unsigned> > & A);

template <typename T>
ostream & operator<<(ostream & os, const set<T> & A);

#endif
