#include "io.hpp"


/**
Returns the boolean value 'true' if the input string contains the string 'true'; else returns 'false'.
*/
bool to_bool(const string & str)
{
  if (str.find("true") != string::npos) return true;
  else return false;
}


/**
Converts a string containing a sequence of numbers separated by white space to a vector of type T.
*/
template <typename T>
vector<T> convert_string_to_vector(string vec)
{ 
  T temp;
  vector<T> data;
  
  istringstream line(vec);
  while (line >> temp) data.push_back(temp);
  
  return data;
  
}


/**
Reads rows of numbers from an input stream and returns a matrix.
*/
template <typename T>
vector<vector<T> > read_matrix_from_disk(ifstream & file)
{
  vector<vector<T> > data;
  
  vector<string> linesOfData;
  string tempstring;
  while (getline(file, tempstring) && (tempstring.find("begin data")) == string::npos) {};
  while (getline(file, tempstring) && (tempstring.find("end data")) == string::npos) 
    linesOfData.push_back(tempstring);

  T temp;
  vector<T> currentData;
  vector<string>::iterator currentLine;
  for (currentLine = linesOfData.begin(); currentLine != linesOfData.end(); ++currentLine) {

    istringstream line(*currentLine);
    currentData.clear();
    while (line >> temp) currentData.push_back(temp);
    data.push_back(currentData);
    
  }

  return data;

}


/**
Reads rows of numbers from a file and returns a matrix.
*/
template <typename T>
vector<vector<T> > read_matrix_from_disk(const string & filename)
{
  ifstream in;
  in.open(filename.c_str());
  if (!in) {
    ostringstream os;
    os << "Error in routine \'read_matrix_from_disk\' : Cannot open file " << filename << "...";
    throw Error(os.str());
  }

  return read_matrix_from_disk<T>(in);
}


/**
Reads an 'offset' number of lines from an input stream and returns a vector of strings.
*/
vector<string> return_next_lines(ifstream & file, unsigned offset)
{
  string currentline;
  vector<string> linesOfData;
  
  unsigned counter = 0;
  while ( (++counter <= offset) && getline(file, currentline) ) {
    linesOfData.push_back(currentline);
  }
  
  return linesOfData;
}


/**
Reads all the lines from an input stream, until a regular expression is matched. Returns a vector of strings.
*/
vector<string> return_next_lines(ifstream & file, const string & regexp)
{
  string currentline;
  vector<string> linesOfData;
  
  bool notcomplete = true;
  while ( getline(file, currentline) and ( notcomplete = (currentline.find(regexp) == string::npos) ) ) {
    linesOfData.push_back(currentline);
  }
  
  if (notcomplete == true) linesOfData.clear();

  return linesOfData;
}


/**
Output operator for a vector of integers.
*/
ostream & operator<<(ostream & os, const intVector & vec) 
{
  for (unsigned i=0; i<vec.size(); i++) {
    os.width(4); os.fill(' ');
    os << vec[i];
  }
 
  return os;

}


/**
Output operator for a vector of rational numbers.
*/
ostream & operator<<(ostream & os, const rationalVector & vec) 
{
  for (unsigned i=0; i<vec.size(); i++) {
    os.width(4); os.fill(' ');
    os << vec[i];
  }
 
  return os;

}


/**
Output operator for a vector of real numbers.
*/
ostream & operator<<(ostream & os, const doubleVector & vec) 
{
  os << setprecision(6);

  for (unsigned i=0; i<vec.size(); i++) {
    os << setw(10) << vec[i];
  }
  
  return os;
}


/**
Output operator for a vector of complex numbers.
*/
ostream & operator<<(ostream & os, const complexVector & vec) 
{
  for (unsigned i=0; i<vec.size(); i++) {
    os.width(6); os.fill(' ');
    os << vec[i];
  }
  
  return os;
}


/**
Output operator for a vector of complex numbers with rational real and imaginary parts.
*/
ostream & operator<<(ostream & os, const vector<complex<rational<int> > > & vec) 
{
  for (unsigned i=0; i<vec.size(); i++) {
    os.width(6); os.fill(' ');

    if (vec[i].imag() == rational<int>(0)) os << vec[i].real();
    else if (vec[i] == complex<rational<int> >(rational<int> (0),rational<int> (1))) os << " i";
    else if (vec[i] == complex<rational<int> >(rational<int> (0),rational<int> (-1))) os << " -i";
    else if (vec[i] == complex<rational<int> >(rational<int> (0),rational<int> (0))) os << " 0";
    else os << vec[i];
  }
  
  return os;
}


/**
Output operator for a matrix of integers.
*/
ostream & operator<<(ostream & os, const intMatrix & A)
{
  os.width(4); os.fill(' ');
  
  for (int i=0; i<(int)A.size(); i++) os << A[i] << endl;
  
  return os;
}


/**
Output operator for a matrix of rational numbers.
*/
ostream & operator<<(ostream & os, const rationalMatrix & A)
{
  os.width(4); os.fill(' ');
  
  for (int i=0; i<(int)A.size(); i++) os << A[i] << endl;
  
  return os;
}


/**
Output operator for a matrix of real numbers.
*/
ostream & operator<<(ostream & os, const doubleMatrix & A)
{
  os << setprecision(4);

  for (unsigned i=0; i<A.size(); i++) {
    for (unsigned j=0; j<A[i].size(); j++) {
      os << setw(10) << A[i][j];
    }
    os << endl;
  }
  return os;
}


/**
Output operator for a vector of matrices of real numbers.
*/
ostream & operator<<(ostream & os, const vector<doubleMatrix> & A)
{
  for (unsigned i=0; i<A.size(); i++) {
    os << "Matrix " << i << "...\n" << A[i] << endl << endl;
    }

  return os;
}


/**
Output operator for a matrix of complex numbers.
*/
ostream & operator<<(ostream & os, const complexMatrix & A)
{
  os.width(4); os.fill(' ');
  
  for (int i=0; i<(int)A.size(); i++) os << A[i] << endl;
  
  return os;
}


/**
Output operator for a matrix of complex numbers with rational real and imaginary parts.
*/
ostream & operator<<(ostream & os, const vector<vector<complex<rational<int> > > > & A)
{
  os.width(4); os.fill(' ');
  
  for (int i=0; i<(int)A.size(); i++) os << A[i] << endl;
  
  return os;
}


/**
Output operator for a set of non-negative integers.
*/
ostream & operator<<(ostream & os, const set<unsigned> & A)
{
  set<unsigned>::const_iterator iter;

  for (iter = A.begin(); iter != A.end(); ++iter) os << *iter << " ";
  
  return os;
}


/**
Output operator for a set of integers.
*/
ostream & operator<<(ostream & os, const set<int> & A)
{
  set<int>::const_iterator iter;

  for (iter = A.begin(); iter != A.end(); ++iter) os << *iter << " ";
  
  return os;
}


/**
Output operator for a set of vectors of rational numbers.
*/
ostream & operator<<(ostream & os, const set<rationalVector> & A)
{
  set<rationalVector>::const_iterator iter;

  for (iter = A.begin(); iter != A.end(); ++iter) os << *iter << endl;
  
  return os;
}


/**
Output operator for a set of vectors of integers.
*/
ostream & operator<<(ostream & os, const set<intVector> & A)
{
  set<intVector>::const_iterator iter;

  for (iter = A.begin(); iter != A.end(); ++iter) os << *iter << endl;
  
  return os;
}


/**
Output operator for a set of vectors of real numbers.
*/
ostream & operator<<(ostream & os, const set<doubleVector> & A)
{
  set<doubleVector>::const_iterator iter;

  for (iter = A.begin(); iter != A.end(); ++iter) os << *iter << endl;
  
  return os;
}


/**
Output operator for a multi-set of vectors of real numbers.
*/
ostream & operator<<(ostream & os, const multiset<doubleVector> & A)
{
  multiset<doubleVector>::const_iterator iter;

  for (iter = A.begin(); iter != A.end(); ++iter) os << *iter << endl;
  
  return os;
}


/**
Output operator for a multi-set of vectors of integers.
*/
ostream & operator<<(ostream & os, const multiset<intVector> & A)
{
  multiset<intVector>::const_iterator iter;

  for (iter = A.begin(); iter != A.end(); ++iter) os << *iter << endl;
  
  return os;
}


/**
Output operator for a multi-set of vectors of rational numbers.
*/
ostream & operator<<(ostream & os, const multiset<rationalVector> & A)
{
  multiset<rationalVector>::const_iterator iter;

  for (iter = A.begin(); iter != A.end(); ++iter) os << *iter << endl;
  
  return os;
}


/**
Output operator for a map between integer vectors and integers.
*/
ostream & operator<<(ostream & os, const map<intVector,int>  & A)
{
  map<intVector,int>::const_iterator iter;

  for (iter = A.begin(); iter != A.end(); ++iter) os << "vec: " << iter->first << "\t\tint: " << iter->second << endl;
  
  return os;
}


/**
Output operator for a map between two integer vectors.
*/
ostream & operator<<(ostream & os, const map<intVector,intVector>  & A)
{
  map<intVector,intVector>::const_iterator iter;

  for (iter = A.begin(); iter != A.end(); ++iter) os << "vec1: " << iter->first 
						     << "\t\tvec2: " << iter->second << endl;
  
  return os;
}


/**
Output operator for a map between two double vectors.
*/
ostream & operator<<(ostream & os, const map<doubleVector,doubleVector>  & A)
{
  map<doubleVector,doubleVector>::const_iterator iter;

  for (iter = A.begin(); iter != A.end(); ++iter) os << "vec1: " << iter->first 
						     << "\t\tvec2: " << iter->second << endl;
  
  return os;
}


/**
Output operator for a vector of strings.
*/
ostream & operator<<(ostream & os, const vector<string> & A)
{
  vector<string>::const_iterator iter;

  for (iter = A.begin(); iter != A.end(); ++iter) os << *iter << endl;
  
  return os;
}


/**
Output operator for a vector of non-negative integers.
*/
ostream & operator<<(ostream & os, const vector<unsigned> & A)
{
  vector<unsigned>::const_iterator iter;

  for (iter = A.begin(); iter != A.end(); ++iter){
    os.width(4); os.fill(' ');
    os << *iter;
  }
  
  return os;
}


/**
Output operator for a set of strings.
*/
ostream & operator<<(ostream & os, const set<string> & A)
{
  set<string>::const_iterator iter;

  for (iter = A.begin(); iter != A.end(); ++iter){
    os << *iter << endl;
  }
  
  return os;
}


/**
Output operator for a matrix of non-negative integers.
*/
ostream & operator<<(ostream & os, const vector<vector<unsigned> > & A)
{
  vector<vector<unsigned> >::const_iterator iter;

  for (iter = A.begin(); iter != A.end(); ++iter) os << *iter << endl;
  
  return os;
}


/**
Output operator for a set of type 'T'.
*/
template <typename T>
ostream & operator<<(ostream & os, const set<T> & A)
{
  typename set<T>::const_iterator iter;

  for (iter = A.begin(); iter != A.end(); ++iter){
    os.width(4); os.fill(' ');
    os << *iter;
  }
  
  return os;
}


/*********** S T A R T   E X P L I C I T   I N S T A N T I A T I O N S ***********/

template
vector<double> convert_string_to_vector(string vec);

template
vector<rational<int> > convert_string_to_vector(string vec);

template
vector<vector<rational<int> > > read_matrix_from_disk(ifstream & file);

template
vector<vector<double> > read_matrix_from_disk(ifstream & file);

template
vector<vector<double> > read_matrix_from_disk(const string & filename);

template
ostream & operator<<(ostream & os, const set<double> & A);


/*********** E N D   E X P L I C I T   I N S T A N T I A T I O N S ***********/
