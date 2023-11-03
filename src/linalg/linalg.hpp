#ifndef LINALG_H
#define LINALG_H

#include "io.hpp"


/************ CONSTANT DEFINITIONS ****************/

const int NR_END = 1;
const double TINY = 1.0e-20;
const double EPSILON = 1e-5;



/************ INLINE FUNCTION DEFINITIONS ****************/

template<class T>
inline void SWAP(T &a, T &b) {
  T dum=a; a=b; b=dum;
}


// Vector - vector addition
template <typename T>
inline vector<T> operator+(const vector<T> & a, const vector<T> & b)
{
  assert(a.size() == b.size());
  vector<T> c(a);
  transform(c.begin(), c.end(), b.begin(), c.begin(), plus<T>());

  return c;
}


// Matrix - matrix addition
template <typename T>
inline vector<vector<T> > operator+(const vector<vector<T> > & A, const vector<vector<T> > & B)
{
  assert( (A.size() == B.size()) && (A.at(0).size() == B.at(0).size()) );

  vector<vector<T> > C(A.size(), vector<T>(A[0].size()));
  for (int i=0; i<C.size(); ++i) 
    for (int j=0; j<C[i].size(); ++j) 
      C[i][j] = A[i][j] + B[i][j];

  return C;

}


// Vector - vector substraction
template <typename T>
inline vector<T> operator-(const vector<T> & a, const vector<T> & b)
{
  assert(a.size() == b.size());
  vector<T> c(a);
  transform(c.begin(), c.end(), b.begin(), c.begin(), minus<T>());

  return c;
}


// Matrix - matrix substraction
template <typename T>
inline vector<vector<T> > operator-(const vector<vector<T> > & A, const vector<vector<T> > & B)
{
  assert( (A.size() == B.size()) && (A.at(0).size() == B.at(0).size()) );

  vector<vector<T> > C(A.size(), vector<T>(A[0].size()));
  for (int i=0; i<C.size(); ++i) 
    for (int j=0; j<C[i].size(); ++j) 
      C[i][j] = A[i][j] - B[i][j];

  return C;

}


// Scalar - vector multiplication
template <typename T>
inline vector<T> operator*(const T k, vector<T> vec)
{
  transform( vec.begin(), vec.end(), vec.begin(), bind2nd(multiplies<T>(),k) );

  return vec;
}


// Scalar - matrix multiplication
template <typename T>
inline vector<vector<T> > operator*(const T k, vector<vector<T> > A)
{
  for (int i=0; i<A.size(); ++i)
    for (int j=0; j<A[i].size(); ++j)
      A[i][j] *= k;

  return A;
}


// Vector - vector multiplication : Scalar product
template <typename T>
inline T operator*(const vector<T> & a, const vector<T> & b)
{
  assert(a.size() == b.size());
  return inner_product(a.begin(), a.end(), b.begin(), T(0));
}

// Vector - vector multiplication : Tensor product
template <typename T>
inline vector<vector<T> > tensor_product(const vector<T> & a, const vector<T> & b)
{
  const int m = a.size();
  const int n = b.size();

  vector<vector<T> > A(m, vector<T>(n));

  for (int i=0; i<m; ++i)
    for (int j=0; j<n; ++j)
      A[i][j] = a[i]*b[j];

  return A;
}


// Matrix - vector multiplication
template <typename T>
inline vector<T> operator*(const vector<vector<T> > & A, const vector<T> & b)
{
  assert(A.at(0).size() == b.size());

  vector<T> v(A.size());

  for (int i=0; i < A.size(); ++i)
    for (int k=0; k < A[0].size(); ++k)
      v[i] += A[i][k]*b[k];
  
  return v;
}


// Matrix - matrix multiplication
template <typename T>
inline vector<vector<T> > operator*(const vector<vector<T> > & A, const vector<vector<T> > & B)
{
  assert(A.at(0).size() == B.size());

  vector<vector<T> > R(A.size(), vector<T>(B[0].size()));

  for (int i=0; i<A.size(); ++i)
    for (int j=0; j<B[0].size(); ++j)
      for (int k=0; k<A[0].size(); ++k)
	R[i][j] += A[i][k]*B[k][j];

  return R;
}


inline bool isZero(int t)
{
  return (t == int(0));
}


inline bool isZero(rational<int> t)
{
  return (t == rational<int>(0));
}


inline bool isZero(double t)
{
  static const double EPSILON = 1e-5;

  return (abs(t) < EPSILON);
  
}


inline bool isZero(const complex<double> & t)
{
  return ( isZero(t.real()) and isZero(t.imag()) );
}


inline bool isZero(const complex<rational<int> > & t)
{
  return ( isZero(t.real()) and isZero(t.imag()) );
}


template <typename T>
inline bool isZero(const vector<T> & t)
{
  for (unsigned i=0; i<t.size(); ++i) if (not(isZero(t[i]))) return false;
  
  return true;
}


template <typename T>
inline bool isZero(const vector<vector<T> > & t)
{
  for (unsigned i=0; i<t.size(); ++i) 
    for (unsigned j=0; j<t[i].size(); ++j)
      if (not(isZero(t[i][j]))) return false;
  
  return true;
}


inline bool is_real(const Complex & t)
{
  return (isZero(t.imag()));
}


inline bool is_real(const complexVector & t)
{
  for (unsigned i=0; i<t.size(); ++i) if ( not(is_real(t[i])) ) return false;
  
  return true;
}


template <typename T>
inline bool is_constant(const vector<T> & v)
{
  if (v.size() == 0) 
    return true;
  else
    return isZero(v - vector<T>(v.size(), v[0]));
}


inline bool isEven(rational<int> t)
{
  if ((t.denominator() ==1) && (t.numerator()%2 == 0)) return true;
  else return false;
  
}


inline bool isEven(double t)
{
  int number = static_cast<int>(floor(t + 0.5));
  if (abs(double(number) - t) > EPSILON) return false;

  if (number % 2 == 0) return true;
  else return false;
}


template <typename T>
inline bool isEqual(vector<T> u, vector<T> v)
{
  assert(u.size() == v.size());

  vector<T> w(u.size());
  for (unsigned i=0; i<u.size(); ++i) w[i] = u[i]-v[i];
  return isZero(w);
  
}


template <typename T>
inline bool is_orthogonal(const vector<T> & v, const vector<T> & w)
{
  return isZero(v * w);
}


template <typename T>
inline bool is_orthogonal(const vector<T> & v, const vector<vector<T> > & A)
{
  for (int i=0; i<A.size(); ++i) 
    if ( isZero(v * A[i]) == false) return false;

  return true;
}


inline bool relativelyPrime(vector<int> a)
{
  a.erase( remove_if(a.begin(), a.end(), bind2nd(equal_to<int>(), 0) ), a.end() );

  if (a.size() == 0) return true;
  
  int divisor = a.at(0);

  for (unsigned i=1; i<a.size(); ++i) {
	divisor = gcd(divisor,a[i]);
  }

  if (divisor > 1) return false;
  else return true;

}


inline bool relativelyPrime(vector<rational<int> > a)
{
  vector<int> b(a.size());

  for (unsigned i=0; i<a.size(); ++i) {
    assert(a[i].denominator() == 1);
    b[i] = a[i].numerator(); 
  }

  return relativelyPrime(b);
}


inline bool isInteger(rational<int> t)
{
  return (t.denominator() == 1);
}


inline int roundx(double y)
{

  return (int)floor(y+0.5);

}


inline bool isInteger(double t)
{
  return (abs(roundx(t)-t) < EPSILON);
}


inline bool isInteger(rationalVector t)
{
  for (unsigned i=0; i<t.size(); ++i)
    if (!isInteger(t[i])) return false; 

  return true;
}


inline bool isInteger(doubleVector t)
{
  for (unsigned i=0; i<t.size(); ++i)
    if (!isInteger(t[i])) return false; 

  return true;
}


template <typename T>
inline bool is_half_integer(const T & t)
{
  return ( (isInteger(t) == false) and (isInteger(T(2)*t) == true) );
}


inline int getDenominator(rational<int> t)
{
  return t.denominator();
}


inline int getDenominator(double t)
{
  for (int i=1; i<100; ++i) {
    double x;
    x = double(i)*t;
    if (isInteger(x)) return i;
  }
  throw Error("Error in routine getDenominator: Index out of range ...");
  return 999;
}


inline bool isGreaterZero(int t)
{
  return (t > 0);
}


inline bool isGreaterZero(rational<int> t)
{
  return (t > rational<int>(0));
}


inline bool isGreaterZero(double t)
{
  if ( (t < 0.0) and (abs(t) > EPSILON) ) return false;
  else return true;
}


inline bool isGreaterEqualZero(int t)
{
  return (t >= 0);
}


inline bool isGreaterEqualZero(rational<int> t)
{
  return (t >= rational<int>(0));
}


inline bool isGreaterEqualZero(double t)
{
  if ( (t < 0.0) and (abs(t) > EPSILON) ) return false;
  else return true;
}


template <typename T>
inline bool isGreaterZero(const vector<T> & t)
{
  bool isgreater = true;
  for (unsigned i=0; (i<t.size()) and (isgreater == true); ++i) {
    if ( !isGreaterZero(t[i]) ) isgreater = false;
  }
  
  return isgreater;
}


template <typename T>
inline bool isGreaterEqualZero(const vector<T> & t)
{
  bool isgreater = true;
  for (unsigned i=0; (i<t.size()) and (isgreater == true); ++i) {
    if ( !isGreaterEqualZero(t[i]) ) isgreater = false;
  }
  
  return isgreater;
}


/************ CLASS DEFINITIONS ****************/


class eigenSystem {
 public:
  unsigned n;
  vector<double> eigenvalues;
  vector<vector<double> > eigenvectors;
  vector<vector<double> > inverse_eigenvectors;
  eigenSystem(vector<vector<double> > mat);
};


/************ CRITERIONS AND FUNCTION OBJECTS ****************/


inline bool sp_is_zero_criterion(vector<double> v, vector<double> w) 
{
  return isZero(v*w);
}


inline bool absolute_value_is_equal_crit(vector<double> v, vector<double> w) 
{
  return isZero(v-w) or isZero(v+w);
}


template <typename T>
class sp_with_shift_is_not_integer_crit {
private:
  vector<T> shift;
public:
  sp_with_shift_is_not_integer_crit (const vector<T> & shift1) : shift(shift1) { }
  bool operator() (const vector<T> & v) const
  {
    return not(isInteger(shift*v));
  }

};


template <typename T>
class sp_with_shift_is_not_zero_crit {
private:
  vector<vector<T> > shift;
public:
  sp_with_shift_is_not_zero_crit (const vector<T> & shift1) : shift(1, shift1) { }
  sp_with_shift_is_not_zero_crit (const vector<vector<T> > & shift1) : shift(shift1) { }
  bool operator() (const vector<T> & v) const
  {
    for (int i=0; i<shift.size(); ++i) if (not(isZero(shift[i]*v))) return true;
    
    return false;
  }

};


template <typename T>
class is_zero_crit {
public:
  bool operator() (const T & v) const
  {
    return isZero(v);
  }

};


template <typename T>
class is_not_zero_crit {
public:
  bool operator() (const T & v) const
  {
    return not(isZero(v));
  }

};


template <typename T>
class weight_wrt_semiordering_is_negative_or_zero_crit {
private:
  map< vector<T>, vector<T> > * rootlabels; 
public:
   weight_wrt_semiordering_is_negative_or_zero_crit (map< vector<T>, vector<T> > & labels) { rootlabels = &labels; }
  bool operator() (const vector<T> & v) const
  {
    vector<T> w = (*rootlabels)[v];

    typename vector<T>::iterator iter = find_if(w.begin(), w.end(), is_not_zero_crit<T>());

    if ( (iter != w.end()) and (isGreaterZero(*iter)) ) return false;
    else return true;
    
  }

};


template <typename T>
class smaller_wrt_root_ordering_crit {
private:
  vector<T> X;
public:
  smaller_wrt_root_ordering_crit (const vector<T> & X1) : X(X1) { } ;
  bool operator() (const vector<T> & v, const vector<T> & w) const
  {
    return X*(v - w) < 0;
  }

};


template <typename T>
class differs_by_lattice_shift_crit : public binary_function<vector<T>, vector<T>, bool>
{
private:
  vector<vector<T> > simpleroots;
public:
  differs_by_lattice_shift_crit(vector<vector<T> > simpleroots1) : simpleroots(simpleroots1) {}

  bool operator() (vector<T> v, vector<T> w) const
  {
    vector<T> u = v - w;
    return inLattice(simpleroots, u);
  }
  
};


template <typename T>
class differs_by_Spin32_lattice_vector_crit : public binary_function<vector<T>, vector<T>, bool>
{
public:
  differs_by_Spin32_lattice_vector_crit() {};

  bool operator() (vector<T> v, vector<T> w) const
  {
    return is_in_Spin32_lattice(v-w);
  }
  
};


/************ FUNCTION DECLARATIONS ****************/

void set_double_zeros_to_zero(doubleVector & A);

void set_double_zeros_to_zero(doubleMatrix & A);

void set_double_zeros_to_zero(complexMatrix & A);

void set_double_zeros_to_zero(rationalMatrix & A);

void set_double_zeros_to_zero(complexrationalMatrix & A);

Rational abs(const complex<Rational> & c);

rational<int> to_rational(const string & str);

rationalVector to_rationalVector(const string & str);

rational<int> to_rational(int d);

rational<int> to_rational(double d, double precision);

rational<int> to_rational(const complex<rational<int> > & t);

template <typename T>
rationalVector to_rational(const vector<T> & d);

template <typename T>
rationalMatrix to_rational(const vector<vector<T> > & d);

int to_integer(const string & str);

intVector to_intVector(const string & str);

int to_integer(const double t);

int to_integer(const rational<int> t);

template <typename T>
intVector to_integer(const vector<T> & t);

template <typename T>
intMatrix to_integer(const vector<vector<T> > & t);

double to_double(const double c);

double to_double(const string & str);

doubleVector to_doubleVector(const string & str);

double to_double(const rational<int> A);

double to_double(const int A);

double to_double(const complex<rational<int> > & A);

double to_double(const complex<double> & A);

template <typename T>
doubleVector to_double(const vector<T> & t);

template <typename T>
doubleMatrix to_double(const vector<vector<T> > & t);

template <typename T>
vector<doubleMatrix> to_double(const vector<vector<vector<T> > > & A);

template <typename T>
vector<vector<T> > findBasis(vector<vector<T> > a);

template <typename T>
void gaussj(vector<vector<T> > & a, vector<vector<T> > & b);

template <typename T>
vector<vector<T> > linSolve(vector<vector<T> > a, int q);

template<typename T>
vector<T> linSolve(const vector<vector<T> > & A, const vector<T> & b);

template<typename T>
vector<vector<T> > linSolve(const vector<vector<T> > & A, const vector<vector<T> > & b);

doubleVector linSolve(const intMatrix & A, const doubleVector & B);

rationalVector linSolve(const vector<vector<int> > & A, const rationalVector & B);

template <typename T>
vector<T> linSolve_return(const intMatrix & A, const intVector & B);

template <typename T>
vector<vector<T> > transpose(const vector<vector<T> > & A);

template <typename T>
T Trace(const vector<vector<T> > A);

doubleMatrix inverse(const doubleMatrix & A);

rationalMatrix inverse(rationalMatrix A);

rationalMatrix inverse(const intMatrix & A);

template <typename T>
vector<vector<T> > identity_matrix(unsigned n);

doubleMatrix diagonalize_matrix(const doubleMatrix & matrix, doubleMatrix & change_of_basis);

template <typename T>
vector<vector<T> > gramSchmidt_orthogonal(vector<vector<T> > A);

int factorial(const int n);

void Comb(vector<intVector> & result, int n, int k, vector<int> & a, int j, int m); 

vector<intVector> combination(int n, int k); 

vector<intVector> combination(int n); 

void permutation(int n, vector<int> & P, int & N, int & t, vector<int> & a, vector<intVector> & result); 

vector<intVector> permutation(const vector<int> & a1);

vector<intVector> permutation(int N);

void partition_greatest_part_k(int n, int k, int t, vector<int> & p, vector<intVector> & result); 

vector<intVector> partition_greatest_part_k(int n, int k);

vector<intVector> partition(int n);

template <typename T>
bool gcd_rossers_algorithm_vector_greater_criterion(vector<T> v, vector<T> w); 

pair<int, intMatrix> gcd_rossers_algorithm(const vector<int> & b);

int gcd(const vector<int> & b);

pair<intVector, intMatrix> solution_of_single_diophantine_equation(const vector<int> & b, int c);

intVector particular_solution_of_single_diophantine_equation(const vector<int> & b, int c);

pair<intVector, intMatrix> solve_system_of_linear_diophantine_equations(const vector<vector<int> > & A, const vector<int> & b);

#endif
