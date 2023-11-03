#include "linalg.hpp"
#include "io.hpp"

// GSL
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>


/**
Sets each entry equal to 0.0, if its absolute value is smaller than 'EPSILON'.
*/
void set_double_zeros_to_zero(doubleVector & A)
{
  for (unsigned i=0; i<A.size(); ++i) {
      if (abs(A[i]) < EPSILON) A[i] = 0.0;
  }

}


/**
Sets each entry equal to 0.0, if its absolute value is smaller than 'EPSILON'.
*/
void set_double_zeros_to_zero(doubleMatrix & A)
{
  for (unsigned i=0; i<A.size(); ++i) {
    for (unsigned j=0; j<A[i].size(); ++j) {
      if (abs(A[i][j]) < EPSILON) A[i][j] = 0.0;
    }
  }

}


/**
Sets each entry equal to 0.0, if its absolute value is smaller than 'EPSILON'.
*/
void set_double_zeros_to_zero(complexMatrix & A)
{
  for (unsigned i=0; i<A.size(); ++i) {
    for (unsigned j=0; j<A[i].size(); ++j) {
      complex<double> cnum1, cnum2;
      if (abs(A[i][j].real()) < EPSILON) cnum1 = complex<double>(double(0), double(0));
      else cnum1 = complex<double>(A[i][j].real(), double(0));
      if (abs(A[i][j].imag()) < EPSILON) cnum2 = complex<double>(double(0), double(0));
      else cnum2 = complex<double>(double(0), A[i][j].imag());
      A[i][j] = cnum1 + cnum2;
    }
  }

}


/**
Dummy function, needed for consistency. No effect.
*/
void set_double_zeros_to_zero(rationalMatrix & A) { }


/**
Dummy function, needed for consistency. No effect.
*/
void set_double_zeros_to_zero(complexrationalMatrix & A) { }


/*
Returns the modulus of a complex number. Note that this *awkward* implementation serves the only purpose that the function 'linsolve' below works for complex numbers.
*/
Rational abs(const complex<Rational> & c)
{
  double c1 = to_double(c.real());
  double c2 = to_double(c.imag());

  return ( to_rational( sqrt(c1*c1 + c2*c2) ) );
}


/**
Converts a string to a rational number.
*/
rational<int> to_rational(const string & str)
{
  istringstream os(str);
  Rational d;
  
  if (os >> d) return d;
  else throw Error("Could not read rational value ...");
}


/**
Converts a string to a rational vector.
*/
rationalVector to_rationalVector(const string & str)
{
  istringstream os(str);
  Rational d;
  rationalVector D;

  while (os >> d) D.push_back(d);

  return D;
}


/**
Converts an integer to a rational number.
*/
rational<int> to_rational(int d)
{
  return rational<int>(d);
}


/**
Converts a double to a rational number such that the difference between the double and the rational number is less than 'precision'.
*/
rational<int> to_rational(double d, double precision)
{
  rational<int> sg;
 
  if (d < 0.0) {
    sg = rational<int>(-1);
    d *= -1.0;
  } else {
    sg = rational<int>(1);
  }

  vector<int> num;
  
  double x,y;

  x = d;
  y = fmod(x, 1.0);
  x -= y;
  num.push_back(static_cast<int>(x));
  
  while (y > precision) {
    x = 1.0/y;
    y = fmod(x, 1.0);
    x -= y;
    num.push_back(static_cast<int>(x));
  }

  rational<int> sum(0);
  for (int i=num.size()-1; i>=1; --i) {
    sum += num[i];
    sum = 1/sum;
  }
  sum += num[0];

  return sg*sum;

}


/**
Converts a complex number with vanishing imaginary part to a rational number. Throws an exception if the imaginary part is not equal to zero.
*/
rational<int> to_rational(const complex<rational<int> > & t)
{

  if (t.imag() != rational<int>(0))
    throw Error("Error : Number has imaginary part ...");
  else 
    return t.real();

}


/**
Converts a vector with entries of type 'T' to a vector of rational numbers.
*/
template <typename T>
rationalVector to_rational(const vector<T> & d)
{
  rationalVector r(d.size());

  for (int i=0; i<r.size(); ++i) r[i] = to_rational(d[i]);

  return r;
}


/**
Converts a matrix with entries of type 'T' to a matrix of rational numbers.
*/
template <typename T>
rationalMatrix to_rational(const vector<vector<T> > & d)
{
  rationalMatrix r(d.size(), rationalVector(d[0].size()));

  for (int i=0; i<r.size(); ++i)
    for (int j=0; j<r[i].size(); ++j) 
      r[i][j] = to_rational(d[i][j]);

  return r;
}


/**
Converts a string to an integer.
*/
int to_integer(const string & str)
{
  istringstream os(str);
  int d;

  if (os >> d) return d;
  else throw Error("Could not read int value ...");
}


/**
Converts a string to a vector of integers.
*/
intVector to_intVector(const string & str)
{
  istringstream os(str);
  int d;
  intVector D;

  while (os >> d) D.push_back(d);

  return D;
}


/**
Converts a double to an integer.
*/
int to_integer(const double t)
{
  int number = static_cast<int>(floor(t + 0.5));

  if (abs(double(number) - t) > EPSILON)
    throw Error("to_integer: Failed ...");

  return number;

}


/**
Converts a rational number to an integer. Throws an exception if the denominator is not equal to 1.
*/
int to_integer(const rational<int> t)
{
  if (t.denominator() != 1)
    throw Error("convert_rational_to_integer: Failed ...");

  return t.numerator();

}


/**
Converts a vector with entries of type 'T' to a vector of integers.
*/
template <typename T>
intVector to_integer(const vector<T> & t)
{
  vector<int> number(t.size());

  for (unsigned i=0; i<t.size(); ++i) {
    number[i] = to_integer(t[i]);
  }

  return number;

}


/**
Converts a string to a vector of integers.
*/
template <typename T>
intMatrix to_integer(const vector<vector<T> > & t)
{
  vector<vector<int> > number(t.size());

  for (unsigned i=0; i<t.size(); ++i) {
    number[i] = to_integer(t[i]);
  }

  return number;

}


/**
Dummy function. Needed for consistency.
*/
double to_double(const double c)
{
  return c;
}


/**
Converts a string to a double.
*/
double to_double(const string & str)
{
  istringstream os(str);
  double d;

  if (os >> d) return d;
  else throw Error("Could not read double value ...");
}


/**
Converts a string to a vector of doubles.
*/
doubleVector to_doubleVector(const string & str)
{
  istringstream os(str);
  double d;
  doubleVector D;

  while (os >> d) D.push_back(d);

  return D;
}


/**
Converts a rational number to a double.
*/
double to_double(const rational<int> A)
{
  return boost::rational_cast<double>(A);
}


/**
Converts an integer to a double.
*/
double to_double(const int A)
{
  return double(A);
}


/**
Converts a rational number with complex numerator and denominator to a double. Throws an exception if the imaginary part is not equal to zero.
*/
double to_double(const complex<rational<int> > & A)
{
  if (A.imag() == rational<int>(0))
    return boost::rational_cast<double>(A.real());
  else 
    throw Error("Error : Number has imaginary part ...");
  
}


/**
Converts a complex number to a double. Throws an exception if the imaginary part is not equal to zero.
*/
double to_double(const complex<double> & A)
{
  if (isZero(A.imag()) == true)
    return A.real();
  else 
    throw Error("Error : Number has imaginary part ...");
  
}


/**
Converts a vector with entries of type 'T' to a vector of doubles.
*/
template <typename T>
doubleVector to_double(const vector<T> & t)
{
  doubleVector number(t.size());

  for (unsigned i=0; i<t.size(); ++i) {
    number[i] = to_double(t[i]);
  }

  return number;

}


/**
Converts a matrix with entries of type 'T' to a matrix of doubles.
*/
template <typename T>
doubleMatrix to_double(const vector<vector<T> > & t)
{
  doubleMatrix number(t.size());

  for (unsigned i=0; i<t.size(); ++i) {
    number[i] = to_double(t[i]);
  }

  return number;

}


/**
Converts a vector of matrices with entries of type 'T' to a vector of matrices of doubles.
*/
template <typename T>
vector<doubleMatrix> to_double(const vector<vector<vector<T> > > & A)
{
  vector<doubleMatrix> result(A.size());

  for (unsigned i=0; i<A.size(); ++i) {
    result[i] = to_double(A[i]);
  }
    
  return result;
}


/**
Given a list of vectors, returns a basis for the vector space they span.
*/
template <typename T>
vector<vector<T> > findBasis(vector<vector<T> > a)
{
  vector<vector<T> > copy_of_a = a;

  unsigned m = a.size();
  unsigned n = a[0].size();

  vector<vector<T> > basis;
  intVector indx(m);

  unsigned i,j,k,h,maxrow = 0, maxcol = 0, itemp;
  T tmp;

  for (unsigned i=1; i<=m; i++) indx[i-1]=i-1;
  
  for (j=1; j<=n; j++) {

    /* Find the row with the largest first value */
    maxrow = j-1; // leads to error 
    maxcol=j-1;

    if (maxrow >= m) {
      for (unsigned i=0; i<m; ++i) basis.push_back(copy_of_a[indx[i]]);
      return basis;
    }
    
    for (i=j;i<=m;i++) // vormals for (i=j+1;i<=m;i++), geaendert am Feb 15, 2005
      for (h=j; h<=n; h++) {
	if (abs(a[i-1][h-1]) > abs(a[maxrow][maxcol])) {
	  maxrow = i-1;
	  maxcol = h-1;
	}
      }

    /* Swap the maxrow and jth row */
    for (k=j;k<=n;k++) {
      tmp=a[j-1][k-1];
      a[j-1][k-1]=a[maxrow][k-1];
      a[maxrow][k-1]=tmp;
    }
    itemp=indx[j-1]; indx[j-1]=indx[maxrow]; indx[maxrow]=itemp;

     /* Swap the maxcol and jth col */
     for (k=1;k<=m;k++) {
       tmp = a[k-1][j-1];
       a[k-1][j-1] = a[k-1][maxcol];
       a[k-1][maxcol] = tmp;
     }

    /* Singular matrix? */
     if (isZero(a[j-1][j-1])) {
       for (unsigned i=0; i<j-1; ++i) basis.push_back(copy_of_a[indx[i]]);
       return basis;
     }
     
    /* Eliminate the ith element of the jth row */

    for (i=j+1;i<=m;i++) {
      for (k=n;k>=j;k--) {
	a[i-1][k-1] = a[i-1][k-1] - a[i-1][j-1] * a[j-1][k-1] / a[j-1][j-1];
      }
    }

  }

  for (unsigned i=0; i<n; ++i) basis.push_back(copy_of_a[indx[i]]);
  return basis;
   
}


/**
Gauss-Jordan elimination. Throws an exception if the matrix is singular.
*/
template <typename T>
void gaussj(vector<vector<T> > & a, vector<vector<T> > & b)
{
  int i, icol, irow, j, k, l, ll;
  T big, dum, pivinv;

  const int n = a.size();
  const int m = b.at(0).size();
  intVector indxc(n), indxr(n), ipiv(n);
  for (j=0; j<n; ++j) ipiv[j]=0;
  for (i=0; i<n; ++i) {
    big=T(0);
    for (j=0;j<n; ++j)
      if (ipiv[j] != 1)
	for (k=0; k<n; ++k) {
	  if (ipiv[k] == 0) {
	    if ( abs(a[j][k]) >= abs(big) ) {
	      big=abs(a[j][k]);
	      irow=j;
	      icol=k;
	    }
	  }
	}
    ++(ipiv[icol]);
    if (irow != icol) {
      for (l=0;l<n;l++) SWAP(a[irow][l],a[icol][l]);
      for (l=0;l<m;l++) SWAP(b[irow][l],b[icol][l]);
    }
    indxr[i]=irow;
    indxc[i]=icol;
    if ( isZero(a[icol][icol]) == true) {
      PRINTLN(icol);
      PRINTLN(a);
      throw Error("gaussj: Singular Matrix");
    }
    pivinv=T(1)/a[icol][icol];
    a[icol][icol]=T(1);
    for (l=0;l<n;l++) a[icol][l] *= pivinv;
    for (l=0;l<m;l++) b[icol][l] *= pivinv;
    for (ll=0;ll<n;ll++)
      if (ll != icol) {
	dum=a[ll][icol];
	a[ll][icol]=T(0);
	for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
	for (l=0;l<m;l++) b[ll][l] -= b[icol][l]*dum;
      }
  }
  for (l=n-1;l>=0;l--) {
    if (indxr[l] != indxc[l])
      for (k=0;k<n;k++)
	SWAP(a[k][indxr[l]],a[k][indxc[l]]);
  }
}


/**
Returns the solutions of the system linear equations Ax=b where A is given by the first n-q columns of 'a', b is given by the last q columns of 'a', and n is the number of columns in 'a'.
*/
template <typename T>
vector<vector<T> > linSolve(vector<vector<T> > a, int q)
{
  const int m = a.size();
  const int n = a[0].size();

  vector<int> indx(m);

  int i,j,k,maxrow, itemp;
  T tmp;

  for (int i=1; i<=m; i++) indx[i-1]=i-1;

  for (j=1;j<=n-q;j++) {

    maxrow = j-1;
    for (i=j+1;i<=m;i++) {
      if (abs(a[i-1][j-1])>abs(a[maxrow][j-1])) maxrow = i-1;
    }
     
    for (k=j;k<=n;k++) {
      tmp=a[j-1][k-1];
      a[j-1][k-1]=a[maxrow][k-1];
      a[maxrow][k-1]=tmp;
    }
    itemp=indx[j-1]; indx[j-1]=indx[maxrow]; indx[maxrow]=itemp;

    assert(!isZero(a[j-1][j-1]));
     
    for (i=j+1;i<=m;i++) {
      for (k=n;k>=j;k--) {
	a[i-1][k-1] = a[i-1][k-1] - a[i-1][j-1] * a[j-1][k-1] / a[j-1][j-1];
      }
    }

  }
 
  vector<T> x;
  vector<vector<T> > solutions;

  for (i=q; i>=1; --i) {

    x = vector<T>(n-q);

    for (j=n-q; j>=1; j--) {
      tmp=T(0);
      for (k=j+1;k<=n-q;k++) {
	tmp = tmp + a[j-1][k-1] * x[k-1];
      }
      x[j-1] = (a[j-1][n-i] - tmp) / a[j-1][j-1];
    }
    solutions.push_back(x);
  }

  return solutions;
   
}


/**
Returns the solution of the system linear equations Ax=b.
*/
template<typename T>
vector<T> linSolve(const vector<vector<T> > & A, const vector<T> & b)
{
  assert(A.size() == b.size());

  vector<vector<T> > A1(A);

  for (unsigned i=0; i<A1.size(); ++i) A1[i].push_back(b[i]);

  vector<T> x = linSolve<T>(A1, 1).at(0);

  return x;
   
}


/**
Returns the solution of the system linear equations Ax=b where b is a list of vectors.
*/
template<typename T>
vector<vector<T> > linSolve(const vector<vector<T> > & A, const vector<vector<T> > & b)
{
  assert(A.size() == b.at(0).size());

  vector<vector<T> > A1(A);

  for (int i=0; i<A1.size(); ++i) 
    for (int j=0; j<b.size(); ++j) 
      A1[i].push_back(b[j][i]);

  return linSolve<T>(A1, b.size());

}



/**
Returns the solution of the system linear equations Ax=b where A is a matrix of integers and b is a vector of doubles.
*/
doubleVector linSolve(const intMatrix & A, const doubleVector & B)
{
  return linSolve<double>(to_double(A), B);
}


/**
Returns the solution of the system linear equations Ax=b where A is a matrix of integers and b is a vector of rational numbers.
*/
rationalVector linSolve(const vector<vector<int> > & A, const rationalVector & B)
{
  return linSolve<Rational>(to_rational<int>(A), B);
}


/**
Returns the solution of the system linear equations Ax=b where A is a matrix of integers and b is a vector of integers.
*/
template <typename T>
vector<T> linSolve_return(const intMatrix & A, const intVector & B)
{
  vector<vector<T> > a(A.size(), vector<T>(A.at(0).size()));
  vector<T> b(B.size());

  for (int i=0; i<A.size(); ++i) {
    b[i] = T(B[i]);
    for (int j=0; j<A[i].size(); ++j)
      a[i][j] = T(A[i][j]);
  }

  return linSolve(a, b);
}


/**
Returns the transpose of a matrix.
*/
template <typename T>
vector<vector<T> > transpose(const vector<vector<T> > & A)
{
  if (A.size() == 0) return A;

  vector<vector<T> > C(A[0].size(), vector<T>(A.size()));

  for (int i=0; i<A[0].size(); ++i)
    for (int j=0; j<A.size(); ++j)
      C[i][j] = A[j][i];

  return C;

}


/**
Returns the trace of a matrix.
*/
template <typename T>
T Trace(const vector<vector<T> > A)
{
  assert(A.size() == A.at(0).size());
  
  T trace(0);
  for (unsigned i=0; i<A.size(); ++i)
    trace += A[i][i];

  return trace;
}


/**
Returns the inverse of a matrix.
*/
doubleMatrix inverse(const doubleMatrix & M)
{
  unsigned n = M.size();
  gsl_matrix * A = gsl_matrix_alloc (n,n);
  for (unsigned i = 0; i < n; ++i)
    for (unsigned j = 0; j < n; ++j)
      gsl_matrix_set (A, i, j, M[i][j]);
     
  gsl_permutation * p = gsl_permutation_alloc (n);
  int s;

  gsl_linalg_LU_decomp (A, p, &s);

  gsl_matrix * invA = gsl_matrix_alloc (n,n);
  gsl_linalg_LU_invert (A, p, invA);

  doubleMatrix C(n,doubleVector(n));
  for (unsigned i = 0; i < n; ++i)
    for (unsigned j = 0; j < n; ++j)
      C[i][j] = gsl_matrix_get (invA, i, j);

  gsl_matrix_free (A);
  gsl_matrix_free (invA);
  gsl_permutation_free (p);

  return C;
  
}


/**
Returns the inverse of a matrix with rational entries.
*/
rationalMatrix inverse(rationalMatrix A)
{
  rationalMatrix b(A.size(), rationalVector(1));

  gaussj(A, b);

  return A;
}


/**
Returns the inverse of a matrix with integer entries. The result is a matrix with rational entries.
*/
rationalMatrix inverse(const intMatrix & A)
{
  return inverse(to_rational(A));
}
  

template <typename T>
vector<vector<T> > identity_matrix(unsigned n)
{
  vector<vector<T> > A(n, vector<T>(n));
  for (unsigned i=0; i<n; ++i)
    A[i][i] = T(1);

  return A;

}


/**
Calculates the eigenvalues and eigenvectors of a normal matrix A (i.e. A*Transpose(A) - Transpose(A)*A =0).
*/
eigenSystem::eigenSystem(vector<vector<double> > mat) 
  : n(mat.size()), eigenvalues(mat.size()), eigenvectors(mat.size(), vector<double>(mat.size()))
{
  const unsigned NMAX = 250;
  
  #ifdef DEBUG
  cout << mat << endl;
  cout << "AAT- ATA" << endl << (mat*transpose(mat)) - (transpose(mat)*mat) << endl;
  #endif

  assert(n < NMAX);

  gsl_matrix * m = gsl_matrix_alloc (n,n);
  for (unsigned i = 0; i < n; ++i)
    for (unsigned j = 0; j < n; ++j)
      gsl_matrix_set (m, i, j, mat[i][j]);
     
  gsl_vector_complex *eval = gsl_vector_complex_alloc (n);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc (n, n);
  
  gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (n);
  gsl_eigen_nonsymmv (m, eval, evec, w);
  gsl_eigen_nonsymmv_free (w);
  gsl_eigen_nonsymmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);

  for (unsigned p=0; p<n; ++p) {

    gsl_complex z = gsl_vector_complex_get(eval, p);

    eigenvalues[p] = GSL_REAL(z);
    if (abs(GSL_IMAG(z)) > EPSILON) 
      throw Error("Error in routine eigensystemOfRealMatrix: Eigenvalues are complex ...");

    gsl_vector_complex_view evec_i = gsl_matrix_complex_column (evec, p);

    for (unsigned q=0; q<n; ++q) {
      z = gsl_vector_complex_get(&evec_i.vector, q);
      eigenvectors[p][q] = GSL_REAL(z);
      if (GSL_IMAG(z)>EPSILON) 
	throw Error("Error in routine eigensystemOfRealMatrix: Eigenvectors are complex ...");
    }

  }

  #ifdef DEBUG
  cout << "Eigenvalues ...\n"
       << "===============\n\n";
  for (unsigned i=0; i<eigenvalues.size(); ++i)
    cout << eigenvalues[i] << endl;
  #endif

  inverse_eigenvectors = inverse(eigenvectors);

}


/**
Diagonalizes a matrix and remembers the change of basis, i.e. if Diag = U*A*Inverse(U), the routine returns Diag and saves U in 'change_of_basis'.
*/
doubleMatrix diagonalize_matrix(const doubleMatrix & matrix, 
				doubleMatrix & change_of_basis)
{

  doubleMatrix eigenvectors;
  try {
    eigenSystem es(matrix);
    eigenvectors = es.eigenvectors;
  } catch (Error) {
    throw Error("Error in routine eigensystemOfRealMatrix: Eigenvectors are complex ...");
  }

  doubleMatrix diagonal = inverse(eigenvectors)*matrix*eigenvectors;
  set_double_zeros_to_zero(diagonal);

  doubleMatrix tempmat = inverse(eigenvectors);
  set_double_zeros_to_zero(tempmat);
  change_of_basis = tempmat;
  return diagonal;
}


/**
Gram-Schmidt orthogonalization.
*/
template <typename T>
vector<vector<T> > gramSchmidt_orthogonal(vector<vector<T> > A)
{
  typedef vector<T> VectorType;
  typedef vector<VectorType> MatrixType;

  assert(A.size() != 0);

  MatrixType C;

  for (unsigned i=0; i<A.size(); ++i) {
    VectorType b = A[i];
    for (unsigned k=0; k<b.size(); ++k) {
      for (unsigned j=0; j<C.size(); ++j) {
	T mu = (A[i]*C[j])/(C[j]*C[j]);
	b[k] -= mu*C[j][k];
      }
    }  
    C.push_back(b);
  }

  return C;  

}


/**
Returns the factorial of an integer.
*/
int factorial(const int n)
{
  assert(n>=0);

  if (n == 0)
    return 1;
  else
    return n * factorial(n-1);
}


/**
C program for distribution from the Combinatorial Object Server.
Generate combinations in lexicographic order. This is the same version
used in the book "Combinatorial Generation."  The program can be
modified, translated to other languages, etc., so long as proper
acknowledgement is given (author and source).  Programmer: Joe Sawada,
1997.  The latest version of this program may be found at the site
http://sue.uvic.ca/~cos/inf/comb/CombinationsInfo.html
*/
void Comb(vector<intVector> & result, int n, int k, vector<int> & a, int j, int m) {

  if (j > n) {
    intVector v(n);
    for (int i=0; i<n; i++) v[i] = a[i+1];
    result.push_back(v);
  }
  else {
    if (k-m < n-j+1) {
      a[j] = 0; Comb(result, n, k, a, j+1, m);
    }
    if (m<k) {
      a[j] = 1; Comb(result, n, k, a, j+1, m+1);
    }
  }
}


/* 
Returns the combinations k out of n elements (i.e. permutations without order).
*/
vector<intVector> combination(int n, int k) 
{
  vector<intVector> result;
  vector<int> a(n+1);

  if (n<=0) exit(1); 

  Comb(result, n, k, a, 1, 0);

  return result;
}


/* 
Returns the combinations of n elements (i.e. permutations without order).
*/
vector<intVector> combination(int n) 
{
  vector<intVector> result;
  
  if (n == 0) return result;

  for (int k=0; k<=n; ++k) {
    vector<intVector> tempresult = combination(n, k);
    for (vector<intVector>::iterator iter = tempresult.begin(); iter != tempresult.end(); ++iter) {
      result.push_back(*iter);
    }
  }
  return result;
}


/**
C program for distribution from the Combinatorial Object Server.
Generate permutations of a multiset in lexicographic order. This is
the same version used in the book "Combinatorial Generation."  The
program can be modified, translated to other languages, etc., so long
as proper acknowledgement is given (author and source).  Programmer:
Joe Sawada, 1997.  The latest version of this program may be found at
the site http://sue.uvic.ca/~cos/inf/mult/Multiset.html
*/
void permutation(int n, vector<int> & P, int & N, int & t, vector<int> & a, vector<intVector> & result) 
{
  if (a[0] == n) {
    intVector temp(N);
    for (int j=N-1; j>=0; j--) temp[N-1-j] = P[j]+1;
    result.push_back(temp);
  }
  else {
    for (int j=0; j<=t; j++)  {
      if (a[j] > 0) {
	P[n-1] = j;
	a[j]--;
	permutation(n-1, P, N, t, a, result);
	a[j]++;
	P[n-1] = 0;
      }
    }
  }
}


/**
Returns the permutations of a multiset. The elements of the
multiset are i = 1,...,N, where the multiplicity of i is a[i-1]. To
get the permutation of N elements, simply choose all multiplicities to
be 1, i.e. a[i] = 1.
*/
vector<intVector> permutation(const vector<int> & a1)
{
  vector<intVector> result;
  
  vector<int> P(100);	/* max size of perm = 100 */

  vector<int> a = a1;

  int n = 0;  
  for (unsigned i=0; i<a.size(); ++i)
    n += a[i];
  int N = n;
  int t = a.size()-1;
  permutation(n, P, N, t, a, result);

  return result;

}


/* 
Returns the permutations N elements.
*/
vector<intVector> permutation(int N)
{
  intVector a(N,1);

  return permutation(a);
}


/**
C program for distribution from Frank Ruskey's Combinatorial Object
Server.  Generates all numerical partitions of n whose largest part is
k.  No input error checking.  Assumes 0 <= k <= n <= MAX.  A simple
modification will generate all partitions of n (see comment at end of
program.)  Algorithm is CAT (Constant Amortized Time).  The program
can be modified, translated to other languages, etc., so long as
proper acknowledgement is given (author and source).  Algorithm and
original Pascal implementation: Frank Ruskey, 1995.  Translation to C:
Joe Sawada, 1997.  The latest version of this program may be found at
the site http://theory.cs.uvic.ca/~cos/inf/nump/NumPartition.html or
http://theory.uvic.ca/~cos/dis/programs.html
*/
void partition_greatest_part_k(int n, int k, int t, vector<int> & p, vector<intVector> & result) 
{

  p[t] = k;
  if (n==k) {
    intVector temp(t);
    for(int i=1; i<=t; ++i) temp[i-1] = p[i];
    result.push_back(temp);
  }
  for (int j=(k < n-k ? k : n-k); j>=1; j--) partition_greatest_part_k(n-k,j,t+1, p, result);
}


/**
Generates all numerical partitions of n whose largest part is k.
*/    
vector<intVector> partition_greatest_part_k(int n, int k)
{
  vector<int> p(100);
  vector<intVector> result;

  partition_greatest_part_k(n, k, 1, p, result);

  return result;
}


/**
Generates all numerical partitions of n.
*/    
vector<intVector> partition(int n)
{
  vector<int> p(100);
  vector<intVector> result;

  partition_greatest_part_k( 2*n, n, 0, p, result);

  return result;
}


/*******************************************************************************************************
   SOLVING DIOPHANTINE EQUATIONS
   =============================
The following routines solve one or more linear Diophantine equations. The general solution is returned.
Based on: "A class of ABS algorithms for Diophantine linear systems", Numer.Math.90:101-115 (2001) by
H. Esmaeili, N. Mahdavi-Amiri, E. Spedicato.
********************************************************************************************************/

/**
Implementation of a criterion used in following routine 'gcd_rossers_algorithm'.
*/
template <typename T>
bool gcd_rossers_algorithm_vector_greater_criterion(vector<T> v, vector<T> w) 
{
  return (v.at(0) > w.at(0));
}


/**
Solves in principle the equation x0*b[0] + ... + xn*b[n] = c.
Solution space is spanned by
--> ( c / gcd ) * basis[0] + any integer * basis[1] + ...
Returns gcd and basis !
*/
pair<int, intMatrix> gcd_rossers_algorithm(const vector<int> & b)
{
  const int n = b.size();

  // Keep track of the signs; algorithm can handle only non-negative numbers
  intVector signindex(n);
  for (int i=0; i<n; ++i) {
    if (b[i] != 0) signindex[i] = b[i]/abs(b[i]);
  }

  // Keep track of the variables
  vector<intVector> ordering(n, intVector(2));
  for (int i=0; i<n; ++i) {
    ordering[i][0] = abs(b[i]);
    ordering[i][1] = i;
  }

  sort(ordering.begin(), ordering.end(), ptr_fun(gcd_rossers_algorithm_vector_greater_criterion<int>));

  intVector a(n), index(n); 
  for (int i=0; i<n; ++i) {
    a[i] = ordering[i][0];
    index[i] = ordering[i][1];
  }

  assert(a[0] > 0);

  intMatrix D_t(n, intVector(n+1));
  for (int i=0; i<n; ++i) {
    D_t[i][0] = a[i];
    D_t[i][i+1] = 1;
  }  
  intMatrix D = transpose(D_t);
  
  intMatrix C, C_cols;

  // Step 1

  C = D;
  C_cols = transpose(C);

  // Step 2

  while (C[0][1] != 0) {

    C_cols[0] = C_cols[0] - to_integer(floor(double(C[0][0])/double(C[0][1])))*C_cols[1];

    sort(C_cols.begin(), C_cols.end(), ptr_fun(gcd_rossers_algorithm_vector_greater_criterion<int>));

    C = transpose(C_cols);

  }

  // Step 3
  
  intVector d = C[0]; // First element is gcd

  intMatrix U(n);
  copy(C.begin()+1, C.end(), U.begin());

  intMatrix basis(n);
  for (int i=0; i<basis.size(); ++i) {
    basis[index[i]] = U[i];
  }

  for (int i=0; i<n; ++i) {
    if (signindex[i] == -1)
      basis[i] = (-1)*basis[i];
  }

  basis = transpose(basis); 

  pair<int, intMatrix> solution;
  solution.first = d[0];
  solution.second = basis;

  return solution;
}


/**
Returns the greatest common divisor of a vector of integers.
*/
int gcd(const vector<int> & b)
{
  return gcd_rossers_algorithm(b).first;
}


/**
Solves a single diophantine equation.
*/
pair<intVector, intMatrix> solution_of_single_diophantine_equation(const vector<int> & b, int c)
{
  pair<int, intMatrix> solution1 = gcd_rossers_algorithm(b);

  if (c % solution1.first != 0) 
    throw Error("In \'solution_of_single_diophantine_equation\' : Equation has no solution ...");

  pair<intVector, intMatrix> solution2;
  solution2.first = (c/solution1.first)*solution1.second.at(0);
  solution2.second.resize(solution1.second.size()-1);
  copy(solution1.second.begin()+1, solution1.second.end(), solution2.second.begin());

  return solution2;
}


/**
Returns a particular solution of a diophantine equation.
*/
intVector particular_solution_of_single_diophantine_equation(const vector<int> & b, int c)
{
  return solution_of_single_diophantine_equation(b, c).first;
}


/**
Solves the system of linear equations Ax = b. The particular solution is stored in pair.first, all solutions are obtained by adding to the particular solutions arbitrary linear combinations of vectors stored in the rows (!) of the matrix pair.second, i.e. first + c0*second[0] + c1*second[1] + ...
*/
pair<intVector, intMatrix> solve_system_of_linear_diophantine_equations(const vector<vector<int> > & A, const vector<int> & b)
{
  int m = A.size();
  int n = A.at(0).size();
  
  int r, t, d, alpha, divisor;
  vector<int> x, s, p, z, w, vec1, vec2;
  vector<vector<int> > H, W, mat;
  
  // step_1

  x = vector<int>(n);
  H = identity_matrix<int>(n);
  r = 0;

 
  for (int i=0; i<m; ++i) {

    // step_2

    t = A[i]*x - b[i];
    s = H*A[i];
  
    // step_3

    if (isZero(s)) {
      if (t == 0) {
	continue;
      } else {
	throw Error("Equations incompatible ...");
      }
    }

    // step_4

    d = gcd(s);

    z = particular_solution_of_single_diophantine_equation(s, d);

    p = transpose(H)*z;

    if (t % d != 0) {
      throw Error("Equations incompatible ...");
    }
    else {
      alpha = t/d;
      x = x - alpha*p;
    }

    // step_5

    w = particular_solution_of_single_diophantine_equation(s, d);

    vec1 = H*A[i];

    vec2 = transpose(H)*w;

    mat = tensor_product(vec1, vec2);

    divisor = w*(H*A[i]);
    for (int k=0; k<mat.size(); ++k) {
      for (int l=0; l<mat[k].size(); ++l) {
	if (mat[k][l] % divisor == 0) {
	  mat[k][l] /= divisor;
	} else {
	  throw Error("Division gives rest ...");
	}
      }
    }

    H = H - mat;

    // step_6

    ++r;

  } // <-- step_7

  pair<intVector, intMatrix> solution;
  solution.first = x;
  solution.second = H;

  return solution;

}


/********************** 
EXPLICIT INSTANTIATIONS
**********************/


template
void gaussj(vector<vector<Rational> > & a, vector<vector<Rational> > & b);

template
void gaussj(vector<vector<double> > & a, vector<vector<double> > & b);

template 
vector<vector<double> > linSolve(vector<vector<double> > a, int q);

template
vector<vector<Rational> > linSolve(vector<vector<Rational> > a, int q);

template
vector<vector<complex<Rational> > > linSolve(vector<vector<complex<Rational> > > a, int q);

template
vector<Rational> linSolve(const vector<vector<Rational> > & A, const vector<Rational> & b);

template
vector<complex<Rational> > linSolve(const vector<vector<complex<Rational> > > & A, const vector<complex<Rational> > & b);

template
vector<double> linSolve(const vector<vector<double> > & A, const vector<double> & b);

template
vector<vector<double> > linSolve(const vector<vector<double> > & A, const vector<vector<double> > & b);

template
vector<vector<Rational> > linSolve(const vector<vector<Rational> > & A, const vector<vector<Rational> > & b);

template
vector<double> linSolve_return(const intMatrix & A, const intVector & B);

template
vector<Rational> linSolve_return(const intMatrix & A, const intVector & B);

template
vector<vector<rational<int> > > findBasis(vector<vector<rational<int> > > a);

template
vector<vector<double> > findBasis(vector<vector<double> > a);

template
vector<vector<rational<int> > > gramSchmidt_orthogonal(vector<vector<rational<int> > >);

template
vector<vector<double> > gramSchmidt_orthogonal(vector<vector<double> >);

template
vector<vector<rational<int> > > identity_matrix(unsigned n);

template
vector<vector<double> > identity_matrix(unsigned n);

template
vector<vector<int> > identity_matrix(unsigned n);

template
vector<vector<int> > transpose(const vector<vector<int> > & A);

template
vector<vector<double> > transpose(const vector<vector<double> > & A);

template
vector<vector<Rational> > transpose(const vector<vector<Rational> > & A);

template
vector<vector<complex<Rational> > > transpose(const vector<vector<complex<Rational> > > & A);

template
rational<int> Trace(const vector<vector<rational<int> > > A);

template
complex<rational<int> > Trace(const vector<vector<complex<rational<int> > > > A);

template
double Trace(const vector<vector<double> > A);

template
complex<double> Trace(const vector<vector<complex<double> > > A);

template
bool isZero(const vector<int> & t);

template
bool isZero(const vector<double> & t);

template
bool isZero(const vector<rational<int> > & t);

template
bool isZero(const vector<complex<double> > & t);

template
bool isZero(const vector<complex<rational<int> > > & t);

template
bool isZero(const vector<vector<int> > & t);

template
bool isZero(const vector<vector<double> > & t);

template
bool isZero(const vector<vector<rational<int> > > & t);

template
bool isZero(const vector<vector<complex<double> > > & t);

template
bool isZero(const vector<vector<complex<rational<int> > > > & t);

template
bool is_constant(const vector<double> & v);

template
bool is_constant(const vector<Rational> & v);

template
bool is_constant(const vector<int> & v);

template
bool isEqual(vector<rational<int> > u, vector<rational<int> > v);

template
bool isEqual(vector<double> u, vector<double> v);

template
bool isEqual(vector<int> u, vector<int> v);

template
bool isEqual(vector<complex<double> > u, vector<complex<double> > v);

template
bool is_orthogonal(const vector<double> & v, const vector<double> & w);

template
bool is_orthogonal(const vector<Rational> & v, const vector<Rational> & w);

template
bool is_orthogonal(const vector<double> & v, const vector<vector<double> > & A);

template
bool is_orthogonal(const vector<Rational> & v, const vector<vector<Rational> > & A);

template
bool is_half_integer(const double & t);

template
bool is_half_integer(const Rational & t);


// Conversions

template
rationalVector to_rational(const vector<int> & d);

template
rationalVector to_rational(const vector<double> & d);

template
rationalVector to_rational(const vector<complex<rational<int> > > & d);

template
rationalMatrix to_rational(const vector<vector<int> > & d);

template
rationalMatrix to_rational(const vector<vector<double> > & d);

template
rationalMatrix to_rational(const vector<vector<complex<rational<int> > > > & d);

template
intVector to_integer(const vector<double> & t);

template
intVector to_integer(const vector<rational<int> > & t);

template
intMatrix to_integer(const vector<vector<double> > & t);

template
intMatrix to_integer(const vector<vector<rational<int> > > & t);


template
doubleVector to_double(const vector<double> & t);

template
doubleVector to_double(const vector<int> & t);

template
doubleVector to_double(const vector<rational<int> > & t);

template
doubleVector to_double(const vector<complex<rational<int> > > & t);

template
doubleVector to_double(const vector<complex<double> > & t);

template
doubleMatrix to_double(const vector<vector<int> > & t);

template
doubleMatrix to_double(const vector<vector<rational<int> > > & t);

template
doubleMatrix to_double(const vector<vector<complex<rational<int> > > > & t);

template
doubleMatrix to_double(const vector<vector<complex<double> > > & t);

template
vector<doubleMatrix> to_double(const vector<vector<vector<complex<rational<int> > > > > & A);

// Scalar - vector multiplication

template
vector<int> operator*(const int k, vector<int> a);

template
vector<double> operator*(const double k, vector<double> a);

template
vector<rational<int> > operator*(const rational<int> k, vector<rational<int> > a);

template
vector<complex<double> > operator*(const complex<double>  k, vector<complex<double> > a);

template
vector<complex<rational<int> > > operator*(const complex<rational<int> >  k,
					   vector<complex<rational<int> > > a);


// Scalar - matrix multiplication

template
vector<vector<double> > operator*(const double k, vector<vector<double> > A);

template
vector<vector<rational<int> > > operator*(const rational<int> k, vector<vector<rational<int> > > A);

template
vector<vector<complex<double> > > operator*(const complex<double>  k, vector<vector<complex<double> > > A);

template
vector<vector<complex<rational<int> > > > operator*(const complex<rational<int> >  k, vector<vector<complex<rational<int> > > > A);


// Matrix - matrix multiplication

template
vector<vector<int> > operator*(const vector<vector<int> > & A, const vector<vector<int> > & B);

template
vector<vector<double> > operator*(const vector<vector<double> > & A, const vector<vector<double> > & B);

template
vector<vector<rational<int> > > operator*(const vector<vector<rational<int> > > & A, const vector<vector<rational<int> > > & B);

template
vector<vector<complex<double> > > operator*(const vector<vector<complex<double> > > & A, const vector<vector<complex<double> > > & B);

template
vector<vector<complex<rational<int> > > > operator*(const vector<vector<complex<rational<int> > > > & A, const vector<vector<complex<rational<int> > > > & B);

// Vector - vector multiplication

template
int operator*(const vector<int> & a, const vector<int> & b);

template
double operator*(const vector<double> & a, const vector<double> & b);

template
rational<int>  operator*(const vector<rational<int> > & a, const vector<rational<int> > & b);

template
complex<double>  operator*(const vector<complex<double> > & a, const vector<complex<double> > & b);

template
complex<rational<int> >  operator*(const vector<complex<rational<int> > > & a, const vector<complex<rational<int> > > & b);

// Matrix - vector multiplication

template
vector<int> operator*(const vector<vector<int> > & A, const vector<int> & b);

template
vector<double> operator*(const vector<vector<double> > & A, const vector<double> & b);

template
vector<rational<int> > operator*(const vector<vector<rational<int> > > & A, const vector<rational<int> > & b);

template
vector<complex<double> > operator*(const vector<vector<complex<double> > > & A, const vector<complex<double> > & b);

template
vector<complex<rational<int> > > operator*(const vector<vector<complex<rational<int> > > > & A, const vector<complex<rational<int> > > & b);

template
vector<vector<int> > tensor_product(const vector<int> & a, const vector<int> & b);

// Vector - vector addition

template
vector<int> operator+(const vector<int> & a, const vector<int> & b);

template
vector<double> operator+(const vector<double> & a, const vector<double> & b);

template
vector<rational<int> > operator+(const vector<rational<int> > & a, const vector<rational<int> > & b);

template
vector<complex<double> > operator+(const vector<complex<double> > & a, const vector<complex<double> > & b);

template
vector<complex<rational<int> > > operator+(const vector<complex<rational<int> > > & a, 
					   const vector<complex<rational<int> > > & b);


// Vector - vector substraction

template
vector<int> operator-(const vector<int> & a, const vector<int> & b);

template
vector<double> operator-(const vector<double> & a, const vector<double> & b);

template
vector<rational<int> > operator-(const vector<rational<int> > & a, const vector<rational<int> > & b);

template
vector<complex<double> > operator-(const vector<complex<double> > & a, const vector<complex<double> > & b);

template
vector<complex<rational<int> > > operator-(const vector<complex<rational<int> > > & a, 
					   const vector<complex<rational<int> > > & b);


// Matrix - matrix addition

template
vector<vector<int> > operator+(const vector<vector<int> > & A, const vector<vector<int> > & B);

template
vector<vector<double> > operator+(const vector<vector<double> > & A, const vector<vector<double> > & B);

template
vector<vector<rational<int> > > operator+(const vector<vector<rational<int> > > & A, const vector<vector<rational<int> > > & B);

template
vector<vector<complex<double> > > operator+(const vector<vector<complex<double> > > & A, const vector<vector<complex<double> > > & B);

template
vector<vector<complex<rational<int> > > > operator+(const vector<vector<complex<rational<int> > > > & A, const vector<vector<complex<rational<int> > > > & B);


// Matrix - matrix subtraction

template
vector<vector<int> > operator-(const vector<vector<int> > & A, const vector<vector<int> > & B);

template
vector<vector<double> > operator-(const vector<vector<double> > & A, const vector<vector<double> > & B);

template
vector<vector<rational<int> > > operator-(const vector<vector<rational<int> > > & A, const vector<vector<rational<int> > > & B);

template
vector<vector<complex<double> > > operator-(const vector<vector<complex<double> > > & A, const vector<vector<complex<double> > > & B);

template
vector<vector<complex<rational<int> > > > operator-(const vector<vector<complex<rational<int> > > > & A, const vector<vector<complex<rational<int> > > > & B);

// Comparisons

template
bool isGreaterZero(const vector<int> & t);

template
bool isGreaterZero(const vector<double> & t);

template
bool isGreaterZero(const vector<rational<int> > & t);

template
bool isGreaterEqualZero(const vector<int> & t);

template
bool isGreaterEqualZero(const vector<double> & t);

template
bool isGreaterEqualZero(const vector<rational<int> > & t);

//////////// CRITERIONS AND FUNCTION OBJECTS //////////////

// template
// SWAP<double>;

// template
// SWAP<Rational>;

template
class sp_with_shift_is_not_integer_crit<double>;

template
class sp_with_shift_is_not_integer_crit<Rational>;

template
class sp_with_shift_is_not_zero_crit<double>;

template
class sp_with_shift_is_not_zero_crit<Rational>;

template
class weight_wrt_semiordering_is_negative_or_zero_crit<double>;

template
class weight_wrt_semiordering_is_negative_or_zero_crit<Rational>;

template
class is_zero_crit<double>;

template
class is_zero_crit<Rational>;

template
class is_not_zero_crit<double>;

template
class is_not_zero_crit<Rational>;

template
class smaller_wrt_root_ordering_crit<double>;

template
class smaller_wrt_root_ordering_crit<Rational>;


