#ifndef CLINALG_H
#define CLINALG_H

#include <vector>
#include <boost/config.hpp>
#include <boost/rational.hpp>
#include "chugeint.h"


using boost::rational;
using std::vector;

template <class T>
class CLinAlg{
public:
    CLinAlg();
    ~CLinAlg();

  bool FindKernel(const vector<vector<rational<T> > > &M, vector<vector<rational<T> > > &BasisOfKernel, bool TestResult = true);
  bool FindPositiveIntegerKernel(const vector<vector<rational<T> > > &M, vector<vector<rational<T> > > &BasisOfKernel);
  bool GetParticularSolution(const vector<vector<rational<T> > > &A, const vector<rational<T> > &b, vector<rational<T> > &ParticularSolution);
  bool GramSchmidt(vector<vector<rational<T> > > Basis, vector<vector<rational<T> > > &OrthogonalBasis, bool TestResult = true);
  bool ReducedRowEchelonForm(vector<vector<rational<T> > > &M);
  bool ScaleToIntegerVector(vector<rational<T> > &Vector);

  bool FindOrdersOfWilsonLines(const vector<vector<rational<T> > > &input, vector<vector<rational<T> > > &output);
  bool SimplifyRelations(vector<vector<rational<T> > > &M);
  bool SimplifyInequality(vector<vector<rational<T> > > &M);
  void SolveSystemOfLinearInequalities(const vector<vector<rational<T> > > &A, vector<vector<rational<T> > > &Cone);

  const rational<T> One;
  const rational<T> Zero;
};

#endif
