#ifndef GROUPTHEORY_H
#define GROUPTHEORY_H

#include "io.hpp"
#include "linalg.hpp"


/** CLASS DYNKIN **
Implements the Dynkin algorithm to calculate the weights of any irreducible representation for any Lie algebra.
*/

class dynkin {
  void dynkin_recursion(intVector lambda);
  unsigned rank;
  intMatrix A;
  set<intVector> weights;
 public:
  dynkin(rationalMatrix simpleroots, intVector highestweight);
  dynkin(doubleMatrix simpleroots, intVector highestweight);
  dynkin(intMatrix simpleroots, rationalMatrix QFM, intVector highestweight);
  intMatrix getA() const;
  set<intVector> getWeights() const;
};



/** CLASS FREUDENTHAL **
Implements the Freudenthal algorithm that calculates the *multiplicity* with which each weight appears in an irreducible representation. Used together with 'class dynkin'.
*/

class freudenthal {
  unsigned rank;
  unsigned dimension;
  intVector highestweight; ///< The highest weight vector
  set<intVector> posweights;
  set<intVector> adjointweights;
  map<intVector, int> weightmap;
  void freudenthal_recursion(intVector lambda);
  intVector add(intVector a, intVector b);
  intVector mult(int k, const intVector a);
  void setPosWeights();
  int multiplicity(intVector w);
  void calculateMultiplicities();
public:
  intVector highestroot;
  rational<int> sp(intVector a, intVector b);
  rational<int> sp(rationalVector a, rationalVector b);
  rational<int> sp(intVector a, rationalVector b);
  rational<int> sp(rationalVector a, intVector b);
  intMatrix A;
  rationalMatrix Q;
  intVector rho;
  set<intVector> weights;
  multiset<intVector> multiweights;
  void print() const;
  freudenthal(doubleMatrix simpleroots, intVector highestroot1, intVector highestweight1);
  freudenthal(rationalMatrix simpleroots, intVector highestroot1, intVector highestweight1);
  freudenthal(intMatrix simpleroots, rationalMatrix QFM, intVector highestroot1, intVector highestweight1);
  set<intVector> getPosWeights() { return posweights; }
  map<intVector, int> getWeightMap();
};


/** CLASS CIRREP **
Encodes basic information of an irreducible representation.
*/

template <typename T>
class cirrep {
 public:
  vector<T> highestweight;
  vector<T> u1charge;
};


/** CLASS GAUGE GROUP **
Encodes basic information on a gauge group like the rank the simple roots, the human readable name according to the ADE classification, the semi-ordering on the weight lattice, etc.
*/

template <typename T>
class gaugeGroupFactor {
 public:
  bool operator<(const gaugeGroupFactor & w) const { return this->algebra < w.algebra; }
  string algebra;
  unsigned rank;
  vector<vector<T> > simpleroots;
};

template <typename T>
class gaugeGroup {
private:
  typedef vector<T> VectorType;
  typedef vector<VectorType> MatrixType;

  vector<T> semiordering;
public:
  string algebra;
  vector<gaugeGroupFactor<T> > factor;
  MatrixType basis;
  MatrixType u1directions;

  MatrixType getSimpleRoots() const 
  {
    MatrixType simpleroots;
    for (unsigned i=0; i<factor.size(); ++i)
      copy(factor[i].simpleroots.begin(), factor[i].simpleroots.end(), back_inserter(simpleroots));
    
    return simpleroots;
  };

  void set_semiordering(const vector<T> & X) { semiordering = X; };
  VectorType get_semiordering(void) { return semiordering; };
};



/** CLASS WEIGHTSYSTEM **
Basic routines concerning weight systems like finding positive roots, simple roots, etc.
*/

template <typename T>
class weightsystem {
private:
  typedef vector<T> VectorType;
  typedef vector<VectorType> MatrixType;

  bool isLess(const VectorType & a, const VectorType & b) {
    return isPositive(b-a);
  };

  MatrixType findPositiveRoots(const MatrixType & weights);
  void findSimpleRoots(const MatrixType & weights);
  void setDynkinLabels(const MatrixType & weights);

  VectorType semiordering; // Defines semi-ordering for root vectors
public:

  bool isElement(const VectorType & a) { 
    return (find(weights.begin(), weights.end(), a) != weights.end()); 
  };

  unsigned pos(const vector<T> & vec) {
    typename MatrixType::iterator iter = find(weights.begin(), weights.end(), vec);
    assert(iter != weights.end());
    return (iter - weights.begin());
  };
  
  static bool isGreater(const VectorType & a, const VectorType & b) {
    return isPositive(a-b);
  };

  static bool isPositive(const VectorType & a);
  MatrixType weights;
  MatrixType basis;
  unsigned rank;
  MatrixType simpleroots;
  map<VectorType, intVector> dl;
  map<VectorType, VectorType> rl;
  weightsystem(const MatrixType & weights1);
  vector<T> get_semiordering(void) { return semiordering; };
};


/** CLASS DYNKIN DIAGRAM **
Encodes the Dynkin diagram of the Lie algebras; routines to identify the unbroken gauge group from its surviving roots.
*/

template <typename T>
class dynkinDiagram {
  typedef vector<T> VectorType;
  typedef vector<VectorType> MatrixType;

  void subAlgebra(set<VectorType> & allroots, MatrixType & subalgebra, VectorType node);
  string determineAlgebra(MatrixType subalgebra);
  int sp(VectorType a, VectorType b);
  unsigned rank;
  string algebra;
  map<VectorType, MatrixType> adjacent;
  VectorType findHead(string alg, const MatrixType & sroots);
  void reorderDiagram(string alg, MatrixType & sroots);
  intMatrix A;
  gaugeGroup<T> group;
public:
  void setBasis(const MatrixType & basis1) { group.basis = basis1; }
  dynkinDiagram(MatrixType & simpleroots);
  dynkinDiagram(MatrixType & simpleroots, const MatrixType & cartanMatrix);
  const gaugeGroup<T> & getGaugeGroup() const { return group; }
};


/******************** 
FUNCTION DECLARATIONS
********************/

template <typename T>
vector<vector<T> > create_simpleroots_E8();

template <typename T>
vector<vector<T> > create_roots_E8();

template <typename T>
vector<vector<T> > create_simpleroots_E8xE8();

template <typename T>
vector<vector<T> > create_basis_Spin32();

rationalMatrix calculateCartanMatrix(const rationalMatrix & sroots);

doubleMatrix calculateCartanMatrix(const doubleMatrix & sroots);

template <typename T>
intMatrix setDynkinLabels(vector< vector<T> > simpleroots, vector< vector<T> > weights);

intMatrix findHighestWeights(const multiset<intVector> & weights, intMatrix newroots, 
				       rationalMatrix Q, intMatrix & original_hw);

multiset<intVector> findHighestWeights(const multiset<intVector> & weights);

vector<vector<int> > findHighestWeights(const vector<vector<int> > & weights);

template <typename T>
gaugeGroup<T> determineAlgebra(const vector<vector<T> > & wsystem);

template <typename T>
std::vector<int> findDynkinLabels(const gaugeGroupFactor<T> & myfactor, const std::vector<T> & mw);

template <typename T>
std::vector< std::vector<int> > findDynkinLabels(const gaugeGroupFactor<T> & myfactor, 
						 const std::vector< std::vector<T> > & mw);

template <typename T>
vector<int> findDynkinLabels(const vector<vector<T> > & myfactor, 
				      const vector<T> & mw);

template <typename T>
vector<vector<int> > findDynkinLabels(const vector<vector<T> > & myfactor, 
				      const vector<vector<T> > & mw);

template <typename T>
set<vector<int> > dynkinAlgorithm(const gaugeGroupFactor<T> & myfactor, 
					     const vector<int> & lambda);

template <typename T>
cirrep<T> findHighestWeight(gaugeGroup<T> myGaugeGroup, vector<vector<T> > myWeights);

intVector highestRoot(string algebra, unsigned rank);

template <typename T>
vector<T> KacLabels(vector<vector<T> > & simpleroots);

template <typename T>
vector<T> extendedKacLabels(vector<vector<T> > & simpleroots);

template <typename T>
vector<T> extendedKacLabels(gaugeGroup<T> & group);

template <typename T>
ostream & operator<<(ostream & os, const class gaugeGroup<T> & A);

#endif
