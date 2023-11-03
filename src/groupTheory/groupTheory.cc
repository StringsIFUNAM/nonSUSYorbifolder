#include "groupTheory.hpp"

extern unsigned SELFDUALLATTICE;


/*************************************** 
IMPLEMENTATION OF THE CLASS WEIGHTSYSTEM 
***************************************/

/**
Constructor. Initializes the weight system with the vectors given in 'weights1'. Defines a semi-ordering on the E8xE8 or Spin(32) lattice depending on the value of the parameter 'SELFDUALLATTICE'.
*/
template <typename T>
weightsystem<T>::weightsystem(const MatrixType & weights1) : weights(weights1), semiordering(weights1.at(0).size())
{
  vector<T> tempX(16);

  if (SELFDUALLATTICE == 1) {

    tempX[0]  = T(23); tempX[1]  = T(6); tempX[2]  = T(5); tempX[3]  = T(4);   
    tempX[4]  = T(3);  tempX[5]  = T(2); tempX[6]  = T(1); tempX[7]  = T(0);

    for (int i=0; i<weights.at(0).size(); ++i) semiordering[i] = tempX[i % 8];

  } else if (SELFDUALLATTICE == 2) {

    tempX[0]  = T(23); tempX[1]  = T(14); tempX[2]  = T(13); tempX[3]  = T(12);   
    tempX[4]  = T(11); tempX[5]  = T(10); tempX[6]  = T(9);  tempX[7]  = T(8);
    tempX[8]  = T(7);  tempX[9]  = T(6);  tempX[10] = T(5);  tempX[11] = T(4);   
    tempX[12] = T(3);  tempX[13] = T(2);  tempX[14] = T(1);  tempX[15] = T(0);

    for (int i=0; i<weights.at(0).size(); ++i) semiordering[i] = tempX[i];

  } else {
    throw Error("Initial algebra must be E8xE8 or SO(32) ...");
  }

  basis = findBasis<T>(weights);

  rank = basis.size();

  if (rank > 0) {
    this->findSimpleRoots(weights);
    this->setDynkinLabels(weights);
  }

}


/**
Returns true if the weight is positive with respect to the semi-ordering of the lattice.
*/
template <typename T>
bool weightsystem<T>::isPositive(const VectorType & a)
{
  typename VectorType::const_iterator iter = a.begin(); 
  while ( (isZero(*iter)) && (iter != a.end()) ) ++iter;
  
  if ( (iter != a.end()) && (isGreaterZero(*iter)) ) return true;
  else return false;
  
}


/**
Returns the positive roots of the roots given in 'weights'.
*/
template <typename T>
typename weightsystem<T>::MatrixType weightsystem<T>::findPositiveRoots(const MatrixType & weights)
{
  MatrixType posRoots;
  for (typename MatrixType::const_iterator iter = weights.begin(); iter != weights.end(); ++iter) {
    if ( semiordering * (*iter) > 0 ) posRoots.push_back(*iter);
  }

  return posRoots;

}


/**
Returns the simple roots of the roots given in 'weights'.
*/
template <typename T>
void weightsystem<T>::findSimpleRoots(const MatrixType & weights)
{
  MatrixType positiveRoots = findPositiveRoots(weights);
  
  set<VectorType> candidates;
  typename MatrixType::iterator iter1;
  typename MatrixType::iterator iter2;
  for (iter1 = positiveRoots.begin(); iter1 != positiveRoots.end(); ++iter1) candidates.insert(*iter1);
  
  for (iter1 = positiveRoots.begin(); iter1 != positiveRoots.end(); ++iter1)
    for (iter2 = iter1; iter2 != positiveRoots.end(); ++iter2)
      candidates.erase((*iter1) + (*iter2));

  assert(rank == candidates.size());

  simpleroots.resize(candidates.size());

  copy(candidates.begin(), candidates.end(), simpleroots.begin());

}


/**
Calculates the Dynkin labels of the weights given in 'weights'.
*/
template <typename T>
void weightsystem<T>::setDynkinLabels(const MatrixType & weights)
{
  intVector x(rank);

  typename MatrixType::const_iterator iter;
  for (iter = weights.begin(); iter != weights.end(); ++iter) {
    for (int i=0; i<rank; ++i) {
      x[i] = to_integer(simpleroots[i]*(*iter));
     }
    dl[*iter] = x;
  }

}


/********************************* 
IMPLEMENTATION OF THE CLASS DYNKIN
*********************************/


/**
Constructor. Sets the simple roots and the highest weight that defines the irreducible representation.
*/
dynkin::dynkin(rationalMatrix simpleroots, intVector highestweight) 
  : A(simpleroots.size(),intVector(simpleroots.size()))
{
  rank = simpleroots.size();

  assert(rank == highestweight.size());

  rationalMatrix transposesimpleroots = transpose(simpleroots);
  rationalMatrix cartanMatrix = simpleroots*transposesimpleroots;
  for (unsigned i=0; i<rank; ++i) 
    for (unsigned j=0; j<rank; ++j) {
      A[i][j] = to_integer(cartanMatrix[i][j]);
    }

  dynkin_recursion(highestweight);

  return;

}


/**
Constructor. Sets the simple roots and the highest weight that defines the irreducible representation.
*/
dynkin::dynkin(doubleMatrix simpleroots, intVector highestweight) 
  : A(simpleroots.size(),intVector(simpleroots.size()))
{
  rank = simpleroots.size();

  assert(rank == highestweight.size());

  doubleMatrix transposesimpleroots = transpose(simpleroots);
  doubleMatrix cartanMatrix = simpleroots*transposesimpleroots;
  for (unsigned i=0; i<rank; ++i) 
    for (unsigned j=0; j<rank; ++j) A[i][j] = roundx(cartanMatrix[i][j]);

  dynkin_recursion(highestweight);

  return;

}


/**
Constructor. Sets the simple roots and the highest weight that defines the irreducible representation.
*/
dynkin::dynkin(intMatrix simpleroots, rationalMatrix QFM, intVector highestweight) 
  : A(simpleroots.size(),intVector(simpleroots.size()))
{
  rank = simpleroots.size();

  assert((rank == highestweight.size()));

  rational<int> sum(0);
  intVector a,b;

  for (unsigned i=0; i<rank; ++i) {
    for (unsigned j=0; j<rank; ++j) { 
      a = simpleroots[i];
      b = simpleroots[j];
      sum = rational<int>(0);
      for (unsigned p=0; p<a.size(); ++p) for (unsigned q=0; q<a.size(); ++q) sum += a[p]*QFM[p][q]*b[q];
      assert(sum.denominator() == 1);
      A[i][j] = sum.numerator();
    }
  }

  dynkin_recursion(highestweight);

  return;

}

/**
Given the highest weight calculates all the weights of an irreducible representation.
*/
void dynkin::dynkin_recursion(intVector lambda) 
{ 
  weights.insert(lambda);
  intVector mu(rank);
  for ( unsigned k=1; k<=rank; k++) {
    if (lambda[k-1]>0) {
      for (int i=1; i<=lambda[k-1]; i++) {
	for (unsigned j=1; j<=rank; j++) mu[j-1]=lambda[j-1]-i*A[k-1][j-1]; 
	if (weights.find(mu) == weights.end()) {
	  dynkin_recursion(mu);
	}
      }
    }
  }
  return;
}


/**
Returns the Cartan matrix.
*/
intMatrix dynkin::getA() const
{
  return A;
}


/**
Returns the weights.
*/
std::set<intVector> dynkin::getWeights() const
{
  return weights;
}


/******************************** 
IMPLEMENTATION OF THE FREUDENTHAL
********************************/


/**
Constructor. Sets the simple roots and the highest weight that defines the irreducible representation. Note that the Freudenthal algorithm not only gives the weights of a representation, but only their multiplicities.
*/
freudenthal::freudenthal(doubleMatrix simpleroots, intVector highestroot1, intVector highestweight1) 
  : rho(simpleroots.size(),1)
{
  rank = simpleroots.size();
  assert( (rank == highestroot1.size()) && (rank == highestweight1.size()) );

  highestroot = highestroot1;
  highestweight = highestweight1;
  
  dimension = 0;

  dynkin a(simpleroots, highestweight);
  weights = a.getWeights();

  dynkin b(simpleroots, highestroot);
  adjointweights = b.getWeights();

  A = a.getA();

  Q = inverse(A);

  setPosWeights();

  calculateMultiplicities();

}


/**
Constructor. Sets the simple roots and the highest weight that defines the irreducible representation. Note that the Freudenthal algorithm not only gives the weights of a representation, but only their multiplicities.
*/
freudenthal::freudenthal(rationalMatrix simpleroots, intVector highestroot1, intVector highestweight1) 
  : rho(simpleroots.size(),1)
{
  rank = simpleroots.size();
  assert( (rank == highestroot1.size()) && (rank == highestweight1.size()) );

  highestroot = highestroot1;
  highestweight = highestweight1;
  
  dimension = 0;

  dynkin a(simpleroots, highestweight);
  weights = a.getWeights();

  dynkin b(simpleroots, highestroot);
  adjointweights = b.getWeights();

  A = a.getA();

  Q = inverse(A);

  setPosWeights();

  calculateMultiplicities();
}


/**
Constructor. Sets the simple roots and the highest weight that defines the irreducible representation. Note that the Freudenthal algorithm not only gives the weights of a representation, but only their multiplicities.
*/
freudenthal::freudenthal(intMatrix simpleroots, rationalMatrix QFM, intVector highestroot1, 
			 intVector highestweight1) : rho(simpleroots.size(),1)
{
  rank = simpleroots.size();
  assert( (rank == highestroot1.size()) && (rank == highestweight1.size()) );

  highestroot = highestroot1;
  highestweight = highestweight1;
  
  dimension = 0;

  dynkin a(simpleroots, QFM, highestweight);
  weights = a.getWeights();

  dynkin b(simpleroots, QFM, highestroot);
  adjointweights = b.getWeights();

  A = a.getA();

  Q = inverse(A);

  setPosWeights();

  calculateMultiplicities();
}


/**
Calculates the multiplicities of all the weights of an irreducible representation given by the highest weight.
*/
void freudenthal::calculateMultiplicities()
{
  set<intVector>::iterator iter;
  for (iter = weights.begin(); iter!= weights.end(); ++iter) {
    dimension += multiplicity(*iter);
    for (int i=0; i<multiplicity(*iter); ++i) multiweights.insert(*iter);
  }

}


/**
Prints the weights of an irreducible representation and the corresponding multiplicities.
*/
void freudenthal::print() const
{
  cout << endl << "Multiplicities ...\n\n";
  
  int count = 0;
  map<intVector, int>::const_iterator pos;
  for (pos = weightmap.begin(); pos!= weightmap.end(); ++pos) {
    cout << "[ " << ++count << " ]" << "\t";
    cout.precision(3); cout.width(6);
    copy((pos->first).begin(),(pos->first).end(),std::ostream_iterator<int>(cout," "));
    cout << "\t" << pos->second << endl;
  }

  cout << endl << "Dimension = " << dimension << endl;

}


/**
Finds the positive weights of an irreducible representation.
*/
void freudenthal::setPosWeights()
{
  rationalVector x;
  rationalVector y(A.size());
  
  rational<int> zero(0);

  rationalMatrix temproots(A.size(), rationalVector(A.size()));
  for (unsigned i=0; i<A.size(); ++i) 
    for (unsigned j=0; j<A[0].size(); ++j) 
      temproots[i][j] = rational<int>(A[i][j]);
  
  rationalMatrix temprootstranspose = transpose(temproots);

  set<intVector>::iterator iter;
  rationalVector::iterator pos;
  for (iter = adjointweights.begin(); iter!= adjointweights.end(); ++iter) {
    for (unsigned i=0; i<(*iter).size(); ++i) y[i] = rational<int>((*iter)[i]);
    x = linSolve(temprootstranspose, y);
    if ( (find_if(x.begin(), x.end(), bind2nd(less<rational<int> >(),zero)) == x.end()) 
	 && (find_if(x.begin(), x.end(), bind2nd(greater<rational<int> >(),zero)) != x.end()))
      posweights.insert(*iter);
  }

}


/**
Gives the scalar product between two weights where the quadratic form matrix is the inverse of the Cartan matrix.
*/
rational<int> freudenthal::sp(intVector a, intVector b)
{
  assert((a.size() == b.size()) && (a.size() == Q.size()));

  rational<int> sum(0);

  for (unsigned i=0; i<a.size(); ++i) for (unsigned j=0; j<a.size(); ++j) sum += a[i]*Q[i][j]*b[j];

  return sum;

}


/**
Gives the scalar product between two weights where the quadratic form matrix is the inverse of the Cartan matrix.
*/
rational<int> freudenthal::sp(rationalVector a, rationalVector b)
{
  assert((a.size() == b.size()) && (a.size() == Q.size()));

  rational<int> sum(0);

  for (unsigned i=0; i<a.size(); ++i) for (unsigned j=0; j<a.size(); ++j) sum += a[i]*Q[i][j]*b[j];

  return sum;

}


/**
Gives the scalar product between two weights where the quadratic form matrix is the inverse of the Cartan matrix.
*/
rational<int> freudenthal::sp(intVector a, rationalVector b)
{
  assert((a.size() == b.size()) && (a.size() == Q.size()));

  rational<int> sum(0);

  for (unsigned i=0; i<a.size(); ++i) for (unsigned j=0; j<a.size(); ++j) sum += a[i]*Q[i][j]*b[j];

  return sum;

}


/**
Gives the scalar product between two weights where the quadratic form matrix is the inverse of the Cartan matrix.
*/
rational<int> freudenthal::sp(rationalVector a, intVector b)
{
  assert((a.size() == b.size()) && (a.size() == Q.size()));

  rational<int> sum(0);

  for (unsigned i=0; i<a.size(); ++i) for (unsigned j=0; j<a.size(); ++j) sum += a[i]*Q[i][j]*b[j];

  return sum;

}

/**
Gives the sum of two weights.
*/
intVector freudenthal::add(intVector a, intVector b)
{
  assert(a.size() == b.size());

  intVector sum(a.size());

  for (unsigned i=0; i<a.size(); ++i) sum[i] = a[i]+b[i];

  return sum;

}


/**
Multiplies the weight 'a' by the number 'k' (multiplication by a scalar).
*/
intVector freudenthal::mult(int k, const intVector a)
{
  intVector b(a.size());

  for (unsigned i=0; i<a.size(); ++i) b[i] = k*a[i];

  return b;

}


/**
Calculates the multiplicity of a weight 'w'.
*/
int freudenthal::multiplicity(intVector w)
{
  if (weightmap.find(w) != weightmap.end()) return weightmap[w];

  rational<int> mul(0);

  if (weights.find(w) == weights.end()) return 0;
  else if (w == highestweight) { weightmap[w] = 1; return 1; }
  else {

    set<intVector>::iterator iter;
    for (iter = posweights.begin(); iter!= posweights.end(); ++iter) {
      for (int k=1; weights.find(add(w,mult(k,*iter))) != weights.end(); ++k) {
	
	mul += sp(add(w,mult(k,*iter)),*iter)*multiplicity(add(w,mult(k,*iter)));


      }
    }

    rational<int> temp =  2*mul/( sp(add(highestweight,rho),add(highestweight,rho)) - sp(add(w,rho),add(w,rho)) );
    assert(temp.denominator() == 1);
    weightmap[w] = temp.numerator();
    return temp.numerator();
    
  }
 
}


/**
Returns a map that associates with each weight its multiplicity.
*/
map<intVector, int> freudenthal::getWeightMap()
{
  
  return weightmap;

}


/********************************** 
IMPLEMENTATION OF THE DYNKINDIAGRAM
**********************************/


/**
Constructor. Creates the Dynkin diagram from the simple roots.
*/
template <typename T>
dynkinDiagram<T>::dynkinDiagram(MatrixType & simpleroots)
  : rank(simpleroots.size()), A(rank,intVector(rank))
{
  MatrixType temp;
  for (unsigned i=0; i<rank; ++i) {
    temp.clear();
    for (unsigned j=0; j<rank; ++j) {
      A[i][j] = sp(simpleroots[i],simpleroots[j]);
      assert(abs(A[i][j])<=2);
      if ((i != j) && (A[i][j] != 0)) temp.push_back(simpleroots[j]);
    }
    adjacent[simpleroots[i]] = temp;
  }

  set<VectorType> allroots;
  typename MatrixType::iterator iter;
  for (iter = simpleroots.begin(); iter != simpleroots.end(); ++iter)
    allroots.insert(*iter);

  gaugeGroupFactor<T> factor;
  while (!allroots.empty()) {
    subAlgebra(allroots, factor.simpleroots, *(allroots.begin()));
    group.factor.push_back(factor);
    factor.simpleroots.clear();
  }

  simpleroots.clear();  // <-- changed recently
  for (unsigned i=0; i<group.factor.size(); ++i) {
    group.factor[i].algebra = determineAlgebra(group.factor[i].simpleroots);
    group.factor[i].rank = group.factor[i].simpleroots.size();
    reorderDiagram(group.factor[i].algebra, group.factor[i].simpleroots);
  }

  // Reorder algebras lexicographically
  sort(group.factor.begin(), group.factor.end());
  reverse(group.factor.begin(), group.factor.end());

  for (unsigned i=0; i<group.factor.size(); ++i) {
    copy(group.factor[i].simpleroots.begin(), group.factor[i].simpleroots.end(), // <-- changed recently
	 back_inserter(simpleroots));
    algebra += ( i == 0 ? (group.factor[i].algebra) :  (" + " + group.factor[i].algebra) );
  }

  group.algebra = algebra;

}


/**
Constructor. Creates the Dynkin diagram from the Cartan matrix.
*/
template <typename T>
dynkinDiagram<T>::dynkinDiagram(MatrixType & simpleroots, const MatrixType & cartanMatrix)
  : rank(simpleroots.size())
{
  A = to_integer<T>(cartanMatrix);

  MatrixType temp;
  for (unsigned i=0; i<rank; ++i) {
    temp.clear();
    for (unsigned j=0; j<rank; ++j) {
      assert(abs(A[i][j])<=2);
      if ((i != j) && (A[i][j]) != 0) temp.push_back(simpleroots[j]);
    }
    adjacent[simpleroots[i]] = temp;
  }
  
  set<VectorType> allroots;
  typename MatrixType::iterator iter;
  for (iter = simpleroots.begin(); iter != simpleroots.end(); ++iter)
    allroots.insert(*iter);

  gaugeGroupFactor<T> factor;
  while (!allroots.empty()) {
    subAlgebra(allroots, factor.simpleroots, *(allroots.begin()));
    group.factor.push_back(factor);
    factor.simpleroots.clear();
  }

  simpleroots.clear();  // <-- changed recently
  for (unsigned i=0; i<group.factor.size(); ++i) {
    group.factor[i].algebra = determineAlgebra(group.factor[i].simpleroots);
    group.factor[i].rank = group.factor[i].simpleroots.size();
    reorderDiagram(group.factor[i].algebra, group.factor[i].simpleroots);
    copy(group.factor[i].simpleroots.begin(), group.factor[i].simpleroots.end(), // <-- changed recently
	 back_inserter(simpleroots));
    algebra += ( i == 0 ? (group.factor[i].algebra) :  (" + " + group.factor[i].algebra) );
  }

  group.algebra = algebra;


}


/**
Reorders the roots of the algebra so that Dynkin diagram has canonical form.
*/
template <typename T>
void dynkinDiagram<T>::reorderDiagram(string alg, MatrixType & sroots)
{

  set<VectorType> candidates;

  typename MatrixType::iterator pos;

  typename MatrixType::iterator iter;
  for (iter = sroots.begin(); iter!= sroots.end(); ++iter) 
    candidates.insert(*iter);

  VectorType head = findHead(alg, sroots);

  iter = find(sroots.begin(), sroots.end(), head);
  swap(sroots[0],sroots[iter-sroots.begin()]);
  candidates.erase(head);
  if (sroots.size() == 1) return;
  else if ( (alg.substr(0,1)=="A") || (alg.substr(0,1)=="D") ) {
    
    for (iter = sroots.begin(); iter!= sroots.end(); ++iter)
      for (unsigned i=0; i<adjacent[*iter].size(); ++i)
	if (candidates.find(adjacent[*iter].at(i)) != candidates.end()) {
	  pos = find(sroots.begin(), sroots.end(), adjacent[*iter].at(i));
	  swap(sroots[iter+1-sroots.begin()],sroots[pos-sroots.begin()]);
	  candidates.erase(adjacent[*iter].at(i));
	}

  } else if (alg.substr(0,1)=="E") {

    typename MatrixType::iterator endpointer = sroots.end(); 
    for (unsigned i=1; i<sroots.size()-2; ++i)
      if (adjacent[sroots[i]].size() == 1) {
	pos = find(sroots.begin(), sroots.end(),sroots[i]);
	candidates.erase(sroots[pos-sroots.begin()]);
	if (adjacent[sroots[--endpointer-sroots.begin()]].size() != 1) 
	  swap(sroots[pos-sroots.begin()],sroots[endpointer-sroots.begin()]);
	else {
	  candidates.erase(sroots[endpointer-sroots.begin()]);
	  swap(sroots[pos-sroots.begin()],sroots[--endpointer-sroots.begin()]);
	}
      }

    for (unsigned i=sroots.size()-2; i<sroots.size(); ++i) 
      if (adjacent[sroots[i]].size() == 1) candidates.erase(sroots[i]);

    endpointer = sroots.end();
    for (unsigned i=0; i<2; ++i) endpointer--;
    for (iter = sroots.begin(); iter != endpointer; ++iter)
      for (unsigned i=0; i<adjacent[*iter].size(); ++i)
	if (candidates.find(adjacent[*iter].at(i)) != candidates.end()) {
	  pos = find(sroots.begin(), sroots.end(), adjacent[*iter].at(i));
	  swap(sroots[(iter-sroots.begin())+1],sroots[pos-sroots.begin()]);
	  candidates.erase(adjacent[*iter].at(i));
	}

    if (adjacent[adjacent[sroots[sroots.size()-1]].at(0)].size() != 3)
      swap(sroots[sroots.size()-1],sroots[sroots.size()-2]);

  }

}


/**
Finds a "corner" of the Dynkin diagram from where eventually the rest is built up by successive addition of roots.
*/
template <typename T>
typename dynkinDiagram<T>::VectorType dynkinDiagram<T>::findHead(string alg, 
								 const MatrixType & sroots)
{
  set<VectorType> candidates;
  typename set<VectorType>::iterator iter;
  for (unsigned i=0; i<sroots.size(); ++i)
    if (adjacent[sroots[i]].size() == 1) candidates.insert(sroots[i]);

  if (sroots.size() == 1) return sroots[0];
  else if (alg.substr(0,1)=="A") {
    assert(candidates.size() <= 2);
    return *(candidates.begin());
  } else if (alg.substr(0,1)=="D") {
    assert(candidates.size() == 3);
    if (sroots.size() == 4) return *(candidates.begin());
    else {
      for (iter = candidates.begin(); iter!= candidates.end(); ++iter)
	if (adjacent[adjacent[*iter].at(0)].size() == 2) return *iter;
    }
  } else if (alg.substr(0,1)=="E") {
    assert(candidates.size() == 3);
    VectorType nextnode;
    for (iter = candidates.begin(); iter!= candidates.end(); ++iter) {
      nextnode = adjacent[*iter].at(0);
      for (unsigned i=0; i<adjacent[nextnode].size(); ++i)
	if (adjacent[adjacent[nextnode].at(i)].size() == 3) return *iter;
    }
  }

  throw Error("findHead() failed ...");

  return *(candidates.begin());

}


/**
Collects roots of subalgebra that contains 'node'.
*/
template <typename T>
void dynkinDiagram<T>::subAlgebra(set<VectorType> & allroots, MatrixType & subalgebra, VectorType node)
{
  allroots.erase(node);
  subalgebra.push_back(node);

  typename MatrixType::iterator iter;
  for (iter = adjacent[node].begin(); iter != adjacent[node].end(); ++iter)
    if (allroots.find(*iter) != allroots.end()) subAlgebra(allroots, subalgebra, *iter);

}


/**
Returns the name of the algebra according to the ADE classification.
*/
template <typename T>
string dynkinDiagram<T>::determineAlgebra(MatrixType subalgebra)
{
  ostringstream oss;
  oss << subalgebra.size();

  set<VectorType> candidates;

  typename MatrixType::iterator iter;
  for (iter = subalgebra.begin(); iter != subalgebra.end(); ++iter)
    if (adjacent[*iter].size() == 1) candidates.insert(*iter);

  assert(candidates.size() <= 3);

  if (candidates.size() <= 2) return "A" + oss.str();

  typename set<VectorType>::iterator first = candidates.begin();
  typename set<VectorType>::iterator second = ++candidates.begin();
  typename set<VectorType>::iterator third = ++++candidates.begin();

  if ( ( adjacent[*first].at(0) == adjacent[*second].at(0) ) ||
       ( adjacent[*first].at(0) == adjacent[*third].at(0) ) ||
       ( adjacent[*second].at(0) == adjacent[*third].at(0) ) )
    return "D" + oss.str();
  else  return "E" + oss.str();
  

}


/**
Returns the scalar product between two vectors.
*/
template <typename T>
int dynkinDiagram<T>::sp(VectorType a, VectorType b)
{
  T sum(0); 

  assert(a.size() == b.size());
  
  for (unsigned i=0; i<a.size(); ++i) sum += a.at(i)*b.at(i);
  
  return to_integer(sum);
  
}


/************************************************** 
IMPLEMENTATION OF ROOT SYSTEMS FOR VARIOUS ALGEBRAS
**************************************************/

/**
Returns the simple roots of E8 in a standard basis.

      o      
      |      
o--o--o--o--o--o--o

*/

template <typename T>
vector<vector<T> > create_simpleroots_E8()
{
  typedef vector<T> VectorType;
  typedef vector<VectorType> MatrixType;

  T half = T(1)/T(2);
  T one(1);

  MatrixType roots(8,VectorType(8));
  
  roots[0].assign(8,-half);
  roots[0][0] = half; roots[0][7] = half;

  roots[1].assign(8,0);
  roots[1][6] = one; roots[1][7] = -one;

  roots[2].assign(8,0);
  roots[2][5] = one; roots[2][6] = -one;

  roots[3].assign(8,0);
  roots[3][4] = one; roots[3][5] = -one;

  roots[4].assign(8,0);
  roots[4][3] = one; roots[4][4] = -one;

  roots[5].assign(8,0);
  roots[5][2] = one; roots[5][3] = -one;

  roots[6].assign(8,0);
  roots[6][1] = one; roots[6][2] = -one;

  roots[7].assign(8,0);
  roots[7][6] = one; roots[7][7] = one;

  return roots;

}


/**
Returns all the roots of E8.
*/
template <typename T>
vector<vector<T> > create_roots_E8()
{
  typedef vector<T> VectorType;
  typedef vector<VectorType> MatrixType;

  const int n=8;

  T zero(0); 
  T half = T(1)/T(2);
  T one(1);

  MatrixType weights;
  VectorType lambda(n);

  for (int i=0; i<n; ++i) {
    for (int j=i+1; j<n; ++j) {
      lambda.assign(n,zero);
      lambda[i] = one;
      lambda[j] = one;
      weights.push_back(lambda);
      lambda[i] = -one;
      lambda[j] = one;
      weights.push_back(lambda);
      lambda[i] = one;
      lambda[j] = -one;
      weights.push_back(lambda);
      lambda[i] = -one;
      lambda[j] = -one;
      weights.push_back(lambda);
    }
  }

  lambda.assign(n,half);
  weights.push_back(lambda);
  for (int k=0; k<n; ++k) lambda[k] *= -one;
  weights.push_back(lambda);

  for (int i=0; i<n; ++i) {
    for (int j=i+1; j<n; ++j) {
      lambda.assign(n,half);
      lambda[i] *= -one;
      lambda[j] *= -one;
      weights.push_back(lambda);
      for (int k=0; k<n; ++k) lambda[k] *= -one;
      weights.push_back(lambda);
    }
  }
  
  for (int i=0; i<n; ++i) {
    for (int j=i+1; j<n; ++j) {
      for (int p=j+1; p<n; ++p) {
	for (int q=p+1; q<n; ++q) {
	  lambda.assign(n,half);
	  lambda[i] *= -one;
	  lambda[j] *= -one;
	  lambda[p] *= -one;
	  lambda[q] *= -one;
	  weights.push_back(lambda);
	}
      }
    }
  }
  
  return weights;
}



/**
Returns the simple roots of E8xE8.
*/
template <typename T>
vector<vector<T> > create_simpleroots_E8xE8()
{
  vector<vector<T> > roots(16, vector<T>(16));

  vector<vector<T> > E8roots = create_simpleroots_E8<T>();

  for (int i=0; i<E8roots.size(); ++i) {
    copy(E8roots[i].begin(), E8roots[i].end(), roots[i].begin());
    copy(E8roots[i].begin(), E8roots[i].end(), roots[i+8].begin()+8);
  }

  return roots;

}


/**
Returns a basis of Spin(32)/Z_2.
*/
template <typename T>
vector<vector<T> > create_basis_Spin32()
{
  vector<vector<T> > roots;

  vector<T> temp;

  temp = vector<T>(16, T(1)/T(2));
  for (int j=1; j<=14; ++j) temp[j] *= T(-1);
  roots.push_back(temp);

  for (int i=14; i>=1; --i) {
    temp = vector<T>(16);
    temp[i] = T(1);
    temp[i+1] = T(-1);
    roots.push_back(temp);
  }

  temp = vector<T>(16);
  temp[14] = T(1);
  temp[15] = T(1);
  roots.push_back(temp);

  return roots;
}


/************************************************** 
VARIOUS FUNCTIONS
**************************************************/


/**
Calculates the Cartan matrix from the simple roots.
*/
rationalMatrix calculateCartanMatrix(const rationalMatrix & sroots)
{
  unsigned n = sroots.size();

  rationalMatrix A(n, rationalVector(n));
  for (unsigned i=0; i<n; ++i)
    for (unsigned j=0; j<n; ++j)
      A[i][j] = 2*(sroots[i]*sroots[j])/(sroots[j]*sroots[j]);
    
  return A;

}


/**
Calculates the Cartan matrix from the simple roots.
*/
doubleMatrix calculateCartanMatrix(const doubleMatrix & sroots)
{
  unsigned n = sroots.size();

  doubleMatrix A(n, doubleVector(n));
  for (unsigned i=0; i<n; ++i)
    for (unsigned j=0; j<n; ++j)
      A[i][j] = 2*(sroots[i]*sroots[j])/(sroots[j]*sroots[j]);
    
  return A;

}


/**
Input is the matrix of simple roots, and a list of weights. Output is the list of weights in Dynkin labels.
*/
template <typename T>
intMatrix setDynkinLabels(vector< vector<T> > simpleroots, vector< vector<T> > weights)
{
  typedef vector<T> VectorType;
  typedef vector<VectorType> MatrixType;
  
  unsigned rank = simpleroots.size();
  intMatrix mapping;
  intVector x(rank);
  T q;

  typename MatrixType::const_iterator iter;
  for (iter = weights.begin(); iter != weights.end(); ++iter) {
    for (unsigned i=0; i<rank; ++i) {
      q = simpleroots[i]*(*iter);
      x[i] = to_integer(q);
     }
    mapping.push_back(x);
  }

  return mapping;

}


/**
Returns the highest weight of a weight system.
*/
intMatrix findHighestWeights(const multiset<intVector> & weights, intMatrix newroots, 
				       rationalMatrix Q, intMatrix & original_hw)
{

  multiset<intVector> newweights;
  intMatrix highestweights;

  rational<int> prod(0);
  intVector x(newroots.size());

  intVector zerovec((*(weights.begin())).size());

  multiset<intVector>::iterator iter;
  for (iter = weights.begin(); iter!= weights.end(); ++iter) {
    for (unsigned i=0; i < x.size(); ++i) {
      prod = to_rational(newroots[i])*(Q*to_rational(*iter));
      assert(prod.denominator() == 1);
      x[i] = prod.numerator();
    }
    
    if ( (find_if(x.begin(), x.end(), bind2nd(less<int>(),0)) == x.end()) ) {
      highestweights.push_back(x);
      original_hw.push_back(*iter);
    }
  }

  return highestweights;
}


/**
Input is a set of weights in Dynkin labels. Output are the highest weights. A weight is a (candidate) for a highest weight if all entries are nonnegative, and at least one entry is positive.
*/
multiset<intVector> findHighestWeights(const multiset<intVector> & weights)
{

  multiset<intVector> highestweights;

  multiset<intVector>::iterator iter;
  for (iter = weights.begin(); iter!= weights.end(); ++iter) {
    if ( find_if(iter->begin(), iter->end(), bind2nd(less<int>(),0)) == iter->end() )
      highestweights.insert(*iter);
  } 

  return highestweights;
  
}


/**
Returns the highest weight of a weight system.
*/
vector<vector<int> > findHighestWeights(const vector<vector<int> > & weights)
{

  vector<vector<int> > highestweights;

  vector<vector<int> >::const_iterator iter;
  for (iter = weights.begin(); iter!= weights.end(); ++iter) {
    if ( find_if(iter->begin(), iter->end(), bind2nd(less<int>(),0)) == iter->end() )
      highestweights.push_back(*iter);
  } 

  return highestweights;
  
}


/**
Returns the algebra corresponding to the weight system 'wsystem'.
*/
template <typename T>
gaugeGroup<T> determineAlgebra(const vector<vector<T> > & wsystem)
{
  if (wsystem.size() == 0) {
    gaugeGroup<T> tempgroup;
    tempgroup.algebra = "U(1)'s";
    return tempgroup;
  }

  weightsystem<T> ws(wsystem);

  dynkinDiagram<T> diagram(ws.simpleroots);

  diagram.setBasis(ws.basis);

  gaugeGroup<T> gruppe = diagram.getGaugeGroup();

  gruppe.set_semiordering(ws.get_semiordering());

  return gruppe;

}


/**
Returns the Dynkin labels of the vector 'mw'.
*/
template <typename T>
std::vector<int> findDynkinLabels(const gaugeGroupFactor<T> & myfactor, const std::vector<T> & mw)
{
  vector<int> result;
  
  T x = 0;
  for (unsigned i=0; i<myfactor.rank; ++i) {
    x = 2*(mw*myfactor.simpleroots[i])/(myfactor.simpleroots[i]*myfactor.simpleroots[i]);
    result.push_back(to_integer(x));
  }
  
  return result;
  
}


/**
Returns the Dynkin labels of the list of vectors 'mw'.
*/
template <typename T>
std::vector< std::vector<int> > findDynkinLabels(const gaugeGroupFactor<T> & myfactor, 
						 const std::vector< std::vector<T> > & mw)
{
  vector< vector<int> > result;
  vector<int> temp;

  T x = 0;

  for (unsigned i=0; i<mw.size(); ++i) {
    temp.clear();
    for (unsigned j=0; j<myfactor.rank; ++j) {
      x = 2*(mw[i]*myfactor.simpleroots[j])/(myfactor.simpleroots[j]*myfactor.simpleroots[j]);
      temp.push_back(to_integer(x));
    }
    result.push_back(temp);
  }
  
  return result;

}


/**
Returns the Dynkin labels of the vector 'mw'.
*/
template <typename T>
vector<int> findDynkinLabels(const vector<vector<T> > & myfactor, 
				      const vector<T> & mw)
{
  vector<int> result;
  
  T x = 0;
  for (unsigned i=0; i<myfactor.size(); ++i) {
    x = 2*(mw*myfactor[i])/(myfactor[i]*myfactor[i]);
    result.push_back(to_integer(x));
  }
  
  return result;
}


/**
Returns the Dynkin labels of the list of vectors 'mw'.
*/
template <typename T>
vector<vector<int> > findDynkinLabels(const vector<vector<T> > & myfactor, 
				      const vector<vector<T> > & mw)
{
  vector<vector<int> > result;
  
  T x = 0;
  for (unsigned k=0; k<mw.size(); ++k) {
    vector<int> temp(myfactor.size());
    for (unsigned i=0; i<myfactor.size(); ++i) {
      x = 2*(mw[k]*myfactor[i])/(myfactor[i]*myfactor[i]);
      temp[i] = to_integer(x);
    }
    result.push_back(temp);
  }

  return result;

}

/**
Returns the roots of an algebra. Note that the zero weight is not a root and is thus removed from the output.
*/
template <typename T>
set<vector<int> > dynkinAlgorithm(const gaugeGroupFactor<T> & myfactor, 
					     const vector<int> & lambda)
{
  vector<vector<T> > roots = myfactor.simpleroots;

  unsigned n = lambda.size();

  assert(n == roots.size());

  dynkin a(roots, lambda);

  return a.getWeights();

}


/**
Returns the highest weight in 'myWeights' with respect to the semi-ordering of 'myGaugeGroup'.
*/
template <typename T>
cirrep<T> findHighestWeight(gaugeGroup<T> myGaugeGroup, vector<vector<T> > myWeights)
{
  assert(myWeights.size() > 0);

  cirrep<T> irrep;

  vector<T> semiordering = myGaugeGroup.get_semiordering();

  vector<T> hw = *max_element(myWeights.begin(), myWeights.end(), 
			      smaller_wrt_root_ordering_crit<T>(semiordering));

  irrep.highestweight = hw;

  return irrep;
}


/**
Input is the algebra (A,B,C,D,E,G, or F), and the rank. Output is the highest weight of the algebra in Dynkin labels.
*/
intVector highestRoot(string algebra, unsigned rank)
{
  
  intVector hr(rank);
  
  if ((algebra == "A") && (rank == 1)) {
    hr[0] = 2; 
  } else if ((algebra == "A") && (rank>1) ) {
    hr[0] = 1; hr[rank-1] = 1;  
  } else if ((algebra == "B") && (rank > 1)) {
    hr[1] = 1;  
  } else if (algebra == "C") {
    hr[0] = 2;  
  } else if ((algebra == "D") && (rank > 1)) {
    hr[1] = 1;  
  } else if ((algebra == "E") && (rank == 6)) {
     hr[5] = 1;  
  } else if ((algebra == "E") && (rank == 7)) {
     hr[0] = 1;  
  } else if ((algebra == "E") && (rank == 8)) {
     hr[6] = 1;  
  } else if ((algebra == "G") && (rank == 2)) {
     hr[0] = 1;  
  } else if ((algebra == "F") && (rank == 4)) {
     hr[0] = 1;  
  } else throw Error("Error in function highestRoot: Algebra does not exist ...");

  return hr;

}


/**
Input is the set of simple roots in the Cartan-Weyl basis. Output are the Kac labels, i.e. the expansion coefficients of the highest root in terms of the simple roots. Note: The simple roots will be reordered to comply with Slansky.
*/
template <typename T>
vector<T> KacLabels(vector<vector<T> > & simpleroots)
{
  intMatrix cartanMatrix = setDynkinLabels(simpleroots, simpleroots);
  
  dynkinDiagram<T> diagram(simpleroots);
  
  assert(diagram.getGaugeGroup().factor.size() == 1); 

  string algebra = diagram.getGaugeGroup().factor.at(0).algebra.substr(0,1);
  unsigned rank = diagram.getGaugeGroup().factor.at(0).rank;

  intVector hr = highestRoot(algebra, rank);

  vector<T> kaclabel = linSolve_return<T>(transpose(cartanMatrix),hr);

  return kaclabel;

}


/**
Input is the set of simple roots in the Cartan-Weyl basis. Output are the Kac labels, i.e. the expansion coefficients of the highest root in terms of the simple roots. The 0th Kac label will be set to 1, corresponding to the highest root. Note: The simple roots will be reordered to comply with Slansky.
*/
template <typename T>
vector<T> extendedKacLabels(vector<vector<T> > & simpleroots)
{
  typedef vector<T> VectorType;

  VectorType kac1 = KacLabels(simpleroots);

  VectorType kac2(kac1.size()+1);
  
  kac2[0] = T(1);

  copy(kac1.begin(), kac1.end(), kac2.begin()+1);

  return kac2;

}


/**
Input is the gauge group data structure. Output are the Kac labels, i.e. the expansion coefficients of the highest root in terms of the simple roots. The 0th Kac label will be set to 1, corresponding to the highest root. Note: The simple roots will be reordered to comply with Slansky.
*/
template <typename T>
vector<T> extendedKacLabels(gaugeGroup<T> & group)
{
  vector<T> kac;

  for (unsigned i=0; i < group.factor.size(); ++i) {
    vector<T> tempkac = extendedKacLabels(group.factor[i].simpleroots);
    copy(tempkac.begin(), tempkac.end(), back_inserter(kac));
  }

  return kac;
}


/**
Prints the human-readable name of the algebra to standard output.
*/
template <typename T>
ostream & operator<<(ostream & os, const class gaugeGroup<T> & A)
{
  cout << "Algebra is " << A.algebra << endl << endl;

  for (unsigned i=0; i<A.factor.size(); ++i) {
    
    cout << A.factor[i].algebra << endl
	 << A.factor[i].simpleroots << endl;
  }

  return os;
}


/********************** 
EXPLICIT INSTANTIATIONS
**********************/


template
vector<vector<double> > create_simpleroots_E8();

template
vector<vector<rational<int> > > create_simpleroots_E8();

template
vector<vector<double> > create_roots_E8();

template
vector<vector<rational<int> > > create_roots_E8();

template
vector<vector<rational<int> > > create_simpleroots_E8xE8();

template
vector<vector<double> > create_simpleroots_E8xE8();

template
class weightsystem<double>;

template
class weightsystem<rational<int> >;

template
class gaugeGroupFactor<double>;

template
class gaugeGroupFactor<rational<int> >;

template
class gaugeGroup<double>;

template
class gaugeGroup<rational<int> >;

template
class dynkinDiagram<double>;

template
class dynkinDiagram<rational<int> >;

template
gaugeGroup<double> determineAlgebra(const vector<vector<double> > & wsystem);

template
gaugeGroup<rational<int> > determineAlgebra(const vector<vector<rational<int> > > & wsystem);

template
vector<int> findDynkinLabels(const gaugeGroupFactor<double> & myfactor, 
			     const std::vector<double> & mw);

template
vector< std::vector<int> > findDynkinLabels(const gaugeGroupFactor<double> & myfactor, 
					    const vector< std::vector<double> > & mw);

template
vector<int> findDynkinLabels(const gaugeGroupFactor<rational<int> > & myfactor, 
			     const std::vector<rational<int> > & mw);

template
vector< std::vector<int> > findDynkinLabels(const gaugeGroupFactor<rational<int> > & myfactor, 
					    const vector< std::vector<rational<int> > > & mw);

template
vector<int> findDynkinLabels(const vector<vector<double> > & myfactor, 
					    const vector<double> & mw);

template
vector<int> findDynkinLabels(const vector<vector<rational<int> > > & myfactor, 
					    const vector<rational<int> > & mw);

template
vector<vector<int> > findDynkinLabels(const vector<vector<rational<int> > > & myfactor, 
				      const vector<vector<rational<int> > > & mw);
template
vector<vector<int> > findDynkinLabels(const vector<vector<double> > & myfactor, 
				      const vector<vector<double> > & mw);

template
set<vector<int> > dynkinAlgorithm(const gaugeGroupFactor<double> & myfactor, 
				  const vector<int> & lambda);

template
cirrep<double> findHighestWeight(const gaugeGroup<double> myGaugeGroup, 
				 vector<vector<double> > myWeights);

template
intMatrix setDynkinLabels(vector< vector<double> > simpleroots, vector< vector<double> > weights);

template
intMatrix setDynkinLabels(vector< vector<rational<int> > > simpleroots, vector< vector<rational<int> > > weights);

template
vector<double> KacLabels(vector<vector<double > > & simpleroots);

template
vector<rational<int> > KacLabels(vector<vector<rational<int> > > & simpleroots);

template
vector<double> extendedKacLabels(vector<vector<double > > & simpleroots);

template
vector<rational<int> > extendedKacLabels(gaugeGroup<rational<int> > & group);

template
vector<double> extendedKacLabels(gaugeGroup<double> & group);

template
vector<rational<int> > extendedKacLabels(vector<vector<rational<int> > > & simpleroots);

template
ostream & operator<<(ostream & os, const class gaugeGroup<rational<int> > & A);

template
ostream & operator<<(ostream & os, const class gaugeGroup<double> & A);

template
vector<vector<double> > create_basis_Spin32();

template
vector<vector<Rational> > create_basis_Spin32();

