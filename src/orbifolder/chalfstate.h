
#ifndef CHALFTSTATE_H
#define CHALFTSTATE_H

#include <vector>
#include <string>

#include "coscillator.h"

class COrbifoldGroupElement;

using std::vector;

#ifndef STRUCT_S_OSCILLATOREXCITATION
#define STRUCT_S_OSCILLATOREXCITATION
struct S_OscillatorExcitation
{
  double           NumberOperator;
  vector<unsigned> OscillatorIndices;
  double           ZeroPointEnergy;
};
#endif

//! CMasslessHalfState.
/*!
A CMasslessHalfState object contains the solutions for massless left- or right-movers (depending on "MoversType Type;") with possible 
oscillator excitations specified by "S_OscillatorExcitation Excitation".
 */

class CMasslessHalfState {
public:
// member functions
  CMasslessHalfState(MoversType Type, const S_OscillatorExcitation &Excitation);
  ~CMasslessHalfState();

  bool SolveMassEquation(const CVector &constructing_Element, const SelfDualLattice &Lattice);

// member variables
  S_OscillatorExcitation Excitation;
  MoversType             Type;
  vector<CVector>        Weights;
};

class CTachyonHalfState {
public:
// member functions
	CTachyonHalfState(MoversType Type, const S_OscillatorExcitation &Excitation);
    ~CTachyonHalfState();

  bool SolveMassEquation(const CVector &constructing_Element, const SelfDualLattice &Lattice);
  bool TachyonSolver(const unsigned &nindex, const unsigned &kindex, const CVector &constructing_Element, const SelfDualLattice &Lattice);

// member variables
  //S_OscillatorExcitation Excitation;
  S_OscillatorExcitation Excitation;
  MoversType             Type;
  vector<CVector>        Weights;
  vector<double>		 tachyonmass;
};


//! CHalfState.
/*!
A CHalfState object descends from a CMasslessHalfState object by sorting the massless solutions 
(i.e. the weights p_sh or q_sh for left- or right-movers) with respect to their centralizer-eigenvalues.
 */

class CHalfState {
public:
// member functions
  CHalfState();
  CHalfState(MoversType Type, unsigned Index);
  ~CHalfState();

  const MoversType &GetType() const {return this->Type;};
  const unsigned   &GetIndex() const {return this->Index;};

// member variables
  vector<double>    Eigenvalues;
  vector<unsigned>  Weights;
  unsigned 		Excited;		// 0 for ground state, 1,2,.. if right-moving oscillators acting, relevant for double-tachyonic and non-susy geometries

private:
  unsigned          Index; 		// refers to the corresponding CMasslessHalfState object
  MoversType        Type;
};

#endif
