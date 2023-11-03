
#ifndef CSHIFTVECTOR_H
#define CSHIFTVECTOR_H

#include "cvector.h"

//! CShiftVector.
/*!
  A CShiftVector object V is a 16 dim. vector acting on the 16 gauge degrees of freedom. It is associated to the rotational transformation of the twist.
 */

using std::ifstream;

class CShiftVector : public CVector{
public:
// member functions
  CShiftVector();
  CShiftVector(SelfDualLattice Lattice);
  CShiftVector(SelfDualLattice Lattice, double V_0, double V_1, double  V_2, double  V_3, double  V_4, double  V_5, double  V_6, double  V_7,
               double V_8, double V_9, double V_10, double V_11, double V_12, double V_13, double V_14, double V_15);
  ~CShiftVector();

  void            operator=(const CVector &Vector);
  void            operator=(const vector<rational<int> > &RationalVector);

  SelfDualLattice GetLattice() const {return this->Lattice;};
  bool            CreateRandom(unsigned TwistOrder, double TwistSqr, const vector<unsigned> &SymmetryBlocks, const vector<vector<vector<double> > > &VectorialBlocks, const vector<vector<vector<double> > > &SpinorialBlocks);
  bool            CreateRandomBrother(const vector<CVector> &LatticeVectors, const double &TwistSqr, const unsigned TwistOrder);
  bool            LoadShiftVector(SelfDualLattice Lattice, ifstream &in);
  unsigned        OrderOfShift() const {return this->Order;};
  void            SetToZero();

// member variables
  SelfDualLattice Lattice;
  unsigned        Order;
private:
// member functions
  bool            UpdateData();
};

#endif
