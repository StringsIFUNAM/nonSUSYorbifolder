
#ifndef CWILSONLINE_H
#define CWILSONLINE_H

#include "cvector.h"

//! CWilsonLine.
/*!
  A CWilsonLine object W is a 16 dim. vector acting on the 16 gauge degrees of freedom. It is associated to a torus translation.
 */

class CWilsonLine : public CVector{
public:
// member functions
  CWilsonLine();
  CWilsonLine(const SelfDualLattice &Lattice);
  CWilsonLine(SelfDualLattice Lattice, double A_0, double A_1, double  A_2, double  A_3, double  A_4, double  A_5, double  A_6, double  A_7,
              double A_8, double A_9, double A_10, double A_11, double A_12, double A_13, double A_14, double A_15);
  ~CWilsonLine();

  void                   operator=(const CVector &Vector);
  void                   operator=(const vector<rational<int> > &RationalVector);

  bool                   CreateRandom(unsigned WLOrder, const vector<unsigned> &SymmetryBlocks, const vector<vector<vector<double> > > &VectorialBlocks, const vector<vector<vector<double> > > &SpinorialBlocks);
  bool                   CreateRandomBrother(const vector<CVector> &LatticeVectors);
  void                   SetLattice(SelfDualLattice Lattice) {this->Lattice = Lattice;};
  void                   SetToZero();

  const unsigned        &GetOrder() const {return this->Order;};
  const bool            &GetIs_Zero() const {return this->Is_Zero;};
  const SelfDualLattice &GetLattice() const {return this->Lattice;};

private:
// member functions
  bool                   UpdateData();

// member variables
  unsigned               Order;
  SelfDualLattice        Lattice;
  bool                   Is_Zero;
};

#endif
