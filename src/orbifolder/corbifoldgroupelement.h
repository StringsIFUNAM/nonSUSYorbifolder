
#ifndef CORBIFOLDGROUPELEMENT_H
#define CORBIFOLDGROUPELEMENT_H

#include "cspacegroupelement.h"
#include "cshiftvector.h"
#include "ctwist.h"

using std::vector;

class COrbifoldGroupElement {
public: 
  COrbifoldGroupElement();
  COrbifoldGroupElement(SelfDualLattice Lattice);
  COrbifoldGroupElement(const CSpaceGroupElement &SGElement, const CShiftVector &Shift, const CTwistVector &Twist);
  ~COrbifoldGroupElement();

  CSpaceGroupElement SGElement;
  CShiftVector       Shift;
  CTwistVector       Twist;
};

#endif
