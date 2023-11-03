
#ifndef CTWIST_H
#define CTWIST_H

#include "cvector.h"

//! CTwistVector.
/*!
  A CTwistVector object v is a 4 dim. vector of phases rotating simultaneously in the 4 (complex) space coordinates (in light cone gauge).
 */

class CTwistVector : public CVector{
public: 
// member functions
  CTwistVector();
  CTwistVector(double v_0, double v_1, double v_2, double v_3);
  ~CTwistVector();

  void            operator=(const CVector &Vector);
  void            operator+=(const CVector &Vector2);

  const double   &Get_a_R() const {return this->a_R;};
  const double   &Get_a_L() const {return this->a_L;};
  const unsigned &OrderOfTwist() const {return this->Order;};
  void            SetToZero();

private:
  bool     UpdateData();

// member variables
  double   a_L;
  double   a_R;
  unsigned Order;
};

#endif
