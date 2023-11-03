#ifndef COSCILLATOR_H
#define COSCILLATOR_H

#include <iostream>

#include "ctwist.h"

//! CModedOscillator.
/*!
A CModedOscillator object describes a left- or right-moving oscillator excitation.
 */

#ifndef ENUM_MOVERSTYPE
#define ENUM_MOVERSTYPE
enum MoversType {NOT_DEF_MOVTYPE = 0, LeftMover, RightMover};
#endif

#ifndef ENUM_LADDEROPERATORTYPE
#define ENUM_LADDEROPERATORTYPE
enum LadderOperatorType {NOT_DEF_LAD=0, CreationOperator, AnnihilationOperator};
#endif

#ifndef ENUM_CHECKSTATUS
#define ENUM_CHECKSTATUS
enum CheckStatus {UNSPECIFIED_CHECKSTATUS = 0, NotChecked, CheckedAndGood, CheckedAndFailed};
#endif

class CModedOscillator {
public:
// member functions
  CModedOscillator(MoversType type, unsigned index, bool indexbar, const CTwistVector &constructing_Twist, int m);
  ~CModedOscillator();

  bool                     operator==(const CModedOscillator &Oscillator2) const;

  const bool               &GetComplex() const;
  const double             &GetFrequency() const;
  const unsigned           &GetIndex() const;
  double                   GetTransformationProperty(const CTwistVector &Twist) const;
  const MoversType         &GetType() const;
  const double             &GetNumberOperator() const;
  const LadderOperatorType &GetLadderOperatorType() const;

  const CheckStatus        &GetOscillator_CheckStatus() const {return this->Oscillator_CheckStatus;};
  
private:
// member variables
  LadderOperatorType LadderOperator;
  double             NumberOperator;
  bool               indexbar;
  double             frequency;
  unsigned           index;
  MoversType         type;

  CheckStatus        Oscillator_CheckStatus;
};

#endif
