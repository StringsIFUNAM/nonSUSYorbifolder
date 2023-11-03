#ifndef CFIXEDPOINT_H
#define CFIXEDPOINT_H

#include "corbifoldgroupelement.h"
#include "coscillator.h"
#include "cstate.h"
#include "chalfstate.h"
#include "groupTheory.hpp"

#include <string>

#ifndef ENUM_CHECKSTATUS
#define ENUM_CHECKSTATUS
enum CheckStatus {UNSPECIFIED_CHECKSTATUS = 0, NotChecked, CheckedAndGood, CheckedAndFailed};
#endif

using std::string;

class COrbifold;
class CState;
class CSector;
class COrbifoldGroup;

//! CFixedBrane.
/*!
A CFixedBrane object contains: the constructing space-group element "SGElement", 
the solutions to the equation for massless left-movers ("MasslessLeftMovers") involving the local shift V_loc corresponding to "SGElement",
the massless left-movers ("LeftMovers") sorted by their eigenvalues with respect to the centralizer elements and finally the centralizer-invariant 
combinations of left- and right-movers ("InvariantStates").
 */

class CFixedBrane {
public:
// member functions
  CFixedBrane();
  CFixedBrane(const CSpaceGroupElement &SGElement, const unsigned &Index_SGElement, const string &FixedBraneLabel);
  ~CFixedBrane();

  bool                           Create(const COrbifoldGroup &OrbifoldGroup, const CSector &Sector, const vector<COrbifoldGroupElement> &Gamma_Centralizer);
  bool                           Reset();

  bool                           GetFieldIndices(const vector<CField> &Fields, const SUSYMultiplet &Multiplet, vector<unsigned> &FieldIndices) const;
  bool                           TachyonicGetFieldIndices(const vector<CField> &Fields, const SUSYMultiplet &Multiplet, vector<unsigned> &FieldIndices) const;
  
  size_t                         GetNumberOfMasslessLeftMovers() const { return this->MasslessLeftMovers.size();};
  size_t                         GetNumberOfLeftMovers() const { return this->LeftMovers.size();};
  size_t                         GetNumberOfInvariantStates() const { return this->InvariantStates.size();};
  size_t                         TachyonicGetNumberOfInvariantStates() const { return this->TachyonicInvariantStates.size();};
  size_t                         ExcitedTachyonicGetNumberOfInvariantStates() const { return this->ExcitedTachyonicInvariantStates.size();};

  CState                        &AccessInvariantState(const unsigned &i);
  CState                        &TachyonicAccessInvariantState(const unsigned &i);
  CState                        &ExcitedTachyonicAccessInvariantState(const unsigned &i);
  string                        &AccessFixedBraneLabel() {return this->FixedBraneLabel;};
  
  const CSpaceGroupElement      &GetSGElement() const;
  const unsigned                &Getconstructing_Element() const;
  const string                  &GetFixedBraneLabel() const;
  const CMasslessHalfState      &GetMasslessLeftMover(const unsigned &i) const;
  const CMasslessHalfState      &TachyonicGetMassLeftMover(const unsigned &i) const;
  const CMasslessHalfState      &ExcitedTachyonicGetMassLeftMover(const unsigned &i) const;
  const CHalfState              &GetLeftMover(const unsigned &i) const;
  const CState                  &GetInvariantState(const unsigned &i) const;
  const CState                  &TachyonicGetInvariantState(const unsigned &i) const;
  const CState                  &ExcitedTachyonicGetInvariantState(const unsigned &i) const;

private:
  bool                           CreateMasslessLeftMover(const COrbifoldGroupElement &constructing_Element, const vector<S_OscillatorExcitation> &Excitations);
  bool                           TachyonicCreateLeftMover(const COrbifoldGroupElement &constructing_Element, const vector<vector<S_OscillatorExcitation> > &Excitations);
//  bool 							 CreateTachyonicLeftMover(const unsigned &windex, const COrbifoldGroup &OrbifoldGroup, const COrbifoldGroupElement &constructing_Element, const vector<vector<S_OscillatorExcitation> > &Excitations, const vector<CModedOscillator> &all_tachyonicOscillators, const CTachyonHalfState &TachyonicRightMovers);
  bool                           CreateStates(const vector<CHalfState> &RightMovers, const vector<COrbifoldGroupElement> &Centralizer, const vector<COrbifoldGroupElement> &Gamma_Centralizer);
  bool                           TachyonicCreateStates(const vector<CHalfState> &RightMovers, const vector<COrbifoldGroupElement> &Centralizer, const vector<COrbifoldGroupElement> &Gamma_Centralizer);
  bool                           FindSUSYMultiplets(const CSector &Sector, const vector<CVector> &InvariantSupercharges);
  bool                           TachyonicFindSUSYMultiplets(const CSector &Sector, const vector<CVector> &InvariantSupercharges);
  bool                           SortByEigenvalue(const COrbifoldGroup &OrbifoldGroup, const vector<CModedOscillator> &all_Oscillators, const vector<COrbifoldGroupElement> &Gamma_Centralizer);
  bool                           TachyonicSortByEigenvalue(const COrbifoldGroup &OrbifoldGroup, const vector<CModedOscillator> &all_Oscillators, const vector<COrbifoldGroupElement> &Gamma_Centralizer);

  // member variables
  CSpaceGroupElement             SGElement;
  unsigned                       Index_SGElement;
  string                         FixedBraneLabel;
  vector<CMasslessHalfState>     MasslessLeftMovers;
  vector<CMasslessHalfState>     TachyonicMassLeftMovers;				//corresponds to lowest R-moving tachyon
  vector<CMasslessHalfState>     ExcitedTachyonicMassLeftMovers;		//corresponds to possible excited R-moving tachyon

  vector<CHalfState>             LeftMovers;
  vector<CHalfState>             TachyonicLeftMovers;					//refers to TachyonicMassLeftMovers
  vector<CHalfState>             ExcitedTachyonicLeftMovers;			//refers to ExcitedTachyonicMassLeftMovers
  vector<CState>                 InvariantStates;
  vector<CState>                 TachyonicInvariantStates;
  vector<CState>                 ExcitedTachyonicInvariantStates;

  CheckStatus                    FixedBrane_CheckStatus;
};

#endif
