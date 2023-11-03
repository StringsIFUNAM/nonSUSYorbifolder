
#ifndef CSTATE_H
#define CSTATE_H
#include <vector>

#include "chalfstate.h"
#include "globalfunctions.h"
#include "groupTheory.hpp"

//! CState.
/*!
A CState object describes a string state consisting of two CHalfState objects (left- and right-movers).
 */

#ifndef ENUM_CHECKSTATUS
#define ENUM_CHECKSTATUS
enum CheckStatus {UNSPECIFIED_CHECKSTATUS = 0, NotChecked, CheckedAndGood, CheckedAndFailed};
#endif

#ifndef ENUM_SUSYMULTIPLET
#define ENUM_SUSYMULTIPLET
enum SUSYMultiplet {NOT_DEF_SUSY = 0, AnyKind, LeftChiral, RightChiral, Vector, VectorCC, Hyper, Halfhyper, Gravity, GravityCC, LCModulus, RCModulus, BOSE_mu, Gauge, bGauge, Scalar, bScalar, moduli, bmoduli, LeftFermi, RightFermi, Tachyon, bTachyon, Excited_Tachyon, Excited_bTachyon};
#endif

using std::vector;

class  CSector;
class  CFixedBrane;
struct SConfig;
struct SDimension;

class CState {
public: 
// member functions
  CState();
  ~CState();

  bool                            Create(const CHalfState &LeftMover, const CHalfState &RightMover, const vector<COrbifoldGroupElement> &Gamma_Centralizer);
  bool                            CreateRepresentations(const CFixedBrane &FixedBrane, SConfig &Vacuum, vector<unsigned> &FieldCounters);
  bool                            TachyonicCreateRepresentations(const CFixedBrane &FixedBrane, SConfig &Vacuum, vector<unsigned> &FieldCounters);
  bool                            ExcitedTachyonicCreateRepresentations(const CFixedBrane &FixedBrane, SConfig &Vacuum, vector<unsigned> &FieldCounters);
  bool                            FindSUSYMultiplets(const CSector &Sector, const vector<CVector> &InvariantSupercharges, const vector<vector<unsigned> > &AllSUSYChargesCombinations);
  bool                            TachyonicFindSUSYMultiplets(const CSector &Sector, const vector<CVector> &InvariantSupercharges, const vector<vector<unsigned> > &AllSUSYChargesCombinations);
  void                            SetInternalIndex(const unsigned &i, const unsigned &j, const unsigned &k);

  bool                            ContainsSUSYMultiplet(const vector<CField> &Fields, const SUSYMultiplet &Multiplet) const;
  bool                            GetFieldIndices(const vector<CField> &Fields, const SUSYMultiplet &Multiplet, vector<unsigned> &FieldIndices) const;

  bool                            RecalculateU1Charges(const vector<CVector> &AllWeightsFromFixedBrane, SConfig &Vacuum) const;

  size_t                          GetNumberOfGammaPhases() const { return this->gamma_phases.size();};
  size_t                          GetNumberOfFieldIndices() const { return this->FieldIndices.size();};

  void                            SetOsciContribution(const CVector &OsciContribution) {this->OsciContribution = OsciContribution;};

  const CHalfState               &GetLeftMover() const;
  const CHalfState               &GetRightMover() const;
  
  const double                   &GetGammaPhase(const unsigned &i) const;
  const unsigned                 &GetFieldIndex(const unsigned &i) const;
  const CVector                  &GetOsciContribution() const {return this->OsciContribution;};
  const vector<vector<CVector> > &GetqCharges() const {return this->qCharges;};

  bool                            DetermineDimension(const gaugeGroupFactor<double> &ggf, const intVector &HighestWeight_DL, SDimension &DimOfRep) const;
  
private:
// member variables
  CHalfState                      LeftMover;
  CHalfState                      RightMover;

  CVector                         OsciContribution;
  vector<vector<CVector> >        qCharges;
  vector<SUSYMultiplet>           Multiplets;

  vector<unsigned>                FieldIndices;

  vector<double>                  gamma_phases;

  CheckStatus                     State_CheckStatus;
  
  vector<unsigned>                internalIndex;
};

#endif

