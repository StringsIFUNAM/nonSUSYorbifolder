
#ifndef CFIELD_H
#define CFIELD_H
#include <vector>

#include "corbifold.h"

//! CField.
/*!
A CField object contains all information to describe a massless field with non-Abelian representation, U(1) charges, discrete charges, a field label and potentially a vev.
 */

using std::vector;

class CSector;
class CFixedBrane;
class CState;
class CRepresentation;
class CModedOscillator;

class CField{
public: 
  CField();
  ~CField();

  rational<int>               GetDiscreteCharge(const SDiscreteSymmetry &DiscreteSymmetry) const;
  rational<int>               GetModularWeight(const SModularSymmetry &ModularSymmetry, const vector<CSector> &Sectors) const;

  void                        SetInternalIndex(const unsigned &i, const unsigned &j, const unsigned &k);

  size_t                      GetNumberOfOscillators(const vector<CSector> &Sectors) const;
  size_t                      GetNumberOfLMWeights() const { return this->WeightIndices.size();};
  size_t                      GetNumberOfRMWeights(const vector<CSector> &Sectors) const;

  const CModedOscillator     &GetOscillator(const unsigned &i, const vector<CSector> &Sectors) const;
  const CVector              &GetLMWeight(const unsigned &i, const vector<CSector> &Sectors) const;
  const CVector              &GetRMWeight(const unsigned &i, const vector<CSector> &Sectors) const;

  const CSector              &GetSector(const vector<CSector> &Sectors) const;
  const CFixedBrane          &GetFixedBrane(const vector<CSector> &Sectors) const;
  const CState               &GetState(const vector<CSector> &Sectors) const;

  const vector<unsigned>     &GetInternalIndex() const {return this->internalIndex;};

  SUSYMultiplet               Multiplet;
  RepVector                   Dimensions;
  CVector                     U1Charges;

  vector<unsigned>            WeightIndices;
  CVector                     q_sh;
  CVector                     OsciContribution;
  CSpaceGroupElement          SGElement;

  vector<string>              Labels;
  vector<unsigned>            Numbers;

  CVector                     VEVs;

  vector<vector<int> >       HighestWeights_DL;

  double                      BmLCharge;
  vector<double>              gamma_phases;

  vector<rational<CHugeInt> > AccU1Charges;

  int                         OriginStandardConfig;

private:
  vector<unsigned>            internalIndex;
};

#endif

