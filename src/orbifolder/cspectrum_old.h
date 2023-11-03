#ifndef CSPECTRUM_H
#define CSPECTRUM_H

#include <vector>
#include "corbifold.h"

using std::vector;

class CSpectrum {
public:
// member functions
  CSpectrum();
  CSpectrum(const SConfig &Vacuum, SUSYMultiplet Multiplet);
  CSpectrum(const SConfig &Vacuum);
  CSpectrum(const SConfig &Vacuum1, const SConfig &Vacuum2);
  ~CSpectrum();

  bool                        operator==(const CSpectrum &Spectrum2) const;
  bool                        GaugeGroupFromBothE8Factors() const;
  unsigned                    CompareE8Factors() const;
  bool                        InterchangeE8xE8(CSpectrum &result) const;
  bool                        IsSpectrumEmpty() const;
  bool 						  IsTachyonFree() const;
  void                        AddIdentifier(const double &NewIdentifier) {this->AdditionalIdentifier.push_back(NewIdentifier);};
  
  const vector<double>       &GetAdditionalIdentifier() const {return this->AdditionalIdentifier;};
  const SUSYMultiplet        &GetSUSYMultiplet() const {return this->Multiplet;};
  const vector<unsigned>     &GetMultiplicity() const {return this->Multiplicity;};
  const vector<vector<unsigned> >     &GetMultiplicity_Matrix() const {return this->Multiplicity_Matrix;};
  const vector<vector<int> > &GetRepresentations() const {return this->Representations;};
  const vector<vector<vector<int> > > &GetRepresentations_Matrix() const {return this->Representations_Matrix;};
  const vector<vector<int> > &GetAlgebra() const {return this->Algebra;};
  const unsigned             &GetPosition_of_and_in_GaugeGroup() const {return this->Position_of_and_in_GaugeGroup;};
  
private:
// member variables
  vector<double>              AdditionalIdentifier;

  SUSYMultiplet               Multiplet;
  vector<unsigned>            Multiplicity;
  vector<vector<int> >        Representations;
  //vector<vector<int> >        HighestWeights_DL;
  
  //vector<SUSYMultiplet>               Multiplet_Matrix;
  vector<vector<unsigned> >           Multiplicity_Matrix;
  vector<vector<vector<int> > >       Representations_Matrix;

  vector<vector<int> >        Algebra;
  unsigned                    Position_of_and_in_GaugeGroup;
  
// member functions
  bool                        ReOrder();
};

#endif
