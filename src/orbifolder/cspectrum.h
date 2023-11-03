#ifndef CSPECTRUM_H
#define CSPECTRUM_H

#include <vector>
#include <list>
#include "corbifold.h"

using std::vector;
using std::list;


#ifndef STRUCT_SREPRESENTATION
#define STRUCT_SREPRESENTATION
struct SRepresentation
{
  unsigned      Multiplicity;
  vector<int>   Representation;
  double        U1Charge;
  
  unsigned      TotalMultiplicity;
  SUSYMultiplet Multiplet;
};
#endif


#ifndef STRUCT_SSPECTRUM
#define STRUCT_SSPECTRUM
struct SSpectrum
{
  vector<unsigned>        GaugeGroup;
  vector<unsigned>        Rank;
  vector<SRepresentation> Spectrum;
};
#endif


#ifndef STRUCT_SREPCOUNT
#define STRUCT_SREPCOUNT
struct SDimCount
{
  unsigned Dim;
  unsigned Counter;
};
#endif


#ifndef STRUCT_SCOLUMNOFSPECTRUM
#define STRUCT_SCOLUMNOFSPECTRUM
struct SColumnOfSpectrum
{
  unsigned    GaugeGroup;
  unsigned    Rank;
  vector<int> Representations;
  
  unsigned          Sum;
  unsigned          Product;
  unsigned          NumberOfOnes;

  rational<int>     AdditionalU1Sum;
  unsigned          AdditionalMultSum;
  vector<SDimCount> NumberOfReps;
  int               AbsNetNumberOfNForSUN;  
};
#endif


class CSpectrum {
public:
// member functions
  CSpectrum();
  CSpectrum(const SConfig &Vacuum, const vector<SUSYMultiplet> &SUSYTypes, string Label = "", const vector<vector<unsigned> > *SectorsForAdditionalIdentifier = NULL);
  ~CSpectrum();

  bool                        operator==(const CSpectrum &Spectrum2) const;
  bool                        IsSpectrumEmpty() const;
  void                        AddIdentifier(const double &NewIdentifier) {this->AdditionalIdentifier.push_back(NewIdentifier);};
  
  const vector<double>       &GetAdditionalIdentifier() const {return this->AdditionalIdentifier;};

  const unsigned             &GetPosition_of_and_in_GaugeGroup() const {return this->Position_of_and_in_GaugeGroup;};
  const bool                  GetUseSpecialU1() const {return this->UseSpecialU1;};
  
  size_t                      GetNumberOfRepresentations() const {return this->Spectrum.Spectrum.size();};
  unsigned                    GetLargestMultiplicity() const;
  size_t                      GetNumberOfGaugeGroups() const {return this->Spectrum.GaugeGroup.size();};
  
  const SSpectrum            &GetSpectrum() const {return this->Spectrum;};
  const SSpectrum            &GetSortedSpectrum(unsigned i) const;
  size_t                      GetNumberOfSortedSpectra() const {return this->SortedSpectra.size();};
  
  string Label;
private:
// member variables
  unsigned                    NumberOfSupersymmetry;

  bool                        UseSpecialU1;
  unsigned                    Position_of_and_in_GaugeGroup;

  SSpectrum                   Spectrum;
  vector<SSpectrum>           SortedSpectra;
  
  vector<double>              AdditionalIdentifier;
  
// member functions
  bool                        ReOrder();
  unsigned                    OrderE8xE8(vector<SColumnOfSpectrum> &ColumnsPart1, vector<SColumnOfSpectrum> &ColumnsPart2);
  bool                        SSpectrumKnown(const SSpectrum &Spectrum) const;
  bool                        InsertSortedColumnsIntoExistingSpectrum(SSpectrum &Spectrum, const vector<vector<SColumnOfSpectrum> > &Columns) const;
};

#endif
