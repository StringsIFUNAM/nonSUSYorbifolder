#include <iostream>
#include <iomanip>
#include <cstdlib>

#include <iterator>
#include <algorithm>

#include "cspectrum.h"
#include "corbifold.h"

#include "cprint.h"

#define CHECKERROR true

using std::cout;
using std::endl;

bool advanced_sort_collums_of_spectrum = true;


bool NumberOfRepsCriterion(const SDimCount &DimCounter1, const SDimCount &DimCounter2);

bool ColumnCriterion1(const SColumnOfSpectrum &Column1, const SColumnOfSpectrum &Column2);
bool ColumnCriterion2(const SColumnOfSpectrum &Column1, const SColumnOfSpectrum &Column2);

bool RowCriterion1(const SRepresentation &rep1, const SRepresentation &rep2);
bool RowCriterion2(const SRepresentation &rep1, const SRepresentation &rep2);


bool operator==(const SDimCount &DimCounter1, const SDimCount &DimCounter2)
{
  return ((DimCounter1.Counter == DimCounter2.Counter) && (DimCounter1.Dim == DimCounter2.Dim));
}


const SSpectrum &CSpectrum::GetSortedSpectrum(unsigned i) const
{
  if (i < this->SortedSpectra.size())
    return this->SortedSpectra[i];

  return this->SortedSpectra[0];
}


unsigned CSpectrum::GetLargestMultiplicity() const
{
  if (this->Spectrum.Spectrum.size() == 0)
    return 0;
  
  return this->Spectrum.Spectrum[0].Multiplicity;
}



/* ########################################################################################
######   CSpectrum()                                                                 ######
######                                                                               ######
######   Version: 31.10.2014                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard constructor of a CSpectrum object. No content is specified.        ######
######   The spectrum contains a list of representations (stored in the private      ######
######   member variable "Representations") with respect to the non-Abelian part of  ######
######   the gauge group, together with its multiplicity (stored in "Multiplicity"). ######
######################################################################################## */
CSpectrum::CSpectrum()
{
  this->UseSpecialU1 = false;
  this->NumberOfSupersymmetry = 0;
}



/* ########################################################################################
######   CSpectrum(const SConfig &Vacuum, ...)                                       ######
######                                                                               ######
######   Version: 05.11.2014                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Vacuum    : collect all fields from this vacuum and store them in this   ######
######                  CSpectrum object                                             ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Constructor of a CSpectrum object. Creates the spectrum of fields from the  ######
######   vacuum "Vacuum".                                                            ######
######################################################################################## */
CSpectrum::CSpectrum(const SConfig &Vacuum,  const vector<SUSYMultiplet> &SUSYTypes, string Label, const vector<vector<unsigned> > *SectorsForAdditionalIdentifier)
{
  this->Label = Label;

  unsigned pos_of_U1Hypercharge = 0;
  this->UseSpecialU1 = false;
  if (Vacuum.SymmetryGroup.observable_sector_U1s.size() == 1)
  {
    pos_of_U1Hypercharge = Vacuum.SymmetryGroup.observable_sector_U1s[0];
    this->UseSpecialU1 = (Vacuum.SymmetryGroup.U1s_AdditionalLabels[pos_of_U1Hypercharge] == "Y");
  }
  
  this->NumberOfSupersymmetry = Vacuum.InvariantSupercharges.size();
  this->Position_of_and_in_GaugeGroup = Vacuum.SymmetryGroup.Position_of_and_in_GaugeGroup;

  const vector<CField> &Fields = Vacuum.Fields;
  const size_t f1 = Fields.size();

  const vector<gaugeGroupFactor<double> > &GaugeGroupFactors = Vacuum.SymmetryGroup.GaugeGroup.factor;
  const size_t NumberOfGaugeGroups = GaugeGroupFactors.size();

  if (this->Position_of_and_in_GaugeGroup == 0)
    this->Position_of_and_in_GaugeGroup = NumberOfGaugeGroups;
  
  unsigned i = 0;

  // begin: gauge group
  this->Spectrum.GaugeGroup.assign(NumberOfGaugeGroups, 0);
  this->Spectrum.Rank.assign(NumberOfGaugeGroups, 0);

  const string A = "A";
  const string D = "D";
  const string E = "E";
  for (i = 0; i < NumberOfGaugeGroups; ++i)
  {
    const string algebra = GaugeGroupFactors[i].algebra;

    if (algebra[0] == A[0])
      this->Spectrum.GaugeGroup[i] = 1;
    else
    if (algebra[0] == D[0])
      this->Spectrum.GaugeGroup[i] = 2;
    else
    if (algebra[0] == E[0])
      this->Spectrum.GaugeGroup[i] = 3;

    this->Spectrum.Rank[i] = GaugeGroupFactors[i].rank;
  }
  // end: gauge group

  
  // begin: create the summary table of the spectrum
  SRepresentation NewRep;
  NewRep.Multiplicity = 0;
  NewRep.TotalMultiplicity = 0;
  NewRep.U1Charge = 0.0;
  NewRep.Representation.assign(NumberOfGaugeGroups, 0);

  vector<SRepresentation> &List = this->Spectrum.Spectrum;
  List.clear();

  size_t s1 = 0;
  size_t s2 = 0;
  unsigned j = 0;
  unsigned k = 0;
  unsigned n = 0;
  bool rep_equal = true;
  bool not_known = true;

  const bool UseSectorsForAdditionalIdentifier = (SectorsForAdditionalIdentifier != NULL);

  vector<double> NumberOfRepsPerSector;
  if (UseSectorsForAdditionalIdentifier)
  {
    s2 = SectorsForAdditionalIdentifier->size();
    NumberOfRepsPerSector.assign(s2, 0.0);
  }
  
  // define the SUSY types that shall be collected into this spectrum depending on the number of supersymmetries
  /*vector<SUSYMultiplet> SUSYTypes;
  if (this->NumberOfSupersymmetry == 0)
  {
    //SUSYTypes.push_back();
    //SUSYTypes.push_back();
    cout << "\nDefine SUSY types in CSpectrum::CSpectrum(const SConfig &Vacuum, string Label, const vector<vector<unsigned> > *SectorsForAdditionalIdentifier)." << endl;
  }
  else
  if (this->NumberOfSupersymmetry == 1)
  {
    SUSYTypes.push_back(LCModulus);
    SUSYTypes.push_back(LeftChiral);
  }
  else
  if (this->NumberOfSupersymmetry == 2)
  {
    SUSYTypes.push_back(Hyper);
    SUSYTypes.push_back(Halfhyper);
  }*/
    
  if (this->UseSpecialU1)
  {
    for (i = 0; i < f1; ++i)
    {
      const CField &Field = Fields[i];
      
      if (find(SUSYTypes.begin(), SUSYTypes.end(), Field.Multiplet) != SUSYTypes.end())
      {
        NewRep.Multiplet = Field.Multiplet;
        
        // begin: Additional Identifier
        if (UseSectorsForAdditionalIdentifier)
        {
          n = Field.SGElement.Get_n();
          k = Field.SGElement.Get_k();
      
          not_known = true;
          for (j = 0; j < s2; ++j)
          {
            if ((SectorsForAdditionalIdentifier->at(j)[0] == n) && (SectorsForAdditionalIdentifier->at(j)[1] == k))
            {
              ++NumberOfRepsPerSector[j];
              not_known = false;
              break;
            }
          }
          if (not_known)
          {
            cout << "Error in CSpectrum::CSpectrum(...) : sector for additional identifier not known." << endl;
            return;
          }
        }
        // end: Additional Identifier

        NewRep.U1Charge = Field.U1Charges[pos_of_U1Hypercharge];
        for (j = 0; j < NumberOfGaugeGroups; ++j)
          NewRep.Representation[j] = Field.Dimensions[j].Dimension;

        not_known = true;
        s1 = List.size();

        for (j = 0; j < s1; ++j)
        {
          SRepresentation &KnownRep = List[j];

          if ((fabs(NewRep.U1Charge - KnownRep.U1Charge) < 0.001) && (NewRep.Multiplet == KnownRep.Multiplet) && (NewRep.Representation == KnownRep.Representation))
          {
            ++KnownRep.Multiplicity;
            not_known = false;
            break;
          }
        }
          
        if (not_known)
        {
          NewRep.Multiplicity = 1;
          List.push_back(NewRep);
        }
      }
    }
  }
  else
  {
    for (i = 0; i < f1; ++i)
    {
      const CField &Field = Fields[i];
      
      if (find(SUSYTypes.begin(), SUSYTypes.end(), Field.Multiplet) != SUSYTypes.end())
      {
        NewRep.Multiplet = Field.Multiplet;

        // begin: Additional Identifier
        if (UseSectorsForAdditionalIdentifier)
        {
          n = Field.SGElement.Get_n();
          k = Field.SGElement.Get_k();
      
          not_known = true;
          for (j = 0; j < s2; ++j)
          {
            if ((SectorsForAdditionalIdentifier->at(j)[0] == n) && (SectorsForAdditionalIdentifier->at(j)[1] == k))
            {
              ++NumberOfRepsPerSector[j];
              not_known = false;
              break;
            }
          }
          if (not_known)
          {
            cout << "Error in CSpectrum::CSpectrum(...) : sector for additional identifier not known." << endl;
            return;
          }
        }
        // end: Additional Identifier

        for (j = 0; j < NumberOfGaugeGroups; ++j)
          NewRep.Representation[j] = Field.Dimensions[j].Dimension;

        not_known = true;
        s1 = List.size();

        for (j = 0; j < s1; ++j)
        {
          SRepresentation &KnownRep = List[j];

          if ((NewRep.Multiplet == KnownRep.Multiplet) && (NewRep.Representation == KnownRep.Representation))
          {
            ++KnownRep.Multiplicity;
            not_known = false;
            break;
          }
        }
          
        if (not_known)
        {
          NewRep.Multiplicity = 1;
          List.push_back(NewRep);
        }
      }
    }
  }
  // end: create the summary table of the spectrum

  s1 = List.size();
  for (i = 0; i < s1; ++i)
  {
    SRepresentation &rep = List[i];

    rep.TotalMultiplicity = rep.Multiplicity;

    for (j = 0; j < NumberOfGaugeGroups; ++j)
      rep.TotalMultiplicity *= abs(rep.Representation[j]);
  }

  // now sort the lines of the spectrum, using RowCriterion1
  stable_sort(List.begin(), List.end(), RowCriterion1);
  
  if (UseSectorsForAdditionalIdentifier)
  {
    // sorting is needed as some orbifolds can have a permutation symmetry along the twisted sectors, e.g. Z_2 x Z_2
    stable_sort(NumberOfRepsPerSector.begin(), NumberOfRepsPerSector.end());
  
    s1 = NumberOfRepsPerSector.size();
    for (i = 0; i < s1; ++i)
    {
      if (fabs(NumberOfRepsPerSector[i]) > 0.0001)
        this->AdditionalIdentifier.push_back(NumberOfRepsPerSector[i]);
    }
  }
  this->ReOrder();
}



/* ########################################################################################
######   ~CSpectrum()                                                                ######
######                                                                               ######
######   Version: 18.04.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard destructor of a CSpectrum object.                                  ######
######################################################################################## */
CSpectrum::~CSpectrum()
{
}



/* ########################################################################################
######   operator==(const CSpectrum &Spectrum2) const                                ######
######                                                                               ######
######   Version: 04.11.2014                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   Spectrum2    : a second spectrum to compare with                            ######
######   output:                                                                     ######
######   return value : are the two spectra identical?                               ######
###########################################################################################
######   description:                                                                ######
######   Compares two CSpectrum obejcts. If this operator returns "false" the        ######
######   spectra are inequivalent, but if it returns "true" they still might be      ######
######   inequivalent because the comparison method was not good enough.             ###### 
######################################################################################## */
bool CSpectrum::operator==(const CSpectrum &Spectrum2) const
{
  if ((Spectrum2.NumberOfSupersymmetry != this->NumberOfSupersymmetry) || (this->UseSpecialU1 != Spectrum2.UseSpecialU1))
  {
    cout << "Warning in bool CSpectrum::operator==(const CSpectrum &Spectrum2) const." << endl;
    return false;
  }

  // begin: compare AdditionalIdentifier
  size_t s1 = this->AdditionalIdentifier.size();
  if (s1 != Spectrum2.AdditionalIdentifier.size())
    return false;

  unsigned i = 0;
  for (i = 0; i < s1; ++i)
  {
    if (fabs(this->AdditionalIdentifier[i] - Spectrum2.AdditionalIdentifier[i]) > 0.001)
      return false;
  }
  // end: compare AdditionalIdentifier

  // begin: compare gauge group
  // gauge groups are different inside E8xE8 factors
  if (this->Position_of_and_in_GaugeGroup != Spectrum2.Position_of_and_in_GaugeGroup)
    return false;
  
  // number of non-Abelian gauge group factors unequal
  const size_t g1 = this->Spectrum.GaugeGroup.size();
  if (g1 != Spectrum2.Spectrum.GaugeGroup.size())
    return false;

  // non-Abelian gauge group factors unequal
  for (i = 0; i < g1; ++i)
  {
    if (this->Spectrum.GaugeGroup[i] != Spectrum2.Spectrum.GaugeGroup[i])
      return false;

    if (this->Spectrum.Rank[i] != Spectrum2.Spectrum.Rank[i])
      return false;
  }
  // end: compare gauge group

  s1 = this->SortedSpectra.size();
  size_t s2 = Spectrum2.SortedSpectra.size();
  if ((s1 == 0) || (s2 == 0))
  {
    cout << "\n  Warning in bool CSpectrum::operator==(const CSpectrum &Spectrum2) const: list of sorted spectra is empty.\n" << endl;
    return false;
  }

  // spectra have different number of representations
  size_t s3 = this->Spectrum.Spectrum.size();
  if (s3 != Spectrum2.Spectrum.Spectrum.size())
    return false;

  // no representations
  if (s3 == 0)
    return true;
  
  unsigned j = 0;
  unsigned k = 0;
  unsigned l = 0;
  
  bool SpectraEqual = true;
  
  // run through the permutations of this spectrum
  for (i = 0; i < s1; ++i)
  {
    const vector<SRepresentation> &List1 = this->SortedSpectra[i].Spectrum;
    
    // run through the permutations of Spectrum2
    for (j = 0; j < s2; ++j)
    {
      const vector<SRepresentation> &List2 = Spectrum2.SortedSpectra[j].Spectrum;
      
      SpectraEqual = true;
      for (k = 0; SpectraEqual && (k < s3); ++k)
      {
        const SRepresentation &Rep1 = List1[k];
        const SRepresentation &Rep2 = List2[k];

        // if SUSY type are different or multiplicities differ or |q_1| != |q_2|
        if ((Rep1.Multiplet != Rep2.Multiplet) || (Rep1.Multiplicity != Rep2.Multiplicity) || (this->UseSpecialU1 && (fabs(fabs(Rep1.U1Charge) - fabs(Rep2.U1Charge)) > 0.001)))
          SpectraEqual = false;
        else
        {
          for (l = 0; l < g1; ++l)
          {
            if (abs(Rep1.Representation[l]) != abs(Rep2.Representation[l]))
            {
              SpectraEqual = false;
              break;
            }
          }
        }
      }
      if (SpectraEqual)
        return true;
    }
  }
  
  return false;
}



unsigned CSpectrum::OrderE8xE8(vector<SColumnOfSpectrum> &ColumnsPart1, vector<SColumnOfSpectrum> &ColumnsPart2)
{
  const size_t c1 = ColumnsPart1.size();
  const size_t c2 = ColumnsPart2.size();

  // begin: compare the number of non-Abelian factors
  if (c1 < c2)  
  {
    this->Position_of_and_in_GaugeGroup = c2;
    ColumnsPart1.swap(ColumnsPart2);
    return 1;
  }
  if (c1 > c2)
    return 0;
  // end: compare the number of non-Abelian factors

  if (c1 == 0)
    return 0;
  
  // now c1 == c2 and Position_of_and_in_GaugeGroup does not need to be changed

  unsigned val1 = 0;
  unsigned val2 = 0;

  // begin: compare the total rank
  for (vector<SColumnOfSpectrum>::iterator Column = ColumnsPart1.begin(); Column != ColumnsPart1.end(); ++Column)
    val1 += Column->Rank;

  for (vector<SColumnOfSpectrum>::iterator Column = ColumnsPart2.begin(); Column != ColumnsPart2.end(); ++Column)
    val2 += Column->Rank;

  if (val1 > val2)
  {
    ColumnsPart1.swap(ColumnsPart2);
    return 1;
  }
  if (val1 < val2)
    return 0;
  // end: compare the total rank

  // begin: compare the sorted gauge groups
  vector<SColumnOfSpectrum>::iterator Column2 = ColumnsPart2.begin();
  for (vector<SColumnOfSpectrum>::iterator Column1 = ColumnsPart1.begin(); Column1 != ColumnsPart1.end(); ++Column1)
  {
    if (Column1->GaugeGroup > Column2->GaugeGroup)
    {
      ColumnsPart1.swap(ColumnsPart2);
      return 1;
    }
    if (Column1->GaugeGroup < Column2->GaugeGroup)
      return 0;
    
    if (Column1->Rank > Column2->Rank)
    {
      ColumnsPart1.swap(ColumnsPart2);
      return 1;
    }
    if (Column1->Rank < Column2->Rank)
      return 0;
    
    ++Column2;
  }
  // end: compare the sorted gauge groups


  // begin: compare the total sum
  val1 = 0;
  val2 = 0;
  for (vector<SColumnOfSpectrum>::iterator Column = ColumnsPart1.begin(); Column != ColumnsPart1.end(); ++Column)
    val1 += Column->Sum;

  for (vector<SColumnOfSpectrum>::iterator Column = ColumnsPart2.begin(); Column != ColumnsPart2.end(); ++Column)
    val2 += Column->Sum;

  if (val1 > val2)
  {
    ColumnsPart1.swap(ColumnsPart2);
    return 1;
  }
  if (val1 < val2)
    return 0;
  // end: compare the total sum

  
  // begin: compare the total product
  val1 = 0;
  val2 = 0;
  for (vector<SColumnOfSpectrum>::iterator Column = ColumnsPart1.begin(); Column != ColumnsPart1.end(); ++Column)
    val1 += Column->Product;

  for (vector<SColumnOfSpectrum>::iterator Column = ColumnsPart2.begin(); Column != ColumnsPart2.end(); ++Column)
    val2 += Column->Product;

  if (val1 > val2)
  {
    ColumnsPart1.swap(ColumnsPart2);
    return 1;
  }
  if (val1 < val2)
    return 0;
  // end: compare the total product

  
  // begin: compare the total number of ones
  val1 = 0;
  val2 = 0;
  for (vector<SColumnOfSpectrum>::iterator Column = ColumnsPart1.begin(); Column != ColumnsPart1.end(); ++Column)
    val1 += Column->NumberOfOnes;

  for (vector<SColumnOfSpectrum>::iterator Column = ColumnsPart2.begin(); Column != ColumnsPart2.end(); ++Column)
    val2 += Column->NumberOfOnes;

  if (val1 > val2)
  {
    ColumnsPart1.swap(ColumnsPart2);
    return 1;
  }
  if (val1 < val2)
    return 0;
  // end: compare the total number of ones

  unsigned i = 0;
  unsigned j = 0;
  const size_t NumberOfReps = ColumnsPart1[0].Representations.size();
  
  // begin: compare sum over dimensions
  val1 = 0;
  val2 = 0;
  for (i = 0; i < c1; ++i)
  {
    const vector<int> &Reps1 = ColumnsPart1[i].Representations;
    const vector<int> &Reps2 = ColumnsPart2[i].Representations;

    for (j = 0; j < NumberOfReps; ++j)
    {
      val1 += abs(Reps1[j]);
      val2 += abs(Reps2[j]);
    }
  }
  if (val1 > val2)  
  {
    ColumnsPart1.swap(ColumnsPart2);
    return 1;
  }
  if (val1 < val2)
    return 0;
  // end: compare sum over dimensions

  // begin: compare greatest dimensions
  unsigned GreatestDim1 = 1;
  unsigned GreatestDim2 = 1;

  for (j = 0; j < NumberOfReps; ++j)
  {
    val1 = 1;
    val2 = 1;

    for (i = 0; i < c1; ++i)
    {
      val1 *= abs(ColumnsPart1[i].Representations[j]);
      val2 *= abs(ColumnsPart2[i].Representations[j]);
    }
    if (val1 > GreatestDim1)
      GreatestDim1 = val1;
    if (val2 > GreatestDim2)
      GreatestDim2 = val2;
  }
  
  if (GreatestDim1 > GreatestDim2)  
  {
    ColumnsPart1.swap(ColumnsPart2);
    return 1;
  }
  if (GreatestDim1 < GreatestDim2)
    return 0;
  // end: compare greatest dimensions
  
  // cannot decide
  return 2;
}



/* ########################################################################################
######   bool CSpectrum::IsSpectrumEmpty() const                                     ######
######                                                                               ######
######   Version: 06.12.2013                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : is this spectrum empty?                                      ######
###########################################################################################
######   description:                                                                ######
######   Checks whether this spectrum contains some representations.                 ######
######################################################################################## */
bool CSpectrum::IsSpectrumEmpty() const
{
  return (this->Spectrum.Spectrum.size() == 0);
}




bool CSpectrum::SSpectrumKnown(const SSpectrum &NewSpectrum) const
{
  unsigned i = 0;
  unsigned j = 0;
  unsigned k = 0;
  
  const vector<SRepresentation> &NewList = NewSpectrum.Spectrum;
  
  const size_t Number_Of_Representations   = NewList.size();
  const size_t Number_Of_GaugeGroupFactors = NewSpectrum.GaugeGroup.size();
  
  bool SpectraEqual = true;
  const size_t s1 = SortedSpectra.size();
  
  for (i = 0; i < s1; ++i)
  {
    const SSpectrum               &KnownSpectrum = this->SortedSpectra[i];
    const vector<SRepresentation> &KnownList     = KnownSpectrum.Spectrum;
    
    // should be the same because the columns are already sorted (partially)
    if ((Number_Of_Representations == KnownList.size()) && (Number_Of_GaugeGroupFactors == KnownSpectrum.GaugeGroup.size()) 
     && (NewSpectrum.GaugeGroup == KnownSpectrum.GaugeGroup) && (NewSpectrum.Rank == KnownSpectrum.Rank))
    {
      SpectraEqual = true;
      for (j = 0; SpectraEqual && (j < Number_Of_Representations); ++j)
      {
        const SRepresentation &LineOfKnownSpectrum = KnownList[j];
        const SRepresentation &LineOfNewSpectrum   = NewList[j];

        // check SUSY type, the multiplicities and the absolute values of the U(1) charges
        if ((LineOfKnownSpectrum.Multiplet != LineOfNewSpectrum.Multiplet)
         || (LineOfKnownSpectrum.Multiplicity != LineOfNewSpectrum.Multiplicity)
         || (fabs(fabs(LineOfKnownSpectrum.U1Charge) - fabs(LineOfNewSpectrum.U1Charge)) > 0.0001))
          SpectraEqual = false;
        else
        {
          const vector<int> &KnownRep = LineOfKnownSpectrum.Representation;
          const vector<int> &NewRep   = LineOfNewSpectrum.Representation;
          
          // check the representations ignoring complex conjugation
          for (k = 0; k < Number_Of_GaugeGroupFactors; ++k)
          {
            if (abs(KnownRep[k]) != abs(NewRep[k]))
            {
              SpectraEqual = false;
              break;
            }
          }
        }
      }
      if (SpectraEqual)
        return true;
    }
  }
  return false;
}


/* ########################################################################################
######   bool CSpectrum::InsertSortedColumnsIntoExistingSpectrum(...) const          ######
######                                                                               ######
######   Version: 22.01.2014                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Insert sorted columns into an existing SSpectrum. Multiplicity, U1Charge    ######
######   and TotalMultiplicity of the existing spectrum are not affected.            ######
######################################################################################## */
bool CSpectrum::InsertSortedColumnsIntoExistingSpectrum(SSpectrum &Spectrum, const vector<vector<SColumnOfSpectrum> > &Columns) const
{
  unsigned j = 0;
  unsigned k = 0;
  unsigned l = 0;
  unsigned m = 0;

  const size_t NumberOfBlocks = Columns.size();
  size_t s1 = 0;

  const size_t Number_Of_Representations = Spectrum.Spectrum.size();

  if ((Spectrum.GaugeGroup.size() == 0) || (Spectrum.Rank.size() == 0))
  {
    unsigned NumberOfGaugeGroups = 0;
    for (j = 0; j < NumberOfBlocks; ++j)
      NumberOfGaugeGroups += Columns[j].size();
      
    Spectrum.GaugeGroup.assign(NumberOfGaugeGroups, 0);
    Spectrum.Rank.assign(NumberOfGaugeGroups, 0);
  }
  
  // begin: insert columns into existing spectrum
  l = 0;  // l runs through the non-Abelian gauge group factors
  for (j = 0; j < NumberOfBlocks; ++j)
  {
    const vector<SColumnOfSpectrum> &ColumnsBlock = Columns[j];
    
    s1 = ColumnsBlock.size();

    #ifdef CHECKERROR
    for (k = 0; k < s1; ++k)
    {
      if (ColumnsBlock[k].Representations.size() != Number_Of_Representations)
      {
        cout << "\n  Warning in bool CSpectrum::InsertSortedColumnsIntoExistingSpectrum(...): number of representations not correct.\n" << endl;
        return false;
      }
    }
    #endif

    for (k = 0; k < s1; ++k)
    {
      Spectrum.GaugeGroup[l] = ColumnsBlock[k].GaugeGroup;
      Spectrum.Rank[l]       = ColumnsBlock[k].Rank;
      ++l;
    }
  }
      
  for (j = 0; j < Number_Of_Representations; ++j)
  {
    SRepresentation &rep  = Spectrum.Spectrum[j];

    m = 0;  // m runs through the non-Abelian gauge group factors
    for (k = 0; k < NumberOfBlocks; ++k)
    {
      const vector<SColumnOfSpectrum> &ColumnsBlock = Columns[k];
      
      s1 = ColumnsBlock.size();
      for (l = 0; l < s1; ++l)
      {
        rep.Representation[m] = ColumnsBlock[l].Representations[j];
        ++m;
      }
    }
  }
  // end: insert columns into existing spectrum
  
  return true;
}



/* ########################################################################################
######   bool CSpectrum::ReOrder()                                                   ######
######                                                                               ######
######   Version: 05.11.2014                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Reorder the rows and collumns of the spectrum so that it can be compared    ######
######   to other spectra.                                                           ######
######################################################################################## */
bool CSpectrum::ReOrder()
{
  const size_t Number_Of_GaugeGroupFactors = this->Spectrum.GaugeGroup.size();
  if (Number_Of_GaugeGroupFactors == 0)
  {
    this->SortedSpectra.push_back(this->Spectrum);
    if (this->UseSpecialU1)
      cout << "\n  Warning in bool CSpectrum::ReOrder(): no non-Abelian gauge group but U(1) used for CSpectrum object.\n" << endl;
    return true;
  }

  const size_t Number_Of_Representations = this->Spectrum.Spectrum.size();
  
  unsigned i = 0;
  unsigned j = 0;
  unsigned k = 0;
  unsigned l = 0;
  unsigned m = 0;
  unsigned p = 0;
  
  int c = 0;
  
  size_t s1 = 0;
  size_t s2 = 0;
  size_t s3 = 0;

  bool DimUnknown = true;
  
  // begin: create the Columns
  SColumnOfSpectrum EmptyColumn;
  EmptyColumn.Sum                   = 0;
  EmptyColumn.AdditionalMultSum     = 0;
  EmptyColumn.AdditionalU1Sum       = 0;
  EmptyColumn.Product               = 1;
  EmptyColumn.AbsNetNumberOfNForSUN = 0;
  EmptyColumn.NumberOfOnes          = 0;
  EmptyColumn.Representations.assign(Number_Of_Representations,1);

  unsigned LengthPart1 = this->Position_of_and_in_GaugeGroup;
  unsigned LengthPart2 = Number_Of_GaugeGroupFactors - this->Position_of_and_in_GaugeGroup;

  const bool TwoParts = ((this->Position_of_and_in_GaugeGroup > 0) && (this->Position_of_and_in_GaugeGroup < Number_Of_GaugeGroupFactors));
  if (!TwoParts)
  {
    LengthPart1 = Number_Of_GaugeGroupFactors;
    LengthPart2 = 0;
  }
  
  vector<SColumnOfSpectrum> ColumnsPart1(LengthPart1, EmptyColumn);
  vector<vector<SColumnOfSpectrum> > TmpColumns(1, ColumnsPart1);

  vector<vector<vector<SColumnOfSpectrum> > > AllColumnPermutations(1, TmpColumns);

  vector<vector<SColumnOfSpectrum> > &Columns = AllColumnPermutations[0];
  
  if (TwoParts)
  {
    vector<SColumnOfSpectrum> ColumnsPart2(LengthPart2, EmptyColumn);
    Columns.push_back(ColumnsPart2);
  }

  unsigned NumberOfBlocks = Columns.size(); // either 1 for SO(32) or 2 for E8 x E8

  j = 0; // j runs through the non-Abelian gauge group factors
  for (i = 0; i < NumberOfBlocks; ++i)
  {
    vector<SColumnOfSpectrum> &ColumnsPart = Columns[i];
    
    for (vector<SColumnOfSpectrum>::iterator Column = ColumnsPart.begin(); Column != ColumnsPart.end(); Column++)
    {
      // GaugeGroup SU : 1
      // GaugeGroup SO : 2
      // GaugeGroup E  : 3
      Column->GaugeGroup = this->Spectrum.GaugeGroup[j];
      Column->Rank       = this->Spectrum.Rank[j];
        
      for (k = 0; k < Number_Of_Representations; ++k)
      {
        int &Dim = this->Spectrum.Spectrum[k].Representation[j];

        if (Dim == 1)
          ++Column->NumberOfOnes;
        else
        {
          Column->Sum     += abs(Dim);
          Column->Product *= abs(Dim);
        }
      
        Column->Representations[k] = Dim;
      }
      ++j;
    }
  }
  // end: create the Columns


  // begin: order E8 x E8 factors
  // this should be very fast, a unique sorting of E8xE8 is not needed
  if (TwoParts)
  {
    // "OrderE8xE8" assumes that the columns in each E8 factor are ordered by their gauge group, number of ones, sum and product
    if (this->OrderE8xE8(Columns[0], Columns[1]) == 2)
    {
      // order of E8 x E8 factors not knwon, hence both orderings are saved
      vector<vector<SColumnOfSpectrum> > SwapedE8Columns = Columns;
      SwapedE8Columns[0].swap(SwapedE8Columns[1]);

      AllColumnPermutations.push_back(SwapedE8Columns);
    }
  }
  // end: order E8 x E8 factors
  
  size_t NumberOfAllColumnPermutations = AllColumnPermutations.size();
  
  // begin: run through the inequivalent orderings of E8 x E8
  for (p = 0; p < NumberOfAllColumnPermutations; ++p)
  {
    vector<vector<SColumnOfSpectrum> > &Columns = AllColumnPermutations[p];
    NumberOfBlocks = Columns.size();

    // begin: sort the columns
    for (i = 0; i < NumberOfBlocks; ++i)
    {
      vector<SColumnOfSpectrum> &ColumnsPart = Columns[i];
      stable_sort(ColumnsPart.begin(), ColumnsPart.end(), ColumnCriterion1);
    }
    // end: sort the columns


    // begin: compare two adjacent columns and split into sub-blocks
    for (c = 0; c < NumberOfBlocks; ++c)
    {
      vector<SColumnOfSpectrum> &ColumnsBlock = Columns[c];
      
      s1 = ColumnsBlock.size();
      for (j = 1; j < s1; ++j)
      {
        SColumnOfSpectrum &Column1 = ColumnsBlock[j-1];
        SColumnOfSpectrum &Column2 = ColumnsBlock[j];
        
        if ((Column1.GaugeGroup != Column2.GaugeGroup) || (Column1.Rank != Column2.Rank)
        || (Column1.Sum != Column2.Sum) || (Column1.Product != Column2.Product) || (Column1.NumberOfOnes != Column2.NumberOfOnes))
        {
          vector<SColumnOfSpectrum> NewBlock(ColumnsBlock.begin() + j, ColumnsBlock.end());
          ColumnsBlock.resize(j);
          Columns.insert(Columns.begin() + c + 1, NewBlock);
          --c;
          ++NumberOfBlocks;
          break;
        }
      }
    }
    // end: compare two adjacent columns and split into sub-blocks

    NumberOfBlocks = Columns.size();
    
    if (advanced_sort_collums_of_spectrum)
    {
      // if any block still contains more than one column additional sorting might help
      if (NumberOfBlocks != Number_Of_GaugeGroupFactors)
      {
        // begin: create additional information needed to compare columns
        for (i = 0; i < NumberOfBlocks; ++i)
        {
          vector<SColumnOfSpectrum> &ColumnsPart = Columns[i];
        
          for (vector<SColumnOfSpectrum>::iterator Column = ColumnsPart.begin(); Column != ColumnsPart.end(); Column++)
          {
            Column->AdditionalU1Sum       = 0;
            Column->AdditionalMultSum     = 0;
            Column->AbsNetNumberOfNForSUN = 0;
            
            vector<SDimCount> &NumberOfReps = Column->NumberOfReps;
            NumberOfReps.clear();

            for (j = 0; j < Number_Of_Representations; ++j)
            {
              int &Dim = Column->Representations[j];

              const SRepresentation &Rep = this->Spectrum.Spectrum[j];
              
              if (this->UseSpecialU1)
                Column->AdditionalU1Sum += abs(rational<int>(Dim * Rep.Multiplicity) * D2Rat(Rep.U1Charge));
            
              if (Dim != 1)
              {
                Column->AdditionalMultSum += abs(Dim) * Rep.Multiplicity;

                DimUnknown = true;
                s1 = NumberOfReps.size();
                for (k = 0; DimUnknown && (k < s1); ++k)
                {
                  if (NumberOfReps[k].Dim == abs(Dim))
                  {
                    DimUnknown = false;
                    ++NumberOfReps[k].Counter;
                  }
                }
                if (DimUnknown)
                {
                  SDimCount tmpCount;
                  tmpCount.Dim     = abs(Dim);
                  tmpCount.Counter = 1;
                  NumberOfReps.push_back(tmpCount);
                }

                // count net number of N of SU(N) for N>2
                if ((Column->GaugeGroup == 1) && (Column->Rank != 1) && (abs(Dim) == Column->Rank+1))
                {
                  if (Dim > 0)
                    ++Column->AbsNetNumberOfNForSUN;
                  else
                    --Column->AbsNetNumberOfNForSUN;
                }
              }
            }
            stable_sort(NumberOfReps.begin(), NumberOfReps.end(), NumberOfRepsCriterion);

            Column->AbsNetNumberOfNForSUN = abs(Column->AbsNetNumberOfNForSUN);
            Column->AdditionalU1Sum       = abs(Column->AdditionalU1Sum);
          }
        }
        // end: create additional information needed to compare columns

        // begin: sort columns
        for (i = 0; i < NumberOfBlocks; ++i)
        {
          vector<SColumnOfSpectrum> &ColumnsBlock = Columns[i];
          
          s1 = ColumnsBlock.size();
          if (s1 != 1)
            stable_sort(ColumnsBlock.begin(), ColumnsBlock.end(), ColumnCriterion2);
        }
        // end: sort columns
            
        // begin: compare two adjacent columns and split the block if these columns are different
        for (c = 0; c < NumberOfBlocks; ++c)
        {
          vector<SColumnOfSpectrum> &ColumnsBlock = Columns[c];
          
          s1 = ColumnsBlock.size();
          if (s1 != 1)
          {
            for (j = 1; j < s1; ++j)
            {
              SColumnOfSpectrum &Column1 = ColumnsBlock[j-1];
              SColumnOfSpectrum &Column2 = ColumnsBlock[j];

              if ((Column1.AbsNetNumberOfNForSUN != Column2.AbsNetNumberOfNForSUN) || (Column1.AdditionalMultSum != Column2.AdditionalMultSum) 
                || (Column1.NumberOfReps != Column2.NumberOfReps) || (Column1.AdditionalU1Sum != Column2.AdditionalU1Sum))
              {
                #ifdef CHECKERROR
                vector<bool> column_error1(4,true);
                column_error1[0] = (Column1.AbsNetNumberOfNForSUN < Column2.AbsNetNumberOfNForSUN);
                column_error1[1] = (Column1.AdditionalMultSum < Column2.AdditionalMultSum);
                column_error1[2] = (Column1.AdditionalU1Sum < Column2.AdditionalU1Sum);
                column_error1[3] = (Column1.NumberOfReps.size() < Column2.NumberOfReps.size());
                vector<bool> column_error2(4,true);
                column_error2[0] = (Column1.AbsNetNumberOfNForSUN > Column2.AbsNetNumberOfNForSUN);
                column_error2[1] = (Column1.AdditionalMultSum > Column2.AdditionalMultSum);
                column_error2[2] = (Column1.AdditionalU1Sum > Column1.AdditionalU1Sum);
                column_error2[3] = (Column1.NumberOfReps.size() > Column2.NumberOfReps.size());
                if (Column1.NumberOfReps.size() == Column2.NumberOfReps.size())
                {
                  for (k = 0; k < Column1.NumberOfReps.size(); ++k)
                  {
                    column_error1.push_back((Column1.NumberOfReps[k].Counter < Column2.NumberOfReps[k].Counter));
                    column_error1.push_back((Column1.NumberOfReps[k].Dim < Column2.NumberOfReps[k].Dim));
                    column_error2.push_back((Column1.NumberOfReps[k].Counter > Column2.NumberOfReps[k].Counter));
                    column_error2.push_back((Column1.NumberOfReps[k].Dim > Column2.NumberOfReps[k].Dim));
                  }
                }
                for (k = 0; k < column_error1.size(); ++k)
                {
                  if (column_error2[k])
                    break;
                  if (column_error1[k])
                  {
                    cout << "\nSplit columns:\n";
                    cout << "  Column1.AbsNetNumberOfNForSUN = " << Column1.AbsNetNumberOfNForSUN << "\n  Column2.AbsNetNumberOfNForSUN = " << Column2.AbsNetNumberOfNForSUN << endl;
                    cout << "  Column1.AdditionalMultSum = " << Column1.AdditionalMultSum << "\n  Column2.AdditionalMultSum = " << Column2.AdditionalMultSum << endl;
                    cout << "  #Column1.NumberOfReps = " << Column1.NumberOfReps.size() << "\n  #Column2.NumberOfReps = " << Column2.NumberOfReps.size() << endl;
                    cout << "  Column1.AdditionalU1Sum = " << Column1.AdditionalU1Sum << "\n  Column2.AdditionalU1Sum = " << Column2.AdditionalU1Sum << endl;
                  }
                }
                #endif

                vector<SColumnOfSpectrum> NewBlock(ColumnsBlock.begin() + j, ColumnsBlock.end());
                ColumnsBlock.resize(j);
                
                Columns.insert(Columns.begin() + c + 1, NewBlock);
                --c;
                ++NumberOfBlocks;
                break;
              }
            }
          }
        }
        // end: compare two adjacent columns and split the block if these columns are different
      }
    }
  }
  // end: run through the inequivalent orderings of E8 x E8

  vector<vector<vector<SColumnOfSpectrum> > > NewAllPermutationsColumns;
  vector<vector<SColumnOfSpectrum> > NewColumns;
  const vector<SColumnOfSpectrum>    EmptyBlock;
  
  SSpectrum NewSpectrum;
  NewSpectrum.GaugeGroup.assign(Number_Of_GaugeGroupFactors, 0);
  NewSpectrum.Rank.assign(Number_Of_GaugeGroupFactors, 0);
  NewSpectrum = this->Spectrum;
  
  this->SortedSpectra.clear();

  // begin: create all permutations of those columns for which it can not be decided whether they are different or not
  bool NewPermutationsCreated = true;
  while (NewPermutationsCreated)
  {
    NewPermutationsCreated = false;
        
    NumberOfAllColumnPermutations = AllColumnPermutations.size();
    for (p = 0; p < NumberOfAllColumnPermutations; ++p)
    {
      const vector<vector<SColumnOfSpectrum> > &Columns = AllColumnPermutations[p];

      NumberOfBlocks = Columns.size();
       
      // begin: if any block still contains more than one column create all permutations
      if (!NewPermutationsCreated && (NumberOfBlocks != Number_Of_GaugeGroupFactors))
      {
        for (i = 0; i < NumberOfBlocks; ++i)
        {
          const vector<SColumnOfSpectrum> &ColumnsBlock = Columns[i];

          // if this block contains more than one column
          s3 = ColumnsBlock.size(); 
          if (s3 != 1)
          {
            vector<unsigned> Indices(s3, 0);
            for (j = 0; j < s3; ++j)
              Indices[j] = j;

            // begin: go through the permutations
            do
            {
              // begin: create permuted columns
              NewColumns = Columns;
              NewColumns[i].clear();
              NewColumns.insert(NewColumns.begin() + i + 1, s3-1, EmptyBlock);

              for (j = 0; j < s3; ++j)
                NewColumns[i+j].push_back(ColumnsBlock[Indices[j]]);

              NewAllPermutationsColumns.push_back(NewColumns);
              // end: create permuted columns

              NewPermutationsCreated = true;              
            } while ( std::next_permutation(Indices.begin(), Indices.end()) );
            // end: go through the permutations
            
            if (NewPermutationsCreated)
              break;
          }
        }
      }
      // end: if any block still contains more than one column create all permutations
      else
        NewAllPermutationsColumns.push_back(Columns);      
    }
    AllColumnPermutations = NewAllPermutationsColumns;
    NewAllPermutationsColumns.clear();
  }
  // end: create all permutations of those columns for which it can not be decided whether they are different or not
    
  // begin: save all permuted spectra
  NumberOfAllColumnPermutations = AllColumnPermutations.size();
  for (p = 0; p < NumberOfAllColumnPermutations; ++p)
  {
    const vector<vector<SColumnOfSpectrum> > &Columns = AllColumnPermutations[p];

    NumberOfBlocks = Columns.size();
       
    if (NumberOfBlocks != Number_Of_GaugeGroupFactors)
    {
      cout << "\n  Warning in bool CSpectrum::ReOrder(): block of columns should contain more than one column.\n" << endl;
      return false;
    }
    else
    // begin: if every block contains exactly one column
    {
      NewSpectrum = this->Spectrum;
      this->InsertSortedColumnsIntoExistingSpectrum(NewSpectrum, Columns);
      
      // sort the lines
      stable_sort(NewSpectrum.Spectrum.begin(), NewSpectrum.Spectrum.end(), RowCriterion2);
        
      // save new permutation
      if (!this->SSpectrumKnown(NewSpectrum))
        this->SortedSpectra.push_back(NewSpectrum);
    }
    // begin: if every block contains exactly one column
  }
  // end: save all permuted spectra
  
  if (this->SortedSpectra.size() == 0)
  {
    cout << "\n  Warning in bool CSpectrum::ReOrder(): no spectrum saved.\n" << endl;
    return false;
  }
  this->Spectrum = this->SortedSpectra[0];

  /*if (this->SortedSpectra.size() > 1)
  {
    cout << "\n Next example: Sorting not unique. #Perm = " << this->SortedSpectra.size() << endl;
    CPrint Print(Tstandard, &cout);
    for (unsigned i = 0; i < this->SortedSpectra.size(); ++i)
    {
      Print.PrintSSpectrum(this->SortedSpectra[i], false, this->Position_of_and_in_GaugeGroup);
      cout << endl;
    }
  }*/   

  return true;
}



/* ###########################################################################
######   bool RowCriterion1(...)                                        ######
######                                                                  ######
######   Defines an easy " > " for the lines of a spectrum.             ######
######   The order is not unique and does not need to be.               ######
########################################################################### */
bool RowCriterion1(const SRepresentation &rep1, const SRepresentation &rep2)
{
  const size_t Number_Of_GaugeGroupFactors = rep1.Representation.size();

  #ifdef CHECKERROR
  if (Number_Of_GaugeGroupFactors != rep2.Representation.size())
  {
    cout << "\n  Warning in bool RowCriterion1(const SRepresentation &rep1, const SRepresentation &rep2): sizes differ:\n  rep1.Representation.size() = " << rep1.Representation.size() << "\n  rep2.Representation.size() = " << rep2.Representation.size() << "\n  Return false." << endl;
    return false;
  }
  #endif

  // compares the multiplicity
  if (rep1.Multiplicity > rep2.Multiplicity)
    return true;
  if (rep1.Multiplicity < rep2.Multiplicity)
    return false;

  // compares the number of particles corresponding to the representation
  // for example 1(16,3,1) => 1 * 16 * 3 * 1 = 48
  if (rep1.TotalMultiplicity > rep2.TotalMultiplicity)
    return true;
  if (rep1.TotalMultiplicity < rep2.TotalMultiplicity)
    return false;

  // compares the extra U(1) charges
  // if |q_1| > |q_2|
  if (fabs(rep1.U1Charge) - fabs(rep2.U1Charge) > 0.0001)
    return true;
  // if |q_2| > |q_1|
  if (fabs(rep2.U1Charge) - fabs(rep1.U1Charge) > 0.0001)
    return false;

  // compares the SUSY type
  if (rep1.Multiplet < rep2.Multiplet)
    return true;
  if (rep1.Multiplet > rep2.Multiplet)
    return false;

  // cannot decide, but perfect order not needed here
  return false;
}



/* ###########################################################################
######   bool RowCriterion2(...)                                        ######
######                                                                  ######
######   Defines " > " for the lines of a spectrum                      ######
########################################################################### */
bool RowCriterion2(const SRepresentation &rep1, const SRepresentation &rep2)
{
  // begin: identical with RowCriterion1
  const size_t Number_Of_GaugeGroupFactors = rep1.Representation.size();

  #ifdef CHECKERROR
  if (Number_Of_GaugeGroupFactors != rep2.Representation.size())
  {
    cout << "\n  Warning in bool RowCriterion1(const SRepresentation &rep1, const SRepresentation &rep2): sizes differ:\n  rep1.Representation.size() = " << rep1.Representation.size() << "\n  rep2.Representation.size() = " << rep2.Representation.size() << "\n  Return false." << endl;
    return false;
  }
  #endif

  // compares the multiplicity
  if (rep1.Multiplicity > rep2.Multiplicity)
    return true;
  if (rep1.Multiplicity < rep2.Multiplicity)
    return false;

  // compares the number of particles corresponding to the representation
  // for example 1(16,3,1) => 1 * 16 * 3 * 1 = 48
  if (rep1.TotalMultiplicity > rep2.TotalMultiplicity)
    return true;
  if (rep1.TotalMultiplicity < rep2.TotalMultiplicity)
    return false;

  // compares the extra U(1) charges
  // if |q_1| > |q_2|
  if (fabs(rep1.U1Charge) - fabs(rep2.U1Charge) > 0.0001)
    return true;
  // if |q_2| > |q_1|
  if (fabs(rep2.U1Charge) - fabs(rep1.U1Charge) > 0.0001)
    return false;
  // begin: identical with RowCriterion1

  // which line of the spectrum has the first non-trivial entry
  // that is larger than the one of the other line
  for (unsigned i = 0; i < Number_Of_GaugeGroupFactors; ++i)
  {
    if (abs(rep1.Representation[i]) > abs(rep2.Representation[i]))
      return true;
    if (abs(rep1.Representation[i]) < abs(rep2.Representation[i]))
      return false;
  }

  // both lines are exactly equal
  if ((rep1.Multiplet == rep2.Multiplet) && (rep1.Multiplicity == rep2.Multiplicity) && (rep1.Representation == rep2.Representation) && (fabs(rep1.U1Charge - rep2.U1Charge) < 0.0001))
    return true;

  // if two lines of the spectrum are the same except for complex conjugation
  // then for example (16,3,1) is bigger than (-16,3,1)
  for (unsigned i = 0; i < Number_Of_GaugeGroupFactors; ++i)
  {
    if (rep1.Representation[i] > rep2.Representation[i])
      return true;
    if (rep1.Representation[i] < rep2.Representation[i])
      return false;
  }

  if (rep1.U1Charge - rep2.U1Charge > 0.0001)
    return true;
  if (rep2.U1Charge - rep1.U1Charge > 0.0001)
    return false;

  // compares the SUSY type
  if (rep1.Multiplet < rep2.Multiplet)
    return true;
  if (rep1.Multiplet > rep2.Multiplet)
    return false;
  
  cout << "\n  Warning in bool RowCriterion2(const SRepresentation &rep1, const SRepresentation &rep2): two lines of a spectrum must be different by construction.\n" << endl;
  return false;
}



/* ###########################################################################
######   bool ColumnCriterion1(...)                                     ######
######                                                                  ######
######   Defines " > " for the columns of a spectrum                    ######
########################################################################### */
bool ColumnCriterion1(const SColumnOfSpectrum &Column1, const SColumnOfSpectrum &Column2)
{
  if (Column1.GaugeGroup > Column2.GaugeGroup)
    return true;
  
  if (Column1.GaugeGroup < Column2.GaugeGroup)
    return false;

  if (Column1.Rank > Column2.Rank)
    return true;
  
  if (Column1.Rank < Column2.Rank)
    return false;
  
  if (Column1.Sum > Column2.Sum)
    return true;
  
  if (Column1.Sum < Column2.Sum)
    return false;

  if (Column1.Product > Column2.Product)
    return true;
  
  if (Column1.Product < Column2.Product)
    return false;

  if (Column1.NumberOfOnes > Column2.NumberOfOnes)
    return true;
  
  if (Column1.NumberOfOnes < Column2.NumberOfOnes)
    return false;

  // cannot decide, but perfect order not needed here
  return false;
}


/* ###########################################################################
######   bool ColumnCriterion2(...)                                     ######
######                                                                  ######
######   Defines " > " for the columns of a sub-block                   ######
########################################################################### */
bool ColumnCriterion2(const SColumnOfSpectrum &Column1, const SColumnOfSpectrum &Column2)
{
  if (Column1.AbsNetNumberOfNForSUN > Column2.AbsNetNumberOfNForSUN)
    return true;
  
  if (Column1.AbsNetNumberOfNForSUN < Column2.AbsNetNumberOfNForSUN)
    return false;

  if (Column1.AdditionalMultSum > Column2.AdditionalMultSum)
    return true;
  
  if (Column1.AdditionalMultSum < Column2.AdditionalMultSum)
    return false;

  if (Column1.AdditionalU1Sum > Column2.AdditionalU1Sum)
    return true;
  
  if (Column1.AdditionalU1Sum < Column2.AdditionalU1Sum)
    return false;
  
  const size_t s1 = Column1.NumberOfReps.size();
  const size_t s2 = Column2.NumberOfReps.size();
  if (s1 > s2)
    return true;
  
  if (s1 < s2)
    return false;
  
  for (unsigned i = 0; i < s1; ++i)
  {
    if (Column1.NumberOfReps[i].Counter > Column2.NumberOfReps[i].Counter)
      return true;

    if (Column1.NumberOfReps[i].Counter < Column2.NumberOfReps[i].Counter)
      return false;
    
    if (Column1.NumberOfReps[i].Dim > Column2.NumberOfReps[i].Dim)
      return true;

    if (Column1.NumberOfReps[i].Dim < Column2.NumberOfReps[i].Dim)
      return false;
  }
  
  // cannot decide, but perfect order not needed here
  return false;
}



/* ###########################################################################
######   bool NumberOfRepsCriterion(...)                                ######
######                                                                  ######
######   Defines " > " for representation counters                      ######
########################################################################### */
bool NumberOfRepsCriterion(const SDimCount &DimCounter1, const SDimCount &DimCounter2)
{
  if (DimCounter1.Counter > DimCounter2.Counter)
    return true;

  if (DimCounter1.Counter < DimCounter2.Counter)
    return false;

  if (DimCounter1.Dim > DimCounter2.Dim)
    return true;

  if (DimCounter1.Dim < DimCounter2.Dim)
    return false;
  
  return false;
}
