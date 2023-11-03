#include "cinequivalentspectra.h"
#include "cprint.h"


bool RowCriterion2(const SRepresentation &rep1, const SRepresentation &rep2);


/* ########################################################################################
######   CInequivalentModels()                                                       ######
######                                                                               ######
######   Version: 16.01.2014                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard constructor of a CInequivalentModels object. The spectra are       ######
######   sorted in the following way:                                                ######
######   "c1"            : highest multiplicity mod 4 is 0                           ######
######   "c2"            : highest multiplicity mod 4 is 1                           ######
######   "c3"            : highest multiplicity mod 4 is 2                           ######
######   "c4"            : highest multiplicity mod 4 is 3                           ######
######   Rank00 - Rank16 : the number of non-Abelian gauge group factors             ######
######   index 0-89      : number of different non-Abelian representations           ######
######################################################################################## */
CInequivalentModels::CInequivalentModels()
{
  vector<CSpectrum> tmp;
  this->InequivalentModels_c1_Rank00.assign(150, tmp);
  this->InequivalentModels_c1_Rank01.assign(150, tmp);
  this->InequivalentModels_c1_Rank02.assign(150, tmp);
  this->InequivalentModels_c1_Rank03.assign(150, tmp);
  this->InequivalentModels_c1_Rank04.assign(150, tmp);
  this->InequivalentModels_c1_Rank05.assign(150, tmp);
  this->InequivalentModels_c1_Rank06.assign(150, tmp);
  this->InequivalentModels_c1_Rank07.assign(150, tmp);
  this->InequivalentModels_c1_Rank08.assign(150, tmp);
  this->InequivalentModels_c1_Rank09.assign(150, tmp);
  this->InequivalentModels_c1_Rank10.assign(150, tmp);
  this->InequivalentModels_c1_Rank11.assign(150, tmp);
  this->InequivalentModels_c1_Rank12.assign(150, tmp);
  this->InequivalentModels_c1_Rank13.assign(150, tmp);
  this->InequivalentModels_c1_Rank14.assign(150, tmp);
  this->InequivalentModels_c1_Rank15.assign(150, tmp);
  this->InequivalentModels_c1_Rank16.assign(150, tmp);

  this->InequivalentModels_c2_Rank00.assign(150, tmp);
  this->InequivalentModels_c2_Rank01.assign(150, tmp);
  this->InequivalentModels_c2_Rank02.assign(150, tmp);
  this->InequivalentModels_c2_Rank03.assign(150, tmp);
  this->InequivalentModels_c2_Rank04.assign(150, tmp);
  this->InequivalentModels_c2_Rank05.assign(150, tmp);
  this->InequivalentModels_c2_Rank06.assign(150, tmp);
  this->InequivalentModels_c2_Rank07.assign(150, tmp);
  this->InequivalentModels_c2_Rank08.assign(150, tmp);
  this->InequivalentModels_c2_Rank09.assign(150, tmp);
  this->InequivalentModels_c2_Rank10.assign(150, tmp);
  this->InequivalentModels_c2_Rank11.assign(150, tmp);
  this->InequivalentModels_c2_Rank12.assign(150, tmp);
  this->InequivalentModels_c2_Rank13.assign(150, tmp);
  this->InequivalentModels_c2_Rank14.assign(150, tmp);
  this->InequivalentModels_c2_Rank15.assign(150, tmp);
  this->InequivalentModels_c2_Rank16.assign(150, tmp);

  this->InequivalentModels_c3_Rank00.assign(150, tmp);
  this->InequivalentModels_c3_Rank01.assign(150, tmp);
  this->InequivalentModels_c3_Rank02.assign(150, tmp);
  this->InequivalentModels_c3_Rank03.assign(150, tmp);
  this->InequivalentModels_c3_Rank04.assign(150, tmp);
  this->InequivalentModels_c3_Rank05.assign(150, tmp);
  this->InequivalentModels_c3_Rank06.assign(150, tmp);
  this->InequivalentModels_c3_Rank07.assign(150, tmp);
  this->InequivalentModels_c3_Rank08.assign(150, tmp);
  this->InequivalentModels_c3_Rank09.assign(150, tmp);
  this->InequivalentModels_c3_Rank10.assign(150, tmp);
  this->InequivalentModels_c3_Rank11.assign(150, tmp);
  this->InequivalentModels_c3_Rank12.assign(150, tmp);
  this->InequivalentModels_c3_Rank13.assign(150, tmp);
  this->InequivalentModels_c3_Rank14.assign(150, tmp);
  this->InequivalentModels_c3_Rank15.assign(150, tmp);
  this->InequivalentModels_c3_Rank16.assign(150, tmp);

  this->InequivalentModels_c4_Rank00.assign(150, tmp);
  this->InequivalentModels_c4_Rank01.assign(150, tmp);
  this->InequivalentModels_c4_Rank02.assign(150, tmp);
  this->InequivalentModels_c4_Rank03.assign(150, tmp);
  this->InequivalentModels_c4_Rank04.assign(150, tmp);
  this->InequivalentModels_c4_Rank05.assign(150, tmp);
  this->InequivalentModels_c4_Rank06.assign(150, tmp);
  this->InequivalentModels_c4_Rank07.assign(150, tmp);
  this->InequivalentModels_c4_Rank08.assign(150, tmp);
  this->InequivalentModels_c4_Rank09.assign(150, tmp);
  this->InequivalentModels_c4_Rank10.assign(150, tmp);
  this->InequivalentModels_c4_Rank11.assign(150, tmp);
  this->InequivalentModels_c4_Rank12.assign(150, tmp);
  this->InequivalentModels_c4_Rank13.assign(150, tmp);
  this->InequivalentModels_c4_Rank14.assign(150, tmp);
  this->InequivalentModels_c4_Rank15.assign(150, tmp);
  this->InequivalentModels_c4_Rank16.assign(150, tmp);

  this->ModelCounter = 0;
}



/* ########################################################################################
######   ~CInequivalentModels()                                                      ######
######                                                                               ######
######   Version: 07.04.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard destructor of a CInequivalentModels object.                        ######
######################################################################################## */
CInequivalentModels::~CInequivalentModels()
{
}


/* ###########################################################################
######   bool InequivalentModelsCriterion(...                           ######
######                                                                  ######
######   Defines " > " for sorted spectra                               ######
######                                                                  ######
########################################################################### */
bool InequivalentModelsCriterion(const CSpectrum &Model1, const CSpectrum &Model2)
{
  // begin: sort by the number of 4D gauge group factors from each E8 factor
  double val1 = Model1.GetPosition_of_and_in_GaugeGroup();
  double val2 = Model2.GetPosition_of_and_in_GaugeGroup();

  if (val1 > val2)
    return true;

  if (val1 < val2)
    return false;
  // end: sort by the number of 4D gauge group factors from each E8 factor

  // begin: sort by the largest multiplicity
  val1 = Model1.GetLargestMultiplicity();
  val2 = Model2.GetLargestMultiplicity();

  if (val1 > val2)
    return true;

  if (val1 < val2)
    return false;
  // end: sort by the largest multiplicity

  size_t s1 = Model1.GetNumberOfRepresentations();
  
  #ifdef CHECKERROR
  if (s1 != Model2.GetNumberOfRepresentations())
  {
    cout << "\n  Error in bool InequivalentModelsCriterion(const CSpectrum &Model1, const CSpectrum &Model2) : Number of representations must be equal. exit." << endl;
    exit(1);
  }
  #endif

  // begin: sort by the other, sorted multiplicities
  unsigned i = 0;
  for (i = 1; i < s1; ++i) // first multiplicity was checked before
  {
    val1 = Model1.GetSpectrum().Spectrum[i].Multiplicity;
    val2 = Model2.GetSpectrum().Spectrum[i].Multiplicity;

    if (val1 > val2)
      return true;

    if (val1 < val2)
      return false;
  }
  // end: sort by the other, sorted multiplicities
  
  
  // begin: sort by the total multiplicities
  for (i = 0; i < s1; ++i)
  {
    val1 = Model1.GetSpectrum().Spectrum[i].TotalMultiplicity;
    val2 = Model2.GetSpectrum().Spectrum[i].TotalMultiplicity;

    if (val1 > val2)
      return true;

    if (val1 < val2)
      return false;
  }
  // end: sort by the total multiplicities

  // begin: sort by the representations
  for (i = 0; i < s1; ++i)
  {
    if (RowCriterion2(Model1.GetSpectrum().Spectrum[i], Model2.GetSpectrum().Spectrum[i]))
      return true;

    if (RowCriterion2(Model2.GetSpectrum().Spectrum[i], Model1.GetSpectrum().Spectrum[i]))
      return false;
  }
  // end: sort by the representations

  // begin: sort by the additional spectrum identifiers
  const vector<double> &AdditionalIdentifier1 = Model1.GetAdditionalIdentifier();
  const vector<double> &AdditionalIdentifier2 = Model2.GetAdditionalIdentifier();
  
  s1 = AdditionalIdentifier1.size();
  size_t s2 = AdditionalIdentifier2.size();
  
  if (s1 > s2)
    return true;

  if (s1 < s2)
    return false;

  for (i = 0; i < s1; ++i)
  {
    if (AdditionalIdentifier1[i] - AdditionalIdentifier2[i] > 0.001)
      return true;

    if (AdditionalIdentifier2[i] - AdditionalIdentifier1[i] > 0.001)
      return false;
  }
  // end: sort by the additional spectrum identifiers
  
  return false;
}



bool CInequivalentModels::SortInequivalentModels()
{
  unsigned i = 0;
  size_t s1 = 0;

  s1 = this->InequivalentModels_c1_Rank00.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c1_Rank00[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c1_Rank01.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c1_Rank01[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c1_Rank02.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c1_Rank02[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c1_Rank03.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c1_Rank03[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c1_Rank04.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c1_Rank04[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c1_Rank05.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c1_Rank05[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c1_Rank06.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c1_Rank06[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c1_Rank07.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c1_Rank07[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c1_Rank08.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c1_Rank08[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c1_Rank09.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c1_Rank09[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c1_Rank10.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c1_Rank10[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c1_Rank11.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c1_Rank11[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c1_Rank12.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c1_Rank12[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c1_Rank13.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c1_Rank13[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c1_Rank14.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c1_Rank14[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c1_Rank15.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c1_Rank15[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c1_Rank16.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c1_Rank16[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  

  s1 = this->InequivalentModels_c2_Rank00.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c2_Rank00[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c2_Rank01.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c2_Rank01[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c2_Rank02.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c2_Rank02[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c2_Rank03.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c2_Rank03[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c2_Rank04.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c2_Rank04[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c2_Rank05.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c2_Rank05[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c2_Rank06.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c2_Rank06[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c2_Rank07.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c2_Rank07[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c2_Rank08.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c2_Rank08[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c2_Rank09.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c2_Rank09[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c2_Rank10.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c2_Rank10[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c2_Rank11.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c2_Rank11[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c2_Rank12.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c2_Rank12[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c2_Rank13.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c2_Rank13[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c2_Rank14.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c2_Rank14[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c2_Rank15.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c2_Rank15[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c2_Rank16.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c2_Rank16[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }


  s1 = this->InequivalentModels_c3_Rank00.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c3_Rank00[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c3_Rank01.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c3_Rank01[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c3_Rank02.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c3_Rank02[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c3_Rank03.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c3_Rank03[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c3_Rank04.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c3_Rank04[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c3_Rank05.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c3_Rank05[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c3_Rank06.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c3_Rank06[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c3_Rank07.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c3_Rank07[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c3_Rank08.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c3_Rank08[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c3_Rank09.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c3_Rank09[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c3_Rank10.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c3_Rank10[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c3_Rank11.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c3_Rank11[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c3_Rank12.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c3_Rank12[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c3_Rank13.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c3_Rank13[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c3_Rank14.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c3_Rank14[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c3_Rank15.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c3_Rank15[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c3_Rank16.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c3_Rank16[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }


  s1 = this->InequivalentModels_c4_Rank00.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c4_Rank00[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c4_Rank01.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c4_Rank01[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c4_Rank02.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c4_Rank02[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c4_Rank03.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c4_Rank03[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c4_Rank04.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c4_Rank04[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c4_Rank05.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c4_Rank05[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c4_Rank06.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c4_Rank06[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c4_Rank07.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c4_Rank07[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c4_Rank08.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c4_Rank08[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c4_Rank09.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c4_Rank09[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c4_Rank10.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c4_Rank10[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c4_Rank11.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c4_Rank11[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c4_Rank12.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c4_Rank12[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c4_Rank13.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c4_Rank13[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c4_Rank14.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c4_Rank14[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c4_Rank15.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c4_Rank15[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }
  s1 = this->InequivalentModels_c4_Rank16.size();
  for (i = 0; i < s1; ++i)
  {
    vector<CSpectrum> &InequivalentModels = this->InequivalentModels_c4_Rank16[i];
    stable_sort(InequivalentModels.begin(), InequivalentModels.end(), InequivalentModelsCriterion);
  }

  return true;
}



/* ########################################################################################
######   IsSpectrumUnknown(const CSpectrum &Spectrum, bool AddNewModel)              ######
######                                                                               ######
######   Version: 22.01.2014                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Spectrum    : check whether "Spectrum" is known                          ######
######   2) AddNewModel : add "Spectrum" if it is new                                ######
######   output:                                                                     ######
######   return value   : is "Spectrum" known?                                       ######
###########################################################################################
######   description:                                                                ######
######   Search in the corresponding list "InequivalentModels_cX_RankXX" for the     ######
######   CSpectrum object "Spectrum" and return true if it is unkown. Add "Spectrum" ######
######   if it is new and "AddNewModel" is true.                                     ######
######################################################################################## */
bool CInequivalentModels::IsSpectrumUnknown(const CSpectrum &Spectrum, bool AddNewModel)
{
  const bool PrintLabelsOfEquivalentModels = false;

  vector<vector<CSpectrum> > *InequivalentModels = NULL;

  const unsigned HighestMultiplicityMod4 = Spectrum.GetLargestMultiplicity() % 4;

  // use "create random orbifold from(*) if(inequivalent) #models(all)" for all 138 space groups with Abelian point group
  // create a statistics on Spectrum.GetNumberOfGaugeGroups
  // result:
  //   Spectrum.GetNumberOfGaugeGroups : 0  1   2   3    4   5   6   7  8  9 10 11 12 13 14 15 16
  //   #models                      : 1 35 184 695 1069 749 338 127 30 13  9  2  2  0  0  0  0
  // hence use the folling order in switch:
  //   Spectrum.GetNumberOfGaugeGroups : 4 5 3 6 2 7 1 8 9 10 11 12 0 13 14 15 16
  switch(Spectrum.GetNumberOfGaugeGroups())
  {
    case 4:
    {
      switch(HighestMultiplicityMod4)
      {
        case 0:
        {
          InequivalentModels = &this->InequivalentModels_c1_Rank04;
          break;
        }
        case 1:
        {
          InequivalentModels = &this->InequivalentModels_c2_Rank04;
          break;
        }
        case 2:
        {
          InequivalentModels = &this->InequivalentModels_c3_Rank04;
          break;
        }
        case 3:
        {
          InequivalentModels = &this->InequivalentModels_c4_Rank04;
          break;
        }
      }
      break;
    }
    case 5:
    {
      switch(HighestMultiplicityMod4)
      {
        case 0:
        {
          InequivalentModels = &this->InequivalentModels_c1_Rank05;
          break;
        }
        case 1:
        {
          InequivalentModels = &this->InequivalentModels_c2_Rank05;
          break;
        }
        case 2:
        {
          InequivalentModels = &this->InequivalentModels_c3_Rank05;
          break;
        }
        case 3:
        {
          InequivalentModels = &this->InequivalentModels_c4_Rank05;
          break;
        }
      }
      break;
    }
    case 3:
    {
      switch(HighestMultiplicityMod4)
      {
        case 0:
        {
          InequivalentModels = &this->InequivalentModels_c1_Rank03;
          break;
        }
        case 1:
        {
          InequivalentModels = &this->InequivalentModels_c2_Rank03;
          break;
        }
        case 2:
        {
          InequivalentModels = &this->InequivalentModels_c3_Rank03;
          break;
        }
        case 3:
        {
          InequivalentModels = &this->InequivalentModels_c4_Rank03;
          break;
        }
      }
      break;
    }
    case 6:
    {
      switch(HighestMultiplicityMod4)
      {
        case 0:
        {
          InequivalentModels = &this->InequivalentModels_c1_Rank06;
          break;
        }
        case 1:
        {
          InequivalentModels = &this->InequivalentModels_c2_Rank06;
          break;
        }
        case 2:
        {
          InequivalentModels = &this->InequivalentModels_c3_Rank06;
          break;
        }
        case 3:
        {
          InequivalentModels = &this->InequivalentModels_c4_Rank06;
          break;
        }
      }
      break;
    }
    case 2:
    {
      switch(HighestMultiplicityMod4)
      {
        case 0:
        {
          InequivalentModels = &this->InequivalentModels_c1_Rank02;
          break;
        }
        case 1:
        {
          InequivalentModels = &this->InequivalentModels_c2_Rank02;
          break;
        }
        case 2:
        {
          InequivalentModels = &this->InequivalentModels_c3_Rank02;
          break;
        }
        case 3:
        {
          InequivalentModels = &this->InequivalentModels_c4_Rank02;
          break;
        }
      }
      break;
    }
    case 7:
    {
      switch(HighestMultiplicityMod4)
      {
        case 0:
        {
          InequivalentModels = &this->InequivalentModels_c1_Rank07;
          break;
        }
        case 1:
        {
          InequivalentModels = &this->InequivalentModels_c2_Rank07;
          break;
        }
        case 2:
        {
          InequivalentModels = &this->InequivalentModels_c3_Rank07;
          break;
        }
        case 3:
        {
          InequivalentModels = &this->InequivalentModels_c4_Rank07;
          break;
        }
      }
      break;
    }
    case 1:
    {
      switch(HighestMultiplicityMod4)
      {
        case 0:
        {
          InequivalentModels = &this->InequivalentModels_c1_Rank01;
          break;
        }
        case 1:
        {
          InequivalentModels = &this->InequivalentModels_c2_Rank01;
          break;
        }
        case 2:
        {
          InequivalentModels = &this->InequivalentModels_c3_Rank01;
          break;
        }
        case 3:
        {
          InequivalentModels = &this->InequivalentModels_c4_Rank01;
          break;
        }
      }
      break;
    }
    case 8:
    {
      switch(HighestMultiplicityMod4)
      {
        case 0:
        {
          InequivalentModels = &this->InequivalentModels_c1_Rank08;
          break;
        }
        case 1:
        {
          InequivalentModels = &this->InequivalentModels_c2_Rank08;
          break;
        }
        case 2:
        {
          InequivalentModels = &this->InequivalentModels_c3_Rank08;
          break;
        }
        case 3:
        {
          InequivalentModels = &this->InequivalentModels_c4_Rank08;
          break;
        }
      }
      break;
    }
    case 9:
    {
      switch(HighestMultiplicityMod4)
      {
        case 0:
        {
          InequivalentModels = &this->InequivalentModels_c1_Rank09;
          break;
        }
        case 1:
        {
          InequivalentModels = &this->InequivalentModels_c2_Rank09;
          break;
        }
        case 2:
        {
          InequivalentModels = &this->InequivalentModels_c3_Rank09;
          break;
        }
        case 3:
        {
          InequivalentModels = &this->InequivalentModels_c4_Rank09;
          break;
        }
      }
      break;
    }
    case 10:
    {
      switch(HighestMultiplicityMod4)
      {
        case 0:
        {
          InequivalentModels = &this->InequivalentModels_c1_Rank10;
          break;
        }
        case 1:
        {
          InequivalentModels = &this->InequivalentModels_c2_Rank10;
          break;
        }
        case 2:
        {
          InequivalentModels = &this->InequivalentModels_c3_Rank10;
          break;
        }
        case 3:
        {
          InequivalentModels = &this->InequivalentModels_c4_Rank10;
          break;
        }
      }
      break;
    }
    case 11:
    {
      switch(HighestMultiplicityMod4)
      {
        case 0:
        {
          InequivalentModels = &this->InequivalentModels_c1_Rank11;
          break;
        }
        case 1:
        {
          InequivalentModels = &this->InequivalentModels_c2_Rank11;
          break;
        }
        case 2:
        {
          InequivalentModels = &this->InequivalentModels_c3_Rank11;
          break;
        }
        case 3:
        {
          InequivalentModels = &this->InequivalentModels_c4_Rank11;
          break;
        }
      }
      break;
    }
    case 12:
    {
      switch(HighestMultiplicityMod4)
      {
        case 0:
        {
          InequivalentModels = &this->InequivalentModels_c1_Rank12;
          break;
        }
        case 1:
        {
          InequivalentModels = &this->InequivalentModels_c2_Rank12;
          break;
        }
        case 2:
        {
          InequivalentModels = &this->InequivalentModels_c3_Rank12;
          break;
        }
        case 3:
        {
          InequivalentModels = &this->InequivalentModels_c4_Rank12;
          break;
        }
      }
      break;
    }
    case 0:
    {
      switch(HighestMultiplicityMod4)
      {
        case 0:
        {
          InequivalentModels = &this->InequivalentModels_c1_Rank00;
          break;
        }
        case 1:
        {
          InequivalentModels = &this->InequivalentModels_c2_Rank00;
          break;
        }
        case 2:
        {
          InequivalentModels = &this->InequivalentModels_c3_Rank00;
          break;
        }
        case 3:
        {
          InequivalentModels = &this->InequivalentModels_c4_Rank00;
          break;
        }
      }
      break;
    }
    case 13:
    {
      switch(HighestMultiplicityMod4)
      {
        case 0:
        {
          InequivalentModels = &this->InequivalentModels_c1_Rank13;
          break;
        }
        case 1:
        {
          InequivalentModels = &this->InequivalentModels_c2_Rank13;
          break;
        }
        case 2:
        {
          InequivalentModels = &this->InequivalentModels_c3_Rank13;
          break;
        }
        case 3:
        {
          InequivalentModels = &this->InequivalentModels_c4_Rank13;
          break;
        }
      }
      break;
    }
    case 14:
    {
      switch(HighestMultiplicityMod4)
      {
        case 0:
        {
          InequivalentModels = &this->InequivalentModels_c1_Rank14;
          break;
        }
        case 1:
        {
          InequivalentModels = &this->InequivalentModels_c2_Rank14;
          break;
        }
        case 2:
        {
          InequivalentModels = &this->InequivalentModels_c3_Rank14;
          break;
        }
        case 3:
        {
          InequivalentModels = &this->InequivalentModels_c4_Rank14;
          break;
        }
      }
      break;
    }
    case 15:
    {
      switch(HighestMultiplicityMod4)
      {
        case 0:
        {
          InequivalentModels = &this->InequivalentModels_c1_Rank15;
          break;
        }
        case 1:
        {
          InequivalentModels = &this->InequivalentModels_c2_Rank15;
          break;
        }
        case 2:
        {
          InequivalentModels = &this->InequivalentModels_c3_Rank15;
          break;
        }
        case 3:
        {
          InequivalentModels = &this->InequivalentModels_c4_Rank15;
          break;
        }
      }
      break;
    }
    case 16:
    {
      switch(HighestMultiplicityMod4)
      {
        case 0:
        {
          InequivalentModels = &this->InequivalentModels_c1_Rank16;
          break;
        }
        case 1:
        {
          InequivalentModels = &this->InequivalentModels_c2_Rank16;
          break;
        }
        case 2:
        {
          InequivalentModels = &this->InequivalentModels_c3_Rank16;
          break;
        }
        case 3:
        {
          InequivalentModels = &this->InequivalentModels_c4_Rank16;
          break;
        }
      }
      break;
    }
    default:
    {
      cout << "Warning in bool CInequivalentModels::IsSpectrumUnknown(...): rank of gauge group out of range. Return false." << endl;
      return false;
    }
  }
  if (InequivalentModels == NULL)
  {
    cout << "Warning in bool CInequivalentModels::IsSpectrumUnknown(...): pointer not set. Return false." << endl;
    return false;
  }

  const size_t reps = Spectrum.GetNumberOfRepresentations();
  if (reps >= InequivalentModels->size())
  {
    cout << "\n  Warning in bool CInequivalentModels::IsSpectrumUnknown(...): \"reps\" = " << reps << " out of range. Return false." << endl;
    return false;
  }

  vector<CSpectrum> &CurrentList = InequivalentModels->at(reps);

  const size_t s1 = CurrentList.size();
  
  if (PrintLabelsOfEquivalentModels)
  {
    for (unsigned i = 0; i < s1; ++i)
    {
      if (Spectrum == CurrentList[i])
      {
        cout << "\n  Equivalent models: " << Spectrum.Label << " and " << CurrentList[i].Label << endl;
        return false;
      }
    }
  }
  else
  {
    for (unsigned i = 0; i < s1; ++i)
    {
      if (Spectrum == CurrentList[i])
        return false;
    }
  }
  
  // spectrum is new
  if (AddNewModel)
  {
    CurrentList.push_back(Spectrum);
    ++this->ModelCounter;
  }
  return true;
}
