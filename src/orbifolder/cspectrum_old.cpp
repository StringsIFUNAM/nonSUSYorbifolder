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

bool sort_spectrum = true;


bool operator==(const SDimCount &DimCounter1, const SDimCount &DimCounter2)
{
  return ((DimCounter1.Counter == DimCounter2.Counter) && (DimCounter1.Dim == DimCounter2.Dim));
}

bool NumberOfRepsCriterion(const SDimCount &DimCounter1, const SDimCount &DimCounter2);

bool ColumnCriterion1(const SColumnOfSpectrum &Column1, const SColumnOfSpectrum &Column2);
bool ColumnCriterion2(const SColumnOfSpectrum &Column1, const SColumnOfSpectrum &Column2);

bool RowCriterion1(const SRepresentation &rep1, const SRepresentation &rep2);
bool RowCriterion2(const SRepresentation &rep1, const SRepresentation &rep2);



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
######   Standard constructor of a CSpectrum object. No content is specified.        ######
######   The spectrum contains a list of representations (stored in the private      ######
######   member variable "Representations") of SUSY type "Multiplet" with respect to ######
######   the non-Abelian part of the gauge group, together with its multiplicity     ######
######   (stored in "Multiplicity").                                                 ###### 
######################################################################################## */
CSpectrum::CSpectrum()
{
  this->UseSpecialU1 = false;
  this->Multiplet = NOT_DEF_SUSY;
}


/* ########################################################################################
######   CSpectrum(const SConfig &Vacuum, SUSYMultiplet Multiplet)                   ######
######                                                                               ######
######   Version: 22.01.2014                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Vacuum    : collect all fields from this vacuum with SUSY type           ######
######                  "Multiplet" and store them in this CSpectrum object          ######
######   2) Multiplet : the SUSY type (e.g. "LeftChiral") to look for                ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Constructor of a CSpectrum object. Creates the spectrum of fields of SUSY   ######
######   type "Multiplet", in the vacuum "Vacuum".                                   ######
######################################################################################## */
CSpectrum::CSpectrum(const SConfig &Vacuum, SUSYMultiplet Multiplet, string Label, const vector<vector<unsigned> > *SectorsForAdditionalIdentifier)
{
  this->Label = Label;

  unsigned pos_of_U1Hypercharge = 0;
  this->UseSpecialU1 = false;
  if (Vacuum.SymmetryGroup.observable_sector_U1s.size() == 1)
  {
    pos_of_U1Hypercharge = Vacuum.SymmetryGroup.observable_sector_U1s[0];
    this->UseSpecialU1 = (Vacuum.SymmetryGroup.U1s_AdditionalLabels[pos_of_U1Hypercharge] == "Y");
  }

  this->Multiplet = Multiplet;
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
  unsigned l = 0;
  bool rep_equal = true;
  bool not_known = true;

  const bool UseSectorsForAdditionalIdentifier = (SectorsForAdditionalIdentifier != NULL);

  vector<double>   NumberOfRepsPerSector;
  if (UseSectorsForAdditionalIdentifier)
  {
    s2 = SectorsForAdditionalIdentifier->size();
    NumberOfRepsPerSector.assign(s2, 0.0);
  }

  if (this->UseSpecialU1)
  {
    for (i = 0; i < f1; ++i)
    {
      const CField &Field = Fields[i];

      if (Field.Multiplet == this->Multiplet)
      {
        // begin: Additional Identifier
        if (UseSectorsForAdditionalIdentifier)
        {
          k = Field.SGElement.Get_k();
          l = Field.SGElement.Get_l();

          not_known = true;
          for (j = 0; j < s2; ++j)
          {
            if ((SectorsForAdditionalIdentifier->at(j)[0] == k) && (SectorsForAdditionalIdentifier->at(j)[1] == l))
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

          if ((fabs(NewRep.U1Charge - KnownRep.U1Charge) < 0.001) && (NewRep.Representation == KnownRep.Representation))
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

      if (Field.Multiplet == this->Multiplet)
      {
        // begin: Additional Identifier
        if (UseSectorsForAdditionalIdentifier)
        {
          k = Field.SGElement.Get_k();
          l = Field.SGElement.Get_l();

          not_known = true;
          for (j = 0; j < s2; ++j)
          {
            if ((SectorsForAdditionalIdentifier->at(j)[0] == k) && (SectorsForAdditionalIdentifier->at(j)[1] == l))
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

          if (NewRep.Representation == KnownRep.Representation)
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
######   CSpectrum(const SConfig &Vacuum, SUSYMultiplet Multiplet)                   ######
######                                                                               ######
######   Version: 07.04.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Vacuum    : collect all fields from this vacuum with SUSY type           ######
######                  "Multiplet" and store them in this CSpectrum object          ######
######   2) Multiplet : the SUSY type (e.g. "LeftChiral") to look for                ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Constructor of a CSpectrum object. Creates the spectrum of fields of SUSY   ######
######   type "Multiplet", in the vacuum "Vacuum".                                   ######
######################################################################################## */
/*CSpectrum::CSpectrum(const SConfig &Vacuum, SUSYMultiplet Multiplet)
{
	this->Multiplet = Multiplet;

	const vector<CField> &Fields = Vacuum.Fields;
	const size_t f1 = Fields.size();

	const vector<gaugeGroupFactor<double> > &GaugeGroupFactors = Vacuum.SymmetryGroup.GaugeGroup.factor;
	const size_t g1 = GaugeGroupFactors.size();

	string A = "A";
	string D = "D";
	string E = "E";

	unsigned i = 0;
	vector<int> tmp_algebra(2,0);
	for (i = 0; i < g1; ++i)
	{
		const string algebra = GaugeGroupFactors[i].algebra;

		if (algebra[0] == A[0])
			tmp_algebra[0] = 1;
		else
			if (algebra[0] == D[0])
				tmp_algebra[0] = 2;
			else
				if (algebra[0] == E[0])
					tmp_algebra[0] = 3;

		tmp_algebra[1] = GaugeGroupFactors[i].rank;
		this->Algebra.push_back(tmp_algebra);
	}
	this->Position_of_and_in_GaugeGroup = Vacuum.SymmetryGroup.Position_of_and_in_GaugeGroup;

	vector<int> Representation(g1,1);

	size_t s1 = 0;
	unsigned j = 0;
	unsigned k = 0;
	bool rep_equal       = true;
	bool field_not_known = true;

	for (i = 0; i < f1; ++i)
	{
		const CField &Field = Fields[i];

		if (Field.Multiplet == this->Multiplet)
		{
			for (j = 0; j < g1; ++j)
				Representation[j] = Field.Dimensions[j].Dimension;

			field_not_known = true;
			s1 = this->Multiplicity.size();
			for (j = 0; field_not_known && (j < s1); ++j)
			{
				const vector<int> &Representation2 = this->Representations[j];

				rep_equal = true;
				for (k = 0; rep_equal && (k < g1); ++k)
				{
					if (Representation[k] != Representation2[k])
						rep_equal = false;
				}

				if (rep_equal)
				{
					++this->Multiplicity[j];
					field_not_known = false;
				}
			}
			if (field_not_known)
			{
				this->Multiplicity.push_back(1);
				this->Representations.push_back(Representation);
			}
		}
	}
	//if (f1 != 0)						//hacking here!!!
	//	this->ReOrder();				//disable reordering since this case is used only to check for tachyons
}*/



/* ########################################################################################
######   CSpectrum(const SConfig &Vacuum, SUSYMultiplet Multiplet)                   ######
######                                                                               ######
######   Version: 22.01.2014                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Vacuum    : collect all fields from this vacuum with SUSY type           ######
######                  "Multiplet" and store them in this CSpectrum object          ######
######   2) Multiplet : the SUSY type (e.g. "LeftChiral") to look for                ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Constructor of a CSpectrum object. Creates the spectrum of fields of SUSY   ######
######   type "Multiplet", in the vacuum "Vacuum".                                   ######
######################################################################################## */
CSpectrum::CSpectrum(string Label, const vector<vector<unsigned> > *SectorsForAdditionalIdentifier)
{
  this->Label = Label;
  this->Multiplet = NOT_DEF_SUSY;

  unsigned pos_of_U1Hypercharge = 0;
  this->UseSpecialU1 = false;
  if (Vacuum.SymmetryGroup.observable_sector_U1s.size() == 1)
  {
    pos_of_U1Hypercharge = Vacuum.SymmetryGroup.observable_sector_U1s[0];
    this->UseSpecialU1 = (Vacuum.SymmetryGroup.U1s_AdditionalLabels[pos_of_U1Hypercharge] == "Y");
  }

  this->Multiplet = Multiplet;
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
  unsigned l = 0;
  bool rep_equal = true;
  bool not_known = true;

  const bool UseSectorsForAdditionalIdentifier = (SectorsForAdditionalIdentifier != NULL);

  vector<double>   NumberOfRepsPerSector;
  if (UseSectorsForAdditionalIdentifier)
  {
    s2 = SectorsForAdditionalIdentifier->size();
    NumberOfRepsPerSector.assign(s2, 0.0);
  }

  if (this->UseSpecialU1)
  {
    for (i = 0; i < f1; ++i)
    {
      const CField &Field = Fields[i];

      if (Field.Multiplet == Scalar or Field.Multiplet == rScalar)
      {
        // begin: Additional Identifier
        if (UseSectorsForAdditionalIdentifier)
        {
          k = Field.SGElement.Get_k();
          l = Field.SGElement.Get_l();

          not_known = true;
          for (j = 0; j < s2; ++j)
          {
            if ((SectorsForAdditionalIdentifier->at(j)[0] == k) && (SectorsForAdditionalIdentifier->at(j)[1] == l))
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

          if ((fabs(NewRep.U1Charge - KnownRep.U1Charge) < 0.001) && (NewRep.Representation == KnownRep.Representation))
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

      if (Field.Multiplet == Scalar or Field.Multiplet == rScalar)
      {
        // begin: Additional Identifier
        if (UseSectorsForAdditionalIdentifier)
        {
          k = Field.SGElement.Get_k();
          l = Field.SGElement.Get_l();

          not_known = true;
          for (j = 0; j < s2; ++j)
          {
            if ((SectorsForAdditionalIdentifier->at(j)[0] == k) && (SectorsForAdditionalIdentifier->at(j)[1] == l))
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

          if (NewRep.Representation == KnownRep.Representation)
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




/* ########################################################################################		//hacking here!!!
######   CSpectrum(const SConfig &Vacuum, SUSYMultiplet Multiplet)                   ######
######                                                                               ######
######   Version: 07.04.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Vacuum    : collect all fields from this vacuum with SUSY type           ######
######                  "Multiplet" and store them in this CSpectrum object          ######
######   2) Multiplet : the SUSY type (e.g. "LeftChiral") to look for                ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Constructor of a CSpectrum object. Creates the spectrum of fields of SUSY   ######
######   type "Multiplet", in the vacuum "Vacuum". Multpilet types should be given in the same order!    ######
######################################################################################## */
CSpectrum::CSpectrum(const SConfig &Vacuum)
{
	this->Multiplet = NOT_DEF_SUSY;

	const vector<CField> &Fields = Vacuum.Fields;
	const size_t f1 = Fields.size();

	const vector<gaugeGroupFactor<double> > &GaugeGroupFactors = Vacuum.SymmetryGroup.GaugeGroup.factor;
	const size_t g1 = GaugeGroupFactors.size();

	string A = "A";
	string D = "D";
	string E = "E";

	unsigned i = 0;
	vector<int> tmp_algebra(2,0);
	for (i = 0; i < g1; ++i)
	{
		const string algebra = GaugeGroupFactors[i].algebra;

		if (algebra[0] == A[0])
			tmp_algebra[0] = 1;
		else
			if (algebra[0] == D[0])
				tmp_algebra[0] = 2;
			else
				if (algebra[0] == E[0])
					tmp_algebra[0] = 3;

		tmp_algebra[1] = GaugeGroupFactors[i].rank;
		this->Algebra.push_back(tmp_algebra);
	}
	this->Position_of_and_in_GaugeGroup = Vacuum.SymmetryGroup.Position_of_and_in_GaugeGroup;

	vector<int> Representation(g1,1);

	size_t s1 = 0;
	unsigned j = 0;
	unsigned k = 0;
	unsigned typecount=0;
	bool rep_equal       = true;
	bool field_not_known = true;

	//Get Scalar Representations
	this->Multiplicity.clear();
	this->Representations.clear();
	for (i = 0; i < f1; ++i)												//Go through fields of SConfig
	{
		const CField &Field = Fields[i];
		//Merge Scalar and rScalar
		if (Field.Multiplet == Scalar or Field.Multiplet == rScalar)
		{
			for (j = 0; j < g1; ++j)
				Representation[j] = Field.Dimensions[j].Dimension;

			field_not_known = true;
			s1 = this->Multiplicity.size();
			for (j = 0; field_not_known && (j < s1); ++j)
			{
				const vector<int> &Representation2 = this->Representations[j];

				rep_equal = true;
				for (k = 0; rep_equal && (k < g1); ++k)
				{
					if (Representation[k] != Representation2[k])			//merge CPTs for scalars disabled
						rep_equal = false;
				}

				if (rep_equal)
				{
					++this->Multiplicity[j];
					field_not_known = false;
				}
			}
			if (field_not_known)
			{
				this->Multiplicity.push_back(1);
				this->Representations.push_back(Representation);
			}
		}
	}
	Multiplicity_Matrix.push_back(Multiplicity);
	Representations_Matrix.push_back(Representations);

	this->Multiplicity.clear();
	this->Representations.clear();
	for (i = 0; i < f1; ++i)												//Go through fields of SConfig
	{
		const CField &Field = Fields[i];
		//Merge LeftFermi_S and LeftFermi_C
		if (Field.Multiplet == LeftFermi_S or Field.Multiplet == LeftFermi_C)
		{
			for (j = 0; j < g1; ++j)
				Representation[j] = Field.Dimensions[j].Dimension;

			field_not_known = true;
			s1 = this->Multiplicity.size();
			for (j = 0; field_not_known && (j < s1); ++j)
			{
				const vector<int> &Representation2 = this->Representations[j];

				rep_equal = true;
				for (k = 0; rep_equal && (k < g1); ++k)
				{
					if (Representation[k] != Representation2[k])
						rep_equal = false;
				}

				if (rep_equal)
				{
					++this->Multiplicity[j];
					field_not_known = false;
				}
			}
			if (field_not_known)
			{
				this->Multiplicity.push_back(1);
				this->Representations.push_back(Representation);
			}
		}
	}
	Multiplicity_Matrix.push_back(Multiplicity);
	Representations_Matrix.push_back(Representations);

	//if (f1 != 0)
		this->ReOrder();
}



/* ########################################################################################		//hacking here!!!
######   CSpectrum(const SConfig &Vacuum, SUSYMultiplet Multiplet)                   ######
######                                                                               ######
######   Version: 07.04.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Vacuum    : collect all fields from this vacuum with SUSY type           ######
######                  "Multiplet" and store them in this CSpectrum object          ######
######   2) Multiplet : the SUSY type (e.g. "LeftChiral") to look for                ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Constructor of a CSpectrum object. Creates the spectrum of fields of SUSY   ######
######   type "Multiplet", in the vacuum "Vacuum". Multpilet types should be given in the same order!    ######
######################################################################################## */
CSpectrum::CSpectrum(const SConfig &Vacuum1, const SConfig &Vacuum2)
{
	this->Multiplet = NOT_DEF_SUSY;

	const vector<CField> &Fields = Vacuum1.Fields;
	const size_t f1 = Fields.size();

	const vector<gaugeGroupFactor<double> > &GaugeGroupFactors = Vacuum1.SymmetryGroup.GaugeGroup.factor;
	const size_t g1 = GaugeGroupFactors.size();

	string A = "A";
	string D = "D";
	string E = "E";

	unsigned i = 0;
	vector<int> tmp_algebra(2,0);
	for (i = 0; i < g1; ++i)
	{
		const string algebra = GaugeGroupFactors[i].algebra;

		if (algebra[0] == A[0])
			tmp_algebra[0] = 1;
		else
			if (algebra[0] == D[0])
				tmp_algebra[0] = 2;
			else
				if (algebra[0] == E[0])
					tmp_algebra[0] = 3;

		tmp_algebra[1] = GaugeGroupFactors[i].rank;
		this->Algebra.push_back(tmp_algebra);
	}
	this->Position_of_and_in_GaugeGroup = Vacuum1.SymmetryGroup.Position_of_and_in_GaugeGroup;

	vector<int> Representation(g1,1);

	size_t s1 = 0;
	unsigned j = 0;
	unsigned k = 0;
	unsigned typecount=0;
	bool rep_equal       = true;
	bool field_not_known = true;

	//Get Scalar Representations
	this->Multiplicity.clear();
	this->Representations.clear();
	for (i = 0; i < f1; ++i)												//Go through fields of SConfig
	{
		const CField &Field = Fields[i];
		//Merge Scalar and rScalar
		if (Field.Multiplet == Scalar or Field.Multiplet == rScalar)
		{
			for (j = 0; j < g1; ++j)
				Representation[j] = Field.Dimensions[j].Dimension;

			field_not_known = true;
			s1 = this->Multiplicity.size();
			for (j = 0; field_not_known && (j < s1); ++j)
			{
				const vector<int> &Representation2 = this->Representations[j];

				rep_equal = true;
				for (k = 0; rep_equal && (k < g1); ++k)
				{
					if (Representation[k] != Representation2[k])			//merge CPTs for scalars disabled
						rep_equal = false;
				}

				if (rep_equal)
				{
					++this->Multiplicity[j];
					field_not_known = false;
				}
			}
			if (field_not_known)
			{
				this->Multiplicity.push_back(1);
				this->Representations.push_back(Representation);
			}
		}
	}
	Multiplicity_Matrix.push_back(Multiplicity);
	Representations_Matrix.push_back(Representations);

	//Get Fermionic Representations
	this->Multiplicity.clear();
	this->Representations.clear();
	for (i = 0; i < f1; ++i)												//Go through fields of SConfig
	{
		const CField &Field = Fields[i];
		//Merge LeftFermi_S and LeftFermi_C
		if (Field.Multiplet == LeftFermi_S or Field.Multiplet == LeftFermi_C)
		{
			for (j = 0; j < g1; ++j)
				Representation[j] = Field.Dimensions[j].Dimension;

			field_not_known = true;
			s1 = this->Multiplicity.size();
			for (j = 0; field_not_known && (j < s1); ++j)
			{
				const vector<int> &Representation2 = this->Representations[j];

				rep_equal = true;
				for (k = 0; rep_equal && (k < g1); ++k)
				{
					if (Representation[k] != Representation2[k])
						rep_equal = false;
				}

				if (rep_equal)
				{
					++this->Multiplicity[j];
					field_not_known = false;
				}
			}
			if (field_not_known)
			{
				this->Multiplicity.push_back(1);
				this->Representations.push_back(Representation);
			}
		}
	}
	Multiplicity_Matrix.push_back(Multiplicity);
	Representations_Matrix.push_back(Representations);

	//Load tachyonic fields of multiplet type rScalar
	const vector<CField> &tachyonicFields = Vacuum2.Fields;
	const size_t f2 = tachyonicFields.size();
	this->Multiplicity.clear();
	this->Representations.clear();

	for (i = 0; i < f2; ++i)												//Go through fields of SConfig
	{
		const CField &Field = tachyonicFields[i];

		for (j = 0; j < g1; ++j)
			Representation[j] = Field.Dimensions[j].Dimension;

		field_not_known = true;
		s1 = this->Multiplicity.size();
		for (j = 0; field_not_known && (j < s1); ++j)
		{
			const vector<int> &Representation2 = this->Representations[j];

			rep_equal = true;
			for (k = 0; rep_equal && (k < g1); ++k)
			{
				if (Representation[k] != Representation2[k])			//merge CPTs for scalars disabled
					rep_equal = false;
			}

			if (rep_equal)
			{
				++this->Multiplicity[j];
				field_not_known = false;
			}
		}
		if (field_not_known)
		{
			this->Multiplicity.push_back(1);
			this->Representations.push_back(Representation);
		}
	}
	this->Multiplicity_Matrix.push_back(Multiplicity);
	this->Representations_Matrix.push_back(Representations);

	//if (f1 != 0)
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



/* ########################################################################################				//hacking here!!!
######   operator==(const CSpectrum &Spectrum2) const                                ######
######                                                                               ######
######   Version: 20.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   Spectrum2    : a second spectrum to compare with                            ######
######   output:                                                                     ######
######   return value : are the two spectra identical?                               ######
###########################################################################################
######   description:                                                                ######
######   Compares two CSPectrum obejcts. If this operator returns "false" the        ######
######   spectra are inequivalent, but if it returns "true" they still might be      ######
######   inequivalent because the comparison method was not good enough.             ###### 
######################################################################################## */
bool CSpectrum::operator==(const CSpectrum &Spectrum2) const
				{
	size_t s1 = this->AdditionalIdentifier.size();

#ifdef CHECKERROR
	if (this->Multiplet != NOT_DEF_SUSY)
	{
		cout << "\n  Warning in bool CSpectrum::operator==(...) const : Spectra contain different SUSY multiplets. Return false." << endl;
		return false;
	}
	if (s1 != Spectrum2.AdditionalIdentifier.size())
	{
		cout << "\n  Warning in bool CSpectrum::operator==(...) const : Number of additional identifiers differ. Return false." << endl;
		return false;
	}
#endif

	unsigned i = 0;
	for (i = 0; i < s1; ++i)
	{
		if (fabs(this->AdditionalIdentifier[i] - Spectrum2.AdditionalIdentifier[i]) > 0.001)
			return false;
	}

	// number of different representations unequal
	//size_t s0=this->Multiplet_Matrix.size();
	const size_t s0 = this->Multiplicity_Matrix.size();            //Consider only Scalars and Fermions separately
	unsigned typecount=0;
	for (typecount=0; typecount<s0; typecount++)							//Run through various multiplet types
	{
		s1 = this->Multiplicity_Matrix[typecount].size();
		if (s1 != Spectrum2.Multiplicity_Matrix[typecount].size())
			return false;
	}

	// gauge group different from E8xE8 factors
	if (this->Position_of_and_in_GaugeGroup != Spectrum2.Position_of_and_in_GaugeGroup)
		return false;

	// number of non-Abelian gauge group factors unequal
	const size_t g1 = this->Algebra.size();
	if (g1 != Spectrum2.Algebra.size())
		return false;

	// non-Abelian gauge group factors unequal
	for (i = 0; i < g1; ++i)
	{
		if (this->Algebra[i] != Spectrum2.Algebra[i])
			return false;
	}

	// same gauge group but no matter representations
	bool flag=true;														//assume empty spectra
	for (typecount=0; typecount<s0; ++typecount)							//Run through various multiplet types
	{
		s1 = this->Multiplicity_Matrix[typecount].size();
		if (s1 != 0)
			flag=false;
	}
	if (flag)
		return true;

	// multiplicities unequal
	flag=true;																//assume  equal multiplicities
	for (typecount=0; typecount<s0; ++typecount)							//Run through various multiplet types
	{
		if (this->Multiplicity_Matrix[typecount] != Spectrum2.Multiplicity_Matrix[typecount])
			flag=false;
	}
	if (!flag)
		return false;

	unsigned j;
	for (typecount=0; typecount<s0; ++typecount)							//Run through various multiplet types
	{
		s1 = this->Multiplicity_Matrix[typecount].size();

		// assume that the two spectra are exactly the same
		for (i = 0; i < s1; ++i)
		{
			const vector<int> &Rep1 = this->Representations_Matrix[typecount][i];
			const vector<int> &Rep2 = Spectrum2.Representations_Matrix[typecount][i];

			for (j = 0; j < g1; ++j)
			{
				if (abs(Rep1[j]) != abs(Rep2[j]))
					return false;
			}
		}
	}

	return true;
				}



/* ########################################################################################
######   GaugeGroupFromBothE8Factors() const                                         ######
######                                                                               ######
######   Version: 18.04.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : does the non-Abelian gauge group originate from both E_8     ######
######                  factors?                                                     ######
###########################################################################################
######   description:                                                                ######
######   In the case of E8xE8, checks whether the non-Abelian part of the gauge      ######
######   group originates from both E_8 factors. For Spin32 it returns "false".      ######
######################################################################################## */
bool CSpectrum::GaugeGroupFromBothE8Factors() const
{
	const size_t s1 = this->Algebra.size();
	if ((s1 <= 1) || (s1 == this->Position_of_and_in_GaugeGroup) || (this->Position_of_and_in_GaugeGroup == 0) || (this->Position_of_and_in_GaugeGroup == 16))
		return false;

	return true;
}



/* ########################################################################################
######   CompareE8Factors() const                                                    ######
######                                                                               ######
######   Version: 18.04.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : 0 if val1 < val2 i.e. first E8 smaller than second           ######
######                  1 if val1 > val2 i.e. first E8 greater than second           ######
######                  2 if val1 == val2 i.e. cannot decide                         ######
###########################################################################################
######   description:                                                                ######
######   Defines an order between the non-Ablian gauge group factors coming from the ######
######   first and second E_8.                                                       ######
######################################################################################## */
unsigned CSpectrum::CompareE8Factors() const
{
	const size_t s1 = this->Algebra.size();

	unsigned val1 = 0;
	unsigned val2 = 0;

	unsigned i = 0;
	for (i = 0; i < this->Position_of_and_in_GaugeGroup; ++i)
		val1 += this->Algebra[i][0] * (this->Algebra[i][1] * this->Algebra[i][1]);

	for (i = this->Position_of_and_in_GaugeGroup; i < s1; ++i)
		val2 += this->Algebra[i][0] * (this->Algebra[i][1] * this->Algebra[i][1]);

	if (val1 == val2)
		return 2;

	return (val1 > val2);
}



/* ########################################################################################
######   InterchangeE8xE8(CSpectrum &result) const                                   ######
######                                                                               ######
######   Version: 27.05.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   1) result    : a CSpectrum object where the result is stored                ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Creates a CSpectrum object where the two E_8 factors have been interchanged.######
######################################################################################## */
bool CSpectrum::InterchangeE8xE8(CSpectrum &result) const
{
	result.AdditionalIdentifier = this->AdditionalIdentifier;
	result.Multiplicity_Matrix         = this->Multiplicity_Matrix;

	if (!this->GaugeGroupFromBothE8Factors())
	{
		result.Representations_Matrix               = this->Representations_Matrix;
		result.Algebra                       = this->Algebra;
		result.Position_of_and_in_GaugeGroup = this->Position_of_and_in_GaugeGroup;
		cout << "\n  Warning in bool CSpectrum::InterchangeE8xE8(...) const : Do not need to interchange E8xE8 in this case. Return false." << endl;
		return false;
	}

	const size_t s1 = this->Algebra.size();
	vector<int> tmp;
	result.Algebra.assign(s1, tmp);
	result.Position_of_and_in_GaugeGroup = s1 - this->Position_of_and_in_GaugeGroup;

	unsigned i = 0;
	unsigned j = 0;

	vector<int> NewRepVector(s1, 1);
	//const size_t s0 = this->Multiplet_Matrix.size();
	const size_t s0 = this->Multiplicity_Matrix.size();            //Consider only Scalars and Fermions separately

	// begin: interchange in Representations
	for (unsigned typecount=0; typecount<s0; ++typecount)					//Run through different multiplet types
	{
		result.Representations.clear();
		const size_t s2 = this->Representations_Matrix[typecount].size();
		for (i = 0; i < s2; ++i)
		{
			const vector<int> &Representation = this->Representations_Matrix[typecount][i];

			for (j = 0; j < result.Position_of_and_in_GaugeGroup; ++j)
				NewRepVector[j] = Representation[j + this->Position_of_and_in_GaugeGroup];

			for (j = 0; j < this->Position_of_and_in_GaugeGroup; ++j)
				NewRepVector[j + result.Position_of_and_in_GaugeGroup] = Representation[j];

			result.Representations.push_back(NewRepVector);
		}
		result.Representations_Matrix.push_back(result.Representations);
	}
	// end: interchange in Representations

	// begin: interchange the algebra
	for (i = 0; i < result.Position_of_and_in_GaugeGroup; ++i)
		result.Algebra[i] = this->Algebra[i + this->Position_of_and_in_GaugeGroup];

	for (i = 0; i < this->Position_of_and_in_GaugeGroup; ++i)
		result.Algebra[i + result.Position_of_and_in_GaugeGroup] = this->Algebra[i];
	// end: interchange the algebra

	result.ReOrder();

	return true;
}



/* ########################################################################################
######   bool CSpectrum::IsSpectrumEmpty() const                                     ######
######                                                                               ######
######   Version: 18.04.2011                                                         ######
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
	bool flag=true;											//assume spectrum has no reps
	for (int i=0; i<this->Representations_Matrix.size(); ++i)
	{
		if (this->Representations_Matrix[i].size() != 0)
			flag=false;
	}
	if (flag)
		return true;
	else
		return false;
}



/* ########################################################################################
######   bool CSpectrum::IsSpectrumEmpty() const                                     ######
######                                                                               ######
######   Version: 18.04.2011                                                         ######
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
bool CSpectrum::IsTachyonFree() const
{
	if (this->Representations.size() == 0)
		return true;
	else
		return false;
}




/* ########################################################################################
######   bool CSpectrum::ReOrder()                                                   ######
######                                                                               ######
######   Version: 07.04.2011                                                         ######
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
	const size_t g1 = this->Algebra.size();
	if (g1 == 0)
		return true;

	// begin: sort the gauge group factors
	vector<unsigned> GG_factor_from;
	vector<unsigned> GG_factor_to;
	if ((this->Position_of_and_in_GaugeGroup > 0) && (this->Position_of_and_in_GaugeGroup < g1))
	{
		GG_factor_from.push_back(0);
		GG_factor_to.push_back(this->Position_of_and_in_GaugeGroup);

		GG_factor_from.push_back(this->Position_of_and_in_GaugeGroup);
		GG_factor_to.push_back(g1);
	}
	else
	{
		GG_factor_from.push_back(0);
		GG_factor_to.push_back(g1);
	}

	//const size_t s0 = this->Multiplet_Matrix.size();
	const size_t s0 = this->Multiplicity_Matrix.size();            //Consider only Scalars and Fermions separately
	for (unsigned typecount=0; typecount<s0; ++typecount)			//Run through mutliplet types
	{
		const size_t s1 = this->Multiplicity_Matrix[typecount].size();
		if (s1 != 0)
		{
			vector<vector<vector<int> > > table1;

			// in the case of E8 x E8' the gauge group is seperated by the word "and"
			const size_t sh = GG_factor_to.size();
			for (unsigned h = 0; h < sh; ++h)
			{
				const unsigned from = GG_factor_from[h];
				const unsigned to   = GG_factor_to[h];
				const unsigned diff = (unsigned)(to - from);

				// The first entry is the multiplcity, the second one is the representation.
				// The last "diff" entries are some additional multiplicities, which
				// are defined for sorting the gauge groups.

				vector<int>                   tmp1(2 + diff,0);
				vector<vector<int> >          tmp2(s1+1,tmp1);
				vector<vector<vector<int> > > table1_tmp(diff, tmp2);

				// table1 looks like:
				//
				// (D, 5, ...)  (8, 16, ...)  (5, 1, ...)
				// (A, 2, ...)  (8,  1, ...)  (5, 3, ...)
				//
				// which is equivalent to:
				//   D5 + A2 + ...
				// 8 (16, 1, ...)
				// 5 ( 1, 3, ...)

				// runs through the gauge group factors
				for (unsigned i = 0; i < diff; ++i)
				{
					table1_tmp.at(i).at(0) = this->Algebra[from + i];

					for (unsigned j = 0; j < s1; ++j)
					{
						table1_tmp.at(i).at(j+1).at(0) = this->Multiplicity_Matrix[typecount].at(j);
						table1_tmp.at(i).at(j+1).at(1) = this->Representations_Matrix[typecount].at(j).at(from + i);

						for (unsigned k = 0; k < diff; ++k)
						{
							unsigned mult = this->Multiplicity_Matrix[typecount].at(j);
							for (unsigned l = 0; l < g1; ++l)
							{
								if (k != l)
									mult *= abs(this->Representations_Matrix[typecount].at(j).at(l));
							}
							table1_tmp.at(i).at(j+1).at(k+2) = mult;
						}
					}
				}
				// sort table1 using criterion1
				stable_sort(table1_tmp.begin(), table1_tmp.end(), criterion1);

				table1.insert(table1.end(), table1_tmp.begin(), table1_tmp.end());
				// it is important to have a unique sorting of the gauge group
				// which also defines a unique ordering of (for example)
				// A2 + A2 + A2
			}
			// end: sort the gauge group factors


			// table2 looks like the "normal" spectrum, for example:
			// 8 (16, 1)
			// 5 ( 1, 3)

			vector<int> tmp3(g1+1,0);
			vector<vector<int> > table2(s1, tmp3);

			this->Algebra.clear();
			vector<int> algebra_vector(2,0);

			for (unsigned i = 0; i < s1; ++i)
			{
				table2.at(i).at(0) = this->Multiplicity_Matrix[typecount].at(i);

				// read the new order of the gauge group
				if (i == 0)
				{
					for (unsigned j = 0; j < g1; ++j)
					{
						algebra_vector[0] = table1.at(j).at(0).at(0);
						algebra_vector[1] = table1.at(j).at(0).at(1);
						this->Algebra.push_back(algebra_vector);
					}
				}

				for (unsigned j = 0; j < g1; ++j)
					table2.at(i).at(j+1) = table1.at(j).at(i+1).at(1);
			}

			// now sort the lines of the spectrum, using criterion2
			sort(table2.begin(), table2.end(), criterion2);

			this->Multiplicity_Matrix[typecount].clear();
			this->Representations_Matrix[typecount].clear();

			// save the sorted spectrum to "result"
			for (unsigned i = 0; i < s1; ++i)
			{
				this->Multiplicity_Matrix[typecount].push_back(table2.at(i).at(0));

				vector<int> Rep;
				for (unsigned j = 1; j < g1+1; ++j)
					Rep.push_back(table2.at(i).at(j));

				this->Representations_Matrix[typecount].push_back(Rep);
			}
		}
	}
	return true;
}



/* ###########################################################################
######   bool criterion1(...)                                           ######
######                                                                  ######
######   Defines "s1 > s2" for gauge groups                             ######
######   (abs of the dimension of the representation only)              ######
########################################################################### */
bool criterion1(const vector<vector<int> > &s1, const vector<vector<int> > &s2)
{
	const size_t s = s1.size();
	const size_t t = s2.size();

#ifdef CHECKERROR
	if (s != t)
	{
		cout << "\n  Warning in bool criterion1(...). Return false." << endl;
		return false;
	}
#endif

	/* ########################################################################### */
	// compares the names of the gauge groups, i.e. E > D > A
	if (s1.at(0).at(0) > s2.at(0).at(0))
		return true;
	if (s1.at(0).at(0) < s2.at(0).at(0))
		return false;

	/* ########################################################################### */
	// compares the rank of the gauge groups
	if (s1.at(0).at(1) > s2.at(0).at(1))
		return true;
	if (s1.at(0).at(1) < s2.at(0).at(1))
		return false;

	/* ########################################################################### */
	// This criterion uses the dimension of the representations
	// and their multiplicities.
	// It does not depend on the "sign" of the representation.
	unsigned sum1 = 0;
	unsigned sum2 = 0;
	for (unsigned i = 1; i < s; ++i)
	{
		sum1 += (s1.at(i).at(0) * abs(s1.at(i).at(1)));
		sum2 += (s2.at(i).at(0) * abs(s2.at(i).at(1)));
	}

	if (sum1 > sum2)
		return true;
	if (sum1 < sum2)
		return false;

	/* ########################################################################### */
	// The gauge group factor, having the largest representation with the
	// highest multiplicity, is " > " than the other.
	// Therefore, the lines of the corresponding spectra are sorted.
	vector<vector<int> > sorted_s1 = s1;
	vector<vector<int> > sorted_s2 = s2;

	stable_sort((++sorted_s1.begin()), sorted_s1.end(), criterion3);
	stable_sort((++sorted_s2.begin()), sorted_s2.end(), criterion3);

	unsigned max_rep1  = 0;
	unsigned max_rep2  = 0;

	unsigned max_mult1 = 0;
	unsigned max_mult2 = 0;

	for (unsigned i = 1; i < s; ++i)
	{
		const unsigned tmp_rep1 = abs(sorted_s1.at(i).at(1));
		const unsigned tmp_rep2 = abs(sorted_s2.at(i).at(1));

		const unsigned tmp_mult1 = sorted_s1.at(i).at(0);
		const unsigned tmp_mult2 = sorted_s2.at(i).at(0);

		if ( ((tmp_mult1 != tmp_mult2) || (tmp_rep1 != tmp_rep2)) )
			//&& (((tmp_rep1 != 1) && (tmp_rep2 == 1)) || ((tmp_rep1 == 1) && (tmp_rep2 != 1))) )
		{
			if (tmp_rep1 > max_rep1)
			{
				max_rep1  = tmp_rep1;
				max_mult1 = tmp_mult1;
			}
			else
				if ((tmp_rep1 == max_rep1) && (tmp_mult1 > max_mult1))
					max_mult1 = tmp_mult1;

			if (tmp_rep2 > max_rep2)
			{
				max_rep2  = tmp_rep2;
				max_mult2 = tmp_mult2;
			}
			else
				if ((tmp_rep2 == max_rep2) && (tmp_mult2 > max_mult2))
					max_mult2 = tmp_mult2;
		}
	}

	if (max_rep1 > max_rep2)
		return true;

	if (max_rep1 < max_rep2)
		return false;

	if (max_mult1 > max_mult2)
		return true;

	if (max_mult1 < max_mult2)
		return false;


	/* ########################################################################### */
	// If both gauge group are exactly the same, the order does not matter.
	// For example:
	//
	//    A2 + A2
	// 3 ( 3, 3)
	// 3 ( 1, 3)
	// 3 ( 3, 1)

	bool are_rows_equal = true;
	for (unsigned i = 0; are_rows_equal && (i < s); ++i)
	{
		if ((sorted_s1.at(i).at(0) != sorted_s2.at(i).at(0)) || (sorted_s1.at(i).at(1) != sorted_s2.at(i).at(1)))
			are_rows_equal = false;
	}
	if (are_rows_equal)
		return true;

	/* ########################################################################### */
	// The following is also allowed:
	//
	//    A2 + A2
	// 3 ( 3,-3)
	// 3 ( 1, 3)
	// 3 (-3, 1)
	// In this case the ordering is again not important.

	vector<bool>                   max_dim_of_rep(221, true);
	vector<vector<bool> >          max_rank(16, max_dim_of_rep);
	vector<vector<vector<bool> > > Is_Representation_Real(4,max_rank);

	for (unsigned i = 2; i < 16; ++i)
		Is_Representation_Real.at(1).at(i).at(i+1) = false;
	// for example           A     A4    5 is complex

	for (unsigned i = 4; i < 16; ++i)
		Is_Representation_Real.at(1).at(i).at((i+1)*i/2) = false;
	//for example            A     A4    10 is complex
	//for example            A     A8    36 is complex

	for (unsigned i = 6; i < 12; ++i)
		Is_Representation_Real.at(1).at(i).at((i+1)*i*(i-1)/(3*2)) = false;
	//for example            A     A7    56 is complex
	//for example            A     A8    84 is complex

	//D5: 16 is complex
	Is_Representation_Real.at(2).at(5).at(16) = false;

	//E6: 27 is complex
	Is_Representation_Real.at(3).at(6).at(27) = false;

	const unsigned     algebra  = sorted_s1.at(0).at(0);
	const unsigned     rank     = sorted_s1.at(0).at(1);
	const vector<bool> &Is_Real = Is_Representation_Real.at(algebra).at(rank);

	vector<bool> reps_found_s1(s,false);
	vector<bool> reps_found_s2(s,false);

	reps_found_s1[0] = true;
	reps_found_s2[0] = true;

	bool found_rep_i = false;

	// first assume that all representations under these two gauge groups are complex conjugate
	for (unsigned i = 1; i < s; ++i)
	{
		found_rep_i = false;

		const int mult_i = sorted_s1.at(i).at(0);
		const int rep_i  = sorted_s1.at(i).at(1);

		// if the representation is complex
		if (!Is_Real.at(abs(rep_i)))
		{
			for (unsigned j = 1; !found_rep_i && (j < s); ++j)
			{
				// and the cc partner is found
				if (!reps_found_s2.at(j) && (rep_i == -sorted_s2.at(j).at(1)) && (mult_i == sorted_s2.at(j).at(0)))
				{
					reps_found_s1.at(i) = true;
					reps_found_s2.at(j) = true;
					found_rep_i = true;
				}
			}
		}
		// if the representation is real
		else
		{
			for (unsigned j = 1; !found_rep_i && (j < s); ++j)
			{
				// and the partner is found
				if (!reps_found_s2.at(j) && (rep_i == sorted_s2.at(j).at(1)) && (mult_i == sorted_s2.at(j).at(0)))
				{
					reps_found_s1.at(i) = true;
					reps_found_s2.at(j) = true;
					found_rep_i = true;
				}
			}
		}
	}

	vector<bool> found_all(s,true);

	if ((reps_found_s1 == found_all) && (reps_found_s2 == found_all))
		return true;

	/* ########################################################################### */
	// search all representations in "sorted_s1" that are complex, e.g. the 3 of SU(3)
	vector<unsigned> complex_reps_in_sorted_s1;

	for (unsigned i = 1; i < s; ++i)
	{
		const unsigned rep = abs(sorted_s1.at(i).at(1));
		if (!Is_Real.at(rep))
		{
			bool this_rep_is_known = false;

			const size_t cs = complex_reps_in_sorted_s1.size();
			for (unsigned j = 0; !this_rep_is_known && (j < cs); ++j)
			{
				if (complex_reps_in_sorted_s1.at(j) == rep)
					this_rep_is_known = true;
			}
			if (!this_rep_is_known)
				complex_reps_in_sorted_s1.push_back(rep);
		}
	}

	if (complex_reps_in_sorted_s1.size() > 1)
	{
		//const size_t crs = complex_reps_in_sorted_s1.size();
		//cout << "new sorting...\n";

		//for (unsigned i = 0; i < crs; ++i)
		//  cout << complex_reps_in_sorted_s1[i] << " ";
		//cout << endl;

		stable_sort(complex_reps_in_sorted_s1.begin(), complex_reps_in_sorted_s1.end());

		//for (unsigned i = 0; i < crs; ++i)
		//  cout << complex_reps_in_sorted_s1[i] << " ";
		//cout << endl;
	}

	// count the number of e.g. 3 and -3 in sorted_s1 and sorted_s2
	const size_t cs = complex_reps_in_sorted_s1.size();

	unsigned plus_counter_in_sorted_s1 = 0;
	unsigned plus_counter_in_sorted_s2 = 0;

	unsigned minus_counter_in_sorted_s1 = 0;
	unsigned minus_counter_in_sorted_s2 = 0;

	for (unsigned i = 0; i < cs; ++i)
	{
		plus_counter_in_sorted_s1 = 0;
		plus_counter_in_sorted_s2 = 0;

		minus_counter_in_sorted_s1 = 0;
		minus_counter_in_sorted_s2 = 0;

		const int rep = complex_reps_in_sorted_s1.at(i);

		for (unsigned j = 1; j < s; ++j)
		{
			if (sorted_s1.at(j).at(1) == rep)
				plus_counter_in_sorted_s1 += sorted_s1.at(j).at(0);
			else
				if (sorted_s1.at(j).at(1) == -rep)
					minus_counter_in_sorted_s1 += sorted_s1.at(j).at(0);

			if (sorted_s2.at(j).at(1) == rep)
				plus_counter_in_sorted_s2 += sorted_s1.at(j).at(0);
			else
				if (sorted_s2.at(j).at(1) == -rep)
					minus_counter_in_sorted_s2 += sorted_s1.at(j).at(0);
		}

		unsigned diff_in_sorted_s1 = abs(int(plus_counter_in_sorted_s1 - minus_counter_in_sorted_s1));
		unsigned diff_in_sorted_s2 = abs(int(plus_counter_in_sorted_s2 - minus_counter_in_sorted_s2));

		if (diff_in_sorted_s1 > diff_in_sorted_s2)
			return true;
		if (diff_in_sorted_s1 < diff_in_sorted_s2)
			return false;
	}

	/* ########################################################################### */
	// same idea, but this time using the additional multiplicities
	const size_t ms = sorted_s1.at(0).size();

	for (unsigned k = 2; k < ms; ++k)
	{
		for (unsigned i = 0; i < cs; ++i)
		{
			unsigned plus_counter_in_sorted_s1 = 0;
			unsigned plus_counter_in_sorted_s2 = 0;

			unsigned minus_counter_in_sorted_s1 = 0;
			unsigned minus_counter_in_sorted_s2 = 0;

			int rep = complex_reps_in_sorted_s1.at(i);

			for (unsigned j = 1; j < s; ++j)
			{
				if (sorted_s1.at(j).at(1) == rep)
					plus_counter_in_sorted_s1 += sorted_s1.at(j).at(k);

				if (sorted_s1.at(j).at(1) == -rep)
					minus_counter_in_sorted_s1 += sorted_s1.at(j).at(k);

				if (sorted_s2.at(j).at(1) == rep)
					plus_counter_in_sorted_s2 += sorted_s1.at(j).at(k);

				if (sorted_s2.at(j).at(1) == -rep)
					minus_counter_in_sorted_s2 += sorted_s1.at(j).at(k);
			}

			unsigned diff_in_sorted_s1 = abs(int(plus_counter_in_sorted_s1 - minus_counter_in_sorted_s1));
			unsigned diff_in_sorted_s2 = abs(int(plus_counter_in_sorted_s2 - minus_counter_in_sorted_s2));

			if (diff_in_sorted_s1 > diff_in_sorted_s2)
				return true;
			if (diff_in_sorted_s1 < diff_in_sorted_s2)
				return false;
		}
	}

	/* ########################################################################### */
	// This criterion uses the dimension of the representations
	// and the additional multiplicities.
	// It does not depend on the "sign" of the representation.
	const size_t qs = s1.at(0).size();

	unsigned sum_q1 = 0;
	unsigned sum_q2 = 0;
	for (unsigned q = 2; q < qs; ++q)
	{
		sum_q1 = 0;
		sum_q2 = 0;

		for (unsigned i = 1; i < s; ++i)
		{
			sum_q1 += (s1.at(i).at(1) * abs(s1.at(i).at(q)));
			sum_q2 += (s2.at(i).at(1) * abs(s2.at(i).at(q)));
		}

		if (sum_q1 > sum_q2)
			return true;
		if (sum_q1 < sum_q2)
			return false;
	}

	cout << "\nWarning 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
	cout << "\ns1 = ";
	for (unsigned i = 0; i < s; ++i)
	{
		cout << "(" << sorted_s1.at(i).at(0) << " x " << sorted_s1.at(i).at(1) << "; ";
		for (unsigned j = 2; j < ms-1; ++j)
			cout << sorted_s1.at(i).at(j) << ", ";
		if ((int)(ms)-1 >= 0)
			cout << sorted_s1.at(i).at(ms-1);
		cout << ")\n";
	}

	cout << endl;

	cout << "\ns2 = ";
	for (unsigned i = 0; i < s; ++i)
	{
		cout << "(" << sorted_s2.at(i).at(0) << " x " << sorted_s2.at(i).at(1) << "; ";
		for (unsigned j = 2; j < ms-1; ++j)
			cout << sorted_s2.at(i).at(j) << ", ";
		if ((int)(ms)-1 >= 0)
			cout << sorted_s2.at(i).at(ms-1);
		cout << ")\n";
	}

	cout << "\n(first entry is the gauge group)\n";

	return false;
}


/* ###########################################################################
######   bool criterion2(...)                                           ######
######                                                                  ######
######   Defines " > " for the lines of a spectrum                      ######
########################################################################### */
bool criterion2(const vector<int> &s1, const vector<int> &s2)
{
	const size_t s = s1.size();
	const size_t t = s2.size();

#ifdef CHECKERROR
	if (s != t)
	{
		cout << "\n  Warning in bool criterion2(const vector<int> &s1, const vector<int> &s2). Return false." << endl;
		return false;
	}
#endif

	// compares the multiplicity
	if (s1.at(0) > s2.at(0))
		return true;
	if (s1.at(0) < s2.at(0))
		return false;

	// compares the number of particles corresponding to the representation
	// for example 1(16,3,1) => 1 * 16 * 3 * 1 = 48
	unsigned mult1 = 1;
	unsigned mult2 = 1;
	for (unsigned i = 0; i < s; ++i)
	{
		mult1 *= abs(s1.at(i));
		mult2 *= abs(s2.at(i));
	}
	if (mult1 > mult2)
		return true;
	if (mult1 < mult2)
		return false;

	// which line of the spectrum has the first non-trivial entry
	// that is larger than the one of the other line
	for (unsigned i = 1; i < s; ++i)
	{
		if (abs(s1.at(i)) > abs(s2.at(i)))
			return true;
		if (abs(s1.at(i)) < abs(s2.at(i)))
			return false;
	}

	if (s1 == s2)
		return true;

	// if two lines of the spectrum are the same except for complex conjugation
	// then for example (16,3,1) is bigger than (-16,3,1)
	bool equal_except_for_complex_conjugation = true;

	for (unsigned i = 1; equal_except_for_complex_conjugation && (i < s); ++i)
	{
		if (abs(s1.at(i)) != abs(s2.at(i)))
			equal_except_for_complex_conjugation = false;
	}
	if (equal_except_for_complex_conjugation)
	{
		for (unsigned i = 1; i < s; ++i)
		{
			if (s1.at(i) > s2.at(i))
				return true;
			if (s1.at(i) < s2.at(i))
				return false;
		}
	}

	cout << "!!!!!!!!!!!!!!!!!!! WARNUNG Crit 2 !!!!!!!!!!!!!!!!!!!!" << endl;
	return false;
}



/* ###########################################################################
######   bool criterion3(...)                                           ######
########################################################################### */
bool criterion3(const vector<int> &s3_1, const vector<int> &s3_2)
{
#ifdef CHECKERROR
	if (s3_1.size() != s3_2.size())
	{
		cout << "\n  Warning in bool criterion3(const vector<int> &s1, const vector<int> &s2).\ns1.size() = " << s3_1.size() << " and s2.size() = " << s3_2.size() << ". Return false." << endl;
		return false;;
	}
#endif

	// compares the multiplicities
	if (s3_1[0] > s3_2[0])
		return true;
	if (s3_1[0] < s3_2[0])
		return false;

	// compares the size of the representations
	if (abs(s3_1[1]) > abs(s3_2[1]))
		return true;
	if (abs(s3_1[1]) < abs(s3_2[1]))
		return false;

	// compares the dimension of the representations
	if (s3_1[1] > s3_2[1])
		return true;
	if (s3_1[1] < s3_2[1])
		return false;

	return true;
}
