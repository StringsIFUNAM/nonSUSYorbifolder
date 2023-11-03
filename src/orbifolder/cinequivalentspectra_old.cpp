#include "cinequivalentspectra.h"
#include "cprint.h"



/* ########################################################################################
######   CInequivalentModels()                                                       ######
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
######   Standard constructor of a CInequivalentModels object. The spectra are       ######
######   sorted in the following way:                                                ######
######   "c1" = case 1   : if the highest multiplicity of the (ordered) spectrum is  ######
######                     less than 100                                             ######
######   "c2" = case 2   : otherwise                                                 ######
######   Rank00 - Rank16 : the number of non-Abelian gauge group factors             ######
######   index 0-89      : number of different non-Abelian representations           ######
######################################################################################## */
CInequivalentModels::CInequivalentModels()
{
	vector<CSpectrum> tmp1;
	vector<vector<CSpectrum> > tmp;
	tmp.assign(70, tmp1);						//this initialization pattern is memory consuming
												//alternatively use as list identifier only either Scalars or Fermions
	this->InequivalentModels_c1_Rank00.assign(70, tmp);
	this->InequivalentModels_c1_Rank01.assign(70, tmp);
	this->InequivalentModels_c1_Rank02.assign(70, tmp);
	this->InequivalentModels_c1_Rank03.assign(70, tmp);
	this->InequivalentModels_c1_Rank04.assign(70, tmp);
	this->InequivalentModels_c1_Rank05.assign(70, tmp);
	this->InequivalentModels_c1_Rank06.assign(70, tmp);
	this->InequivalentModels_c1_Rank07.assign(70, tmp);
	this->InequivalentModels_c1_Rank08.assign(70, tmp);
	this->InequivalentModels_c1_Rank09.assign(70, tmp);
	this->InequivalentModels_c1_Rank10.assign(70, tmp);
	this->InequivalentModels_c1_Rank11.assign(70, tmp);
	this->InequivalentModels_c1_Rank12.assign(70, tmp);
	this->InequivalentModels_c1_Rank13.assign(70, tmp);
	this->InequivalentModels_c1_Rank14.assign(70, tmp);
	this->InequivalentModels_c1_Rank15.assign(70, tmp);
	this->InequivalentModels_c1_Rank16.assign(70, tmp);

	this->InequivalentModels_c2_Rank00.assign(70, tmp);
	this->InequivalentModels_c2_Rank01.assign(70, tmp);
	this->InequivalentModels_c2_Rank02.assign(70, tmp);
	this->InequivalentModels_c2_Rank03.assign(70, tmp);
	this->InequivalentModels_c2_Rank04.assign(70, tmp);
	this->InequivalentModels_c2_Rank05.assign(70, tmp);
	this->InequivalentModels_c2_Rank06.assign(70, tmp);
	this->InequivalentModels_c2_Rank07.assign(70, tmp);
	this->InequivalentModels_c2_Rank08.assign(70, tmp);
	this->InequivalentModels_c2_Rank09.assign(70, tmp);
	this->InequivalentModels_c2_Rank10.assign(70, tmp);
	this->InequivalentModels_c2_Rank11.assign(70, tmp);
	this->InequivalentModels_c2_Rank12.assign(70, tmp);
	this->InequivalentModels_c2_Rank13.assign(70, tmp);
	this->InequivalentModels_c2_Rank14.assign(70, tmp);
	this->InequivalentModels_c2_Rank15.assign(70, tmp);
	this->InequivalentModels_c2_Rank16.assign(70, tmp);

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



/* ########################################################################################
######   IsSpectrumUnknown(const CSpectrum &Spectrum, bool AddNewModel)              ######
######                                                                               ######
######   Version: 18.04.2011                                                         ######
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
######   CSpectrum object "Spectrum" and return true if it is unknown. Add "Spectrum" ######
######   if it is new and "AddNewModel" is true.                                     ######
######################################################################################## */
bool CInequivalentModels::IsSpectrumUnknown(const CSpectrum &Spectrum, bool AddNewModel)
{
	//vector<vector<CSpectrum> > *InequivalentModels = NULL;

	//in InequivalentModels_(i,j,k): i classifies according to number of Scalar reps and j according to number of rScalar reps
	//i.e. in the (i,j) component of InequivalentModels the list(index k) has models with i Scalar-reps and j rScalar-reps
	vector<vector<vector<CSpectrum> > >  *InequivalentModels = NULL;

	const unsigned HighestMultiplicity = Spectrum.GetMultiplicity_Matrix()[0][0];				//use Fermi_L to get their highest multiplicities (choice arbitrary)

	switch(Spectrum.GetAlgebra().size())
	{
	case 0:
	{
		if (HighestMultiplicity < 100)
			InequivalentModels = &this->InequivalentModels_c1_Rank00;
		else
			InequivalentModels = &this->InequivalentModels_c2_Rank00;
		break;
	}
	case 1:
	{
		if (HighestMultiplicity < 100)
			InequivalentModels = &this->InequivalentModels_c1_Rank01;
		else
			InequivalentModels = &this->InequivalentModels_c2_Rank01;
		break;
	}
	case 2:
	{
		if (HighestMultiplicity < 100)
			InequivalentModels = &this->InequivalentModels_c1_Rank02;
		else
			InequivalentModels = &this->InequivalentModels_c2_Rank02;
		break;
	}
	case 3:
	{
		if (HighestMultiplicity < 100)
			InequivalentModels = &this->InequivalentModels_c1_Rank03;
		else
			InequivalentModels = &this->InequivalentModels_c2_Rank03;
		break;
	}
	case 4:
	{
		if (HighestMultiplicity < 100)
			InequivalentModels = &this->InequivalentModels_c1_Rank04;
		else
			InequivalentModels = &this->InequivalentModels_c2_Rank04;
		break;
	}
	case 5:
	{
		if (HighestMultiplicity < 100)
			InequivalentModels = &this->InequivalentModels_c1_Rank05;
		else
			InequivalentModels = &this->InequivalentModels_c2_Rank05;
		break;
	}
	case 6:
	{
		if (HighestMultiplicity < 100)
			InequivalentModels = &this->InequivalentModels_c1_Rank06;
		else
			InequivalentModels = &this->InequivalentModels_c2_Rank06;
		break;
	}
	case 7:
	{
		if (HighestMultiplicity < 100)
			InequivalentModels = &this->InequivalentModels_c1_Rank07;
		else
			InequivalentModels = &this->InequivalentModels_c2_Rank07;
		break;
	}
	case 8:
	{
		if (HighestMultiplicity < 100)
			InequivalentModels = &this->InequivalentModels_c1_Rank08;
		else
			InequivalentModels = &this->InequivalentModels_c2_Rank08;
		break;
	}
	case 9:
	{
		if (HighestMultiplicity < 100)
			InequivalentModels = &this->InequivalentModels_c1_Rank09;
		else
			InequivalentModels = &this->InequivalentModels_c2_Rank09;
		break;
	}
	case 10:
	{
		if (HighestMultiplicity < 100)
			InequivalentModels = &this->InequivalentModels_c1_Rank10;
		else
			InequivalentModels = &this->InequivalentModels_c2_Rank10;
		break;
	}
	case 11:
	{
		if (HighestMultiplicity < 100)
			InequivalentModels = &this->InequivalentModels_c1_Rank11;
		else
			InequivalentModels = &this->InequivalentModels_c2_Rank11;
		break;
	}
	case 12:
	{
		if (HighestMultiplicity < 100)
			InequivalentModels = &this->InequivalentModels_c1_Rank12;
		else
			InequivalentModels = &this->InequivalentModels_c2_Rank12;
		break;
	}
	case 13:
	{
		if (HighestMultiplicity < 100)
			InequivalentModels = &this->InequivalentModels_c1_Rank13;
		else
			InequivalentModels = &this->InequivalentModels_c2_Rank13;
		break;
	}
	case 14:
	{
		if (HighestMultiplicity < 100)
			InequivalentModels = &this->InequivalentModels_c1_Rank14;
		else
			InequivalentModels = &this->InequivalentModels_c2_Rank14;
		break;
	}
	case 15:
	{
		if (HighestMultiplicity < 100)
			InequivalentModels = &this->InequivalentModels_c1_Rank15;
		else
			InequivalentModels = &this->InequivalentModels_c2_Rank15;
		break;
	}
	case 16:
	{
		if (HighestMultiplicity < 100)
			InequivalentModels = &this->InequivalentModels_c1_Rank16;
		else
			InequivalentModels = &this->InequivalentModels_c2_Rank16;
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

	const size_t s0 = Spectrum.GetRepresentations_Matrix().size();
	for (unsigned typecounter=0; typecounter<s0; ++typecounter)
	{
		const size_t reps = Spectrum.GetRepresentations_Matrix()[typecounter].size();
		if (reps > 69)
		{
			cout << "Warning in bool CInequivalentModels::IsSpectrumUnknown(...): \"reps\" = " << reps << " out of range. Return false." << endl;
			return false;
		}

	}

	//vector<CSpectrum> &CurrentList = InequivalentModels->at(Spectrum.GetRepresentations_Matrix()[s0-1].size());		//use the Scalar_r as list classification

	//using the Scalars and Fermions as list classification flags, see above
	vector<vector<CSpectrum> > &List_tmp = InequivalentModels->at(Spectrum.GetRepresentations_Matrix()[0].size());
	vector<CSpectrum> &CurrentList=List_tmp.at( Spectrum.GetRepresentations_Matrix()[1].size() );

	const size_t s1 = CurrentList.size();

	unsigned i = 0;

	if (Spectrum.GaugeGroupFromBothE8Factors())
	{
		switch(Spectrum.CompareE8Factors())
		{
		// first E_8 smaller than second one
		case 0:
		{
			CSpectrum SpectrumE8sInterchanged;
			Spectrum.InterchangeE8xE8(SpectrumE8sInterchanged);

			for (i = 0; i < s1; ++i)
			{
				if (SpectrumE8sInterchanged == CurrentList[i])
					return false;
			}

			// spectrum is new
			if (AddNewModel)
			{
				//cout << "interchange E_8 factors" << endl;

				CurrentList.push_back(SpectrumE8sInterchanged);
				++this->ModelCounter;
			}
			return true;
		}
		// first E_8 greater than second one
		case 1:
		{
			break;
		}
		// cannot decide
		case 2:
		{
			CSpectrum SpectrumE8sInterchanged;
			Spectrum.InterchangeE8xE8(SpectrumE8sInterchanged);

			for (i = 0; i < s1; ++i)
			{
				if ((Spectrum == CurrentList[i]) || (SpectrumE8sInterchanged == CurrentList[i]))
					return false;
			}

			// spectrum is new
			if (AddNewModel)
			{
				//cout << "both E_8 factors too similar" << endl;

				CurrentList.push_back(Spectrum);
				++this->ModelCounter;
			}
			return true;
		}
		}
	}

	// if first E_8 greater than second one or spectrum not from both E_8 factors
	for (i = 0; i < s1; ++i)
	{
		if (Spectrum == CurrentList[i])
			return false;
	}
	// spectrum is new
	if (AddNewModel)
	{
		//cout << "do not interchange E_8 factors " << endl;

		CurrentList.push_back(Spectrum);
		++this->ModelCounter;
	}

	return true;
}
