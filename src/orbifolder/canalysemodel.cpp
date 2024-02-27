#include "canalysemodel.h" 
#include "corbifold.h"
#include "clinalg.h"
#include "globalfunctions.h"
#include "cprint.h"
#include "cfield.h"


CAnalyseModel::CAnalyseModel()
{
}


CAnalyseModel::~CAnalyseModel()
{
}


/* ########################################################################################
######   AreSomeFixedBranesEmpty(const COrbifold &Orbifold, ...) const               ######
######                                                                               ######
######   Version: 19.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Orbifold  : the orbifold of the vev-config "VEVConfig"                   ######
######   2) VEVConfig : contains the massless fields                                 ######
######   output:                                                                     ######
######   return value : are there some fixed points without massless twisted strings?######
###########################################################################################
######   description:                                                                ######
######   Goes through the list of fixed points of "Orbifold" and checks whether one  ######
######   of them is empty, i.e. has no massless twisted strings.                     ######
######################################################################################## */
bool CAnalyseModel::AreSomeFixedBranesEmpty(const COrbifold &Orbifold, const SConfig &VEVConfig) const
{
	if (VEVConfig.InvariantSupercharges.size() != 1)
	{
		cout << "\n  Warning: AreSomeFixedBranesEmpty(...) only works for N=1 SUSY." << endl;
		return false;
	}

	const vector<CField> &Fields = VEVConfig.Fields;

	unsigned i = 0;
	unsigned j = 0;
	size_t s2 = 0;

	vector<unsigned> FieldIndices;

	const size_t s1 = Orbifold.GetNumberOfSectors();
	for(i = 0; i < s1; ++i)
	{
		const CSector &Sector = Orbifold.GetSector(i);

		if (Sector.SectorHasLeftchiralRightmover())
		{
			s2 = Sector.GetNumberOfFixedBranes();
			for(j = 0; j < s2; ++j)
			{
				FieldIndices.clear();
				Sector.GetFixedBrane(j).GetFieldIndices(Fields, LeftChiral, FieldIndices);
				if (FieldIndices.size() == 0)
					return true;
			}
		}
	}
	return false;
}



/* ########################################################################################
######   ComputeN2BetaFunctionCoefficient(const COrbifold &Orbifold, ...) const      ######
######                                                                               ######
######   Version: 29.02.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Orbifold                 : find N=2 subsectors of this orbifold and      ######
######                                 compute their beta function coefficients      ######
######   2) VEVConfig                : contains the gauge group factors whose beta   ######
######                                 function coefficients shall be computed       ######
######   3) GaugeIndices             : contains the quadratic and cubic indices of   ######
######                                 non-Abelian representations                   ######
######   4) N2_Sector_Exists         : vector of 3 booleans, specifying whether the  ######
######                                 i-th torus yields a N=2 subsector             ######
######   5) BetaFunctionCoefficients : contains three (one per T^2 torus) vectors of ######
######                                 beta function coefficients (one per gauge     ######
######                                 group factor)                                 ######
######   output:                                                                     ######
######   return value                : finished succesfully?                         ######
###########################################################################################
######   description:                                                                ######
######   Find N=2 subsectors of "Orbifold" and compute their beta function           ######
######   coefficients (multiplied by |P_j|/|P|) for all gauge group factors of       ######
######   "VEVConfig".                                                                ######
######################################################################################## */
bool CAnalyseModel::ComputeN2BetaFunctionCoefficient(const COrbifold *Orbifold, const SConfig &VEVConfig, const CGaugeIndices &GaugeIndices, vector<bool> &N2_Sector_Exists, vector<vector<rational<int> > > &BetaFunctionCoefficients) const
{
	N2_Sector_Exists.assign(3, false);
	BetaFunctionCoefficients.clear();

	const COrbifoldGroup  &OrbifoldGroup = Orbifold->OrbifoldGroup;
	const vector<CSector> &Sectors       = Orbifold->GetSectors();
	const bool            &ZMxZN         = OrbifoldGroup.GetSpaceGroup().IsZMxZN();

	const CGaugeGroup &GaugeGroup       = VEVConfig.SymmetryGroup.GaugeGroup;
	const size_t      number_of_factors = GaugeGroup.factor.size();

	CPrint Print(Tstandard, &cout);

	const int P = Sectors.size();
	int Pi = 0;

	unsigned j = 0;
	unsigned k = 0;

	vector<unsigned> first_Sector_with_Ti_fixed;
	vector<unsigned> TwistFactors2;

	vector<rational<int> > pre_bs;
	vector<rational<int> > ith_bs;

	for (unsigned ith_torus = 1; ith_torus < 4; ++ith_torus)
	{
		ith_bs.assign(number_of_factors, 0);
		first_Sector_with_Ti_fixed.clear();

		for (j = 1; j < P; ++j)
		{
			const CTwistVector &Twist = Sectors[j].GetTwist();
			if (is_integer(Twist[ith_torus]))
			{
				first_Sector_with_Ti_fixed.push_back(Sectors[j].Get_m());			//hacking here!!! not good for ZMxZNxZK
				if (ZMxZN)
					first_Sector_with_Ti_fixed.push_back(Sectors[j].Get_n());
				break;
			}
		}

		// create the suborbifolds
		if (first_Sector_with_Ti_fixed.size() != 0)
		{
			N2_Sector_Exists[ith_torus-1] = true;

			COrbifoldGroup SubOrb_OG_Ti_fixed;
			if (!OrbifoldGroup.CreateSubGroup(Print, first_Sector_with_Ti_fixed, TwistFactors2, true, SubOrb_OG_Ti_fixed))
			{
				cout << "\n  Warning in bool CAnalyseModel::ComputeN2BetaFunctionCoefficient(...) const : Cannot compute suborbifold. Return false." << endl;
				return false;
			}
			COrbifold SubOrbifold_Ti_fixed(SubOrb_OG_Ti_fixed);

			Pi = SubOrb_OG_Ti_fixed.GetOrderZM();
			if ((SubOrb_OG_Ti_fixed.GetOrderZN() != 0) && (SubOrb_OG_Ti_fixed.GetOrderZN() != 1))
				Pi *= SubOrb_OG_Ti_fixed.GetOrderZN();

			pre_bs.clear();

			const size_t number_of_factors_Ti_fixed = SubOrbifold_Ti_fixed.StandardConfig.SymmetryGroup.GaugeGroup.factor.size();

			for (j = 0; j < number_of_factors_Ti_fixed; ++j)
			{
				rational<int> tmp = this->ComputeBetaFunctionCoefficient(GaugeIndices, SubOrbifold_Ti_fixed.StandardConfig, j);
				pre_bs.push_back(tmp * rational<int>(Pi,P));
			}

			for (j = 0; j < number_of_factors; ++j)
			{
				const gaugeGroupFactor<double> &ggfA = GaugeGroup.factor[j];
				for (k = 0; k < number_of_factors_Ti_fixed; ++k)
				{
					if (this->GroupAFromB(ggfA, SubOrbifold_Ti_fixed.StandardConfig.SymmetryGroup.GaugeGroup.factor[k]))
						ith_bs[j] = pre_bs[k];
				}
			}
		}
		BetaFunctionCoefficients.push_back(ith_bs);
	}
	return true;
}



/* ########################################################################################
######   ComputeBetaFunctionCoefficient(const SConfig &VEVConfig, ...) const         ######
######                                                                               ######
######   Version: 29.02.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) GaugeIndices : contains the quadratic and cubic indices of non-Abelian   ######
######                     representations                                           ######
######   2) VEVConfig    : contains the gauge group and the fields                   ######
######   3) factor       : specifies which gauge group factor to use                 ######
######   output:                                                                     ######
######   return value    : the beta function coefficient                             ######
###########################################################################################
######   description:                                                                ######
######   Compute the beta function coefficient of the "factor"-th gauge group factor ######
######   and return the result.                                                      ######
######################################################################################## */
rational<int> CAnalyseModel::ComputeBetaFunctionCoefficient(const CGaugeIndices &GaugeIndices, const SConfig &VEVConfig, const unsigned factor) const
{
	unsigned Index = 0;
	rational<int> beta = 0;
	rational<int> Factor = 0;
    rational<int> Factors = 0;
	
	unsigned k = 0;
	const size_t number_of_factors = VEVConfig.SymmetryGroup.GaugeGroup.factor.size();
	if (factor >= number_of_factors)
	{
		cout << "\n  Warning in bool CAnalyseModel::ComputeBetaFunctionCoefficient(...) : \"factor\" out of range. Return 0." << endl;
		return 0;
	}
	
    rational<int> adj_factor(11,3);
    rational<int> facf(2,3);
    rational<int> facs(1,3);
    	
    vector<SUSYMultiplet> Multiplets(2);				    //Particle types to be printed
	Multiplets[0]=Scalar;									//and to be given for equivalence check
	Multiplets[1]=LeftFermi;

	const vector<CField>           &Fields = VEVConfig.Fields;
	const gaugeGroupFactor<double> &ggf    = VEVConfig.SymmetryGroup.GaugeGroup.factor[factor];

	if (!GaugeIndices.GetQuadraticIndexAdj(ggf, Index))
		return 0;

	beta = -1 * (adj_factor * rational<int>(Index,2));

	// begin: run through all left-chiral fields
	for( vector<CField>::const_iterator it_field = Fields.begin(); it_field != Fields.end(); ++it_field)
	{
		// left-chiral multiplet
		if (it_field->Multiplet == Multiplets[1])
		{
			const RepVector &Dimensions = it_field->Dimensions;

			// begin: compute the charges under the modular symmetries
			if (!GaugeIndices.GetQuadraticIndex(ggf, abs(Dimensions[factor].Dimension), Index))
				return 0;

			if (Index != 0)
			{
				Factor = Index;

				for (k = 0; k < number_of_factors; ++k)
				{
					if (k != factor)
						Factor *= abs(Dimensions[k].Dimension);
				}
				Factor = Factor*rational<int>(1,2);
				beta = beta + facf*Factor;
			}	
			// end: compute the charges under the modular symmetries
		}
		
		// left-chiral multiplet
		if (it_field->Multiplet == Multiplets[0])
		{
			const RepVector &Dimensions = it_field->Dimensions;

			// begin: compute the charges under the modular symmetries
			if (!GaugeIndices.GetQuadraticIndex(ggf, abs(Dimensions[factor].Dimension), Index))
				return 0;

			if (Index != 0)
			{
				Factors = Index;

				for (k = 0; k < number_of_factors; ++k)
				{
					if (k != factor)
						Factors *= abs(Dimensions[k].Dimension);
				}
				Factors = Factors*rational<int>(1,2);
				beta = beta + facs*Factors;
			}	
			// end: compute the charges under the modular symmetries
		}
						
	} 
	// end: run through all left-chiral fields

  
    return beta;  

}

bool CAnalyseModel::GroupAFromB(const gaugeGroupFactor<double> &ggfA, const gaugeGroupFactor<double> &ggfB) const
{
	const vector<vector<double> > &ggf_SimpleRoots_A = ggfA.simpleroots;
	const vector<vector<double> > &ggf_SimpleRoots_B = ggfB.simpleroots;

	const size_t sA = ggf_SimpleRoots_A.size();
	const size_t sB = ggf_SimpleRoots_B.size();

	vector<bool> AfromB(sA, false);

	unsigned i = 0;
	unsigned j = 0;

	for (i = 0; i < sA; ++i)
	{
		const vector<double> &ggf_SimpleRoot_A = ggf_SimpleRoots_A[i];

		for (j = 0; j < sB; ++j)
		{
			const vector<double> &ggf_SimpleRoot_B = ggf_SimpleRoots_B[j];
			if (fabs(ggf_SimpleRoot_A * ggf_SimpleRoot_B) >= 0.0001)
				AfromB[i] = true;
		}
	}
	if (find(AfromB.begin(), AfromB.end(), false) == AfromB.end())
		return true;

	if (find(AfromB.begin(), AfromB.end(), true) == AfromB.end())
		return false;

	cout << "\n  Warning in bool CAnalyseModel::GroupAFromB(...) : Some roots of A originate from B, some do not. Return false." << endl;
	return false;
}


/* ########################################################################################
######   AccidentalU1Charges_Load(SConfig &VEVConfig, ifstream &in)                  ######
######                                                                               ######
######   Version: 18.05.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) VEVConfig : contains the fields                                          ######
######   2) in        : ifstream object containing the accidental U(1) charges       ######
######   output:                                                                     ######
######   return value : accidental U(1) charges loaded succesfully?                  ######
###########################################################################################
######   description:                                                                ######
######   Load accidental U(1) charges from the file associated to "in" and store     ######
######   them in the vev-config "VEVConfig".                                         ######
######################################################################################## */
bool CAnalyseModel::AccidentalU1Charges_Load(SConfig &VEVConfig, ifstream &in)
{
	const size_t f1 = VEVConfig.Fields.size();
	if ((f1 == 0) || (!in.is_open()) || (!in.good()))
		return false;

	string currentline = "";
	string str1 = "-1234567890";

	string tmp_string1 = "";

	string::size_type n1 = 0;
	string::size_type n2 = 0;

	string::size_type loc1 = 0;
	string::size_type loc2 = 0;

	unsigned j = 0;

	bool AccU1ChargesLoad = false;

	vector<rational<CHugeInt> > AccU1Charges;
	while (getline(in, currentline))
	{
		if (currentline.find("end") == string::npos)
		{
			AccU1Charges.clear();

			loc1 = 0;
			loc2 = 0;
			while (loc1 != string::npos)
			{
				loc1 = currentline.find(" ", loc2);
				if (loc1 != string::npos)
				{
					tmp_string1 = currentline.substr(loc2, loc1 - loc2);
					loc2 = loc1 + 1;

					n1 = tmp_string1.find_first_of(str1, 0);
					n2 = tmp_string1.find("/", n1);

					if ((n1 != string::npos) && (n2 != string::npos))
					{
						CHugeInt num(tmp_string1.substr(n1, n2-n1));
						CHugeInt den(tmp_string1.substr(n2-n1+1, string::npos));
						if (!num.Check() || !den.Check())
						{
							cout << "\n  Warning in bool CAnalyseModel::AccidentalU1Charges_Load(...) : Problem with CHugeInt. Return false." << endl;
							return false;
						}
						AccU1Charges.push_back(rational<CHugeInt>(num,den));
					}
					else
					{
						cout << "\n  Warning in bool CAnalyseModel::AccidentalU1Charges_Load(...) : CHugeInt object not found. Return false." << endl;
						return false;
					}
				}
			}

			if (AccU1Charges.size() != f1)
			{
				cout << "\n  Warning in bool CAnalyseModel::AccidentalU1Charges_Load(...): number of charges differs from number of fields." << endl;
				return false;
			}

			AccU1ChargesLoad = true;
			for (j = 0; j < f1; ++j)
				VEVConfig.Fields[j].AccU1Charges.push_back(AccU1Charges[j]);
		}
	}
	return AccU1ChargesLoad;
}



/* ########################################################################################
######   AccidentalU1Charges_Save(SConfig &VEVConfig, ostream &out)                  ######
######                                                                               ######
######   Version: 16.04.2009                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) VEVConfig : contains the fields                                          ######
######   2) in        : ostream object containing the file where the charges are     ######
######                  saved                                                        ######
######   output:                                                                     ######
######   return value : accidental U(1) charges saved succesfully?                   ######
###########################################################################################
######   description:                                                                ######
######   Save accidental U(1) charges from the vev-config "VEVConfig" to the file    ######
######   associated to "out".                                                        ######
######################################################################################## */
bool CAnalyseModel::AccidentalU1Charges_Save(const SConfig &VEVConfig, ostream &out)
{
	const size_t f1 = VEVConfig.Fields.size();
	if (f1 == 0)
		return true;

	unsigned i = 0;
	unsigned j = 0;

	unsigned Number_Of_AccSymmetries = 0;
	for (j = 0; j < f1; ++j)
	{
		const CField &Field = VEVConfig.Fields[j];
		if (Field.Multiplet == LeftChiral)
		{
			Number_Of_AccSymmetries = Field.AccU1Charges.size();
			break;
		}
	}

	if (Number_Of_AccSymmetries == 0)
		return false;

	for (i = 0; i < Number_Of_AccSymmetries; ++i)
	{
		for (j = 0; j < f1; ++j)
		{
			const CField &Field = VEVConfig.Fields[j];
			if (Field.Multiplet == LeftChiral)
			{
				if (!Field.AccU1Charges[i].numerator().Check() || !Field.AccU1Charges[i].denominator().Check())
				{
					cout << "\n  Warning in bool CAnalyseModel::AccidentalU1Charges_Save(...) : Problem with CHugeInt. Return false." << endl;
					return false;
				}
				out << Field.AccU1Charges[i] << " ";
			}
			else
				out << "0/1 ";
		}
		out << endl;
	}
	return true;
}



/* ########################################################################################
######   CreatePhenoScheme(const string &SchemeLabel, PhenoScheme &Scheme) const     ######
######                                                                               ######
######   Version: 15.04.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) SchemeLabel : can be "SU5", "SO10", "SM" for Standard Model or "PS" for  ######
######                    Pati-Salam                                                 ######
######   2) Scheme      : store the PhenoScheme here                                 ######
######   output:                                                                     ######
######   return value   : PhenoScheme created succesfully?                           ######
###########################################################################################
######   description:                                                                ######
######   A PhenoScheme contains the gauge group, matter representations and their    ######
######   labels. For example for SU(5), there can be:                                ######
######     10-plets called "T",                                                      ######
######      5-plets called "F",                                                      ######
######     \bar{10}-plets called "bT",                                               ######
######      \bar{5}-plets called "bF" and                                            ######
######      singlets called "N".                                                      ######
######################################################################################## */
bool CAnalyseModel::CreatePhenoScheme(const string &SchemeLabel, PhenoScheme &Scheme) const
{
	Scheme.GaugeGroupFactors.clear();
	Scheme.SetOfDimensions.clear();
	Scheme.SetOfU1Charges.clear();
	Scheme.Labels.clear();

	SDimension tmp;
	tmp.Dimension = 1;

	if (SchemeLabel == "SU5")
	{
		Scheme.SchemeLabel = SchemeLabel;
		Scheme.IndexOfSetToNormalizeCharges = 1;

		Scheme.GaugeGroupFactors.push_back("A4");
		Scheme.GGs_AdditionalLabels.push_back("");
		CVector NoU1Charges;

		tmp.AdditionalLabel = "";
		RepVector Dimensions(1, tmp);

		Dimensions[0].Dimension = 1;
		Scheme.SetOfDimensions.push_back(Dimensions);
		Scheme.SetOfU1Charges.push_back(NoU1Charges);
		Scheme.Labels.push_back("N");
		Scheme.Multiplets.push_back(Scalar);

		Dimensions[0].Dimension = 10;
		Scheme.SetOfDimensions.push_back(Dimensions);
		Scheme.SetOfU1Charges.push_back(NoU1Charges);
		Scheme.Labels.push_back("T");
		Scheme.Multiplets.push_back(LeftFermi);

		Dimensions[0].Dimension = -10;
		Scheme.SetOfDimensions.push_back(Dimensions);
		Scheme.SetOfU1Charges.push_back(NoU1Charges);
		Scheme.Labels.push_back("bT");
		Scheme.Multiplets.push_back(LeftFermi); 

		Dimensions[0].Dimension = -5;
		Scheme.SetOfDimensions.push_back(Dimensions);
		Scheme.SetOfU1Charges.push_back(NoU1Charges);
		Scheme.Labels.push_back("F");
		Scheme.Multiplets.push_back(LeftFermi);  
		Scheme.SetOfDimensions.push_back(Dimensions); 
		Scheme.SetOfU1Charges.push_back(NoU1Charges);
		Scheme.Labels.push_back("H"); 
		Scheme.Multiplets.push_back(Scalar);  

		Dimensions[0].Dimension = 5;
		Scheme.SetOfDimensions.push_back(Dimensions);
		Scheme.SetOfU1Charges.push_back(NoU1Charges);
		Scheme.Labels.push_back("bF");
		Scheme.Multiplets.push_back(LeftFermi);  
		Scheme.SetOfDimensions.push_back(Dimensions);
		Scheme.SetOfU1Charges.push_back(NoU1Charges);
		Scheme.Labels.push_back("bH");
		Scheme.Multiplets.push_back(Scalar);  


		return true;
	}
	if (SchemeLabel == "SO10")
	{
		Scheme.SchemeLabel = SchemeLabel;
		Scheme.IndexOfSetToNormalizeCharges = 1;

		Scheme.GaugeGroupFactors.push_back("D5");
		Scheme.GGs_AdditionalLabels.push_back("");
		CVector NoU1Charges;

		tmp.AdditionalLabel = "";
		RepVector Dimensions(1, tmp);

		Dimensions[0].Dimension = 1;
		Scheme.SetOfDimensions.push_back(Dimensions);
		Scheme.SetOfU1Charges.push_back(NoU1Charges);
		Scheme.Labels.push_back("N");
		Scheme.Multiplets.push_back(Scalar); 

		Dimensions[0].Dimension = 16;
		Scheme.SetOfDimensions.push_back(Dimensions);
		Scheme.SetOfU1Charges.push_back(NoU1Charges);
		Scheme.Labels.push_back("F");
		Scheme.Multiplets.push_back(LeftFermi); 

		Dimensions[0].Dimension = -16;
		Scheme.SetOfDimensions.push_back(Dimensions);
		Scheme.SetOfU1Charges.push_back(NoU1Charges);
		Scheme.Labels.push_back("bF");
		Scheme.Multiplets.push_back(LeftFermi);

		Dimensions[0].Dimension = 10;
		Scheme.SetOfDimensions.push_back(Dimensions);
		Scheme.SetOfU1Charges.push_back(NoU1Charges);
		Scheme.Labels.push_back("H");
		Scheme.Multiplets.push_back(Scalar);

		return true;
	}
	if (SchemeLabel == "SM")
	{
		Scheme.SchemeLabel = SchemeLabel;
		Scheme.IndexOfSetToNormalizeCharges = 0;

		Scheme.NormalizationOfLengthOfU1Generators.push_back(5.0/6.0);

		Scheme.GaugeGroupFactors.push_back("A2");
		Scheme.GaugeGroupFactors.push_back("A1");
		Scheme.GGs_AdditionalLabels.push_back("C");
		Scheme.GGs_AdditionalLabels.push_back("L");
		CVector YCharges(1);

		tmp.AdditionalLabel = "";
		RepVector Dimensions(2, tmp);

		Dimensions[0].Dimension =  3;
		Dimensions[1].Dimension =  2;
		YCharges[0] =  1.0/6.0;
		Scheme.SetOfDimensions.push_back(Dimensions);
		Scheme.SetOfU1Charges.push_back(YCharges);
		Scheme.Labels.push_back("q");
		Scheme.Multiplets.push_back(LeftFermi);  
		Scheme.SetOfDimensions.push_back(Dimensions); 
		Scheme.SetOfU1Charges.push_back(YCharges); 
		Scheme.Labels.push_back("sq");            
		Scheme.Multiplets.push_back(Scalar);      

		Dimensions[0].Dimension = -3;
		Dimensions[1].Dimension =  2;
		YCharges[0] = -1.0/6.0;
		Scheme.SetOfDimensions.push_back(Dimensions);
		Scheme.SetOfU1Charges.push_back(YCharges);
		Scheme.Labels.push_back("bq");
		Scheme.Multiplets.push_back(LeftFermi);  
		Scheme.SetOfDimensions.push_back(Dimensions); 
		Scheme.SetOfU1Charges.push_back(YCharges); 
		Scheme.Labels.push_back("bsq");            
		Scheme.Multiplets.push_back(Scalar);      

		Dimensions[0].Dimension = -3;
		Dimensions[1].Dimension =  1;
		YCharges[0] = -2.0/3.0;
		Scheme.SetOfDimensions.push_back(Dimensions);
		Scheme.SetOfU1Charges.push_back(YCharges);
		Scheme.Labels.push_back("bu");
		Scheme.Multiplets.push_back(LeftFermi);  
		Scheme.SetOfDimensions.push_back(Dimensions); 
		Scheme.SetOfU1Charges.push_back(YCharges); 
		Scheme.Labels.push_back("bsu");            
		Scheme.Multiplets.push_back(Scalar);      

		Dimensions[0].Dimension =  3;
		Dimensions[1].Dimension =  1;
		YCharges[0] =  2.0/3.0;
		Scheme.SetOfDimensions.push_back(Dimensions);
		Scheme.SetOfU1Charges.push_back(YCharges);
		Scheme.Labels.push_back("u");
		Scheme.Multiplets.push_back(LeftFermi);  
		Scheme.SetOfDimensions.push_back(Dimensions);
		Scheme.SetOfU1Charges.push_back(YCharges);
		Scheme.Labels.push_back("su");            
		Scheme.Multiplets.push_back(Scalar);      

		Dimensions[0].Dimension = -3;
		Dimensions[1].Dimension =  1;
		YCharges[0] =  1.0/3.0;
		Scheme.SetOfDimensions.push_back(Dimensions);
		Scheme.SetOfU1Charges.push_back(YCharges);
		Scheme.Labels.push_back("bd");
		Scheme.Multiplets.push_back(LeftFermi);  
		Scheme.SetOfDimensions.push_back(Dimensions); 
		Scheme.SetOfU1Charges.push_back(YCharges); 
		Scheme.Labels.push_back("bsd");            
		Scheme.Multiplets.push_back(Scalar);      


		Dimensions[0].Dimension =  3;
		Dimensions[1].Dimension =  1;
		YCharges[0] = -1.0/3.0;
		Scheme.SetOfDimensions.push_back(Dimensions);
		Scheme.SetOfU1Charges.push_back(YCharges);
		Scheme.Labels.push_back("d");
		Scheme.Multiplets.push_back(LeftFermi);  
		Scheme.SetOfDimensions.push_back(Dimensions); 
		Scheme.SetOfU1Charges.push_back(YCharges); 
		Scheme.Labels.push_back("sd");            
		Scheme.Multiplets.push_back(Scalar);      

		Dimensions[0].Dimension =  1;
		Dimensions[1].Dimension =  2;
		YCharges[0] = -1.0/2.0;
		Scheme.SetOfDimensions.push_back(Dimensions);
		Scheme.SetOfU1Charges.push_back(YCharges);
		Scheme.Labels.push_back("l");
		Scheme.Multiplets.push_back(LeftFermi);  
		Scheme.SetOfDimensions.push_back(Dimensions); 
		Scheme.SetOfU1Charges.push_back(YCharges); 
		Scheme.Labels.push_back("h");            
		Scheme.Multiplets.push_back(Scalar);     

		Dimensions[0].Dimension =  1;
		Dimensions[1].Dimension =  2;
		YCharges[0] =  1.0/2.0;
		Scheme.SetOfDimensions.push_back(Dimensions);
		Scheme.SetOfU1Charges.push_back(YCharges);
		Scheme.Labels.push_back("bl");
		Scheme.Multiplets.push_back(LeftFermi);  
		Scheme.SetOfDimensions.push_back(Dimensions); 
		Scheme.SetOfU1Charges.push_back(YCharges); 
		Scheme.Labels.push_back("bh");      
		Scheme.Multiplets.push_back(Scalar);  

		Dimensions[0].Dimension =  1;
		Dimensions[1].Dimension =  1;
		YCharges[0] =  1.0;
		Scheme.SetOfDimensions.push_back(Dimensions);
		Scheme.SetOfU1Charges.push_back(YCharges);
		Scheme.Labels.push_back("be");
		Scheme.Multiplets.push_back(LeftFermi);  
		Scheme.SetOfDimensions.push_back(Dimensions); 
		Scheme.SetOfU1Charges.push_back(YCharges); 
		Scheme.Labels.push_back("bse");            
		Scheme.Multiplets.push_back(Scalar);      

		Dimensions[0].Dimension =  1;
		Dimensions[1].Dimension =  1;
		YCharges[0] = -1.0;
		Scheme.SetOfDimensions.push_back(Dimensions);
		Scheme.SetOfU1Charges.push_back(YCharges);
		Scheme.Labels.push_back("e");
		Scheme.Multiplets.push_back(LeftFermi);  
		Scheme.SetOfDimensions.push_back(Dimensions); 
		Scheme.SetOfU1Charges.push_back(YCharges); 
		Scheme.Labels.push_back("se");            
		Scheme.Multiplets.push_back(Scalar);      

		Dimensions[0].Dimension =  1;
		Dimensions[1].Dimension =  1;
		YCharges[0] =  0.0;
		Scheme.SetOfDimensions.push_back(Dimensions);
		Scheme.SetOfU1Charges.push_back(YCharges);
		Scheme.Labels.push_back("n");
		Scheme.Multiplets.push_back(LeftFermi);  
		Scheme.SetOfDimensions.push_back(Dimensions); 
		Scheme.SetOfU1Charges.push_back(YCharges); 
		Scheme.Labels.push_back("sn");            
		Scheme.Multiplets.push_back(Scalar);      

		return true;
	}
	if (SchemeLabel == "PS")
	{
		Scheme.SchemeLabel = SchemeLabel;
		Scheme.IndexOfSetToNormalizeCharges = 0;

		Scheme.GaugeGroupFactors.push_back("A3");
		Scheme.GaugeGroupFactors.push_back("A1");
		Scheme.GaugeGroupFactors.push_back("A1");
		Scheme.GGs_AdditionalLabels.push_back("C");
		Scheme.GGs_AdditionalLabels.push_back("L");
		Scheme.GGs_AdditionalLabels.push_back("R");
		CVector Charges;

		tmp.AdditionalLabel = "";
		RepVector Dimensions(3, tmp);

		Dimensions[0].Dimension =  4;
		Dimensions[1].Dimension =  2;
		Dimensions[2].Dimension =  1;
		Scheme.SetOfDimensions.push_back(Dimensions);
		Scheme.SetOfU1Charges.push_back(Charges);
		Scheme.Labels.push_back("f");
		Scheme.Multiplets.push_back(LeftFermi);

		Dimensions[0].Dimension = -4;
		Dimensions[1].Dimension =  2;
		Dimensions[2].Dimension =  1;
		Scheme.SetOfDimensions.push_back(Dimensions);
		Scheme.SetOfU1Charges.push_back(Charges);
		Scheme.Labels.push_back("bf");
		Scheme.Multiplets.push_back(LeftFermi); 

		Dimensions[0].Dimension = -4;
		Dimensions[1].Dimension =  1;
		Dimensions[2].Dimension =  2;
		Scheme.SetOfDimensions.push_back(Dimensions);
		Scheme.SetOfU1Charges.push_back(Charges);
		Scheme.Labels.push_back("fc");
		Scheme.Multiplets.push_back(LeftFermi);

		Dimensions[0].Dimension =  4;
		Dimensions[1].Dimension =  1;
		Dimensions[2].Dimension =  2;
		Scheme.SetOfDimensions.push_back(Dimensions);
		Scheme.SetOfU1Charges.push_back(Charges);
		Scheme.Labels.push_back("bfc");
		Scheme.Multiplets.push_back(LeftFermi);

		Dimensions[0].Dimension =  4;
		Dimensions[1].Dimension =  1;
		Dimensions[2].Dimension =  1;
		Scheme.SetOfDimensions.push_back(Dimensions);
		Scheme.SetOfU1Charges.push_back(Charges);
		Scheme.Labels.push_back("q");
		Scheme.Multiplets.push_back(LeftFermi); 

		Dimensions[0].Dimension = -4;
		Dimensions[1].Dimension =  1;
		Dimensions[2].Dimension =  1;
		Scheme.SetOfDimensions.push_back(Dimensions);
		Scheme.SetOfU1Charges.push_back(Charges);
		Scheme.Labels.push_back("bq");
		Scheme.Multiplets.push_back(LeftFermi);

		Dimensions[0].Dimension =  6;
		Dimensions[1].Dimension =  1;
		Dimensions[2].Dimension =  1;
		Scheme.SetOfDimensions.push_back(Dimensions);
		Scheme.SetOfU1Charges.push_back(Charges);
		Scheme.Labels.push_back("c");
		Scheme.Multiplets.push_back(LeftFermi); 

		Dimensions[0].Dimension =  1;
		Dimensions[1].Dimension =  2;
		Dimensions[2].Dimension =  2;
		Scheme.SetOfDimensions.push_back(Dimensions);
		Scheme.SetOfU1Charges.push_back(Charges);
		Scheme.Labels.push_back("h");
		Scheme.Multiplets.push_back(LeftFermi); 

		Dimensions[0].Dimension =  1;
		Dimensions[1].Dimension =  2;
		Dimensions[2].Dimension =  1;
		Scheme.SetOfDimensions.push_back(Dimensions);
		Scheme.SetOfU1Charges.push_back(Charges);
		Scheme.Labels.push_back("dl");
		Scheme.Multiplets.push_back(LeftFermi); 

		Dimensions[0].Dimension =  1;
		Dimensions[1].Dimension =  1;
		Dimensions[2].Dimension =  2;
		Scheme.SetOfDimensions.push_back(Dimensions);
		Scheme.SetOfU1Charges.push_back(Charges);
		Scheme.Labels.push_back("dr");
		Scheme.Multiplets.push_back(LeftFermi);

		Dimensions[0].Dimension =  1;
		Dimensions[1].Dimension =  1;
		Dimensions[2].Dimension =  1;
		Scheme.SetOfDimensions.push_back(Dimensions);
		Scheme.SetOfU1Charges.push_back(Charges);
		Scheme.Labels.push_back("s");
		Scheme.Multiplets.push_back(Scalar); 

		return true;
	}

	return false;
}



/* ########################################################################################
######   FindPositionOfGaugeGroup(const PhenoScheme &Scheme, ...) const              ######
######                                                                               ######
######   Version: 17.05.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Scheme            : a PhenoScheme, e.g. "SM" for Standard Model          ######
######   2) OriginalVEVConfig : vev-config to analyze                                ######
######   3) GoodVEVConfigs    : resulting vev-configs                                ######
######   output:                                                                     ######
######   return value         : does the gauge group contain the one of "Scheme"?    ######
###########################################################################################
######   description:                                                                ######
######   Identifies whether the scheme gauge group of "Scheme" is inside the one of  ######
######   "OriginalVEVConfig". The scheme gauge group must be inside one of the E8.   ######
######################################################################################## */
bool CAnalyseModel::FindPositionOfGaugeGroup(const PhenoScheme &Scheme, const SConfig &OriginalVEVConfig, vector<SConfig> &GoodVEVConfigs) const
{
	const CGaugeGroup                       &GaugeGroup     = OriginalVEVConfig.SymmetryGroup.GaugeGroup;
	const vector<gaugeGroupFactor<double> > &factors        = GaugeGroup.factor;

	const size_t number_of_factors = factors.size();

	size_t s1 = OriginalVEVConfig.SymmetryGroup.Position_of_and_in_GaugeGroup;
	if (s1 > number_of_factors)
	{
		cout << "\n  Warning in bool CAnalyseModel::FindPositionOfGaugeGroup(...) const: Position_of_and_in_GaugeGroup out of range. Return false." << endl;
		return false;
	}

	const vector<string> &SchemeGaugeGroupFactors = Scheme.GaugeGroupFactors;
	const size_t s2 = SchemeGaugeGroupFactors.size();

	unsigned i = 0;
	unsigned j = 0;

	vector<vector<unsigned> > positions_of_one_factor_part1(s2);
	vector<vector<unsigned> > positions_of_one_factor_part2(s2);

	vector<unsigned>          PositionOfUnequalDigits(2,0);
	vector<vector<unsigned> > PositionsOfUnequalDigits;

	// begin: find the gauge group of the scheme inside E8 x E8 / SO(32)
	for (i = 0; i < s2; ++i)
	{
		const string &SchemeFactor = SchemeGaugeGroupFactors[i];

		for (j = 0; j < i; ++j)
		{
			const string &SchemeFactor2 = SchemeGaugeGroupFactors[j];
			if (SchemeFactor == SchemeFactor2)
			{
				PositionOfUnequalDigits[0] = j;
				PositionOfUnequalDigits[1] = i;
				PositionsOfUnequalDigits.push_back(PositionOfUnequalDigits);
			}
		}

		vector<unsigned> &positions_part1 = positions_of_one_factor_part1[i];
		vector<unsigned> &positions_part2 = positions_of_one_factor_part2[i];

		for (j = 0; j < number_of_factors; ++j)
		{
			const gaugeGroupFactor<double> &ggf = factors[j];

			if (j < s1)
			{
				if (ggf.algebra == SchemeFactor)
					positions_part1.push_back(j);
			}
			else
			{
				if (ggf.algebra == SchemeFactor)
					positions_part2.push_back(j);
			}
		}
	}
	// end: find the gauge group of the scheme inside E8 x E8 / SO(32)

	/*cout << "\n\nbegin PositionsOfUnequalDigits\n";
  for (i = 0; i < PositionsOfUnequalDigits.size(); ++i)
    cout << PositionsOfUnequalDigits[i][0] << " and " << PositionsOfUnequalDigits[i][1] << endl;
  cout << "end PositionsOfUnequalDigits\n" << endl;

  cout << "begin positions_of_one_factor_part1\n";
  for (i = 0; i < s2; ++i)
  {
    vector<unsigned> &positions_part1 = positions_of_one_factor_part1[i];
    for (j = 0; j < positions_part1.size(); ++j)
      cout << positions_part1[j];
    cout << endl;
  }
  cout << "end positions_of_one_factor_part1\n" << endl;

  cout << "begin positions_of_one_factor_part2\n";
  for (i = 0; i < s2; ++i)
  {
    vector<unsigned> &positions_part2 = positions_of_one_factor_part2[i];
    for (j = 0; j < positions_part2.size(); ++j)
      cout << positions_part2[j];
    cout << endl;
  }
  cout << "end positions_of_one_factor_part2\n" << endl;

  cout << "can it work?" << endl;*/
	bool part1 = true;
	bool part2 = true;
	for (i = 0; part1 && (i < s2); ++i)
	{
		if (positions_of_one_factor_part1[i].size() == 0)
			part1 = false;
	}
	for (i = 0; part2 && (i < s2); ++i)
	{
		if (positions_of_one_factor_part2[i].size() == 0)
			part2 = false;
	}
	if (!part1 && !part2)
	{
		//cout << "no\n" << endl;
		return false;
	}
	//cout << "yes" << endl;

	unsigned counter = 1;

	SConfig NewVEVConfig;
	NewVEVConfig = OriginalVEVConfig;
	NewVEVConfig.SymmetryGroup.observable_sector_GGs.clear();
	NewVEVConfig.SymmetryGroup.observable_sector_U1s.clear();
	NewVEVConfig.SymmetryGroup.GGs_AdditionalLabels.assign(number_of_factors, "");
	NewVEVConfig.SymmetryGroup.U1s_AdditionalLabels.assign(GaugeGroup.u1directions.size(), "");

	vector<unsigned> Number(s2, 0);
	vector<unsigned> MaxDigits(s2, 0);

	bool positions_ok = true;
	const size_t s3 = PositionsOfUnequalDigits.size();
	if (part1)
	{
		for (i = 0; i < s2; ++i)
			MaxDigits[i] = positions_of_one_factor_part1[i].size();

		//cout << "MaxDigits " << MaxDigits << endl;
		Number.assign(s2, 0);
		do
		{
			//cout << "Number " << Number << endl;
			positions_ok = true;
			for (i = 0; positions_ok && (i < s3); ++i)
			{
				const vector<unsigned> &Positions = PositionsOfUnequalDigits[i];
				if (Number[Positions[0]] >= Number[Positions[1]])
					positions_ok = false;
			}
			if (positions_ok)
			{
				//cout << "ok ... " << flush;
				ostringstream os;
				os << counter;
				++counter;

				NewVEVConfig.ConfigLabel = Scheme.SchemeLabel + " VEVConfig " + os.str();
				NewVEVConfig.SymmetryGroup.observable_sector_GGs.clear();
				for (i = 0; i < s2; ++i)
				{
					const vector<unsigned> &positions_part1 = positions_of_one_factor_part1[i];
					NewVEVConfig.SymmetryGroup.observable_sector_GGs.push_back(positions_part1[Number[i]]);
					NewVEVConfig.SymmetryGroup.GGs_AdditionalLabels[positions_part1[Number[i]]] = Scheme.GGs_AdditionalLabels[i];
				}
				GoodVEVConfigs.push_back(NewVEVConfig);
				//cout << "good vauum" << endl;
			}
		} while (NextNumberNoTwins(Number, MaxDigits, PositionsOfUnequalDigits, 0));
	}

	if (part2)
	{
		for (i = 0; i < s2; ++i)
			MaxDigits[i] = positions_of_one_factor_part2[i].size();

		//cout << "MaxDigits " << MaxDigits << endl;
		Number.assign(s2, 0);
		do
		{
			//cout << "Number " << Number << endl;
			positions_ok = true;
			for (i = 0; positions_ok && (i < s3); ++i)
			{
				const vector<unsigned> &Positions = PositionsOfUnequalDigits[i];
				if (Number[Positions[0]] >= Number[Positions[1]])
					positions_ok = false;
			}
			if (positions_ok)
			{
				//cout << "ok ... " << flush;
				ostringstream os;
				os << counter;
				++counter;

				NewVEVConfig.ConfigLabel = Scheme.SchemeLabel + " VEVConfig " + os.str();
				NewVEVConfig.SymmetryGroup.observable_sector_GGs.clear();
				for (i = 0; i < s2; ++i)
				{
					const vector<unsigned> &positions_part2 = positions_of_one_factor_part2[i];
					NewVEVConfig.SymmetryGroup.observable_sector_GGs.push_back(positions_part2[Number[i]]);
					NewVEVConfig.SymmetryGroup.GGs_AdditionalLabels[positions_part2[Number[i]]] = Scheme.GGs_AdditionalLabels[i];
				}
				GoodVEVConfigs.push_back(NewVEVConfig);
				//cout << "good vauum" << endl;
			}
		} while (NextNumberNoTwins(Number, MaxDigits, PositionsOfUnequalDigits, 0));
	}

	return true;
}


/* ########################################################################################
######   DTerms_FindDMonomials(const COrbifold &Orbifold, ...) const                 ######
######                                                                               ######
######   Version: 11.05.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Orbifold    : the orbifold of the vev-config "VEVConfig"                 ######
######   2) VEVConfig   : list of D-monomials will be saved here                     ######
######   3) FieldIndices: list of field indices that may be contained in D-monomial  ######
######   4) D0_withFI   : shall a Fayet Iliopoulos D-term for D_0 be used ?          ######
######                                                                               ######
######   output:                                                                     ######
######   5) GaugeEquivalentFields: field indices with the same charges               ######
######   6) AllFieldsInMonomial  : summary of field indices collected from all       ######
######                             D-monomials                                       ######
###########################################################################################
######   description:                                                                ######
######   Identifies gauge invaraint monomials in the fields of "VEVConfig" with      ######
######   indices "FieldIndices". These monomials correspond to D=0.                  ######
######################################################################################## */
bool CAnalyseModel::DTerms_FindDMonomials(const COrbifold &Orbifold, SConfig &VEVConfig, const vector<unsigned> &FieldIndices, bool D0_withFI, vector<vector<unsigned> > &GaugeEquivalentFields, vector<unsigned> &AllFieldsInMonomial) const
{
	if ((AllFieldsInMonomial.size() != 0) || (GaugeEquivalentFields.size() != 0))
	{
		cout << "\n  Warning in bool CAnalyseModel::DTerms_FindDMonomials(..): one of the paramters is not empty - now cleared." << endl;
		AllFieldsInMonomial.clear();
		GaugeEquivalentFields.clear();
	}

	if (D0_withFI && !VEVConfig.SymmetryGroup.IsFirstU1Anomalous)
	{
		cout << "\n  Warning in bool CAnalyseModel::DTerms_FindDMonomials(..): anomalous U(1) not defined. Return false." << endl;
		return false;
	}

	unsigned i = 0;
	unsigned j = 0;

	const size_t number_of_U1s       = VEVConfig.SymmetryGroup.GaugeGroup.u1directions.size();
	size_t       number_of_used_U1s  = number_of_U1s;
	unsigned     shift               = 0;
	if (D0_withFI)
	{
		shift = 1;
		number_of_used_U1s -= 1;
	}

	const rational<CHugeInt> Zero;

	vector<rational<CHugeInt> >          Charges(number_of_used_U1s, Zero);
	vector<vector<rational<CHugeInt> > > ChargeMatrixTranspose;

	unsigned FieldIndex = 0;

	// begin: collect the U(1) charges of those fields that shall be involved in D-flat monomials
	//        but only use gauge-inequivalent fields
	{
		vector<RepVector> Known_Representations;
		vector<unsigned>  tmp_GaugeEquivalentFields;
		bool field_is_unknown = true;

		// run through the field indices
		const size_t s0 = FieldIndices.size();

		for (i = 0; i < s0; ++i)
		{
			FieldIndex = FieldIndices[i];
			const CField &tmp_Field = VEVConfig.Fields[FieldIndex];

			const RepVector &Representation = tmp_Field.Dimensions;
			const CVector   &U1Charges      = tmp_Field.U1Charges;

			// convert the U(1) charges to ratioanl of CHugeInt
			for (j = 0; j < number_of_used_U1s; ++j)
				Charges[j] = D2HugeInt(U1Charges[j+shift]);

			// begin: are the charges of the current field known?
			field_is_unknown = true;
			const size_t s1 = ChargeMatrixTranspose.size();
			for (j = 0; field_is_unknown && (j < s1); ++j)
			{
				if ((ChargeMatrixTranspose[j] == Charges) && AreRepVectorsEqual(Known_Representations[j], Representation))
				{
					field_is_unknown = false;
					GaugeEquivalentFields[j].push_back(FieldIndex);
				}
			}
			if (field_is_unknown)
			{
				// save the charges
				Known_Representations.push_back(Representation);
				ChargeMatrixTranspose.push_back(Charges);

				// save the field index
				tmp_GaugeEquivalentFields.assign(1,FieldIndex);
				GaugeEquivalentFields.push_back(tmp_GaugeEquivalentFields);
			}
			// end: are the charges of the current field known?
		}
	}
	// end: collect the U(1) charges of those fields that shall be involved in D-flat monomials
	//      but only use gauge-inequivalent fields

	const size_t s1 = ChargeMatrixTranspose.size();

	// begin: transpose the charge matrix
	vector<rational<CHugeInt> > lineA(s1, Zero);
	vector<vector<rational<CHugeInt> > > ChargeMatrix(number_of_used_U1s,lineA);

	for (i = 0; i < s1; ++i)
	{
		for (j = 0; j < number_of_used_U1s; ++j)
			ChargeMatrix[j][i] = ChargeMatrixTranspose[i][j];
	}
	// end: transpose the charge matrix

	// begin: find the (positive integer) kernel of the charge matrix
	CLinAlg<CHugeInt> LA;
	vector<vector<rational<CHugeInt> > > BasisOfKernel;

	if (!LA.FindPositiveIntegerKernel(ChargeMatrix, BasisOfKernel))
	{
		cout << "\n  Warning in bool CAnalyseModel::DTerms_FindDMonomials(...): Could not create the positive, integer kernel. Return false." << endl;
		return false;
	}
	// end: find the (positive integer) kernel of the charge matrix

	const size_t s2 = BasisOfKernel.size();

	// begin: check and save the result
	{
		unsigned k = 0;
		CVector U1Sum(number_of_U1s);
		long long int num = 0;

		for (i = 0; i < s2; ++i)
		{
			const vector<rational<CHugeInt> > &Vector_of_Exponents = BasisOfKernel[i];

			// begin: save monomials
			CMonomial Monomial(number_of_U1s);
			U1Sum.Assign(number_of_U1s);
			num = 0;

			for (j = 0; j < Vector_of_Exponents.size(); ++j)
			{
				if (Vector_of_Exponents[j].denominator().ToLongLongInt() != 1)
				{
					cout << "\n  Warning in bool CAnalyseModel::DTerms_FindDMonomials(...): Exponent is not integer. Return false." << endl;
					return false;
				}

				num = Vector_of_Exponents[j].numerator().ToLongLongInt();
				if (num != 0)
				{
					Monomial.GaugeEquivalentFields.push_back(GaugeEquivalentFields[j]);
					Monomial.Exponents.push_back(num);
					U1Sum += VEVConfig.Fields[GaugeEquivalentFields[j][0]].U1Charges * (double)num;
				}
			}

			Monomial.U1Charges = U1Sum;
			if (Monomial.CheckGaugeInvariance(Orbifold, VEVConfig, false, true))
				VEVConfig.SetOfMonomials.push_back(Monomial);
			// end: save monomials

			// begin: collect all fields that are involved in the D-flat directions
			for (j = 0; j < Vector_of_Exponents.size(); ++j)
			{
				if (Vector_of_Exponents[j] != Zero)
				{
					const vector<unsigned> &GaugeEquivalentFields_j = GaugeEquivalentFields[j];
					for (k = 0; k < GaugeEquivalentFields_j.size(); ++k)
					{
						unsigned tmp = GaugeEquivalentFields_j[k];

						if (find(AllFieldsInMonomial.begin(), AllFieldsInMonomial.end(), tmp) == AllFieldsInMonomial.end())
							AllFieldsInMonomial.push_back(tmp);
					}
				}
			}
			// end: collect all fields that are involved in the D-flat directions
		}
	}
	// end: check and save the result

	return true;
}



/* ########################################################################################
######   AnalyseModel(const COrbifold &Orbifold, ...) const                          ######
######                                                                               ######
######   Version: 21.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Orbifold               : the orbifold of the vev-config                  ######
######                               "OriginalVEVConfig"                             ######
######   2) OriginalVEVConfig      : analyze this vev-config for SM, PS and/or SU5   ######
######   3) SM                     : call function with SM=true to look for Standard ######
######                               Model; function sets SM=true/false depending    ######
######                               whether "OriginalVEVConfig" contains a Standard ######
######                               Model                                           ######
######   4) PS                     : call function with PS=true to look for Pati-    ######
######                               Salam, see 3)                                   ######
######   5) SU5                    : call function with SU5=true to look for SU(5),  ######
######                               see 3)                                          ######
######   6) AllVEVConfigs          : the resulting vev-configs are stored here       ######
######   7) Print                  : print the output to this CPrint object          ######
######   8) NumberOfGenerations    : look for models with this number of generations ######
######   9) SM_PrintSU5SimpleRoots : for SM=true, print the simple roots of the      ######
######                               SU(5) used to identify hypercharge              ######
######   output:                                                                     ######
######   return value              : phenomenologically interesting vev-config found?######
###########################################################################################
######   description:                                                                ######
######   Analyzes the vev-config "OriginalVEVConfig" for "SM", "PS" and "SU5" and    ######
######   store the result in "AllVEVConfigs"                                         ######
######################################################################################## */
bool CAnalyseModel::AnalyseModel(const COrbifold &Orbifold, const SConfig &OriginalVEVConfig, bool &SM, bool &PS, bool &SU5, vector<SConfig> &AllVEVConfigs, CPrint &Print, unsigned NumberOfGenerations, bool SM_PrintSU5SimpleRoots) const
{
	const SelfDualLattice Lattice = Orbifold.OrbifoldGroup.GetLattice();

	unsigned i = 0;
	unsigned j = 0;

	size_t s1 = 0;
	size_t s2 = 0;

	vector<CVector> Hypercharges;

	// begin: create the E8 x E8' or SO(32) gauge group
	const CVector Null_Shift(16);
	S_OscillatorExcitation Excitation;
	Excitation.NumberOperator = 0.0;
	Excitation.ZeroPointEnergy = 2.0;
	CMasslessHalfState Roots10D(LeftMover, Excitation);
	Roots10D.SolveMassEquation(Null_Shift, Lattice);

	if (Roots10D.Weights.size() != 480)
	{
		cout << "\n  Warning in bool CAnalyseModel::AnalyseModel(...) const : check the 10D gauge group. Return false." << endl;
		return false;
	}
	// end: create the E8 x E8' or SO(32) gauge group

	vector<SConfig> TestVEVConfigs;
	bool Good_VEVConfig_Found = false;

	unsigned pos_of_U1 = 0;
	if (OriginalVEVConfig.SymmetryGroup.IsFirstU1Anomalous)
		pos_of_U1 = 1;

	if (SM)
	{
		SM = false;

		vector<vector<CVector> > PossibleRootsOfSU5;

		PhenoScheme SMScheme;
		this->CreatePhenoScheme("SM", SMScheme);

		unsigned SMcounter = 1;

		if (this->FindPositionOfGaugeGroup(SMScheme, OriginalVEVConfig, TestVEVConfigs))
		{
			//cout<<"SM gauge group found!"<<endl;
			s1 = TestVEVConfigs.size();

			for (i = 0; i < s1; ++i)
			{
				SConfig &SMVEVConfig = TestVEVConfigs[i];

				if (abs(this->SM_NetNumberOfGenerations(SMVEVConfig)) == NumberOfGenerations)
				{
					//cout<<"Three families found!"<<endl;
					Hypercharges.clear();
					PossibleRootsOfSU5.clear();
					if (this->SM_GetHypercharges(SMVEVConfig, Hypercharges, PossibleRootsOfSU5, Roots10D.Weights))
					{
						s2 = Hypercharges.size();
						for (j = 0; j < s2; ++j)
						{
							Orbifold.Config_SetU1Direction(SMVEVConfig, Hypercharges[j], pos_of_U1);

							SMVEVConfig.SymmetryGroup.observable_sector_U1s.clear();
							SMVEVConfig.SymmetryGroup.observable_sector_U1s.push_back(pos_of_U1);
							SMVEVConfig.SymmetryGroup.U1s_AdditionalLabels[pos_of_U1] = "Y";

							if (this->SM_CheckVectorlikeness(SMVEVConfig, NumberOfGenerations, true, false, SMVEVConfig.HiggsNo))
							{

								SConfig NewSMVEVConfig = SMVEVConfig;
								this->AutoCreateLabels(SMScheme, NewSMVEVConfig);			//hacking here!!! disable labeling

								NewSMVEVConfig.ConfigLabel  = "SMConfig";
								NewSMVEVConfig.ConfigNumber = SMcounter;
								++SMcounter;

								SM = true;
								AllVEVConfigs.push_back(NewSMVEVConfig);
								Good_VEVConfig_Found = true;

								if (SM_PrintSU5SimpleRoots)
								{
									(*Print.out) << "\n  " << Print.cbegin << "Hypercharge generator:" << Print.cend << "\n" << setiosflags(ios::fixed);
									Print.PrintRational(Hypercharges[j], Lattice);
									(*Print.out) << "\n";
									(*Print.out) << "  " << Print.cbegin << "SU(5) origin:" << Print.cend << "\n" << setiosflags(ios::fixed);
									Print.PrintRational(PossibleRootsOfSU5[j][0], Lattice);
									(*Print.out) << "\n";
									Print.PrintRational(PossibleRootsOfSU5[j][1], Lattice);
									(*Print.out) << "\n";
									Print.PrintRational(PossibleRootsOfSU5[j][2], Lattice);
									(*Print.out) << "\n";
									Print.PrintRational(PossibleRootsOfSU5[j][3], Lattice);
									(*Print.out) << endl;
								}
							}
						}
					}
				}
			}
		}
	}
	if (SU5)
	{
		vector<CVector> FlippedU1s;

		SU5 = false;
		PhenoScheme SU5Scheme;
		this->CreatePhenoScheme("SU5", SU5Scheme);

		unsigned SU5counter = 1;

		TestVEVConfigs.clear();
		if (this->FindPositionOfGaugeGroup(SU5Scheme, OriginalVEVConfig, TestVEVConfigs))
		{
			s1 = TestVEVConfigs.size();

			for (i = 0; i < s1; ++i)
			{
				SConfig &SU5VEVConfig = TestVEVConfigs[i];

				if (abs(this->SU5_NetNumberOfGenerations(SU5VEVConfig)) == NumberOfGenerations)
				{
					SConfig NewSU5VEVConfig = SU5VEVConfig;
					this->AutoCreateLabels(SU5Scheme, NewSU5VEVConfig);

					NewSU5VEVConfig.ConfigLabel  = "SU5Config";
					NewSU5VEVConfig.ConfigNumber = SU5counter;
					++SU5counter;

					SU5 = true;
					AllVEVConfigs.push_back(NewSU5VEVConfig);
					Good_VEVConfig_Found = true;

					/*if (this->SU5_GetFlippedU1(SU5VEVConfig, Roots10D.Weights, FlippedU1s))
          {
            s2 = FlippedU1s.size();
            for (j = 0; j < s2; ++j)
            {
              Orbifold.Config_SetU1Direction(SU5VEVConfig, FlippedU1s[j], pos_of_U1);

              SU5VEVConfig.SymmetryGroup.observable_sector_U1s.clear();
              SU5VEVConfig.SymmetryGroup.observable_sector_U1s.push_back(pos_of_U1);
              SU5VEVConfig.SymmetryGroup.U1s_AdditionalLabels[pos_of_U1] = "fl";

              if (this->SU5_CheckVectorlikeness(Orbifold, SU5VEVConfig, pos_of_U1, Print, false))
              {
                SConfig NewSU5VEVConfig = SU5VEVConfig;
                this->AutoCreateLabels(SU5Scheme, NewSU5VEVConfig);

                NewSU5VEVConfig.ConfigLabel  = "SU5flConfig";
                NewSU5VEVConfig.ConfigNumber = SU5counter;
                ++SU5counter;

                SU5 = true;
                AllVEVConfigs.push_back(NewSU5VEVConfig);
                Good_VEVConfig_Found = true;
              }
            }
          }*/
				}
			}
		}
	}

	if (PS)
	{
		PS = false;
		PhenoScheme PSScheme;
		this->CreatePhenoScheme("PS", PSScheme);

		unsigned PScounter = 1;

		TestVEVConfigs.clear();
		if (this->FindPositionOfGaugeGroup(PSScheme, OriginalVEVConfig, TestVEVConfigs))
		{
			s1 = TestVEVConfigs.size();

			for (i = 0; i < s1; ++i)
			{
				SConfig &PSVEVConfig = TestVEVConfigs[i];

				if (abs(this->PS_NetNumberOfGenerations(PSVEVConfig)) == NumberOfGenerations)
				{
					if (this->PS_CheckVectorlikeness(PSVEVConfig, Print, NumberOfGenerations, false))
					{
						SConfig NewPSVEVConfig = PSVEVConfig;
						this->AutoCreateLabels(PSScheme, NewPSVEVConfig);

						NewPSVEVConfig.ConfigLabel  = "PSConfig";
						NewPSVEVConfig.ConfigNumber = PScounter;
						++PScounter;

						PS = true;
						AllVEVConfigs.push_back(NewPSVEVConfig);
						Good_VEVConfig_Found = true;
					}
				}
			}
		}
	}
	return Good_VEVConfig_Found;
}



/* ########################################################################################
######   FindEmptyFixedBranes(const COrbifold &Orbifold, ...) const                  ######
######                                                                               ######
######   Version: 19.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Orbifold         : the orbifold of the vev-config "OriginalVEVConfig"    ######
######   2) VEVConfig        : contains the massless fields                          ######
######   3) EmptyFixedBranes : the space group elements of those fixed points that   ######
######                         do not have massless twisted strings are stored here  ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Goes through the list of fixed points of "Orbifold" and identifies those    ######
######   who are empty, i.e. have no massless twisted strings. The result is stored  ######
######   in "EmptyFixedBranes".                                                      ######
######################################################################################## */
void CAnalyseModel::FindEmptyFixedBranes(const COrbifold &Orbifold, const SConfig &VEVConfig, vector<CSpaceGroupElement> &EmptyFixedBranes) const
{
	if (VEVConfig.InvariantSupercharges.size() != 1)
	{
		cout << "\n  Warning: FindEmptyFixedBranes(...) only works for N=1 SUSY." << endl;
		return;
	}

	if (EmptyFixedBranes.size() != 0)
	{
		cout << "\n  Warning in void CAnalyseModel::FindEmptyFixedBranes(...) const : EmptyFixedBranes not empty - now cleared!" << endl;
		EmptyFixedBranes.clear();
	}

	unsigned i = 0;
	unsigned j = 0;
	size_t s2 = 0;

	const vector<CField> &Fields = VEVConfig.Fields;

	vector<unsigned> FieldIndices;

	const size_t s1 = Orbifold.GetNumberOfSectors();
	for(i = 0; i < s1; ++i)
	{
		const CSector &Sector = Orbifold.GetSector(i);

		if (Sector.SectorHasLeftchiralRightmover())
		{
			s2 = Sector.GetNumberOfFixedBranes();
			for(j = 0; j < s2; ++j)
			{
				FieldIndices.clear();
				const CFixedBrane &FixedBrane = Sector.GetFixedBrane(j);

				FixedBrane.GetFieldIndices(Fields, LeftChiral, FieldIndices);
				if (FieldIndices.size() == 0)
					EmptyFixedBranes.push_back(FixedBrane.GetSGElement());
			}
		}
	}
}



/* ########################################################################################
######   AutoCreateLabels(const PhenoScheme &Scheme, SConfig &VEVConfig) const       ######
######                                                                               ######
######   Version: 06.07.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Scheme    : specifies the labels for "SU5", "SO10", "SM" or "PS"         ######
######                  phenomenology                                                ######
######   2) VEVConfig : contains the massless fields                                 ######
######   output:                                                                     ######
######   return value : field labels created succesfully?                            ######
###########################################################################################
######   description:                                                                ######
######   Create labels for the fields of "VEVConfig" using the PhenoScheme "Scheme". ######
######################################################################################## */
bool CAnalyseModel::AutoCreateLabels(const PhenoScheme &Scheme, SConfig &VEVConfig) const
{
	// Set the precision
	const double prec = 0.0001;

	vector<CField> &Fields = VEVConfig.Fields;
	const size_t f1 = Fields.size();

	const SSymmetryGroup &SymmetryGroup = VEVConfig.SymmetryGroup;

	const size_t number_of_usedU1s     = SymmetryGroup.observable_sector_U1s.size();
	const size_t number_of_usedfactors = SymmetryGroup.observable_sector_GGs.size();

	const vector<gaugeGroupFactor<double> > &factors      = SymmetryGroup.GaugeGroup.factor;
	const vector<vector<double> >           &U1Generators = SymmetryGroup.GaugeGroup.u1directions;
	const size_t number_of_factors = factors.size();
	const size_t number_of_U1s     = U1Generators.size();

	unsigned i = 0;
	unsigned j = 0;
	unsigned k = 0;
	unsigned Index = 0;

	const size_t s1 = Scheme.SetOfDimensions.size();

	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// begin: Error
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if (f1 == 0)
	{
		cout << "\n  Warning in bool CAnalyseModel::AutoCreateLabels(...) const: Spectrum is empty. Return false." << endl;
		return false;
	}
	if ((s1 == 0) || (Scheme.SetOfU1Charges.size() != s1) || (Scheme.Labels.size() != s1))
	{
		cout << "\n  Warning in bool CAnalyseModel::AutoCreateLabels(...) const: \"Scheme\" is ill-defined. Return false." << endl;
		return false;
	}

	// begin: check "SetOfDimensions", "observable_sector_GGs" and "GaugeGroupFactors"
	if ((Scheme.SetOfDimensions[0].size() != number_of_usedfactors) || (Scheme.GaugeGroupFactors.size() != number_of_usedfactors))
	{
		cout << "\n  Warning in bool CAnalyseModel::AutoCreateLabels(...) const: Cannot find the gauge group factors specified by \"Scheme\". Return false." << endl;
		return false;
	}

	for (i = 0; i < number_of_usedfactors; ++i)
	{
		Index = SymmetryGroup.observable_sector_GGs[i];
		if (Index >= number_of_factors)
		{
			cout << "\n  Warning in bool CAnalyseModel::AutoCreateLabels(...) const: \"observable_sector_GGs\" out of range. Return false." << endl;
			return false;
		}
	}
	// end: check "SetOfDimensions", "observable_sector_GGs" and "GaugeGroupFactors"

	// begin: check "observable_sector_U1s" and "SetOfU1Charges"
	if (Scheme.SetOfU1Charges[0].size() != number_of_usedU1s)
	{
		cout << "\n  Warning in bool CAnalyseModel::AutoCreateLabels(...) const: Cannot find the U(1) factor specified by \"Scheme\". Return false." << endl;
		return false;
	}

	for (i = 0; i < number_of_usedU1s; ++i)
	{
		Index = SymmetryGroup.observable_sector_U1s[i];
		if (Index >= number_of_U1s)
		{
			cout << "\n  Warning in bool CAnalyseModel::AutoCreateLabels(...) const: \"observable_sector_U1s\" out of range. Return false." << endl;
			return false;
		}
	}
	if (SymmetryGroup.IsFirstU1Anomalous)
	{
		for (i = 0; i < number_of_usedU1s; ++i)
		{
			Index = SymmetryGroup.observable_sector_U1s[i];
			if (SymmetryGroup.observable_sector_U1s[i] == 0)
				cout << "\n  Warning in bool CAnalyseModel::AutoCreateLabels(...) const: one factor of \"observable_sector_U1s\" is the anomalous U(1). Continue..." << endl;
		}
	}
	// end: check "observable_sector_U1s" and "SetOfU1Charges"
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// end: Error
	////////////////////////////////////////////////////////////////////////////////////////////////////////////

	vector<double> UseNormalizationFactors(number_of_usedU1s, 1.0);
	// begin: find the normalization factors to be used in the schemes U(1) charges
	if (Scheme.NormalizationOfLengthOfU1Generators.size() != 0)
	{
		double sqr = 0.0;
		if (Scheme.NormalizationOfLengthOfU1Generators.size() != number_of_usedU1s)
		{
			cout << "\n  Warning in bool CAnalyseModel::AutoCreateLabels(...) const: Normalization of U(1) Generators failed. Return false." << endl;
			return false;
		}
		for (i = 0; i < number_of_usedU1s; ++i)
		{
			const vector<double> &U1Generator = U1Generators[SymmetryGroup.observable_sector_U1s[i]];

			sqr = 0.0;
			for (j = 0; j < 16; ++j)
				sqr += U1Generator[j] * U1Generator[j];

			if (fabs(sqr - Scheme.NormalizationOfLengthOfU1Generators[i]) > prec)
				UseNormalizationFactors[i] = sqrt(sqr/Scheme.NormalizationOfLengthOfU1Generators[i]);
		}
	}
	// end: find the normalization factors to be used in the schemes U(1) charges

	// begin: find the reordering of the scheme
	// the ordering of the gauge group factors in "SymmetryGroup.observable_sector_GGs" and "Scheme.GaugeGroupFactors" might differ
	vector<unsigned> Scheme2observable_sector_GGs(number_of_usedfactors, 0);
	{
		bool not_found = true;
		vector<bool> PosUsed(number_of_usedfactors, false);
		for (i = 0; i < number_of_usedfactors; ++i)
		{
			const string &ggf = Scheme.GaugeGroupFactors[i];

			not_found = true;
			for (j = 0; not_found && (j < number_of_usedfactors); ++j)
			{
				if (!PosUsed[j] && (factors[SymmetryGroup.observable_sector_GGs[j]].algebra == ggf))
				{
					Scheme2observable_sector_GGs[i] = j;
					PosUsed[j] = true;
					not_found = false;
				}
			}
		}
		for (i = 0; i < number_of_usedfactors; ++i)
		{
			if (!PosUsed[i])
			{
				cout << "\n  Warning in bool CAnalyseModel::AutoCreateLabels(...) const: Cannot find all gauge group factor. Return false." << endl;
				return false;
			}
		}
	}
	// end: find the reordering of the scheme

	if (Scheme.IndexOfSetToNormalizeCharges >= s1)
	{
		cout << "\n  Warning in bool CAnalyseModel::AutoCreateLabels(...) const: \"IndexOfSetToNormalizeCharges\" is ill-defined. Return false." << endl;
		return false;
	}

	const RepVector &SpecialDimensions = Scheme.SetOfDimensions[Scheme.IndexOfSetToNormalizeCharges];
	const CVector   &SpecialU1Charges  = Scheme.SetOfU1Charges[Scheme.IndexOfSetToNormalizeCharges];

	size_t s2 = 0;

	SDimension tmp;
	tmp.Dimension = 1;
	tmp.AdditionalLabel = "";
	RepVector newDimensions(number_of_usedfactors, tmp);
	CVector   newU1Charges(number_of_usedU1s);

	RepVector newDimensionsCC(number_of_usedfactors, tmp);
	CVector   newU1ChargesCC(number_of_usedU1s);

	// begin: find the RepVectors "SpecialDimensions" and "SpecialDimensionsCC"
	//        where the number of "SpecialDimensions" representations must be larger than the one of "SpecialDimensionsCC"
	vector<int> U1Flip(number_of_usedU1s, 1);
	vector<int> RepFlip(number_of_usedfactors, 1);
	bool equal = true;

	{
		unsigned mult = 1;
		vector<RepVector> KnownRepsOfSpecialType;
		vector<CVector>   KnownU1ChargesOfSpecialType;
		vector<int>       Counter;

		for (i = 0; i < f1; ++i)
		{
			const CField &Field = Fields[i];

			// only left-chiral fields
			if (Field.Multiplet == LeftFermi || Field.Multiplet == Scalar)
			{
				// get the relevant dimensions and U(1) charges
				for (j = 0; j < number_of_usedfactors; ++j)
					newDimensions[j] = Field.Dimensions[SymmetryGroup.observable_sector_GGs[j]];

				for (j = 0; j < number_of_usedU1s; ++j)
					newU1Charges[j] = Field.U1Charges[SymmetryGroup.observable_sector_U1s[j]];

				// check equality up to the sign
				equal = true;
				for (j = 0; equal && (j < number_of_usedfactors); ++j)
				{
					if (abs(newDimensions[j].Dimension) != abs(SpecialDimensions[Scheme2observable_sector_GGs[j]].Dimension))
						equal = false;
				}
				for (j = 0; equal && (j < number_of_usedU1s); ++j)
				{
					if (fabs(fabs(newU1Charges[j]) - fabs(UseNormalizationFactors[j] * SpecialU1Charges[j])) > prec)
						equal = false;
				}
				// current "Dimensions" is of special type
				if (equal)
				{
					mult = Multiplicity(SymmetryGroup, Field.Dimensions);

					// begin: charge conjugation
					for (j = 0; j < number_of_usedU1s; ++j)
						newU1ChargesCC[j]  = -newU1Charges[j];

					CC(newDimensions, newDimensionsCC);
					// end: charge conjugation

					equal = false;

					s2 = KnownRepsOfSpecialType.size();
					for (j = 0; !equal && (j < s2); ++j)
					{
						if (AreRepVectorsEqual(newDimensions, KnownRepsOfSpecialType[j]) && (newU1Charges == KnownU1ChargesOfSpecialType[j]))
						{
							equal = true;
							Counter[j] += mult;
						}
						else
							if (AreRepVectorsEqual(newDimensionsCC, KnownRepsOfSpecialType[j]) && (newU1ChargesCC == KnownU1ChargesOfSpecialType[j]))
							{
								equal = true;
								Counter[j] -= mult;
							}
					}
					if (!equal)
					{
						KnownRepsOfSpecialType.push_back(newDimensions);
						KnownU1ChargesOfSpecialType.push_back(newU1Charges);
						Counter.push_back(mult);
					}
				}
			}
		}

		unsigned number_of_chiral_reps = 0;
		s2 = KnownRepsOfSpecialType.size();
		for (i = 0; i < s2; ++i)
		{
			if (Counter[i] != 0)
				++number_of_chiral_reps;
		}
		if (number_of_chiral_reps != 1)
		{
			cout << "\n  Warning in bool CAnalyseModel::AutoCreateLabels(...) const: number of chiral representations of special type not equal 1. Return false." << endl;
			return false;
		}

		int extra_sign = 1;
		for (i = 0; i < s2; ++i)
		{
			if (Counter[i] != 0)
			{
				const RepVector &KnownRep = KnownRepsOfSpecialType[i];

				RepVector FoundSpecialDimensions(number_of_usedfactors, tmp);
				if (Counter[i] < 0)
				{
					extra_sign = -1;
					CC(KnownRep, FoundSpecialDimensions);
				}
				else
					FoundSpecialDimensions = KnownRep;

				const CVector &FoundU1ChargesOfSpecialType = KnownU1ChargesOfSpecialType[i];

				for (j = 0; j < number_of_usedfactors; ++j)
				{
					if (FoundSpecialDimensions[j].Dimension != SpecialDimensions[Scheme2observable_sector_GGs[j]].Dimension)
						RepFlip[j] = -1;
				}
				for (j = 0; j < number_of_usedU1s; ++j)
				{
					if (fabs((extra_sign * FoundU1ChargesOfSpecialType[j]) - (UseNormalizationFactors[j] * SpecialU1Charges[j])) > prec)
						U1Flip[j] = -1;
				}

				break;
			}
		}
	}
	// end: find the RepVectors "SpecialDimensions" and "SpecialDimensionsCC"

	vector<RepVector> new_SetOfDimensions;
	vector<CVector>   new_SetOfU1Charges;
	vector<SUSYMultiplet> new_SetOfMultiplets;
	vector<string>    new_Labels;
	vector<unsigned>  new_counter;

	string labelAdd = "";
	string label = "";
	vector<string> XLabels;
	XLabels.push_back("v");
	XLabels.push_back("w");
	XLabels.push_back("x");
	XLabels.push_back("y");
	XLabels.push_back("z");
	unsigned Xcounter = 0;
	unsigned XcounterAdd = 2;

	vector<unsigned>  counter(s1 ,1);

	bool label_not_found = true;
	for (i = 0; i < f1; ++i)
	{
		CField &Field = Fields[i];
		if (Field.Multiplet == LeftFermi || Field.Multiplet == Scalar)
		{
			const RepVector &Dimensions = Field.Dimensions;
			const CVector   &U1Charges  = Field.U1Charges;

			label_not_found = true;

			// search the current "Dimensions" and "U1Charges" in the label scheme
			for (j = 0; label_not_found && (j < s1); ++j)
			{
				const RepVector &SchemeDimensions = Scheme.SetOfDimensions[j];
				const CVector   &SchemeU1Charges  = Scheme.SetOfU1Charges[j];
				const SUSYMultiplet &SchemeMultiplet = Scheme.Multiplets[j];

				equal = true;
				for (k = 0; equal && (k < number_of_usedfactors); ++k)
				{
					const SDimension &dim1 = Dimensions[SymmetryGroup.observable_sector_GGs[k]];
					const SDimension &dim2 = SchemeDimensions[Scheme2observable_sector_GGs[k]];

					// if the dimension is complex
					if ((dim1.Dimension != 1) && (dim1.Dimension != 2))
					{
						if ((dim1.AdditionalLabel != dim2.AdditionalLabel) || (dim1.Dimension != RepFlip[k] * dim2.Dimension))
							equal = false;
					}
					else
					{
						if ((dim1.AdditionalLabel != dim2.AdditionalLabel) || (dim1.Dimension != dim2.Dimension))
							equal = false;
					}
				}
				for (k = 0; equal && (k < number_of_usedU1s); ++k)
				{
					if (fabs(U1Charges[SymmetryGroup.observable_sector_U1s[k]] - (UseNormalizationFactors[k] * U1Flip[k] * SchemeU1Charges[k])) > prec)
						equal = false;
				}
				if (equal && Field.Multiplet == SchemeMultiplet)
				{
					Field.Labels.push_back(Scheme.Labels[j]);
					Field.Numbers.push_back(counter[j]);
					++counter[j];

					label_not_found = false;
				}
			}
			// if the current "Dimensions" and "U1Charges" are not part of the label scheme
			if (label_not_found)
			{
				for (k = 0; k < number_of_usedfactors; ++k)
					newDimensions[k] = Dimensions[SymmetryGroup.observable_sector_GGs[k]];

				for (k = 0; k < number_of_usedU1s; ++k)
					newU1Charges[k] = U1Charges[SymmetryGroup.observable_sector_U1s[k]];

				s2 = new_SetOfDimensions.size();

				equal = false;
				for (k = 0; !equal && (k < s2); ++k)
				{
					if ((AreRepVectorsEqual(new_SetOfDimensions[k], newDimensions)) && (new_SetOfU1Charges[k] == newU1Charges) && (new_SetOfMultiplets[k] == Field.Multiplet))
					{
						Field.Labels.push_back(new_Labels[k]);
						Field.Numbers.push_back(new_counter[k]);
						++new_counter[k];

						equal = true;
						break;
					}
				}
				if (!equal)
				{
					// begin: charge conjugation
					for (k = 0; k < number_of_usedU1s; ++k)
						newU1ChargesCC[k]  = -newU1Charges[k];

					CC(newDimensions, newDimensionsCC);
					// end: charge conjugation

					// if the current "Dimensions" and "U1Charges" are the charge conjugate of a known one add a "b" in front of the label
					for (k = 0; !equal && (k < s2); ++k)
					{
						if ((AreRepVectorsEqual(new_SetOfDimensions[k], newDimensionsCC)) && (new_SetOfU1Charges[k] == newU1ChargesCC) && (new_SetOfMultiplets[k] == Field.Multiplet))
						{
							label = "b" + new_Labels[k];
							Field.Labels.push_back(label);
							Field.Numbers.push_back(1);

							new_SetOfDimensions.push_back(newDimensions);
							new_SetOfU1Charges.push_back(newU1Charges);
							new_SetOfMultiplets.push_back(Field.Multiplet);
							new_Labels.push_back(label);
							new_counter.push_back(2);

							equal = true;
							break;
						}
					}
					if (!equal)
					{
						if (Xcounter == 5)
						{
							ostringstream os;
							os << XcounterAdd;
							labelAdd = os.str();

							Xcounter = 0;
							++XcounterAdd;
						}
						if (Field.Multiplet == Scalar || Field.Multiplet == bScalar) 
						label = "s" + XLabels[Xcounter] + labelAdd;   
						else                                          
						label = XLabels[Xcounter] + labelAdd;
						++Xcounter;

						Field.Labels.push_back(label);
						Field.Numbers.push_back(1);

						new_SetOfDimensions.push_back(newDimensions);
						new_SetOfU1Charges.push_back(newU1Charges);
						new_SetOfMultiplets.push_back(Field.Multiplet);
						new_Labels.push_back(label);
						new_counter.push_back(2);
					}
				}
			}
		}
		else
		{
			Field.Labels.push_back(Field.Labels[Field.Labels.size()-1]);
			Field.Numbers.push_back(Field.Numbers[Field.Numbers.size()-1]);
		}
	}
	// use the newly created labels
	++VEVConfig.use_Labels;

	return true;
}


/* ########################################################################################
######   SM_FindPositionOfGaugeGroupFromGG(...) const                                ######
######                                                                               ######
######   Version: 19.05.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) OriginalVEVConfig : vev-config to analyze                                ######
######   2) GG_SimpleRoots    : simple roots of some larger gauge group, like SU(5)  ######
######   3) SMVEVConfigs      : resulting vev-configs                                ######
######   output:                                                                     ######
######   return value         : new vev-config identified?                           ######
###########################################################################################
######   description:                                                                ######
######   Identifies whether the Standard Model gauge group is inside the one spanned ######
######   by "GG_SimpleRoots".                                                        ######
######################################################################################## */
bool CAnalyseModel::SM_FindPositionOfGaugeGroupFromGG(const SConfig &OriginalVEVConfig, const vector<vector<double> > &GG_SimpleRoots, vector<SConfig> &SMVEVConfigs) const
{
	// Set the precision
	const double prec = 0.0001;

	const vector<gaugeGroupFactor<double> > &factors = OriginalVEVConfig.SymmetryGroup.GaugeGroup.factor;

	const size_t number_of_GG_simpleroots = GG_SimpleRoots.size();
	const size_t number_of_factors        = factors.size();

	if (number_of_GG_simpleroots == 0)
	{
		cout << "\n  Warning in bool CAnalyseModel::SM_FindPositionOfGaugeGroupFromGG(...) const: \"GG_SimpleRoots\" is empty. Return false." << endl;
		return false;
	}

	const string SU3 = "A2";
	const string SU2 = "A1";

	unsigned i = 0;
	unsigned j = 0;
	unsigned k = 0;
	unsigned l = 0;

	double tmp_sp = 0.0;

	bool orthogonal = true;

	unsigned counter = 1;
	SConfig NewVEVConfig;
	NewVEVConfig = OriginalVEVConfig;
	NewVEVConfig.SymmetryGroup.observable_sector_GGs.clear();
	NewVEVConfig.SymmetryGroup.observable_sector_U1s.clear();
	NewVEVConfig.SymmetryGroup.GGs_AdditionalLabels.assign(number_of_factors, "");
	NewVEVConfig.SymmetryGroup.U1s_AdditionalLabels.assign(OriginalVEVConfig.SymmetryGroup.GaugeGroup.u1directions.size(), "");
	bool new_VEVConfig_found = false;

	// run through all gauge group factors
	for (i = 0; i < number_of_factors; ++i)
	{
		const gaugeGroupFactor<double> &ggf = factors[i];

		// if current gauge group factor is SU(3)
		if (ggf.algebra == SU3)
		{
			const vector<vector<double> > &SU3_SimpleRoots = ggf.simpleroots;

			// assume that the two simple roots of SU(3) are orthogonal to all GG_SimpleRoots
			orthogonal = true;
			for (j = 0; orthogonal && (j < 2); ++j)
			{
				const vector<double> &SU3_SimpleRoot = SU3_SimpleRoots[j];

				for (k = 0; orthogonal && (k < number_of_GG_simpleroots); ++k)
				{
					const vector<double> &GG_SimpleRoot = GG_SimpleRoots[k];

					tmp_sp = 0.0;
					for (l = 0; l < 16; ++l)
						tmp_sp += SU3_SimpleRoot[l] * GG_SimpleRoot[l];

					if (fabs(tmp_sp) > prec)
						orthogonal = false;
				}
			}
			// if SU(3) is inside GG
			if (!orthogonal)
			{
				// run through all gauge group factors
				for (j = 0; j < number_of_factors; ++j)
				{
					const gaugeGroupFactor<double> &ggf = factors[j];

					// if current gauge group factor is SU(2)
					if (ggf.algebra == SU2)
					{
						const vector<double> &SU2_SimpleRoot = ggf.simpleroots[0];

						// assume that the simple root of SU(2) is orthogonal to all GG_SimpleRoots
						orthogonal = true;
						for (k = 0; orthogonal && (k < number_of_GG_simpleroots); ++k)
						{
							const vector<double> &GG_SimpleRoot = GG_SimpleRoots[k];

							tmp_sp = 0.0;
							for (l = 0; l < 16; ++l)
								tmp_sp += SU2_SimpleRoot[l] * GG_SimpleRoot[l];

							if (fabs(tmp_sp) > prec)
								orthogonal = false;
						}
						// if SU(2) is inside GG
						if (!orthogonal)
						{
							new_VEVConfig_found = true;

							ostringstream os;
							os << counter;
							++counter;

							NewVEVConfig.ConfigLabel = "SM VEVConfig " + os.str();
							NewVEVConfig.SymmetryGroup.observable_sector_GGs.clear();
							NewVEVConfig.SymmetryGroup.observable_sector_GGs.push_back(i);
							NewVEVConfig.SymmetryGroup.observable_sector_GGs.push_back(j);
							NewVEVConfig.SymmetryGroup.GGs_AdditionalLabels[i] = "C";
							NewVEVConfig.SymmetryGroup.GGs_AdditionalLabels[j] = "L";

							SMVEVConfigs.push_back(NewVEVConfig);
						}
					}
				}
			}
		}
	}
	return new_VEVConfig_found;
}



/* ########################################################################################
######   int CAnalyseModel::SM_NetNumberOfGenerations(...) const                     ######
######                                                                               ######
######   Version: 10.12.2010                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) SMVEVConfig : a Standard Model vev-config for which the net number of    ######
######                    (3,2)-plets is counted                                     ######
######   output:                                                                     ######
######   return value   : net number of Standard Model (3,2)-plets                   ######
###########################################################################################
######   description:                                                                ######
######   Counts the net number of (3,2)-plets in the vev-config "SMVEVConfig".       ######
######################################################################################## */
int CAnalyseModel::SM_NetNumberOfGenerations(const SConfig &SMVEVConfig) const
{
	int net_no_families = 0;
	int mult            = 0;

	unsigned i = 0;
	unsigned j = 0;

	const vector<CField> &Fields = SMVEVConfig.Fields;
	const size_t f1 = Fields.size();

	if (f1 == 0)
	{
		cout << "\n  Warning in int CAnalyseModel::SM_NetNumberOfGenerations(...) const: Spectrum is empty. Return 0." << endl;
		return 0;
	}

	const vector<gaugeGroupFactor<double> > &factors = SMVEVConfig.SymmetryGroup.GaugeGroup.factor;

	if (SMVEVConfig.SymmetryGroup.observable_sector_GGs.size() < 2)
	{
		cout << "\n  Warning in int CAnalyseModel::SM_NetNumberOfGenerations(...) const: Cannot find the Standard Model. Return 0." << endl;
		return 0;
	}

	const unsigned pos_of_SU3 = SMVEVConfig.SymmetryGroup.observable_sector_GGs[0];
	const unsigned pos_of_SU2 = SMVEVConfig.SymmetryGroup.observable_sector_GGs[1];

	const size_t number_of_factors = factors.size();
	if ((pos_of_SU3 >= number_of_factors) || (pos_of_SU2 >= number_of_factors))
	{
		cout << "\n  Warning in int CAnalyseModel::SM_NetNumberOfGenerations(...) const: Cannot find the Standard Model group. Return 0." << endl;
		return 0;
	}

	if ((factors[pos_of_SU3].algebra != "A2") || (factors[pos_of_SU2].algebra != "A1"))
	{
		cout << "\n  Warning in int CAnalyseModel::SM_NetNumberOfGenerations(...) const: Cannot find the Standard Model group. Return 0." << endl;
		return 0;
	}

	for (i = 0; i < f1; ++i)
	{
		const CField &Field = Fields[i];
		if (Field.Multiplet == LeftFermi)			//hacking here!!!
		{
			const RepVector &Dimensions = Field.Dimensions;
			if (abs(Dimensions[pos_of_SU3].Dimension) == 3)
			{
				// count the net number of (3,2)-plets
				if (Dimensions[pos_of_SU2].Dimension == 2)
				{
					mult = 1;

					for (j = 0; j < number_of_factors; ++j)
					{
						if (j == pos_of_SU3)
							mult *= Dimensions[j].Dimension;
						else
							mult *= abs(Dimensions[j].Dimension);
					}
					net_no_families += mult/6;
				}
			}
		}
	}
	return net_no_families;
}



/* ########################################################################################
######   SM_GetHypercharges(...) const                                               ######
######                                                                               ######
######   Version: 07.03.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) SMVEVConfig        : a Standard Model vev-config where to look for       ######
######                           hypercharge                                         ######
######   2) Hypercharges       : the resulting hypercharge generators                ######
######   3) PossibleRootsOfSU5 : roots of an intermediate SU(5) used to define the   ######
######                           hypercharge                                         ######
######   4) Roots10D           : 480 roots of E_8 x E_8 or SO(32) used to identify   ######
######                           an intermediate SU(5)                               ######
######   output:                                                                     ######
######   return value          : hypercharge identified?                             ######
###########################################################################################
######   description:                                                                ######
######   Analyze the vev-config "SMVEVConfig" to identify possible hypercharge       ######
######   generators.                                                                 ######
######################################################################################## */
bool CAnalyseModel::SM_GetHypercharges(const SConfig &SMVEVConfig, vector<CVector> &Hypercharges, vector<vector<CVector> > &PossibleRootsOfSU5, const vector<CVector> &Roots10D) const
{
	const double prec = 0.0001;
	unsigned i = 0, j = 0, k = 0;

	if (Hypercharges.size() != 0)
	{
		cout << "\n  Warning in bool CAnalyseModel::SM_GetHypercharges(...) const: Hypercharges not empty - now cleared." << endl;
		Hypercharges.clear();
	}

	if (PossibleRootsOfSU5.size() != 0)
	{
		cout << "\n  Warning in bool CAnalyseModel::SM_GetHypercharges(...) const: PossibleRootsOfSU5 not empty - now cleared." << endl;
		PossibleRootsOfSU5.clear();
	}

	const vector<gaugeGroupFactor<double> > &factors = SMVEVConfig.SymmetryGroup.GaugeGroup.factor;

	if (SMVEVConfig.SymmetryGroup.observable_sector_GGs.size() < 2)
	{
		cout << "\n  Warning in bool CAnalyseModel::SM_GetHypercharges(...): Cannot find the Standard Model. Return false." << endl;
		return false;
	}

	const unsigned pos_of_SU3 = SMVEVConfig.SymmetryGroup.observable_sector_GGs[0];
	const unsigned pos_of_SU2 = SMVEVConfig.SymmetryGroup.observable_sector_GGs[1];

	const size_t number_of_factors = factors.size();
	if ((pos_of_SU3 >= number_of_factors) || (pos_of_SU2 >= number_of_factors))
	{
		cout << "\n  Warning in bool CAnalyseModel::SM_GetHypercharges(...): Cannot find the Standard Model group. Return false." << endl;
		return false;
	}

	if ((factors[pos_of_SU3].algebra != "A2") || (factors[pos_of_SU2].algebra != "A1"))
	{
		cout << "\n  Warning in bool CAnalyseModel::SM_GetHypercharges(...): Cannot find the Standard Model group. Return false." << endl;
		return false;
	}

	if (Roots10D.size() != 480)
	{
		cout << "\n  Warning in bool CAnalyseModel::SM_GetHypercharges(...) const: \"Roots10D\" does not contain 480 roots. Return false." << endl;
		return false;
	}

	const vector<vector<double> > &SU3_Simpleroots = factors[pos_of_SU3].simpleroots;
	const vector<vector<double> > &SU2_Simpleroots = factors[pos_of_SU2].simpleroots;

	if ((SU3_Simpleroots.size() != 2) || (SU2_Simpleroots.size() != 1))
	{
		cout << "\n  Warning in bool CAnalyseModel::SM_GetHypercharges(...) const: Simple roots of SU(3) or SU(2) not correct. Return false." << endl;
		return false;
	}

	CVector SU3_Simpleroot1, SU3_Simpleroot2, SU2_Simpleroot1, SU5_Simpleroot4;
	CVector tmp_root(16);

	SU3_Simpleroot1 = SU3_Simpleroots[0];
	SU3_Simpleroot2 = SU3_Simpleroots[1];
	SU2_Simpleroot1 = SU2_Simpleroots[0];

	vector<CVector> RootsOfSU5;
	RootsOfSU5.push_back(SU3_Simpleroot1);
	RootsOfSU5.push_back(SU3_Simpleroot2);
	RootsOfSU5.push_back(SU2_Simpleroot1);
	RootsOfSU5.push_back(SU5_Simpleroot4);

	CVector TestHypercharge(16);
	vector<CVector> TestHypercharges;
	vector<vector<CVector> > TestRootsOfSU5;

	vector<vector<double> > Orig_Y_Orthogonal;
	vector<vector<double> > Orig_all_SimpleRoots;
	vector<vector<double> > Orig_SU5_Orthogonal;

	vector<vector<double> > Y_Orthogonal;
	vector<vector<double> > all_SimpleRoots;
	vector<vector<double> > SU5_Orthogonal;

	for (j = 0; j < number_of_factors; ++j)
	{
		const vector<vector<double> > &tmp_SimpleRoots = factors[j].simpleroots;

		Orig_all_SimpleRoots.insert(Orig_all_SimpleRoots.end(), tmp_SimpleRoots.begin(), tmp_SimpleRoots.end());
		if ((j == pos_of_SU3) || (j == pos_of_SU2))
			Orig_Y_Orthogonal.insert(Orig_Y_Orthogonal.end(), tmp_SimpleRoots.begin(), tmp_SimpleRoots.end());
		else
			Orig_SU5_Orthogonal.insert(Orig_SU5_Orthogonal.end(), tmp_SimpleRoots.begin(), tmp_SimpleRoots.end());
	}

	const size_t s3 = Orig_SU5_Orthogonal.size();

	// run through the 480 roots of E8 x E8
	for (i = 0; i < 480; ++i)
	{
		SU5_Simpleroot4 = Roots10D[i];

		Y_Orthogonal    = Orig_Y_Orthogonal;
		all_SimpleRoots = Orig_all_SimpleRoots;
		SU5_Orthogonal  = Orig_SU5_Orthogonal;

		// begin: find a candidate for the fourth SU(5) simple root
		bool is_SU5_Simpleroot4 = true;

		for (j = 0; (j < s3) && is_SU5_Simpleroot4; ++j)
		{
			tmp_root = SU5_Orthogonal[j];

			if (fabs(SU5_Simpleroot4 * tmp_root) > prec)
				is_SU5_Simpleroot4 = false;
		}

		if (is_SU5_Simpleroot4)
		{
			is_SU5_Simpleroot4 = false;
			if (fabs(SU5_Simpleroot4 * SU2_Simpleroot1 + 1.0) < prec)
			{
				const double SU3Root1xSU5Root4 = SU5_Simpleroot4 * SU3_Simpleroot1;
				const double SU3Root2xSU5Root4 = SU5_Simpleroot4 * SU3_Simpleroot2;

				if (((fabs(SU3Root1xSU5Root4 + 1.0) < prec) && (fabs(SU3Root2xSU5Root4) < prec)) || ((fabs(SU3Root1xSU5Root4) < prec) && (fabs(SU3Root2xSU5Root4 + 1.0) < prec)))
					is_SU5_Simpleroot4 = true;
			}
		}
		// end: find a candidate for the fourth SU(5) simple root

		if (is_SU5_Simpleroot4)
		{
			all_SimpleRoots.push_back(SU5_Simpleroot4);

			// begin: find the U(1)'s that are orthogonal to SU(5)
			vector<vector<double> > tmp_U1s;
			Find_Basis_Of_Orthogonal_Space(all_SimpleRoots, UNSPECIFIED_LATTICE, 16, tmp_U1s);

			if (tmp_U1s.size() != SMVEVConfig.SymmetryGroup.GaugeGroup.u1directions.size() - 1)
			{
				cout << "\n  Warning in bool CAnalyseModel::SM_GetHypercharges(...) const: Excesive number of U(1)'s found. Return false." << endl;
				return false;
			}
			// end: find the U(1)'s that are orthogonal to SU(5)

			// begin: find the Hypercharge
			SU5_Orthogonal.insert(SU5_Orthogonal.end(), tmp_U1s.begin(), tmp_U1s.end());
			Y_Orthogonal.insert(Y_Orthogonal.end(), SU5_Orthogonal.begin(), SU5_Orthogonal.end());

			vector<vector<double> > tmp_Y;
			Find_Basis_Of_Orthogonal_Space(Y_Orthogonal, UNSPECIFIED_LATTICE, 16, tmp_Y);

			if (tmp_Y.size() != 1)
			{
				cout << "\n  Warning in bool CAnalyseModel::SM_GetHypercharges(...) const:" << tmp_Y.size() << " Y's were found, but only one is expected. Return false." << endl;
				return false;
			}
			// end: find the Hypercharge

			// begin: check whether the Hypercharge can be spanned by the gauge U(1)'s
			vector<vector<double> > U1directions = SMVEVConfig.SymmetryGroup.GaugeGroup.u1directions;
			const size_t gauge_U1s = U1directions.size();
			U1directions.push_back(tmp_Y[0]);
			U1directions = findBasis<double>(U1directions);

			if (gauge_U1s != U1directions.size())
			{
				cout << "\n  Warning in bool CAnalyseModel::SM_GetHypercharges(...) const: Hypercharge is NOT gauge U1. Return false." << endl;
				return false;
			}
			// end: check whether the Hypercharge can be spanned by the gauge U(1)'s

			TestHypercharge = tmp_Y[0];
			TestHypercharge = TestHypercharge * sqrt((5.0/6.0)/(TestHypercharge.GetSqrTo(16)));

			TestHypercharges.push_back(TestHypercharge);

			RootsOfSU5[3] = SU5_Simpleroot4;
			TestRootsOfSU5.push_back(RootsOfSU5);
		}
	}

	size_t s5 = 0;
	bool orthogonal = true;

	CVector U1DirectionA(16);
	U1DirectionA = SMVEVConfig.SymmetryGroup.GaugeGroup.u1directions[0];

	CVector Simpleroot(16);

	const bool AnomalousU1 = SMVEVConfig.SymmetryGroup.IsFirstU1Anomalous;

	const size_t s4 = TestHypercharges.size();
	for (i = 0; i < s4; ++i)
	{
		const CVector &Test = TestHypercharges[i];

		orthogonal = true;

		for (j = 0; orthogonal && (j < number_of_factors); ++j)
		{
			const vector<vector<double> > &Simpleroots = factors[j].simpleroots;
			s5 = Simpleroots.size();
			for (k = 0; orthogonal && (k < s5); ++k)
			{
				Simpleroot = Simpleroots[k];
				if (fabs(Test * Simpleroot) > prec)
					orthogonal = false;
			}
		}

		if (orthogonal && AnomalousU1 && (fabs(U1DirectionA * Test) > prec))
			orthogonal = false;

		if (orthogonal)
		{
			s5 = Hypercharges.size();
			bool known = false;
			for (j = 0; !known && j < s5; ++j)
			{
				CVector YTmp = Hypercharges[j];
				if((Test == YTmp) || (Test == YTmp * (-1)))
					known = true;
			}
			if (!known)
			{
				PossibleRootsOfSU5.push_back(TestRootsOfSU5[i]);
				Hypercharges.push_back(Test);
			}
		}
	}

	//cout << Hypercharges.size() << " Hypercharges :" << endl;
	//for (i=0; i<Hypercharges.size(); ++i)
	//{
	//  Hypercharges[i].PrintRational(Lattice, cout);
	//  cout << endl;
	//}

	if (Hypercharges.size() != 0)
		return true;
	else
		return false;
}



/* ########################################################################################
######   bool CAnalyseModel::SM_CheckVectorlikeness(...)                             ######
######                                                                               ######
######   Version: 15.07.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) SMVEVConfig         : vev-config with SM gauge group to analyze          ######
######   2) NumberOfGenerations : number of SM generations (e.g. 3)                  ######
######   3) normalizeY          : normalize hypercharge generator to 5/6             ######
######   4) info                : print info                                         ######
######   output:                                                                     ######
######   return value           : spectrum is SM generations + vectorlike exotics?   ######
###########################################################################################
######   description:                                                                ######
######   Checks whether the vev-config "SMVEVConfig" has a net number of             ######
######   "NumberOfGenerations" SM generations plus vectorlike exotics.               ######
######################################################################################## */
bool CAnalyseModel::SM_CheckVectorlikeness(const SConfig &SMVEVConfig, unsigned NumberOfGenerations, bool normalizeY, bool info, unsigned &HiggsNo) const
{
	// Set the precision
	const double prec = 0.0001;

	double Y_Q  = 1.0/6.0;
	double Y_ub =-2.0/3.0;
	double Y_eb = 1.0;
	double Y_db = 1.0/3.0;
	double Y_L  =-1.0/2.0;
	double Y_n  = 0.0;
	double Y_Hu = 1.0/2.0;
	double Y_Hd =-1.0/2.0;

	unsigned i = 0;

	const vector<CField> &Fields = SMVEVConfig.Fields;
	const size_t f1 = Fields.size();
	if (f1 == 0)
	{
		cout << "\n  Warning in bool CAnalyseModel::SM_CheckVectorlikeness(...): Spectrum is empty. Return false." << endl;
		return false;
	}

	const vector<gaugeGroupFactor<double> > &factors = SMVEVConfig.SymmetryGroup.GaugeGroup.factor;
	const size_t number_of_factors = factors.size();
	const size_t number_of_U1s     = SMVEVConfig.SymmetryGroup.GaugeGroup.u1directions.size();

	if (SMVEVConfig.SymmetryGroup.observable_sector_GGs.size() < 2)
	{
		cout << "\n  Warning in bool CAnalyseModel::SM_CheckVectorlikeness(...): Cannot find the Standard Model. Return false." << endl;
		return false;
	}

	const unsigned pos_of_SU3 = SMVEVConfig.SymmetryGroup.observable_sector_GGs[0];
	const unsigned pos_of_SU2 = SMVEVConfig.SymmetryGroup.observable_sector_GGs[1];

	if ((pos_of_SU3 >= number_of_factors) || (pos_of_SU2 >= number_of_factors))
	{
		cout << "\n  Warning in bool CAnalyseModel::SM_CheckVectorlikeness(...): Cannot find the Standard Model group. Return false." << endl;
		return false;
	}

	if ((factors[pos_of_SU3].algebra != "A2") || (factors[pos_of_SU2].algebra != "A1"))
	{
		cout << "\n  Warning in bool CAnalyseModel::SM_CheckVectorlikeness(...): Cannot find the Standard Model group. Return false." << endl;
		return false;
	}

	if (SMVEVConfig.SymmetryGroup.observable_sector_U1s.size() < 1)
	{
		cout << "\n  Warning in bool CAnalyseModel::SM_CheckVectorlikeness(...): Cannot find the Hypercharge. Return false." << endl;
		return false;
	}

	const unsigned pos_of_Y = SMVEVConfig.SymmetryGroup.observable_sector_U1s[0];
	if (pos_of_Y >= number_of_U1s)
	{
		cout << "\n  Warning in bool CAnalyseModel::SM_CheckVectorlikeness(...): Cannot find the Hypercharge. Return false." << endl;
		return false;
	}

	if (SMVEVConfig.SymmetryGroup.IsFirstU1Anomalous && (pos_of_Y == 0))
	{
		cout << "\n  Warning in bool CAnalyseModel::SM_CheckVectorlikeness(...): Hypercharge cannot be the anomalous U(1). Return false." << endl;
		return false;
	}

	if (normalizeY)
	{
		const vector<double> &hypercharge = SMVEVConfig.SymmetryGroup.GaugeGroup.u1directions[pos_of_Y];

		double Ysqr = 0.0;
		for (i = 0; i < 16; ++i)
			Ysqr += hypercharge[i] * hypercharge[i];

		if (fabs(Ysqr - 5.0/6.0) > prec)
		{
			const double factor = sqrt(6.0/5.0 * Ysqr);

			Y_Q  *= factor;
			Y_ub *= factor;
			Y_eb *= factor;
			Y_db *= factor;
			Y_L  *= factor;
			Y_Hu *= factor;
			Y_Hd *= factor;
		}
	}

	// first, it contains all different SM representations;
	// later, after pairing up, it contains the net numbers of all SM representations
	//        (except for the neutrinos) plus possible chiral exotics
	vector<unsigned>     SM_Mults;
	vector<vector<int> > SM_Reps;
	vector<double>       SM_Ys;

	// contains the net number of all representations, charged under the SU(3) x SU(2) x U(1)_y gauge group
	// and the net number of SM singlets (1,1)_0
	vector<int>          NetNumber_SM_Mults;
	vector<vector<int> > NetNumber_SM_Reps;
	vector<double>       NetNumber_SM_Ys;

	size_t SM_size = 0;
	size_t NN_size = 0;

	vector<int> tmp_SM_Rep(2,1); // contains SU(3) x SU(2) reps.

	vector<int> HiggsRep(2,1);
	HiggsRep[1] = 2;

	const vector<int> neutrino(2,1);
	const vector<int> exotic_lepton = HiggsRep;

	unsigned j = 0;
	unsigned k = 0;

	//unsigned Mult    = 0;
	unsigned SM_Mult = 0;
	double   Y              = 0.0;
	double   SM_Y           = 0.0;
	double   NetNumber_SM_Y = 0.0;

	//hacking here!!! to search for SM higgs doublet
	//get scalars
	// begin: find all different SU(3) x SU(2) x U(1)_y representations + singlets
	{
		unsigned tmp_mult     = 1;
		bool     SM_Rep_known = false;

		for (i = 0; i < f1; ++i)
		{
			const CField &Field = Fields[i];
			if (Field.Multiplet == Scalar)
			{
				const RepVector &Dimensions = Field.Dimensions;
				Y = Field.U1Charges[pos_of_Y];

				tmp_SM_Rep[0] = Dimensions[pos_of_SU3].Dimension;
				tmp_SM_Rep[1] = Dimensions[pos_of_SU2].Dimension;

				tmp_mult = 1;
				for (k = 0; k < number_of_factors; ++k)
				{
					if ((k != pos_of_SU3) && (k != pos_of_SU2))
						tmp_mult *= abs(Dimensions[k].Dimension);
				}

				SM_Rep_known = false;

				SM_size = SM_Reps.size();
				for (j = 0; !SM_Rep_known && (j < SM_size); ++j)
				{
					// if SU(3) x SU(2) x U(1)_y rep is known, increase the counter
					if ((SM_Reps[j] == tmp_SM_Rep) && (fabs(SM_Ys[j] - Y) < prec))
					{
						SM_Rep_known = true;
						SM_Mults[j] += tmp_mult;
					}
				}
				// if SU(3) x SU(2) x U(1)_y rep is not known, save it
				if (!SM_Rep_known)
				{
					SM_Mults.push_back(tmp_mult);
					SM_Reps.push_back(tmp_SM_Rep);
					SM_Ys.push_back(Y);
				}
			}
		}
	}
	SM_size = SM_Mults.size();
	// end: find all different SU(3) x SU(2) x U(1)_y representations + singlets

	//begin: find absolute number of Higgs doublets
	{
		HiggsNo=0;
		for ( i=0; i<SM_Mults.size(); i++ ) {
			if ((SM_Reps[i] == HiggsRep) && ((fabs(SM_Ys[i] - Y_Hu) < prec) or (fabs(SM_Ys[i] - Y_Hd) < prec)))
				HiggsNo += SM_Mults[i];
		}
	}
	if (HiggsNo==0)
	{
		if (info)
			cout << "Info SM_CheckVectorlikeness: Number of Higgs doublet is zero.\n";
		return false;
	}
	//end: find absolute number of Higgs doublets

	/*
  // begin: print
  cout << "SU(3) x SU(2) x U(1)_y representations:\n";
  for (i = 0; i < SM_size; ++i)
  cout << SM_Mults[i] << " (" << SM_Reps[i][0] << ", " << SM_Reps[i][1] << ")_" << SM_Ys[i] << "\n";
  cout << endl;
  cout<<"Absolute number of Higgs doublets = "<<HiggsNo<<endl;
  // end: print

	// begin: compute the net number of representations
	{
		bool NetNumber_SM_Rep_known = false;
		int Net_Number_Higgs=0;

		for (i = 0; i < SM_size; ++i)
		{
			const vector<int> &SM_Rep = SM_Reps[i];
			SM_Mult = SM_Mults[i];
			SM_Y    = SM_Ys[i];

			// contains the cc-partner of "Known_SM_Rep"
			tmp_SM_Rep = SM_Rep;

			if (abs(tmp_SM_Rep[0]) == 3)
				tmp_SM_Rep[0] *= -1;

			NetNumber_SM_Rep_known = false;

			NN_size = NetNumber_SM_Mults.size();
			for (j = 0; !NetNumber_SM_Rep_known && (j < NN_size); ++j)
			{
				const vector<int> &NetNumber_SM_Rep = NetNumber_SM_Reps[j];
				NetNumber_SM_Y = NetNumber_SM_Ys[j];

				if ((NetNumber_SM_Rep == SM_Rep) && (fabs(NetNumber_SM_Y - SM_Y) < prec))
				{
					NetNumber_SM_Rep_known = true;
					NetNumber_SM_Mults[j] += SM_Mult;
				}
				else
					if ((NetNumber_SM_Rep == tmp_SM_Rep) && (fabs(NetNumber_SM_Y + SM_Y) < prec))
					{

						NetNumber_SM_Rep_known = true;
						NetNumber_SM_Mults[j] -= SM_Mult;
					}
			}
			if (!NetNumber_SM_Rep_known)
			{
				NetNumber_SM_Mults.push_back(SM_Mult);
				NetNumber_SM_Reps.push_back(SM_Rep);
				NetNumber_SM_Ys.push_back(SM_Y);
			}
		}
	}
	NN_size = NetNumber_SM_Mults.size();
	// end: compute the net number of representations of scalars

	// begin: create the net numbers of all scalar representations
	{
		SM_Mults.clear();
		SM_Reps.clear();
		SM_Ys.clear();

		int tmp_mult = 1;
		for (i = 0; i < NN_size; ++i)
		{
			tmp_mult = NetNumber_SM_Mults[i];
			//if (tmp_mult != 0)							//hacking here!!!
			//{
				tmp_SM_Rep = NetNumber_SM_Reps[i];
				if (tmp_mult < 0)
				{
					if (abs(tmp_SM_Rep[0]) == 3)
						tmp_SM_Rep[0] *= -1;

					SM_Mults.push_back(abs(tmp_mult));
					SM_Reps.push_back(tmp_SM_Rep);
					SM_Ys.push_back(((-1) * NetNumber_SM_Ys[i]));
				}
				else
				{
					SM_Mults.push_back(tmp_mult);
					SM_Reps.push_back(tmp_SM_Rep);
					SM_Ys.push_back(NetNumber_SM_Ys[i]);
				}
			//}
		}
	}
	// end: create the net numbers of all scalar representations
	SM_size = SM_Mults.size();

	// begin: print
	cout << "NN scalar Spectrum:\n";
	for (i = 0; i < SM_size; ++i)
		cout << SM_Mults[i] << " (" << SM_Reps[i][0] << ", " << SM_Reps[i][1] << ")_" << SM_Ys[i] << "\n";
	cout << endl;
	// end: print

	//begin: determine net number of Higgs doublet
	{
		for (i = 0; i < SM_size; ++i)
		{
			const vector<int> &SM_Rep = SM_Reps[i];
			SM_Mult = SM_Mults[i];
			SM_Y=SM_Ys[i];
			// higgs-doublet (1,2)_-1/2 or (1,2)_+1/2
			if (SM_Rep == HiggsRep)
			{
				//cout<<SM_Y<<endl;
				if ( (fabs(SM_Y - Y_Hu) < prec) or (fabs(SM_Y - Y_Hd) < prec) )
				{
					//cout<<SM_Mult<<endl;
					if (SM_Mult == 0) {
						if (info)
							cout << "Info SM_CheckVectorlikeness: Net Number of Higgs doublet is zero.\n";
						return false;
					}
					else if (SM_Mult == 2)
					{
						doubleHiggs = true;
					}
					else if (SM_Mult > 2){
						if (info)
							cout << "Info SM_CheckVectorlikeness: Net Number of Higgs doublet is greater than two.\n";
						return false;
					}
				}
			}
		}
	}
	//end: determine net number of Higgs doublet
	 */


	SM_Mults.clear();
	SM_Reps.clear();
	SM_Ys.clear();
	NetNumber_SM_Mults.clear();
	NetNumber_SM_Reps.clear();
	NetNumber_SM_Ys.clear();


	//get fermions
	// begin: find all different SU(3) x SU(2) x U(1)_y representations + singlets
	{
		unsigned tmp_mult     = 1;
		bool     SM_Rep_known = false;

		for (i = 0; i < f1; ++i)
		{
			const CField &Field = Fields[i];
			if (Field.Multiplet == LeftFermi)
			{
				const RepVector &Dimensions = Field.Dimensions;
				Y = Field.U1Charges[pos_of_Y];

				tmp_SM_Rep[0] = Dimensions[pos_of_SU3].Dimension;
				tmp_SM_Rep[1] = Dimensions[pos_of_SU2].Dimension;

				tmp_mult = 1;
				for (k = 0; k < number_of_factors; ++k)
				{
					if ((k != pos_of_SU3) && (k != pos_of_SU2))
						tmp_mult *= abs(Dimensions[k].Dimension);
				}

				SM_Rep_known = false;

				SM_size = SM_Reps.size();
				for (j = 0; !SM_Rep_known && (j < SM_size); ++j)
				{
					// if SU(3) x SU(2) x U(1)_y rep is known, increase the counter
					if ((SM_Reps[j] == tmp_SM_Rep) && (fabs(SM_Ys[j] - Y) < prec))
					{
						SM_Rep_known = true;
						SM_Mults[j] += tmp_mult;
					}
				}
				// if SU(3) x SU(2) x U(1)_y rep is not known, save it
				if (!SM_Rep_known)
				{
					SM_Mults.push_back(tmp_mult);
					SM_Reps.push_back(tmp_SM_Rep);
					SM_Ys.push_back(Y);
				}
			}
		}
	}
	SM_size = SM_Mults.size();
	// end: find all different SU(3) x SU(2) x U(1)_y representations + singlets

	/*
  // begin: print
  cout << "SU(3) x SU(2) x U(1)_y representations:\n";
  for (i = 0; i < SM_size; ++i)
  cout << SM_Mults[i] << " (" << SM_Reps[i][0] << ", " << SM_Reps[i][1] << ")_" << SM_Ys[i] << "\n";
  cout << endl;
  // end: print
	 */

	// begin: compute the net number of representations
	{
		bool NetNumber_SM_Rep_known = false;
		bool u_d_Higgs              = false;

		for (i = 0; i < SM_size; ++i)
		{
			const vector<int> &SM_Rep = SM_Reps[i];
			SM_Mult = SM_Mults[i];
			SM_Y    = SM_Ys[i];

			// contains the cc-partner of "Known_SM_Rep"
			tmp_SM_Rep = SM_Rep;

			if (abs(tmp_SM_Rep[0]) == 3)
				tmp_SM_Rep[0] *= -1;

			NetNumber_SM_Rep_known = false;

			NN_size = NetNumber_SM_Mults.size();
			for (j = 0; !NetNumber_SM_Rep_known && (j < NN_size); ++j)
			{
				const vector<int> &NetNumber_SM_Rep = NetNumber_SM_Reps[j];
				NetNumber_SM_Y = NetNumber_SM_Ys[j];

				if ((NetNumber_SM_Rep == SM_Rep) && (fabs(NetNumber_SM_Y - SM_Y) < prec))
				{
					NetNumber_SM_Rep_known = true;
					NetNumber_SM_Mults[j] += SM_Mult;
				}
				else
					if ((NetNumber_SM_Rep == tmp_SM_Rep) && (fabs(NetNumber_SM_Y + SM_Y) < prec))
					{
						// look for the up- and down-Higgs (up and down Higgs will pair up here)   //hacking here!!! higgs will be checked in the scalars section

						NetNumber_SM_Rep_known = true;
						NetNumber_SM_Mults[j] -= SM_Mult;
					}
			}
			if (!NetNumber_SM_Rep_known)
			{
				NetNumber_SM_Mults.push_back(SM_Mult);
				NetNumber_SM_Reps.push_back(SM_Rep);
				NetNumber_SM_Ys.push_back(SM_Y);
			}
		}
	}
	NN_size = NetNumber_SM_Mults.size();
	// end: compute the net number of representations

	/*
  // begin: print
  cout << "Net number of SU(3) x SU(2) x U(1)_y representations\n";
  for (i = 0; i < NN_size; ++i)
  cout << NetNumber_SM_Mults[i] << " (" << NetNumber_SM_Reps[i][0] << ", " << NetNumber_SM_Reps[i][1] << ")_" << NetNumber_SM_Ys[i] << "\n";
  cout << endl;
  // end: print
	 */

	// begin: the net number of neutrinos (1,1)_0 must be larger or equal NumberOfGenerations
	//        and the net number of exotic lepton-doublets (1,2)_0 must be even
	//        (due to Witten's anomaly, if net number of exotic lepton-doublets (1,2)_0 is odd,
	//         there must be an odd number of exotic lepton-doublets with non-vanishing hypercharge
	//         -> chiral exotic)
	{
		bool exotic_lepton_found = false;
		bool neutrino_found      = false;

		for (i = 0; !(neutrino_found && exotic_lepton_found) && (i < NN_size); ++i)
		{
			// neutrino
			if ((NetNumber_SM_Reps[i] == neutrino) && (fabs(NetNumber_SM_Ys[i] - Y_n) < prec))
			{
				if (NetNumber_SM_Mults[i] < NumberOfGenerations)
				{
					// if (info)
					cout << "Info SM_CheckVectorlikeness: Net number of neutrinos is not larger than " << NumberOfGenerations << "!\n";
					return false;
				}
				neutrino_found = true;
				NetNumber_SM_Mults[i] = 0;
			}
			else
				// exotic lepton doublet
				if ((NetNumber_SM_Reps[i] == exotic_lepton) && (fabs(NetNumber_SM_Ys[i]) < prec))
				{
					if (NetNumber_SM_Mults[i] % 2 != 0)
					{
						if (info)
							cout << "Info SM_CheckVectorlikeness: Net number of exotic lepton-doublets (1,2)_0 is not even!\n";
						return false;		//hacking here!!! disable if needed
					}
					exotic_lepton_found = true;
					NetNumber_SM_Mults[i] = 0;
				}
		}
		if (!neutrino_found)
		{
			if (info)
				cout << "Info SM_CheckVectorlikeness: Neutrinos (1,1)_0 not found!\n";
			return false;
		}
	}
	// end: the net number of neutrinos (1,1)_0 must be larger or equal NumberOfGenerations
	//      and the net number of exotic lepton-doublets (1,2)_0 must be even

	/*
  // begin: print
  cout << "Net number after checking the neutrinos and the exotic lepton-doublets:\n";
  for (i = 0; i < NN_size; ++i)
  cout << NetNumber_SM_Mults[i] << " (" << NetNumber_SM_Reps[i][0] << ", " << NetNumber_SM_Reps[i][1] << ")_" << NetNumber_SM_Ys[i] << "\n";
  cout << endl;
  //end: print
	 */

	// begin: create the net numbers of all SM representations
	//        (except for the neutrinos) plus possible chiral exotics
	{
		SM_Mults.clear();
		SM_Reps.clear();
		SM_Ys.clear();

		int tmp_mult = 1;
		for (i = 0; i < NN_size; ++i)
		{
			tmp_mult = NetNumber_SM_Mults[i];
			if (tmp_mult != 0)
			{
				tmp_SM_Rep = NetNumber_SM_Reps[i];
				if (tmp_mult < 0)
				{
					if (abs(tmp_SM_Rep[0]) == 3)
						tmp_SM_Rep[0] *= -1;

					SM_Mults.push_back(abs(tmp_mult));
					SM_Reps.push_back(tmp_SM_Rep);
					SM_Ys.push_back(((-1) * NetNumber_SM_Ys[i]));
				}
				else
				{
					SM_Mults.push_back(tmp_mult);
					SM_Reps.push_back(tmp_SM_Rep);
					SM_Ys.push_back(NetNumber_SM_Ys[i]);
				}
			}
		}
		SM_size = SM_Mults.size();

		// begin: print
		/*cout << "NN Spectrum:\n";
		for (i = 0; i < SM_size; ++i)
			cout << SM_Mults[i] << " (" << SM_Reps[i][0] << ", " << SM_Reps[i][1] << ")_" << SM_Ys[i] << "\n";
		cout << endl;*/
		// end: print

		if (SM_size != 5) //hacking here!!!
		{
			if (info)
			cout << "Info SM_CheckVectorlikeness: Model has chiral exotics (paired-up spectrum must contain exactly 5 SM-representations + neutrinos).\n";
			return false;		//hacking here!!! disable if needed
		}
	}
	// end: create the net numbers of all SM representations
	//      (except for the neutrinos) plus possible chiral exotics


	// begin: find the quark doublet
	//        and check wether it is a (3,2) rep. or a (-3,2) with charge 1/6 or -1/6
	int U1sign  = 1;
	int RepSign = 1;
	{
		bool find_3_2 = false;

		for(i = 0; !find_3_2 && (i < SM_size); ++i)
		{
			const vector<int> &SM_Rep = SM_Reps[i];

			if ((abs(SM_Rep[0]) == 3) && (SM_Rep[1] == 2))
			{
				if (SM_Rep[0] == -3)
					RepSign =-1;
				else
					RepSign = 1;

				SM_Y = SM_Ys[i];
				if (fabs(fabs(SM_Y) - Y_Q) < prec)
				{
					find_3_2 = true;

					if (fabs(SM_Y - Y_Q) < prec)
						U1sign = 1;
					else
						U1sign =-1;

					if (SM_Mults[i] != NumberOfGenerations)
					{
						if (info)
							cout << "Info SM_CheckVectorlikeness: Model has not a net number of " << NumberOfGenerations << " generations of quark-doublets.\n";
						return false;
					}
				}
			}
		}
		if (!find_3_2)
		{
			if (info)
				cout << "Info SM_CheckVectorlikeness: No quark-doublet found.\n";
			return false;
		}
	}
	// end: find the quark doublet
	//      and check wether it is a (3,2) rep. or a (-3,2) with charge 1/6 or -1/6

	// begin: final check for vectorlikeness
	bool ub_found = false;
	bool db_found = false;
	bool eb_found = false;
	bool Q_found  = false;
	bool L_found  = false;

	for (i = 0; i < SM_size; ++i)
	{
		const vector<int> &SM_Rep = SM_Reps[i];
		SM_Mult = SM_Mults[i];
		SM_Y    = U1sign * SM_Ys[i];

		// (x,2)
		if (SM_Rep[1] == 2)
		{
			// quark-doublet (3,2)_1/6
			if ((RepSign * SM_Rep[0]) == 3)
			{
				if (fabs(SM_Y - Y_Q) < prec)
				{
					if (SM_Mult == NumberOfGenerations)
						Q_found = true;
					else
						Q_found = false;
				}
			}
			else
				// lepton-doublet (1,2)_-1/2
				if (SM_Rep[0] == 1)
				{
					if (fabs(SM_Y - Y_L) < prec)
					{
						if (SM_Mult == NumberOfGenerations)
							L_found = true;
						else
							L_found = false;
					}
				}
		}
		else
			// (x,1)
			if (SM_Rep[1] == 1)
			{
				// (1,1)
				if (SM_Rep[0] == 1)
				{
					// electron (1,1)_1
					if (fabs(SM_Y - Y_eb) < prec)
					{
						if (SM_Mult == NumberOfGenerations)
							eb_found = true;
						else
							eb_found = false;
					}
				}
				else
					// (-3,1)
					if ((RepSign * SM_Rep[0]) == -3)
					{
						// up quark (-3,1)_-2/3
						if (!ub_found && (fabs(SM_Y - Y_ub) < prec))
						{
							if (SM_Mult == NumberOfGenerations)
								ub_found = true;
							else
								ub_found = false;
						}

						// down quark (-3,1)_1/3
						if (!db_found && (fabs(SM_Y - Y_db) < prec))
						{
							if (SM_Mult == NumberOfGenerations)
								db_found = true;
							else
								db_found = false;
						}
					}
			}
	}
	if (!(Q_found && L_found && eb_found && ub_found && db_found))
	{
		if (info)
			cout << "Info SM_CheckVectorlikeness: Model has chiral exotics.\n";
		return false;
	}
	else
	{
		if (info)
			cout << "Model is Standard Model plus vectorlike exotics!" << endl;
		return true;
	}
	// end: final check for vectorlikeness

	return true;
}



/* ########################################################################################
######   SM_GetProtonHexalityFromSO12(...) const                                     ######
######                                                                               ######
######   Version: 13.04.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Orbifold         : the orbifold of the vev-config "SMVEVConfig"          ######
######   2) SMVEVConfig      : a Standard Model vev-config                           ######
######   3) Roots10D         : 480 roots of E_8 x E_8 or SO(32) used to identify an  ######
######                         intermediate SO(12)                                   ######
######   4) ProtonHexalities : possible generators for Proton hexality               ######
######   output:                                                                     ######
######   return value        : proton hexality identified?                           ######
###########################################################################################
######   description:                                                                ######
######   Analyze the vev-config "SMVEVConfig" to identify possible proton hexality   ######
######   generators.                                                                 ######
######################################################################################## */
bool CAnalyseModel::SM_GetProtonHexalityFromSO12(const COrbifold &Orbifold, const SConfig &SMVEVConfig, const vector<CVector> &Roots10D, vector<CVector> &ProtonHexalities) const
{
	const double prec = 0.0001;
	int i = 0, j = 0, k = 0 , l = 0, m = 0, n = 0;

	if (ProtonHexalities.size() != 0)
	{
		cout << "\n  Warning in bool CAnalyseModel::SM_GetProtonHexalityFromSO12(...) const: ProtonHexalities not empty - now cleared." << endl;
		ProtonHexalities.clear();
	}

	const vector<gaugeGroupFactor<double> > &factors  = SMVEVConfig.SymmetryGroup.GaugeGroup.factor;
	const size_t number_of_factors = factors.size();

	if ((SMVEVConfig.SymmetryGroup.observable_sector_GGs.size() < 2) || (SMVEVConfig.SymmetryGroup.observable_sector_U1s.size() == 0))
	{
		cout << "\n  Warning in bool CAnalyseModel::SM_GetProtonHexalityFromSO12(...) const: Cannot find the Standard Model. Return false." << endl;
		return false;
	}

	const unsigned pos_of_SU3 = SMVEVConfig.SymmetryGroup.observable_sector_GGs[0];
	const unsigned pos_of_SU2 = SMVEVConfig.SymmetryGroup.observable_sector_GGs[1];

	const unsigned pos_of_Y = SMVEVConfig.SymmetryGroup.observable_sector_U1s[0];

	if ((pos_of_SU3 >= number_of_factors) || (pos_of_SU2 >= number_of_factors))
	{
		cout << "\n  Warning in bool CAnalyseModel::SM_GetProtonHexalityFromSO12(...) const: Position of SU(3) or SU(2) out of range. Return false." << endl;
		return false;
	}

	if (Roots10D.size() != 480)
	{
		cout << "\n  Warning in bool CAnalyseModel::SM_GetProtonHexalityFromSO12(...) const: \"Roots10D\" does not contain 480 roots. Return false." << endl;
		return false;
	}

	const vector<vector<double> > &SU3_Simpleroots = factors[pos_of_SU3].simpleroots;
	const vector<vector<double> > &SU2_Simpleroots = factors[pos_of_SU2].simpleroots;

	if ((SU3_Simpleroots.size() != 2) || (SU2_Simpleroots.size() != 1))
	{
		cout << "\n  Warning in bool CAnalyseModel::SM_GetProtonHexalityFromSO12(...) const: Simple roots of SU(3) or SU(2) not correct. Return false." << endl;
		return false;
	}

	CVector U1DirectionA(16);
	U1DirectionA = SMVEVConfig.SymmetryGroup.GaugeGroup.u1directions[0];
	const vector<double> &U1DirectionY = SMVEVConfig.SymmetryGroup.GaugeGroup.u1directions[pos_of_Y];
	const bool AnomalousU1 = SMVEVConfig.SymmetryGroup.IsFirstU1Anomalous;

	CVector SU3_Simpleroot1, SU3_Simpleroot2, SU2_Simpleroot1;
	SU3_Simpleroot1 = SU3_Simpleroots[0];
	SU3_Simpleroot2 = SU3_Simpleroots[1];
	SU2_Simpleroot1 = SU2_Simpleroots[0];

	//cout << "SU(3) simple root 1" << endl;
	//SU3_Simpleroot1.Print(cout, E8xE8);
	//cout << endl;
	//cout << "SU(3) simple root 2" << endl;
	//SU3_Simpleroot2.Print(cout, E8xE8);
	//cout << endl;
	//cout << "SU(2) simple root" << endl;
	//SU2_Simpleroot1.Print(cout, E8xE8);
	//cout << endl;
	//cout << "U(1)_Y" << endl;
	//for (unsigned q = 0; q < 16; ++q)
	//  cout << U1DirectionY[q] << " ";
	//cout << endl;

	vector<doubleVector> UnbrokenRoots;
	vector<CVector> Copy_Roots10D;
	CGaugeGroup localGaugeGroup;
	size_t g1 = 0;
	size_t s2 = 0;

	CVector tmpVector;
	CVector ProHex;
	CVector tmpVectorX;
	CVector tmpVectorZ;

	const COrbifoldGroup &OrbifoldGroup = Orbifold.OrbifoldGroup;
	const CShiftVector &ZM_Shift = OrbifoldGroup.GetShift(0);
	const CShiftVector &ZN_Shift = OrbifoldGroup.GetShift(1);

	unsigned M = OrbifoldGroup.GetSpaceGroup().GetM();
	unsigned N = 1;

	if (OrbifoldGroup.GetSpaceGroup().IsZMxZN())
		N = OrbifoldGroup.GetSpaceGroup().GetN();

	CVector localShift(16);
	for (i = 0; i < M; ++i)
	{
		for (j = 0; j < N; ++j)
		{
			// create the local shift of the (i,j) twisted sector
			localShift = (ZM_Shift * i) + (ZN_Shift * j);
			Copy_Roots10D = Roots10D;

			UnbrokenRoots.clear();
			for (k = 0; k < 480; ++k)
			{
				if (fabs(localShift * Roots10D[k]) < prec)
					UnbrokenRoots.push_back(Roots10D[k]);
			}

			// create the gauge group of the (i,j) twisted sector
			localGaugeGroup = determineAlgebra(UnbrokenRoots);
			g1 = localGaugeGroup.factor.size();

			// search for SO(12) in the gauge group of the (i,j) twisted sector
			for (k = 0; k < g1; ++k)
			{
				if (localGaugeGroup.factor[k].algebra == "D6")
				{
					const vector<vector<double> > &SO12SimpleRoots = localGaugeGroup.factor[k].simpleroots;

					if (SO12SimpleRoots.size() != 6)
					{
						cout << "\n  Warning in bool CAnalyseModel::SM_GetProtonHexalityFromSO12(...) const: Simple roots of SO(12) not correct. Return false." << endl;
						return false;
					}

					bool SM_inside = false;

					// begin: is SU(3) x SU(2) inside this SO(12)?
					bool SP_nontrivial = false;

					// SU(3) simple root 1
					for (l = 0; l < 6; ++l)
					{
						if (fabs(SO12SimpleRoots[l] * SU3_Simpleroot1) > prec)
							SP_nontrivial = true;
					}
					if (SP_nontrivial)
					{
						SP_nontrivial = false;

						// SU(3) simple root 2
						for (l = 0; l < 6; ++l)
						{
							if (fabs(SO12SimpleRoots[l] * SU3_Simpleroot2) > prec)
								SP_nontrivial = true;
						}

						if (SP_nontrivial)
						{
							SP_nontrivial = false;

							// SU(2) simple root
							for (l = 0; l < 6; ++l)
							{
								if (fabs(SO12SimpleRoots[l] * SU2_Simpleroot1) > prec)
									SP_nontrivial = true;
							}
							if (SP_nontrivial)
								SM_inside = true;
						}
					}
					// end: is SU(3) x SU(2) inside this SO(12)?

					if (SM_inside)
					{
						/*cout << "SO(12) simple roots" << endl;
            for (unsigned q = 0; q < 6; ++q)
            {
              CVector tmp;
              tmp = SO12SimpleRoots[q];
              tmp.Print(cout, E8xE8);
              cout << endl;
            }*/

						// begin: find the directions inside SO(12) that are orthogonal to SU(3) x SU(2) x U(1)_Y
						vector<doubleVector> SO12_Orthogonal;
						Find_Basis_Of_Orthogonal_Space(SO12SimpleRoots, UNSPECIFIED_LATTICE, 16, SO12_Orthogonal);

						SO12_Orthogonal.insert(SO12_Orthogonal.end(), SU3_Simpleroots.begin(), SU3_Simpleroots.end());
						SO12_Orthogonal.insert(SO12_Orthogonal.end(), SU2_Simpleroots.begin(), SU2_Simpleroots.end());
						SO12_Orthogonal.push_back(U1DirectionY);

						vector<doubleVector> SM_Orthogonal;
						Find_Basis_Of_Orthogonal_Space(SO12_Orthogonal, UNSPECIFIED_LATTICE, 16, SM_Orthogonal);

						if (SM_Orthogonal.size() != 2)
						{
							cout << "\n  Warning in bool CAnalyseModel::SM_GetProtonHexalityFromSO12(...) const: Number of U(1)s from SO(12) incorrect. Return false." << endl;
							return false;
						}
						// end: find the directions inside SO(12) that are orthogonal to SU(3) x SU(2) x U(1)_Y

						tmpVectorX = SM_Orthogonal[0];
						tmpVectorZ = SM_Orthogonal[1];

						/*cout << "U(1)_X" << endl;
            tmpVectorX.Print(cout, E8xE8);
            cout << endl;
            cout << "U(1)_Z" << endl;
            tmpVectorZ.Print(cout, E8xE8);
            cout << endl;*/

						vector<doubleVector> coeff;
						if (AnomalousU1)
						{
							doubleVector a;
							a.push_back(tmpVectorX * U1DirectionA);
							a.push_back(tmpVectorZ * U1DirectionA);
							vector<doubleVector> A(1,a);

							if ((fabs(a[0]) < prec) && (fabs(a[1]) < prec))
							{
								/*cout << "Both ProHex-U(1)s are non-anomalous - Using 4 linear combinations of the two U(1)s:" << endl;
                CPrint Print(Tstandard, &cout);
                Print.PrintVector(tmpVectorX, E8xE8);
                cout << endl;
                Print.PrintVector(tmpVectorZ, E8xE8);
                cout << "\n" << endl;*/

								doubleVector c1(2,1.0);
								coeff.push_back(c1);
								c1[1] = -1.0;
								coeff.push_back(c1);

								c1[0] = 1.0;
								c1[1] = 0.0;
								coeff.push_back(c1);
								c1[0] = 0.0;
								c1[1] = 1.0;
								coeff.push_back(c1);
							}
							else
							{
								Find_Basis_Of_Orthogonal_Space(A, UNSPECIFIED_LATTICE, 2, coeff);

								if (coeff.size() != 1)
								{
									cout << "\n  Warning in bool CAnalyseModel::SM_GetProtonHexalityFromSO12(...) const: Check number of ProHex candidates. Return false." << endl;
									return false;
								}
							}
						}
						else
						{
							/*cout << "No anomalous U(1). Using 4 linear combinations of the two U(1)s:" << endl;
              CPrint Print(Tstandard, &cout);
              Print.PrintVector(tmpVectorX, E8xE8);
              cout << endl;
              Print.PrintVector(tmpVectorZ, E8xE8);
              cout << "\n" << endl;*/

							doubleVector c1(2,1.0);
							coeff.push_back(c1);
							c1[1] = -1.0;
							coeff.push_back(c1);

							c1[0] = 1.0;
							c1[1] = 0.0;
							coeff.push_back(c1);
							c1[0] = 0.0;
							c1[1] = 1.0;
							coeff.push_back(c1);
						}

						for (n = 0; n < coeff.size(); ++n)
						{
							ProHex = ((tmpVectorX * coeff[n][0]) + (tmpVectorZ * coeff[n][1])) * 0.5;

							bool orthogonal = true;

							for (l = 0; orthogonal && (l < number_of_factors); ++l)
							{
								const vector<vector<double> > &Simpleroots = factors[l].simpleroots;
								s2 = Simpleroots.size();
								for (m = 0; orthogonal && (m < s2); ++m)
								{
									tmpVector = Simpleroots[m];
									if (fabs(ProHex * tmpVector) > prec)
										orthogonal = false;
								}
							}
							if (!orthogonal)
							{
								cout << "\n  Warning in bool CAnalyseModel::SM_GetProtonHexalityFromSO12(...) const: ProHex candidate is not orthogonal to all simple roots. Return false." << endl;
								return false;
							}

							if (orthogonal && AnomalousU1 && (fabs(U1DirectionA * ProHex) > prec))
							{
								cout << "\n  Warning in bool CAnalyseModel::SM_GetProtonHexalityFromSO12(...) const: ProHex candidate is not orthogonal to the anomalous U(1). Return false." << endl;
								return false;
							}

							if (orthogonal && (find(ProtonHexalities.begin(), ProtonHexalities.end(), ProHex) == ProtonHexalities.end()))
							{
								/*cout << "ProHex candidate" << endl;
                CPrint Print(Tstandard, &cout);
                Print.PrintVector(ProHex, E8xE8);
                cout << endl;*/

								ProtonHexalities.push_back(ProHex);
							}
						}
					}
				}
			}
		}
	}
	return (ProtonHexalities.size() != 0);
}



/* ########################################################################################
######   PS_NetNumberOfGenerations(const SConfig &PSVEVConfig) const                 ######
######                                                                               ######
######   Version: 28.01.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) PSVEVConfig : a Pati-Salam vev-config for which the net numbers of       ######
######                    (4,2,1)- and (bar{4},1,2) plets are counted                ######
######   output:                                                                     ######
######   return value   : net number of Pati-Salam generations                       ######
###########################################################################################
######   description:                                                                ######
######   Counts the net number of Pati-Salam generations in the vev-config           ######
######   "PSVEVConfig".                                                              ######
######################################################################################## */
int CAnalyseModel::PS_NetNumberOfGenerations(const SConfig &PSVEVConfig) const
{
	int net_no_familiesL = 0;
	int net_no_familiesR = 0;
	int mult             = 0;

	unsigned i = 0;
	unsigned j = 0;

	const vector<CField> &Fields = PSVEVConfig.Fields;
	const size_t f1 = Fields.size();
	if (f1 == 0)
	{
		cout << "\n  Warning in int CAnalyseModel::PS_NetNumberOfGenerations(...): Spectrum is empty. Return false." << endl;
		return false;
	}

	const vector<gaugeGroupFactor<double> > &factors = PSVEVConfig.SymmetryGroup.GaugeGroup.factor;
	const size_t number_of_factors = factors.size();

	if (PSVEVConfig.SymmetryGroup.observable_sector_GGs.size() < 3)
	{
		cout << "\n  Warning in int CAnalyseModel::PS_NetNumberOfGenerations(...): Cannot find the Pati-Salam group. Return false." << endl;
		return false;
	}

	const unsigned pos_of_SU4  = PSVEVConfig.SymmetryGroup.observable_sector_GGs[0];
	const unsigned pos_of_SU2L = PSVEVConfig.SymmetryGroup.observable_sector_GGs[1];
	const unsigned pos_of_SU2R = PSVEVConfig.SymmetryGroup.observable_sector_GGs[2];

	if ((pos_of_SU4 >= number_of_factors) || (pos_of_SU2L >= number_of_factors) || (pos_of_SU2R >= number_of_factors))
	{
		cout << "\n  Warning in int CAnalyseModel::PS_NetNumberOfGenerations(...): Cannot find the Pati-Salam group. Return false." << endl;
		return false;
	}

	if ((factors[pos_of_SU4].algebra != "A3") || (factors[pos_of_SU2L].algebra != "A1") || (factors[pos_of_SU2R].algebra != "A1"))
	{
		cout << "\n  Warning in int CAnalyseModel::PS_NetNumberOfGenerations(...): Cannot find the Pati-Salam group. Return false." << endl;
		return false;
	}

	for (i = 0; i < f1; ++i)
	{
		const CField &Field = Fields[i];
		if (Field.Multiplet == LeftChiral)
		{
			const RepVector &Dimensions = Field.Dimensions;

			if (abs(Dimensions[pos_of_SU4].Dimension) == 4)
			{
				// count the net number of (4,2,1)-plets
				if ((Dimensions[pos_of_SU2L].Dimension == 2) && (Dimensions[pos_of_SU2R].Dimension == 1))
				{
					mult = 1;
					for (j = 0; j < number_of_factors; ++j)
					{
						if (j == pos_of_SU4)
							mult *= Dimensions[j].Dimension;
						else
							mult *= abs(Dimensions[j].Dimension);
					}
					net_no_familiesL += mult/8;
				}
				// count the net number of (-4,1,2)-plets
				if ((Dimensions[pos_of_SU2L].Dimension == 1) && (Dimensions[pos_of_SU2R].Dimension == 2))
				{
					mult = 1;
					for (j = 0; j < number_of_factors; ++j)
					{
						if (j == pos_of_SU4)
							mult *= Dimensions[j].Dimension;
						else
							mult *= abs(Dimensions[j].Dimension);
					}
					net_no_familiesR += mult/8;
				}
			}
		}
	}
	if (((net_no_familiesL * net_no_familiesR) > 0 ) || (abs(net_no_familiesL) != abs(net_no_familiesR)))
		return 0;

	return net_no_familiesL;
}



/* ########################################################################################
######   PS_CheckVectorlikeness(...) const                                           ######
######                                                                               ######
######   Version: 17.05.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) PSVEVConfig         : vev-config with Pati-Salam gauge group to analyze  ######
######   3) Print               : if "info" is true print the output to this CPrint  ######
######                            object                                             ######
######   3) NumberOfGenerations : number of Pati-Salam generations                   ######
######   4) info                : print info?                                        ######
######   output:                                                                     ######
######   return value           : spectrum is PS generations + vectorlike exotics?   ######
###########################################################################################
######   description:                                                                ######
######   Checks whether the vev-config "PSVEVConfig" has a net number of             ######
######   "NumberOfGenerations" PS generations plus vectorlike exotics.               ######
######################################################################################## */
bool CAnalyseModel::PS_CheckVectorlikeness(SConfig &PSVEVConfig, CPrint &Print, unsigned NumberOfGenerations, bool info) const
{
	unsigned i = 0;
	unsigned j = 0;

	const vector<CField> &Fields = PSVEVConfig.Fields;
	const size_t f1 = Fields.size();
	if (f1 == 0)
	{
		cout << "\n  Warning in bool CAnalyseModel::PS_CheckVectorlikeness(...): Spectrum is empty. Return false." << endl;
		return false;
	}

	const vector<gaugeGroupFactor<double> > &factors = PSVEVConfig.SymmetryGroup.GaugeGroup.factor;
	const size_t number_of_factors = factors.size();

	if (PSVEVConfig.SymmetryGroup.observable_sector_GGs.size() < 3)
	{
		cout << "\n  Warning in bool CAnalyseModel::PS_CheckVectorlikeness(...): Cannot find the Pati-Salam group. Return false." << endl;
		return false;
	}

	const unsigned pos_of_SU4  = PSVEVConfig.SymmetryGroup.observable_sector_GGs[0];
	const unsigned pos_of_SU2L = PSVEVConfig.SymmetryGroup.observable_sector_GGs[1];
	const unsigned pos_of_SU2R = PSVEVConfig.SymmetryGroup.observable_sector_GGs[2];

	if ((pos_of_SU4 >= number_of_factors) || (pos_of_SU2L >= number_of_factors) || (pos_of_SU2R >= number_of_factors))
	{
		cout << "\n  Warning in bool CAnalyseModel::PS_CheckVectorlikeness(...): Cannot find the Pati-Salam group. Return false." << endl;
		return false;
	}

	if ((factors[pos_of_SU4].algebra != "A3") || (factors[pos_of_SU2L].algebra != "A1") || (factors[pos_of_SU2R].algebra != "A1"))
	{
		cout << "\n  Warning in bool CAnalyseModel::PS_CheckVectorlikeness(...): Cannot find the Pati-Salam group. Return false." << endl;
		return false;
	}

	int number_of_421 = 0;
	int number_of_412 = 0;

	int number_of_m421 = 0;
	int number_of_m412 = 0;

	int number_of_411 = 0;
	int number_of_422 = 0;
	int number_of_m411 = 0; // must be equal to number_of_411 for getting mass
	int number_of_m422 = 0; // must be equal to number_of_422 for getting mass

	int number_of_121 = 0; // must be 0 or larger than 1 for getting mass
	int number_of_112 = 0; // must be 0 or larger than 1 for getting mass
	int number_of_122 = 0; // must be larger than 1 for the SM Higgs

	int number_of_621 = 0; // must be 0 or larger than 1 for getting mass
	int number_of_612 = 0; // must be 0 or larger than 1 for getting mass

	int mult = 0;

	for (i = 0; i < f1; ++i)
	{
		const CField &Field = Fields[i];
		if (Field.Multiplet == LeftChiral)
		{
			const RepVector &Dimensions = Field.Dimensions;

			mult = 1;
			for (j = 0; j < number_of_factors; ++j)
				mult *= abs(Dimensions[j].Dimension);

			// count number of (4,x,y)
			if (Dimensions[pos_of_SU4].Dimension == 4)
			{
				// count number of (4,1,y)
				if (Dimensions[pos_of_SU2L].Dimension == 1)
				{
					// count number of (4,1,1)
					if (Dimensions[pos_of_SU2R].Dimension == 1)
						number_of_411 += mult/4;

					// count number of (4,1,2)
					if (Dimensions[pos_of_SU2R].Dimension == 2)
						number_of_412 += mult/8;
				}
				else
					// count number of (4,2,y)
					if (Dimensions[pos_of_SU2L].Dimension == 2)
					{
						// count number of (4,2,1)
						if (Dimensions[pos_of_SU2R].Dimension == 1)
							number_of_421 += mult/8;

						// count number of (4,2,2)
						if (Dimensions[pos_of_SU2R].Dimension == 2)
							number_of_422 += mult/16;
					}
			}
			else
				// count number of (-4,x,y)
				if (Dimensions[pos_of_SU4].Dimension == -4)
				{
					// count number of (-4,1,y)
					if (Dimensions[pos_of_SU2L].Dimension == 1)
					{
						// count number of (-4,1,1)
						if (Dimensions[pos_of_SU2R].Dimension == 1)
							number_of_m411 += mult/4;

						// count number of (-4,1,2)
						if (Dimensions[pos_of_SU2R].Dimension == 2)
							number_of_m412 += mult/8;
					}
					else
						// count number of (-4,2,y)
						if (Dimensions[pos_of_SU2L].Dimension == 2)
						{
							// count number of (-4,2,1)
							if (Dimensions[pos_of_SU2R].Dimension == 1)
								number_of_m421 += mult/8;

							// count number of (-4,2,2)
							if (Dimensions[pos_of_SU2R].Dimension == 2)
								number_of_m422 += mult/16;
						}
				}
				else
					// count number of (1,x,y)
					if (Dimensions[pos_of_SU4].Dimension == 1)
					{
						// count number of (1,1,y)
						if (Dimensions[pos_of_SU2L].Dimension == 1)
						{
							// count number of (1,1,2)
							if (Dimensions[pos_of_SU2R].Dimension == 2)
								number_of_112 += mult/2;
						}
						else
							// count number of (1,2,y)
							if (Dimensions[pos_of_SU2L].Dimension == 2)
							{
								// count number of (1,2,1)
								if (Dimensions[pos_of_SU2R].Dimension == 1)
									number_of_121 += mult/2;
								else
									// count number of (1,2,2)
									if (Dimensions[pos_of_SU2R].Dimension == 2)
										number_of_122 += mult/4;
							}
					}
					else
						// count number of (6,x,y)
						if (Dimensions[pos_of_SU4].Dimension == 6)
						{
							// count number of (6,1,y)
							if (Dimensions[pos_of_SU2L].Dimension == 1)
							{
								// count number of (6,1,2)
								if (Dimensions[pos_of_SU2R].Dimension == 2)
									number_of_612 += mult/12;
							}
							else
								// count number of (6,2,y)
								if (Dimensions[pos_of_SU2L].Dimension == 2)
								{
									// count number of (6,2,1)
									if (Dimensions[pos_of_SU2R].Dimension == 1)
										number_of_621 += mult/12;
								}
						}
		}
	}

	if (number_of_411 != number_of_m411)
		return false;

	if (number_of_422 != number_of_m422)
		return false;

	if (number_of_122 == 0) // need at least one up- and one down-Higgs for electroweak symmetry breaking
	{
		if (info)
			(*Print.out) << "Info PS_CheckVectorlikeness: SM up- and down type Higgs not found!\n";
		return false;
	}
	if ((number_of_121 == 1) || (number_of_112 == 1) || (number_of_621 == 1) || (number_of_612 == 1))
	{
		if (info)
			(*Print.out) << "Info PS_CheckVectorlikeness: Vectorlike exotics (1,2,1), (1,1,2), (6,2,1) or (6,1,2) cannot be paired up due to antisymmetry of 2x2 = 1_a!\n";
		return false;
	}

	if (abs(number_of_421 - number_of_m421) != NumberOfGenerations) // net number of NumberOfGenerations generations
	{
		if (info)
			(*Print.out) << "Info PS_CheckVectorlikeness: Net number of (4,2,1) is " << abs(number_of_421 - number_of_m421) << "!\n";
		return false;
	}

	if (abs(number_of_412 - number_of_m412) != NumberOfGenerations) // net number of NumberOfGenerations generations
	{
		if (info)
			(*Print.out) << "Info PS_CheckVectorlikeness: Net number of (4,1,2) is " << abs(number_of_412 - number_of_m412) << "!\n";
		return false;
	}

	// need at least one PS breaking Higgs
	if (((number_of_421 == 0) || (number_of_m421 == 0)) && ((number_of_412 == 0) || (number_of_m412 == 0)))
	{
		if (info)
			(*Print.out) << "Info PS_CheckVectorlikeness: PS breaking Higgs is missing!\n";
		return false;
	}

	return true;
}



/* ########################################################################################
######   PS_GetProtonHexality(...)                                                   ######
######                                                                               ######
######   Version: 14.04.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) PSVEVConfig      : a Pati-Salam vev-config                               ######
######   2) Roots10D         : 480 roots of E_8 x E_8 or SO(32) used to identify an  ######
######                         intermediate SO(12)                                   ######
######   3) ProtonHexalities : possible generators for Proton hexality               ######
######   output:                                                                     ######
######   return value        : proton hexality identified?                           ######
###########################################################################################
######   description:                                                                ######
######   Analyze the vev-config "PSVEVConfig" to identify possible proton hexality   ######
######   generators.                                                                 ######
######################################################################################## */
bool CAnalyseModel::PS_GetProtonHexality(const SConfig &PSVEVConfig, const vector<CVector> &Roots10D, vector<CVector> &ProtonHexalities)
{
	const unsigned pos_of_SU4  = PSVEVConfig.SymmetryGroup.observable_sector_GGs[0];
	const unsigned pos_of_SU2L = PSVEVConfig.SymmetryGroup.observable_sector_GGs[1];
	const unsigned pos_of_SU2R = PSVEVConfig.SymmetryGroup.observable_sector_GGs[2];

	const double prec = 0.0001;
	int i = 0, j = 0, k = 0;

	if (ProtonHexalities.size() != 0)
	{
		cout << "\n  Warning in bool CAnalyseModel::PS_GetProtonHexality(...): ProtonHexalities not empty - now cleared." << endl;
		ProtonHexalities.clear();
	}

	const vector<gaugeGroupFactor<double> > &factors = PSVEVConfig.SymmetryGroup.GaugeGroup.factor;
	const size_t s1 = factors.size();

	if ((pos_of_SU4 >= s1) || (pos_of_SU2L >= s1) || (pos_of_SU2R >= s1))
	{
		cout << "\n  Warning in bool CAnalyseModel::PS_GetProtonHexality(...) : Position of Pati-Salam gauge group out of range. Return false." << endl;
		return false;
	}

	if (Roots10D.size() != 480)
	{
		cout << "\n  Warning in bool CAnalyseModel::PS_GetProtonHexality(...) : \"Roots10D\" does not contain 480 roots. Return false." << endl;
		return false;
	}

	const vector<vector<double> > &SU4_Simpleroots = factors[pos_of_SU4].simpleroots;
	const vector<vector<double> > &SU2L_Simpleroots = factors[pos_of_SU2L].simpleroots;
	const vector<vector<double> > &SU2R_Simpleroots = factors[pos_of_SU2R].simpleroots;

	if ((SU4_Simpleroots.size() != 3) || (SU2L_Simpleroots.size() != 1) || (SU2R_Simpleroots.size() != 1))
	{
		cout << "\n  Warning in bool CAnalyseModel::PS_GetProtonHexality(...): Simple roots of Pati-Salam gauge group not correct. Return false." << endl;
		return false;
	}

	CVector SU4_Simpleroot1, SU4_Simpleroot2, SU4_Simpleroot3, SU2L_Simpleroot1, SU2R_Simpleroot1, SO12_Simpleroot;
	CVector tmp_root(16);

	SU4_Simpleroot1  = SU4_Simpleroots[0];
	SU4_Simpleroot2  = SU4_Simpleroots[1];
	SU4_Simpleroot3  = SU4_Simpleroots[2];
	SU2L_Simpleroot1 = SU2L_Simpleroots[0];
	SU2R_Simpleroot1 = SU2R_Simpleroots[0];

	CVector TestProtonHexality(16);
	vector<CVector> TestProtonHexalities;

	vector<vector<double> > Orig_X_Orthogonal;
	vector<vector<double> > Orig_all_SimpleRoots;
	vector<vector<double> > Orig_SO10_Orthogonal;

	vector<vector<double> > X_Orthogonal;
	vector<vector<double> > all_SimpleRoots;
	vector<vector<double> > SO10_Orthogonal;

	// begin: run through the gauge group factors and collect the simple roots
	for (j = 0; j < s1; ++j)
	{
		const vector<vector<double> > &tmp_SimpleRoots = factors[j].simpleroots;

		Orig_all_SimpleRoots.insert(Orig_all_SimpleRoots.end(), tmp_SimpleRoots.begin(), tmp_SimpleRoots.end());
		if ((j == pos_of_SU4) || (j == pos_of_SU2L) || (j == pos_of_SU2R))
			Orig_X_Orthogonal.insert(Orig_X_Orthogonal.end(), tmp_SimpleRoots.begin(), tmp_SimpleRoots.end());
		else
			Orig_SO10_Orthogonal.insert(Orig_SO10_Orthogonal.end(), tmp_SimpleRoots.begin(), tmp_SimpleRoots.end());
	}
	// end: run through the gauge group factors and collect the simple roots

	const size_t s3 = SO10_Orthogonal.size();

	// run through the roots of E8 x E8
	for (i = 0; i < 480; ++i)
	{
		SO12_Simpleroot = Roots10D[i];

		X_Orthogonal    = Orig_X_Orthogonal;
		all_SimpleRoots = Orig_all_SimpleRoots;
		SO10_Orthogonal = Orig_SO10_Orthogonal;

		// begin: find a candidate for the sixth SO(12) simple root
		bool is_SO12_Simpleroot = true;

		// begin: check whether the current root is orthogonal to Pati-Salam / SO(10)
		for (j = 0; (j < s3) && is_SO12_Simpleroot; ++j)
		{
			tmp_root = SO10_Orthogonal[j];

			if (fabs(SO12_Simpleroot * tmp_root) > prec)
				is_SO12_Simpleroot = false;
		}
		// end: check whether the current root is orthogonal to Pati-Salam / SO(10)

		// then check the Dynkin diagram, i.e. the scalar products with the simple roots of Pati-Salam
		if (is_SO12_Simpleroot)
		{
			is_SO12_Simpleroot = false;
			if ((fabs(SO12_Simpleroot * SU2L_Simpleroot1 + 1.0) < prec) && (fabs(SO12_Simpleroot * SU2R_Simpleroot1 + 1.0) < prec))
			{
				const double sp1 = SO12_Simpleroot * SU4_Simpleroot1;
				const double sp2 = SO12_Simpleroot * SU4_Simpleroot2;
				const double sp3 = SO12_Simpleroot * SU4_Simpleroot3;

				const double sp12 = SU4_Simpleroot1 * SU4_Simpleroot2;
				const double sp13 = SU4_Simpleroot1 * SU4_Simpleroot3;
				const double sp23 = SU4_Simpleroot2 * SU4_Simpleroot3;

				if (((fabs(sp1 + 1.0) < prec) && (fabs(sp2) < prec) && (fabs(sp3) < prec) && ((fabs(sp12) < prec) || (fabs(sp13) < prec))) || ((fabs(sp1) < prec) && (fabs(sp2 + 1.0) < prec) && (fabs(sp3) < prec) && ((fabs(sp12) < prec) || (fabs(sp23) < prec))) || ((fabs(sp1) < prec) && (fabs(sp2) < prec) && (fabs(sp3 + 1.0) < prec) && ((fabs(sp13) < prec) || (fabs(sp23) < prec))))
					is_SO12_Simpleroot = true;
			}
		}
		// end: find a candidate for the sixth SO(12) simple root

		// if the new root enhances the Dynkin diagram from Pati-Salam to SO(12)
		if (is_SO12_Simpleroot)
		{
			all_SimpleRoots.push_back(SO12_Simpleroot);

			// begin: find the U(1)'s that are orthogonal to SO(12)
			vector<vector<double> > tmp_U1s;
			Find_Basis_Of_Orthogonal_Space(all_SimpleRoots, UNSPECIFIED_LATTICE, 16, tmp_U1s);

			if (tmp_U1s.size() != PSVEVConfig.SymmetryGroup.GaugeGroup.u1directions.size() - 1)
			{
				cout << "\n  Warning in bool CAnalyseModel::PS_GetProtonHexality(...) : Excesive number of U(1)'s found! Return false." << endl;
				return false;
			}
			// end: find the U(1)'s that are orthogonal to SO(12)

			// begin: find the ProHex
			SO10_Orthogonal.insert(SO10_Orthogonal.end(), tmp_U1s.begin(), tmp_U1s.end());
			X_Orthogonal.insert(X_Orthogonal.end(), SO10_Orthogonal.begin(), SO10_Orthogonal.end());

			vector<vector<double> > tmp_X;
			Find_Basis_Of_Orthogonal_Space(X_Orthogonal, UNSPECIFIED_LATTICE, 16, tmp_X);

			if (tmp_X.size() != 1)
			{
				cout << "\n  Warning in bool CAnalyseModel::PS_GetProtonHexality(...) :" << tmp_X.size() << " X's were found, but only one is expected! Return false." << endl;
				return false;
			}
			// end: find the ProHex

			bool X_is_gauge_U1 = true;
			// begin: check whether the ProHex can be spanned by the gauge U(1)'s
			vector<vector<double> > U1directions = PSVEVConfig.SymmetryGroup.GaugeGroup.u1directions;
			const size_t gauge_U1s = U1directions.size();
			U1directions.push_back(tmp_X[0]);
			U1directions = findBasis<double>(U1directions);

			if (gauge_U1s != U1directions.size())
			{
				cout << "\n  Warning in bool CAnalyseModel::PS_GetProtonHexality(...): Proton Hexality is NOT gauge U1!!" << endl;
				X_is_gauge_U1 = false;
			}
			// end: check whether the ProHex can be spanned by the gauge U(1)'s

			if (X_is_gauge_U1)
			{
				TestProtonHexality = tmp_X[0];

				TestProtonHexalities.push_back(TestProtonHexality);
			}
		}
	}

	size_t s5 = 0;
	bool orthogonal = true;

	CVector U1DirectionA(16);
	U1DirectionA = PSVEVConfig.SymmetryGroup.GaugeGroup.u1directions[0];

	CVector Simpleroot(16);

	const bool AnomalousU1 = PSVEVConfig.SymmetryGroup.IsFirstU1Anomalous;

	const size_t s4 = TestProtonHexalities.size();
	for (i = 0; i < s4; ++i)
	{
		const CVector &Test = TestProtonHexalities[i];

		orthogonal = true;

		for (j = 0; orthogonal && (j < s1); ++j)
		{
			const vector<vector<double> > &Simpleroots = factors[j].simpleroots;
			s5 = Simpleroots.size();
			for (k = 0; orthogonal && (k < s5); ++k)
			{
				Simpleroot = Simpleroots[k];
				if (fabs(Test * Simpleroot) > prec)
					orthogonal = false;
			}
		}

		if (orthogonal && AnomalousU1 && (fabs(U1DirectionA * Test) > prec))
			orthogonal = false;

		if (orthogonal)
		{
			s5 = ProtonHexalities.size();
			bool known = false;
			for (j = 0; !known && j < s5; ++j)
			{
				CVector XTmp = ProtonHexalities[j];
				if((Test == XTmp) || (Test == XTmp * (-1)))
					known = true;
			}
			if (!known)
				ProtonHexalities.push_back(Test);
		}
	}

	//cout << ProtonHexalities.size() << " ProtonHexalities :" << endl;
	//for (i=0; i<ProtonHexalities.size(); ++i)
	//{
	//  ProtonHexalities[i].PrintRational(E8xE8, cout);
	//  cout << endl;
	//}

	if (ProtonHexalities.size() != 0)
		return true;
	else
		return false;
}



/* ########################################################################################
######   SU5_NetNumberOfGenerations(const SConfig &SU5VEVConfig) const               ######
######                                                                               ######
######   Version: 28.01.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) SU5VEVConfig : a SU(5) vev-config for which the net numbers of 10- and   ######
######                    5-plets are counted                                        ######
######   output:                                                                     ######
######   return value   : net number of SU(5) generations                            ######
###########################################################################################
######   description:                                                                ######
######   Counts the net number of SU(5) generations in the vev-config "SU5VEVConfig".######
######################################################################################## */
int CAnalyseModel::SU5_NetNumberOfGenerations(const SConfig &SU5VEVConfig) const
{
	int net_no_10 = 0;
	int net_no_5  = 0;
	int mult      = 0;

	unsigned i = 0;
	unsigned j = 0;

	const vector<CField> &Fields = SU5VEVConfig.Fields;
	const size_t f1 = Fields.size();
	if (f1 == 0)
	{
		cout << "\n  Warning in int CAnalyseModel::SU5_NetNumberOfGenerations(...) const: Spectrum is empty. Return false." << endl;
		return false;
	}
	if (SU5VEVConfig.SymmetryGroup.observable_sector_GGs.size() < 1)
	{
		cout << "\n  Warning in int CAnalyseModel::SU5_NetNumberOfGenerations(...) const: Cannot find the SU(5) group. Return false." << endl;
		return false;
	}

	const unsigned pos_of_SU5  = SU5VEVConfig.SymmetryGroup.observable_sector_GGs[0];

	const size_t number_of_factors = SU5VEVConfig.SymmetryGroup.GaugeGroup.factor.size();

	for (i = 0; i < f1; ++i)
	{
		const CField &Field = Fields[i];
		if (Field.Multiplet == LeftChiral)
		{
			const RepVector &Dimensions = Field.Dimensions;

			if (abs(Dimensions[pos_of_SU5].Dimension) == 10)
			{
				mult = 1;

				for (j = 0; j < number_of_factors; ++j)
				{
					if (j == pos_of_SU5)
						mult *= Dimensions[j].Dimension;
					else
						mult *= abs(Dimensions[j].Dimension);
				}
				net_no_10 += mult/10;
			}
			else
				if (abs(Dimensions[pos_of_SU5].Dimension) == 5)
				{
					mult = 1;

					for (j = 0; j < number_of_factors; ++j)
					{
						if (j == pos_of_SU5)
							mult *= Dimensions[j].Dimension;
						else
							mult *= abs(Dimensions[j].Dimension);
					}
					net_no_5 += mult/5;
				}
		}
	}
	if (((net_no_10 * net_no_5) > 0 ) || (abs(net_no_10) != abs(net_no_5)))
		return 0;

	return net_no_10;
}



/* ########################################################################################
######   SU5_GetFlippedU1(const SConfig &SU5VEVConfig, ...) const                    ######
######                                                                               ######
######   Version: 28.01.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) SU5VEVConfig : a SU(5) vev-config                                        ######
######   2) Roots10D     : 480 roots of E_8 x E_8 or SO(32) used to identify an      ######
######                     intermediate SO(10)                                       ######
######   3) FlippedU1s   : possible generators for U(1)_fl                           ######
######   output:                                                                     ######
######   return value    : U(1)_fl identified?                                       ######
###########################################################################################
######   description:                                                                ######
######   Analyze the vev-config "SU5VEVConfig" to identify possible U(1)_fl          ######
######   generators.                                                                 ######
######################################################################################## */
bool CAnalyseModel::SU5_GetFlippedU1(const SConfig &SU5VEVConfig, const vector<CVector> &Roots10D, vector<CVector> &FlippedU1s) const
{
	const double prec = 0.0001;
	int i = 0, j = 0, k = 0;

	if (FlippedU1s.size() != 0)
	{
		cout << "\n  Warning in bool CAnalyseModel::SU5_GetFlippedU1(...) const: FlippedU1s not empty - now cleared." << endl;
		FlippedU1s.clear();
	}

	const CGaugeGroup                       &GaugeGroup = SU5VEVConfig.SymmetryGroup.GaugeGroup;
	const vector<gaugeGroupFactor<double> > &factors    = GaugeGroup.factor;
	const size_t number_of_factors = factors.size();

	if (SU5VEVConfig.SymmetryGroup.observable_sector_GGs.size() < 1)
	{
		cout << "\n  Warning in bool CAnalyseModel::SU5_GetFlippedU1(...) const: Cannot find the SU(5) group. Return false." << endl;
		return false;
	}

	const unsigned pos_of_SU5  = SU5VEVConfig.SymmetryGroup.observable_sector_GGs[0];
	const string SU5 = "A4";
	if ((pos_of_SU5 >= number_of_factors) || (factors[pos_of_SU5].algebra != SU5))
	{
		cout << "\n  Warning in bool CAnalyseModel::SU5_GetFlippedU1(...) const: Position of SU(5) out of range. Return false." << endl;
		return false;
	}

	const size_t s1 = Roots10D.size();

	if (s1 != 480)
	{
		cout << "\n  Warning in bool CAnalyseModel::SU5_GetFlippedU1(...) const: \"Roots10D\" does not contain 480 roots. Return false." << endl;
		return false;
	}

	const vector<vector<double> > &SU5_Simpleroots = factors[pos_of_SU5].simpleroots;

	if (SU5_Simpleroots.size() != 4)
	{
		cout << "\n  Warning in bool CAnalyseModel::SU5_GetFlippedU1(...) const: Simple roots of SU(5) not correct. Return false." << endl;
		return false;
	}

	CVector SU5_Simpleroot1, SU5_Simpleroot2, SU5_Simpleroot3, SU5_Simpleroot4, SO10_Simpleroot5;
	CVector tmp_root(16);

	SU5_Simpleroot1 = SU5_Simpleroots[0];
	SU5_Simpleroot2 = SU5_Simpleroots[1];
	SU5_Simpleroot3 = SU5_Simpleroots[2];
	SU5_Simpleroot4 = SU5_Simpleroots[3];

	if (((SU5_Simpleroot1 * SU5_Simpleroot2) != -1) || ((SU5_Simpleroot2 * SU5_Simpleroot3) != -1) || ((SU5_Simpleroot3 * SU5_Simpleroot4) != -1) || ((SU5_Simpleroot1 * SU5_Simpleroot3) != 0) || ((SU5_Simpleroot1 * SU5_Simpleroot4) != 0) || ((SU5_Simpleroot2 * SU5_Simpleroot4) != 0))
	{
		cout << "\n  Warning in bool CAnalyseModel::SU5_GetFlippedU1(...) const: Dynkin diagram of SU(5) not correct. Return false." << endl;
		return false;
	}

	CVector TestFlippedU1(16);
	vector<CVector> TestFlippedU1s;

	vector<vector<double> > Orig_FlippedU1_Orthogonal;
	vector<vector<double> > Orig_all_SimpleRoots;
	vector<vector<double> > Orig_SO10_Orthogonal;

	vector<vector<double> > FlippedU1_Orthogonal;
	vector<vector<double> > all_SimpleRoots;
	vector<vector<double> > SO10_Orthogonal;

	for (j = 0; j < number_of_factors; ++j)
	{
		const vector<vector<double> > &tmp_SimpleRoots = factors[j].simpleroots;

		Orig_all_SimpleRoots.insert(Orig_all_SimpleRoots.end(), tmp_SimpleRoots.begin(), tmp_SimpleRoots.end());
		if (j == pos_of_SU5)
			Orig_FlippedU1_Orthogonal.insert(Orig_FlippedU1_Orthogonal.end(), tmp_SimpleRoots.begin(), tmp_SimpleRoots.end());
		else
			Orig_SO10_Orthogonal.insert(Orig_SO10_Orthogonal.end(), tmp_SimpleRoots.begin(), tmp_SimpleRoots.end());
	}

	const size_t s3 = Orig_SO10_Orthogonal.size();

	for (i = 0; i < s1; ++i)
	{
		FlippedU1_Orthogonal = Orig_FlippedU1_Orthogonal;
		all_SimpleRoots      = Orig_all_SimpleRoots;
		SO10_Orthogonal      = Orig_SO10_Orthogonal;

		// begin: find a candidate for the fifth SO(10) simple root
		bool is_SO10_Simpleroot5 = true;

		SO10_Simpleroot5 = Roots10D[i];

		// missing simple root of SO(10) must be orthogonal to the other gauge group factors
		for (j = 0; (j < s3) && is_SO10_Simpleroot5; ++j)
		{
			tmp_root = SO10_Orthogonal[j];

			if (fabs(SO10_Simpleroot5 * tmp_root) > prec)
				is_SO10_Simpleroot5 = false;
		}

		// missing simple root of SO(10) must complete the Dynkin diagram of SU(5) to the one of SO(10)
		if (is_SO10_Simpleroot5)
		{
			is_SO10_Simpleroot5 = false;
			if ((fabs(SO10_Simpleroot5 * SU5_Simpleroot1) < prec) && (fabs(SO10_Simpleroot5 * SU5_Simpleroot4) < prec))
			{
				const double SU5Root2xSO10Root5 = SO10_Simpleroot5 * SU5_Simpleroot2;
				const double SU5Root3xSO10Root5 = SO10_Simpleroot5 * SU5_Simpleroot3;

				if (((fabs(SU5Root2xSO10Root5 + 1.0) < prec) && (fabs(SU5Root3xSO10Root5) < prec)) || ((fabs(SU5Root2xSO10Root5) < prec) && (fabs(SU5Root3xSO10Root5 + 1.0) < prec)))
					is_SO10_Simpleroot5 = true;
			}
		}
		// end: find a candidate for the fifth SO(10) simple root

		if (is_SO10_Simpleroot5)
		{
			all_SimpleRoots.push_back(SO10_Simpleroot5);

			// begin: find the U(1)'s that are orthogonal to SO(10)
			vector<vector<double> > tmp_U1s;
			Find_Basis_Of_Orthogonal_Space(all_SimpleRoots, UNSPECIFIED_LATTICE, 16, tmp_U1s);

			if (tmp_U1s.size() != GaugeGroup.u1directions.size() - 1)
			{
				cout << "\n  Warning in bool CAnalyseModel::SU5_GetFlippedU1(...) const: Excesive number of U(1)s found. Return false." << endl;
				return false;
			}
			// end: find the U(1)'s that are orthogonal to SO(10)

			// begin: find the flipped U(1)
			SO10_Orthogonal.insert(SO10_Orthogonal.end(), tmp_U1s.begin(), tmp_U1s.end());
			FlippedU1_Orthogonal.insert(FlippedU1_Orthogonal.end(), SO10_Orthogonal.begin(), SO10_Orthogonal.end());

			vector<vector<double> > tmp_Fl;
			Find_Basis_Of_Orthogonal_Space(FlippedU1_Orthogonal, UNSPECIFIED_LATTICE, 16, tmp_Fl);

			if (tmp_Fl.size() != 1)
			{
				cout << "\n  Warning in bool CAnalyseModel::SU5_GetFlippedU1(...) const:" << tmp_Fl.size() << " flipped U(1)s were found, but only one is expected. Return false." << endl;
				return false;
			}
			// end: find the flipped U(1)

			bool flU1_is_gauge_U1 = true;
			// begin: check whether the flipped U(1) can be spanned by the gauge U(1)'s
			vector<vector<double> > U1directions = GaugeGroup.u1directions;
			const size_t gauge_U1s = U1directions.size();
			U1directions.push_back(tmp_Fl[0]);
			U1directions = findBasis<double>(U1directions);

			if (gauge_U1s != U1directions.size())
			{
				cout << "\n  Warning in bool CAnalyseModel::SU5_GetFlippedU1(...) const: Flipped U(1) is NOT gauge U1." << endl;
				flU1_is_gauge_U1 = false;
			}
			// end: check whether the flipped U(1) can be spanned by the gauge U(1)'s

			if (flU1_is_gauge_U1)
			{
				TestFlippedU1 = tmp_Fl[0];
				//TestFlippedU1.PrintRational(E8xE8, (*Print.out));
				//(*Print.out) << endl;
				//TestFlippedU1 = TestFlippedU1 * sqrt((5.0/6.0)/(TestFlippedU1.GetSqrTo(16)));

				TestFlippedU1s.push_back(TestFlippedU1);
			}
		}
	}

	size_t s5 = 0;
	bool orthogonal = true;

	CVector U1DirectionA(16);
	U1DirectionA = GaugeGroup.u1directions[0];

	CVector Simpleroot(16);

	const bool AnomalousU1 = SU5VEVConfig.SymmetryGroup.IsFirstU1Anomalous;

	const size_t s4 = TestFlippedU1s.size();
	for (i = 0; i < s4; ++i)
	{
		const CVector &Test = TestFlippedU1s[i];

		orthogonal = true;

		for (j = 0; orthogonal && (j < number_of_factors); ++j)
		{
			const vector<vector<double> > &Simpleroots = factors[j].simpleroots;
			s5 = Simpleroots.size();
			for (k = 0; orthogonal && (k < s5); ++k)
			{
				Simpleroot = Simpleroots[k];
				if (fabs(Test * Simpleroot) > prec)
				{
					cout << "\n  Warning in bool CAnalyseModel::SU5_GetFlippedU1(...) const: U(1)_fl direction is not orthogonal to the simple roots." << endl;
					orthogonal = false;
				}
			}
		}

		if (orthogonal && AnomalousU1 && (fabs(U1DirectionA * Test) > prec))
			orthogonal = false;

		if (orthogonal)
		{
			s5 = FlippedU1s.size();
			bool known = false;
			for (j = 0; !known && j < s5; ++j)
			{
				CVector YTmp = FlippedU1s[j];
				if((Test == YTmp) || (Test == YTmp * (-1)))
					known = true;
			}
			if (!known)
				FlippedU1s.push_back(Test);
		}
	}

	//(*Print.out) << FlippedU1s.size() << " FlippedU1s :" << endl;
	//for (i=0; i<FlippedU1s.size(); ++i)
	//{
	//  FlippedU1s[i].PrintRational(E8xE8, (*Print.out));
	//  (*Print.out) << endl;
	//}

	if (FlippedU1s.size() != 0)
		return true;
	else
		return false;
}



/* ########################################################################################
######   SU5_CheckVectorlikeness(const COrbifold &Orbifold, ...) const               ######
######                                                                               ######
######   Version: 19.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Orbifold         : the orbifold of the vev-config "SU5VEVConfig"         ######
######   2) SU5VEVConfig     : vev-config with SU(5) gauge group to analyze          ######
######   3) pos_of_flippedU1 : index of U(1)_fl in the list of all U(1) factors      ######
######   4) Print            : if "info" is true print the output to this CPrint     ######
######                         object                                                ######
######   5) info             : print info?                                           ######
######   output:                                                                     ######
######   return value        : spectrum is SU(5)xU(1)_fl three generations plus      ######
######                         vectorlike exotics?                                   ######
###########################################################################################
######   description:                                                                ######
######   Checks whether the vev-config "SU5VEVConfig" has a net number of three      ######
######   generations of SU(5)xU(1)_fl plus vectorlike exotics.                       ######
######################################################################################## */
bool CAnalyseModel::SU5_CheckVectorlikeness(const COrbifold &Orbifold, SConfig &SU5VEVConfig, unsigned pos_of_flippedU1, CPrint &Print, bool info) const
{
	// Set the precision
	const double prec = 0.0001;

	const double U1_Generation10 =  1.0;
	const double U1_Generationm5 = -3.0;
	const double U1_Generation1  =  5.0;
	const double U1_SMHiggs5     = -2.0;
	const double U1_SMHiggsm5    =  2.0;
	const double U1_SU5Higgs10   =  1.0;
	const double U1_SU5Higgsm10  = -1.0;

	int Rep10 = 10;
	int Repm5 = -5;
	int Rep1  =  1;
	int SMHiggs5    =   5;
	int SMHiggsm5   =  -5;
	int SU5Higgs10  =  10;
	int SU5Higgsm10 = -10;

	unsigned i = 0;
	unsigned j = 0;
	unsigned k = 0;

	const vector<CField> &Fields = SU5VEVConfig.Fields;
	const size_t f1 = Fields.size();
	if (f1 == 0)
	{
		cout << "\n  Warning in bool CAnalyseModel::SU5_CheckVectorlikeness(...) const: Spectrum is empty. Return false." << endl;
		return false;
	}

	const vector<gaugeGroupFactor<double> > &factors = SU5VEVConfig.SymmetryGroup.GaugeGroup.factor;
	const size_t number_of_factors = factors.size();
	const size_t number_of_U1s     = SU5VEVConfig.SymmetryGroup.GaugeGroup.u1directions.size();

	if (pos_of_flippedU1 >= number_of_U1s)
	{
		cout << "\n  Warning in bool CAnalyseModel::SU5_CheckVectorlikeness(...) const: U(1)_" << pos_of_flippedU1 << " does not exist. Return false." << endl;
		return false;
	}

	if (SU5VEVConfig.SymmetryGroup.IsFirstU1Anomalous && (pos_of_flippedU1 == 0))
	{
		cout << "\n  Warning in bool CAnalyseModel::SU5_CheckVectorlikeness(...) const: Flipped U(1) cannot be the anomalous U(1). Return false." << endl;
		return false;
	}

	if (SU5VEVConfig.SymmetryGroup.observable_sector_GGs.size() == 0)
	{
		cout << "\n  Warning in bool CAnalyseModel::SU5_CheckVectorlikeness(...) const: \"observable_sector_GGs\" is empty. Return false." << endl;
		return false;
	}

	const unsigned pos_of_SU5 = SU5VEVConfig.SymmetryGroup.observable_sector_GGs[0];

	if (pos_of_SU5 >= number_of_factors)
	{
		cout << "\n  Warning in bool CAnalyseModel::SU5_CheckVectorlikeness(...) const: Check the position of SU(5). Return false." << endl;
		return false;
	}
	{
		const string SU5 = "A4";

		if (factors[pos_of_SU5].algebra != SU5)
		{
			cout << "\n  Warning in bool CAnalyseModel::SU5_CheckVectorlikeness(...) const: Cannot find the SU(5). Return false." << endl;
			return false;
		}
	}

	// first, it contains all different SU(5) representations;
	// later, after pairing up, it contains the net numbers of all SU(5) representations
	//        (except for the neutrinos) plus possible chiral exotics
	vector<unsigned>   SU5_Mults;
	vector<SDimension> SU5_Reps;
	vector<double>     SU5_U1s;

	// contains the net number of all representations, charged under the SU(5) x U(1)_fl gauge group
	// and the net number of SU(5) singlets (1)_0
	vector<int>        NetNumber_SU5_Mults;
	vector<SDimension> NetNumber_SU5_Reps;
	vector<double>     NetNumber_SU5_U1s;

	size_t SU5_size = 0;
	size_t NN_size = 0;

	SDimension tmp_SU5_Rep; // contains SU(5) reps.
	tmp_SU5_Rep.Dimension = 1;
	tmp_SU5_Rep.AdditionalLabel = "";

	unsigned SU5_Mult = 0;
	double   U1             = 0.0;
	double   SU5_U1         = 0.0;
	double   NetNumber_SU5_U1 = 0.0;

	// begin: find all different SU(5) x U(1)_fl representations + singlets
	{
		unsigned tmp_mult     = 1;
		bool     SU5_Rep_known = false;

		for (i = 0; i < f1; ++i)
		{
			const CField &Field = Fields[i];
			if (Field.Multiplet == LeftChiral)
			{
				const RepVector &Dimensions = Field.Dimensions;
				const CVector   &U1Charges  = Field.U1Charges;

				tmp_SU5_Rep = Dimensions[pos_of_SU5];
				U1          = U1Charges[pos_of_flippedU1];

				tmp_mult = 1;
				for (k = 0; k < number_of_factors; ++k)
				{
					if (k != pos_of_SU5)
						tmp_mult *= abs(Dimensions[k].Dimension);
				}

				SU5_Rep_known = false;

				SU5_size = SU5_Reps.size();
				for (j = 0; !SU5_Rep_known && (j < SU5_size); ++j)
				{
					// if SU(5) x U(1)_fl rep is known, increase the counter
					if ((SU5_Reps[j].Dimension == tmp_SU5_Rep.Dimension) && (fabs(SU5_U1s[j] - U1) < prec))
					{
						SU5_Rep_known = true;
						SU5_Mults[j] += tmp_mult;
					}
				}
				// if SU(5) x U(1)_fl rep is not known, save it
				if (!SU5_Rep_known)
				{
					SU5_Mults.push_back(tmp_mult);
					SU5_Reps.push_back(tmp_SU5_Rep);
					SU5_U1s.push_back(U1);
				}
			}
		}
	}
	SU5_size = SU5_Mults.size();
	// end: find all different SU(5) x U(1)_fl representations + singlets


	// begin: print
	(*Print.out) << "SU(5) x U(1)_fl representations:\n";
	for (i = 0; i < SU5_size; ++i)
		(*Print.out) << SU5_Mults[i] << " (" << SU5_Reps[i].Dimension << ")_" << SU5_U1s[i] << "\n";
	(*Print.out) << endl;
	// end: print

	vector<unsigned>   PairedUp_SU5_Mults;
	vector<SDimension> PairedUp_SU5_Reps;
	vector<double>     PairedUp_SU5_U1s;
	// begin: compute the net number of representations
	{
		bool NetNumber_SU5_Rep_known = false;

		for (i = 0; i < SU5_size; ++i)
		{
			const SDimension SU5_Rep = SU5_Reps[i];
			SU5_Mult = SU5_Mults[i];
			SU5_U1   = SU5_U1s[i];

			// "tmp_SU5_Rep" contains the cc-partner of "SU5_Rep"
			tmp_SU5_Rep = SU5_Rep;

			if ((abs(tmp_SU5_Rep.Dimension) == 5) || (abs(tmp_SU5_Rep.Dimension) == 10))
				tmp_SU5_Rep.Dimension *= -1;

			NetNumber_SU5_Rep_known = false;

			NN_size = NetNumber_SU5_Mults.size();
			for (j = 0; !NetNumber_SU5_Rep_known && (j < NN_size); ++j)
			{
				const int NetNumber_SU5_Rep = NetNumber_SU5_Reps[j].Dimension;
				NetNumber_SU5_U1 = NetNumber_SU5_U1s[j];

				if ((NetNumber_SU5_Rep == SU5_Rep.Dimension) && (fabs(NetNumber_SU5_U1 - SU5_U1) < prec))
				{
					NetNumber_SU5_Rep_known = true;
					NetNumber_SU5_Mults[j] += SU5_Mult;
				}
				else
					if ((NetNumber_SU5_Rep == tmp_SU5_Rep.Dimension) && (fabs(NetNumber_SU5_U1 + SU5_U1) < prec))
					{
						// save the paired-up representations in order to find the Higgses later
						if (SU5_Mult < NetNumber_SU5_Mults[j])
							PairedUp_SU5_Mults.push_back(SU5_Mult);
						else
							PairedUp_SU5_Mults.push_back(NetNumber_SU5_Mults[j]);

						PairedUp_SU5_Reps.push_back(SU5_Rep);
						PairedUp_SU5_U1s.push_back(SU5_U1);

						NetNumber_SU5_Rep_known = true;
						NetNumber_SU5_Mults[j] -= SU5_Mult;
					}
			}
			if (!NetNumber_SU5_Rep_known)
			{
				NetNumber_SU5_Mults.push_back(SU5_Mult);
				NetNumber_SU5_Reps.push_back(SU5_Rep);
				NetNumber_SU5_U1s.push_back(SU5_U1);
			}
		}
	}
	NN_size = NetNumber_SU5_Mults.size();
	// end: compute the net number of representations

	// the list of paired-up representations should contain at least the SM Higgs and the GUT Higgs
	if (PairedUp_SU5_Mults.size() < 2)
	{
		if (info)
			(*Print.out) << "Info SU5_CheckVectorlikeness: Paired-up representations cannot give rise to Higgses!\n";
		return false;
	}

	// begin: make all multiplicities positive
	for (i = 0; i < NN_size; ++i)
	{
		if (NetNumber_SU5_Mults[i] < 0)
		{
			NetNumber_SU5_Mults[i] *= (-1);
			NetNumber_SU5_U1s[i]   *= (-1);
			const SDimension SU5_Rep = NetNumber_SU5_Reps[i];

			if ((abs(SU5_Rep.Dimension) == 5) || (abs(SU5_Rep.Dimension) == 10))
				NetNumber_SU5_Reps[i].Dimension *= (-1);
		}
	}
	// end: make all multiplicities positive

	(*Print.out) << "Net number of SU(5) x U(1)_fl representations\n";
	for (i = 0; i < NN_size; ++i)
		(*Print.out) << NetNumber_SU5_Mults[i] << " (" << NetNumber_SU5_Reps[i].Dimension << ")_" << NetNumber_SU5_U1s[i] << "\n";
	(*Print.out) << endl;

	// begin: remove representations with zero multiplicity
	SU5_Mults.clear();
	SU5_Reps.clear();
	SU5_U1s.clear();

	for (i = 0; i < NN_size; ++i)
	{
		if (NetNumber_SU5_Mults[i] != 0)
		{
			SU5_Mults.push_back(NetNumber_SU5_Mults[i]);
			SU5_Reps.push_back(NetNumber_SU5_Reps[i]);
			SU5_U1s.push_back(NetNumber_SU5_U1s[i]);
		}
	}
	// end: remove representations with zero multiplicity

	NN_size = SU5_Mults.size();

	// begin: print
	(*Print.out) << "Non zero number of SU(5) x U(1)_fl representations\n";
	for (i = 0; i < NN_size; ++i)
		(*Print.out) << SU5_Mults[i] << " (" << SU5_Reps[i].Dimension << ")_" << SU5_U1s[i] << "\n";
	(*Print.out) << endl;

	(*Print.out) << "Paired-up SU(5) x U(1)_fl representations\n";
	for (i = 0; i < PairedUp_SU5_Mults.size(); ++i)
		(*Print.out) << PairedUp_SU5_Mults[i] << " (" << PairedUp_SU5_Reps[i].Dimension << ")_" << PairedUp_SU5_U1s[i] << "\n";
	(*Print.out) << endl;
	// end: print

	bool   change_sign = false;
	double ScaleFactor = 1.0;
	// begin: find three \bar{5}-plets and three 10-plets with the same U(1)_fl charge
	{
		bool stop_10 = false;
		bool stop_m5 = false;

		for (i = 0; i < NN_size; ++i)
		{
			if (!stop_10 && (abs(SU5_Reps[i].Dimension) == 10))
			{
				if (SU5_Mults[i] != 3)
				{
					if (info)
						(*Print.out) << "Info SU5_CheckVectorlikeness: Either not all extra 10-plets are paired up or the three 10-plets have different U(1)_fl charges!\n";
					stop_10 = true;
				}
				else
				{
					if (SU5_Reps[i].Dimension < 0)
						change_sign = true;

					// the 10-plet of a generation should have U(1)_fl charge 1
					ScaleFactor = U1_Generation10/SU5_U1s[i];
				}
			}
			if (!stop_m5 && (abs(SU5_Reps[i].Dimension) == 5))
			{
				if (SU5_Mults[i] != 3)
				{
					if (info)
						(*Print.out) << "Info SU5_CheckVectorlikeness: Either not all extra 5-plets are paired up or the three -5-plets have different U(1)_fl charges!\n";
					stop_m5 = true;
				}
			}
		}
		if (stop_10 || stop_m5)
			return false;
	}
	// end: find three \bar{5}-plets and three 10-plets with the same U(1)_fl charge

	// only the following representations should have remained in the list: (10)_1 + (-5)_-3  and maybe (1)_0 and (1)_5
	if ((NN_size < 2) || (NN_size > 4))
	{
		if (info)
			(*Print.out) << "Info SU5_CheckVectorlikeness: Too many different representations for a flipped SU(5) model!\n";
		return false;
	}

	// begin: rescale U(1)_fl if neccessary
	if (fabs(ScaleFactor - 1.0) > prec)
	{
		CVector tmp;
		tmp = SU5VEVConfig.SymmetryGroup.GaugeGroup.u1directions[pos_of_flippedU1];
		tmp = tmp * ScaleFactor;

		Orbifold.Config_SetU1Direction(SU5VEVConfig, tmp, pos_of_flippedU1);

		if (info)
			(*Print.out) << "Info SU5_CheckVectorlikeness: Rescaling flipped U(1).\n";

		return this->SU5_CheckVectorlikeness(Orbifold, SU5VEVConfig, pos_of_flippedU1, Print, info);
	}
	// end: rescale U(1)_fl if neccessary

	if (change_sign)
	{
		Rep10 = -10;
		Repm5 =  5;
		SMHiggs5    =  -5;
		SMHiggsm5   =   5;
		SU5Higgs10  = -10;
		SU5Higgsm10 =  10;
	}

	// begin: search for the Higgse
	bool SM_Higgs = false;
	for (i = 0; !SM_Higgs && (i < PairedUp_SU5_Mults.size()); ++i)
	{
		if (((PairedUp_SU5_Reps[i].Dimension == SMHiggs5) && (PairedUp_SU5_U1s[i] == U1_SMHiggs5)) || ((PairedUp_SU5_Reps[i].Dimension == SMHiggsm5) && (PairedUp_SU5_U1s[i] == U1_SMHiggsm5)))
			SM_Higgs = true;
	}
	if (!SM_Higgs)
	{
		if (info)
			(*Print.out) << "Info SU5_CheckVectorlikeness: SM Higgs (5)_-2 + (-5)_2 not found!\n";
		return false;
	}

	bool SU5_Higgs = false;
	for (i = 0; !SU5_Higgs && (i < PairedUp_SU5_Mults.size()); ++i)
	{
		if (((PairedUp_SU5_Reps[i].Dimension == SU5Higgs10) && (PairedUp_SU5_U1s[i] == U1_SU5Higgs10)) || ((PairedUp_SU5_Reps[i].Dimension == SU5Higgsm10) && (PairedUp_SU5_U1s[i] == U1_SU5Higgsm10)))
			SU5_Higgs = true;
	}
	if (!SU5_Higgs)
	{
		if (info)
			(*Print.out) << "Info SU5_CheckVectorlikeness: SU5 Higgs (10)_1 + (-10)_-1 not found!\n";
		return false;
	}
	// end: search for the Higgse

	// begin: search for 3 x [ (10)_1 + (-5)_-3 ] and maybe 3 x (1)_5
	bool Rep10_found = false;
	bool Repm5_found = false;
	bool Rep1_found  = false;
	for (i = 0; i < NN_size; ++i)
	{
		if (SU5_Mults[i] == 3)
		{
			if ((SU5_Reps[i].Dimension == Rep10) && (SU5_U1s[i] == U1_Generation10))
				Rep10_found = true;
			if ((SU5_Reps[i].Dimension == Repm5) && (SU5_U1s[i] == U1_Generationm5))
				Repm5_found = true;
			if ((SU5_Reps[i].Dimension == Rep1) && (SU5_U1s[i] == U1_Generation1))
				Rep1_found = true;
		}
	}
	// end: search for 3 x [ (10)_1 + (-5)_-3 ] and maybe 3 x (1)_5

	if (!Rep10_found)
	{
		if (info)
			(*Print.out) << "Info SU5_CheckVectorlikeness: 3 (10)_1 not found.\n";
		return false;
	}
	if (!Repm5_found)
	{
		if (info)
			(*Print.out) << "Info SU5_CheckVectorlikeness: 3 (-5)_-3 not found.\n";
		return false;
	}
	if (!Repm5_found)
	{
		if (info)
			cout << "\n  Warning in SU5_CheckVectorlikeness: (1)_5 not found, but needed from anomaly reasons!\n";
	}

	if (info)
		(*Print.out) << "Info SU5_CheckVectorlikeness: Spectrum has 3 generations plus vectorlike exotics!\n";

	return true;
}



/* ########################################################################################
######   GetAllBlowUpModePerFixedPoint(const COrbifold &Orbifold, ...) const         ######
######                                                                               ######
######   Version: 19.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Orbifold     : the orbifold of the vev-config "VEVConfig"                ######
######   2) VEVConfig    : contains the massless fields                              ######
######   3) FieldIndices : list of fields that may be blowup modes                   ######
######                                                                               ######
######   output:                                                                     ######
######   4) BlowupModes : list of blowup modes per fixed point                       ######
######   return value   : is there at least one blowup mode per                      ######
######                    fixed point?                                               ######
###########################################################################################
######   description:                                                                ######
######   For each fixed point of "Orbifold" collect all field indices (from the set  ######
######   "FieldIndices" in the vev-config "VEVConfig") of fields localized there.    ######
######################################################################################## */
bool CAnalyseModel::GetAllBlowUpModePerFixedPoint(const COrbifold &Orbifold, const SConfig &VEVConfig, const vector<unsigned> &FieldIndices, vector<vector<unsigned> > &BlowupModes) const
{
	if (VEVConfig.InvariantSupercharges.size() != 1)
	{
		cout << "\n  Warning: GetAllBlowUpModePerFixedPoint(...) only works for N=1 SUSY." << endl;
		return false;
	}

	bool OneBlowupModePerFixedPoint = true;

	if (BlowupModes.size() != 0)
	{
		cout << "\n  Warning in bool CAnalyseModel::GetAllBlowUpModePerFixedPoint(...) const: \"BlowupModes\" cleared." << endl;
		BlowupModes.clear();
	}

	vector<unsigned> AllLocalModes;
	vector<unsigned> LocalBlowUpModes;
	size_t s2 = 0;
	size_t s3 = 0;
	unsigned i = 0;
	unsigned j = 0;
	unsigned k = 0;

	string label1 = "";
	string label2 = "";

	const vector<CField> &Fields = VEVConfig.Fields;

	// run through the sectors
	const size_t s1 = Orbifold.GetNumberOfSectors();
	for(i = 0; i < s1; ++i)
	{
		const CSector &Sector = Orbifold.GetSector(i);

		// search in the twisted sectors only
		if (Sector.SectorHasLeftchiralRightmover() && ((Sector.Get_m() != 0) || (Sector.Get_n() != 0)))			//hacking here!!! not refined for ZMxZNxZK
		{
			s2 = Sector.GetNumberOfFixedBranes();

			// run through the fixed points of the current twisted sector
			for(j = 0; j < s2; ++j)
			{
				AllLocalModes.clear();
				Sector.GetFixedBrane(j).GetFieldIndices(Fields, LeftChiral, AllLocalModes);
				s3 = AllLocalModes.size();

				LocalBlowUpModes.clear();
				for (k = 0; k < s3; ++k)
				{
					if (find(FieldIndices.begin(), FieldIndices.end(), AllLocalModes[k]) != FieldIndices.end())
						LocalBlowUpModes.push_back(AllLocalModes[k]);
				}

				BlowupModes.push_back(LocalBlowUpModes);
				if (LocalBlowUpModes.size() == 0)
					OneBlowupModePerFixedPoint = false;
			}
		}
	}
	return OneBlowupModePerFixedPoint;
}



/* ########################################################################################
######   GetOneBlowUpModePerFixedPoint(const COrbifold &Orbifold, ...) const         ######
######                                                                               ######
######   Version: 01.06.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Orbifold     : the orbifold of the vev-config "VEVConfig"                ######
######   2) VEVConfig    : contains the massless fields                              ######
######   3) FieldIndices : list of indices of fields that may be blowup modes        ######
######   output:                                                                     ######
######   4) BlowupModes  : list of blowup modes per fixed point                      ######
######   return value    : is there at least one blowup mode per                     ######
######                     fixed point?                                              ######
###########################################################################################
######   description:                                                                ######
######   For each fixed point of "Orbifold" randomly choose one field index (from    ######
######   the set "FieldIndices" in the vev-config "VEVConfig") of a field localized  ######
######   there.                                                                      ######
######################################################################################## */
bool CAnalyseModel::GetOneBlowUpModePerFixedPoint(const COrbifold &Orbifold, const SConfig &VEVConfig, const vector<unsigned> &FieldIndices, vector<unsigned> &BlowupModes) const
{
	if (VEVConfig.InvariantSupercharges.size() != 1)
	{
		cout << "\n  Warning: GetOneBlowUpModePerFixedPoint(...) only works for N=1 SUSY." << endl;
		return false;
	}

	if (BlowupModes.size() != 0)
	{
		cout << "\n  Warning in bool CAnalyseModel::GetOneBlowUpModePerFixedPoint(...) const: FieldIndices cleared!" << endl;
		BlowupModes.clear();
	}

	srand ( time(NULL) );
	const double tmp1 = (double)RAND_MAX + 0.00001;
	unsigned i = 0;
	unsigned j = 0;

	vector<vector<unsigned> > tmp_BlowupModes;
	const bool OneBlowupModePerFixedPoint = this->GetAllBlowUpModePerFixedPoint(Orbifold, VEVConfig, FieldIndices, tmp_BlowupModes);

	const size_t s1 = tmp_BlowupModes.size();
	for (i = 0; i < s1; ++i)
	{
		const vector<unsigned> &LocalBlowUpModes = tmp_BlowupModes[i];
		if (LocalBlowUpModes.size() != 0)
		{
			j = (unsigned)(floor((double)LocalBlowUpModes.size() * rand()/tmp1));
			if (j >= LocalBlowUpModes.size())
			{
				cout << "\n  Warning in bool CAnalyseModel::GetOneBlowUpModePerFixedPoint(...): random number outside the allowed region!" << endl;
				j = 0;
			}
			BlowupModes.push_back(LocalBlowUpModes[j]);
		}
	}

	return OneBlowupModePerFixedPoint;
}


/* ########################################################################################
######   Labels_Create(istream &in, SConfig &VEVConfig, ...) const                   ######
######                                                                               ######
######   Version: 16.11.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) in        : read field labels from here                                  ######
######   2) VEVConfig : assign new labels to the fields in "VEVConfig"               ######
######   3) Print     : if "info" is true print the output to this CPrint object     ######
######   4) Multiplet : Assign labels to fields of this SUSY type (e.g. LeftChiral)  ######
######   5) info      : print info?                                                  ######
######   output:                                                                     ######
######   return value : finished succesfully?                                        ######
###########################################################################################
######   description:                                                                ######
######   Read new labels for the fields of SUSY type "Multiplet" in the vev-config   ######
######   "VEVConfig" from the input "in".                                            ######
######################################################################################## */
bool CAnalyseModel::Labels_Create(istream &in, SConfig &VEVConfig, CPrint &Print, const vector<SUSYMultiplet>  &Multiplet, bool info) const
{
	if(!in.good())
	{
		(*Print.out) << "Warning! Could not find the input for creating the labels." << endl;
		return false;
	}

	const SSymmetryGroup &SymmetryGroup = VEVConfig.SymmetryGroup;

	vector<CField> NewFields = VEVConfig.Fields;
	const size_t f1 = NewFields.size();
	if (f1 == 0)
		return false;

	unsigned i = 0;
	unsigned j = 0;

	vector<unsigned> tmp_FieldIndices;

	vector<unsigned>          Spec_Multiplicities;
	vector<RepVector>         Spec_Dimensions;
	vector<CVector>           Spec_U1Charges;
	vector<vector<unsigned> > Spec_FieldIndices;

	size_t s1 = 0;
	size_t s2 = 0;
	bool field_not_known = true;

          vector<SUSYMultiplet> Multiplets(2);						
	      Multiplets[0]=Scalar;									
	      Multiplets[1]=LeftFermi;

	for (i = 0; i < f1; ++i)
	{

      for (int k=0; k<Multiplet.size(); k++)
	  { 
		const CField &Field = NewFields[i];


       if ( Field.Multiplet == Multiplet[k] )     

		{
			field_not_known = true;
			s1 = Spec_Multiplicities.size();
			for (j = 0; field_not_known && (j < s1); ++j)
			{
				if (AreRepVectorsEqual(SymmetryGroup, Spec_Dimensions[j], Field.Dimensions) && AreU1ChargesEqual(SymmetryGroup, Spec_U1Charges[j], Field.U1Charges))
				{
					++Spec_Multiplicities[j];
					Spec_FieldIndices[j].push_back(i);

					field_not_known = false;
				}
			}
			if (field_not_known)
			{
				Spec_Multiplicities.push_back(1);
				Spec_Dimensions.push_back(Field.Dimensions);
				Spec_U1Charges.push_back(Field.U1Charges);

				tmp_FieldIndices.clear();
				tmp_FieldIndices.push_back(i);
				Spec_FieldIndices.push_back(tmp_FieldIndices);
			}
		}
      } 
         
	}

	s1 = Spec_Multiplicities.size();

	if (info)
	{
		(*Print.out) << "  Massless spectrum:\n\n";
		for (i = 0; i < s1; ++i)
		{
			(*Print.out) << setw(3) << Spec_Multiplicities[i] << " ";
			Print.PrintRep(Spec_Dimensions[i], SymmetryGroup);
			(*Print.out) << "  ";
			Print.PrintU1Charges(Spec_U1Charges[i], SymmetryGroup);
			(*Print.out) << "\n";
		}
		(*Print.out) << "\n  Please create the labels:" << endl;
	}

	// begin: read the input
	string tmp_string = "";
	const string AllowedCharacters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
	vector<string> corresponding_Labels;
	size_t NewNumberOfLabels = 0;

	for (i = 0; i < s1; ++i)
	{
		if (info)
		{
			(*Print.out) << setw(3) << Spec_Multiplicities[i] << " ";
			Print.PrintRep(Spec_Dimensions[i], SymmetryGroup);
			(*Print.out) << "  ";
			Print.PrintU1Charges(Spec_U1Charges[i], SymmetryGroup);
			(*Print.out) << " = ";
		}

		// read the label
		if (!getline(in, tmp_string))
			return false;

		// the new label shall not contain any of the forbidden characters
		if ((tmp_string.size() == 0) || (tmp_string.find_first_not_of(AllowedCharacters) != string::npos))
		{
			(*Print.out) << "  " << Print.cbegin << "Label Error: The label is only allowed to contain characters and numbers." << Print.cend << endl;
			return false;
		}

		// the new label shall not be in use already for a different representation
		if (find(corresponding_Labels.begin(), corresponding_Labels.end(), tmp_string) != corresponding_Labels.end())
			return false;

		corresponding_Labels.push_back(tmp_string);

		const vector<unsigned> FieldIndices = Spec_FieldIndices[i];
		s2 = FieldIndices.size();
		for (j = 0; j < s2; ++j)
		{
			CField &Field = NewFields[FieldIndices[j]];
			Field.Labels.push_back(tmp_string);
			Field.Numbers.push_back(j+1);

			NewNumberOfLabels = Field.Labels.size();
		}
	}
	// end: read the input

	unsigned counter = 1;
	for (j = 0; j < f1; ++j)
	{
		CField &Field = NewFields[j];
		if (Field.Labels.size() != NewNumberOfLabels)
		{
			Field.Labels.push_back("X");
			Field.Numbers.push_back(counter);
			++counter;
		}
	}
	VEVConfig.Fields = NewFields;

	return true;
}



/* ########################################################################################
######   Labels_Load(string ifilename, SConfig &VEVConfig) const                     ######
######                                                                               ######
######   Version: 17.05.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) ifilename : filename of file to load                                     ######
######   2) VEVConfig : the loaded labels are storerd in this vev-config             ######
######   output:                                                                     ######
######   return value : labels loaded succesfully from file?                         ######
###########################################################################################
######   description:                                                                ######
######   Load field labels for the fields of "VEVConfig" from file named "ifilename".######
######################################################################################## */
bool CAnalyseModel::Labels_Load(string ifilename, SConfig &VEVConfig) const
{
	vector<CField> &Fields = VEVConfig.Fields;

	const size_t f1 = Fields.size();

	if (f1 == 0)
		return false;

	ifstream in;
	in.open(ifilename.data());
	if(!in.is_open() || !in.good())
		return false;

	string::size_type loc1 = 0;

	unsigned index       = 0;
	string   currentline = "";
	string   tmp_string  = "";

	vector<bool>     IndicesFound(f1, false);
	vector<string>   NewLabels(f1, "");
	vector<unsigned> NewNumbers(f1, 1);

	const string ForbiddenCharacters = "^+- ()_\\";

	unsigned counter = 0;
	while (getline(in, currentline))
	{
		// analyse only non-empty lines of the file
		if (currentline.find_first_not_of(" ") != string::npos)
		{
			++counter;

			// begin: find the field index
			loc1 = currentline.find(" ", 0);
			if (loc1 == string::npos)
				return false;

			tmp_string = currentline.substr(0, loc1);
			if (tmp_string.find_first_not_of("0123456789") != string::npos)
				return false;

			index = (unsigned)atoi(tmp_string.c_str());
			if (index >= f1)
				return false;
			// end: find the field index

			// begin: find the label
			currentline = currentline.substr(loc1 + 1, string::npos);
			loc1 = currentline.find("_", 0);
			if (loc1 == string::npos)
				return false;

			NewLabels[index] = currentline.substr(0, loc1);
			// end: find the label

			// begin: find the number
			tmp_string = currentline.substr(loc1 + 1, string::npos);
			if (tmp_string.find_first_not_of("0123456789") != string::npos)
				return false;

			NewNumbers[index]   = (unsigned)atoi(tmp_string.c_str());
			IndicesFound[index] = true;

			if ((NewNumbers[index] == 0) || (NewLabels[index].find_first_of(ForbiddenCharacters) != string::npos))
				return false;
			// end: find the number
		}
	}
	in.close();

	// check that there is only one new label per field
	if (counter != f1)
		return false;

	// check if all fields have a new label
	if (find(IndicesFound.begin(), IndicesFound.end(), false) != IndicesFound.end())
		return false;

	// only after all labels are checked and ok, they are stored now
	for (unsigned i = 0; i < f1; ++i)
	{
		Fields[i].Labels.push_back(NewLabels[i]);
		Fields[i].Numbers.push_back(NewNumbers[i]);
	}
	++VEVConfig.use_Labels;

	return true;
}




/* ########################################################################################
######   SetBmLGenerator(const COrbifold &Orbifold,SConfig &VEVConfig, ...) const    ######
######                                                                               ######
######   Version: 16.09.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Orbifold              : the orbifold of the vev-config "VEVConfig"       ######
######   2) VEVConfig             : contains the massless fields                     ######
######   3) BmLGenerator          : 16-dim vector corresponding to the B-L generator ######
######   4) NeedsToBeNonAnomalous : U(1)_B-L generator shall be orthogonal to the    ######
######                              anomalous one                                    ######
######   output:                                                                     ######
######   return value             : finished succesfully?                            ######
###########################################################################################
######   description:                                                                ######
######   Computes the B-L charges with respect to the generator "BmLGenerator" for   ######
######   the fields of "VEVConfig".                                                  ######
######################################################################################## */
bool CAnalyseModel::SetBmLGenerator(const COrbifold &Orbifold, SConfig &VEVConfig, const CVector &BmLGenerator, bool NeedsToBeNonAnomalous) const
{
	// Set the precision
	const double prec = 0.0001;

	if (BmLGenerator.GetSize() != 16)
	{
		cout << "\n  Warning in bool CAnalyseModel::SetBmLGenerator(...) : the B-L generator must be a vector with 16 entries." << endl;
		return false;
	}

	const CGaugeGroup &GaugeGroup = VEVConfig.SymmetryGroup.GaugeGroup;

	const unsigned number_of_U1 = GaugeGroup.u1directions.size();

	if (number_of_U1 == 0)
	{
		cout << "\n  Warning in bool CAnalyseModel::SetBmLGenerator(...) : Cannot set B-L generator, because there are no U(1) factors." << endl;
		return false;
	}

	unsigned i = 0;
	unsigned j = 0;
	unsigned k = 0;

	bool orthogonal = true;

	if (NeedsToBeNonAnomalous && VEVConfig.SymmetryGroup.IsFirstU1Anomalous)
	{
		if (number_of_U1 == 1)
		{
			cout << "\n  Warning in bool CAnalyseModel::SetBmLGenerator(...) : Cannot set B-L generator, because the only U(1) is anomalous." << endl;
			return false;
		}

		const vector<double> &U1Direction_A = GaugeGroup.u1directions[0];
		CVector U1DirectionA(16);
		U1DirectionA = U1Direction_A;

		// begin: check if B-L direction is orthogonal to the simple roots and to the anomalous U1 direction
		orthogonal = true;
		if (fabs(U1DirectionA * BmLGenerator) > prec)
			orthogonal = false;

		if (!orthogonal && NeedsToBeNonAnomalous)
		{
			cout << "\n  Warning in bool CAnalyseModel::SetBmLGenerator(...) : B-L generator is not orthogonal to the one of the anomalous U(1)." << endl;
			return false;
		}
	}

	long double sp = 0;
	const size_t s1 = GaugeGroup.factor.size();
	size_t s2 = 0;

	// run through all simple roots
	for (i = 0; orthogonal && (i < s1); ++i)
	{
		const vector<vector<double> > &ggf_SimpleRoots = GaugeGroup.factor[i].simpleroots;
		s2 = ggf_SimpleRoots.size();

		// check if the B-L generator is orthogonal to the simple roots
		for (j = 0; orthogonal && (j < s2); ++j)
		{
			const vector<double> &SimpleRoot = ggf_SimpleRoots[j];

			sp = 0;
			for (k = 0; k < 16; ++k)
				sp += SimpleRoot[k] * BmLGenerator[k];

			if (fabs(sp) > prec)
				orthogonal = false;
		}
	}

	if (!orthogonal)
	{
		cout << "\n  Warning in bool CAnalyseModel::SetBmLGenerator(...) : B-L generator is not orthogonal to the simple roots." << endl;
		return false;
	}
	// end: check if B-L direction is orthogonal to the simple roots and to the anomalous U1 direction

	const vector<CSector> &Sectors = Orbifold.GetSectors();

	// begin: compute the B-L charge for all fields
	const size_t f1 = VEVConfig.Fields.size();
	for (i = 0; i < f1; ++i)
	{
		CField &Field = VEVConfig.Fields[i];
		const CVector Weight = Field.GetLMWeight(0, Sectors);

		Field.BmLCharge = Weight * BmLGenerator;
		if (fabs(Field.BmLCharge) < prec)
			Field.BmLCharge = 0.0;
	}
	// end: compute the B-L charge for all fields

	VEVConfig.SymmetryGroup.BmLGenerator = BmLGenerator;

	return true;
}



/* ########################################################################################
######   SetVev(SConfig &VEVConfig, const vector<unsigned> &FieldIndices, ...)       ######
######                                                                               ######
######   Version: 25.02.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) VEVConfig            : contains the massless fields                      ######
######   2) FieldIndices         : indices of fields whose vev shall be changed      ######
######   3) FieldComponents      : for each field of "FieldIndices" the component    ######
######                             that shall obtain a vev                           ######
######   4) VEV                  : the value of the vev                              ######
######   5) random               : random=true for random vevs                       ######
######   6) U1ChargesOfVevFields : the U(1) charges for all fields from              ######
######                             "FieldIndices" with vev                           ######
######   output:                                                                     ######
######   return value            : finished succesfully?                             ######
###########################################################################################
######   description:                                                                ######
######   Assign a vev (either "VEV" or random value) to the fields with indices      ######
######   "FieldIndices" from the vev-config "VEVConfig".                             ######
######################################################################################## */
bool CAnalyseModel::SetVev(SConfig &VEVConfig, const vector<unsigned> &FieldIndices, const vector<unsigned> &FieldComponents, double VEV, bool random, vector<vector<double> > *U1ChargesOfVevFields)
{
	bool vev_activ = true;

	if (random)
		srand ( time(NULL) );
	else
	{
		if (fabs(VEV) < 0.0001)
			vev_activ = false;
	}

	bool SaveU1s = false;
	if (U1ChargesOfVevFields != NULL)
		SaveU1s = true;

	vector<CField> &Fields = VEVConfig.Fields;
	const size_t f1 = Fields.size();

	const size_t s0 = FieldIndices.size();
	if (s0 != FieldComponents.size())
	{
		cout << "\n  Warning in bool CAnalyseModel::SetVev(...): Component indices of fields are ill-defined. Return false." << endl;
		return false;
	}
	for(unsigned i = 0; i < s0; ++i)
	{
		if (FieldIndices[i] >= f1)
		{
			cout << "\n  Warning in bool CAnalyseModel::SetVev(...): Field index is out of range. Return false." << endl;
			return false;
		}
		CField &tmp_Field = Fields[FieldIndices[i]];

		if (vev_activ)
		{
			if (SaveU1s)
				U1ChargesOfVevFields->push_back(tmp_Field.U1Charges);

			if (random)
				tmp_Field.VEVs[FieldComponents[i]] = 1.0 - 1.576453 * ((double)(rand() + 1.0))/((double)(RAND_MAX + 1.0));
			else
				tmp_Field.VEVs[FieldComponents[i]] = VEV;
		}
		else
			tmp_Field.VEVs.SetToZero();
	}
	return true;
}



/* ########################################################################################
######   FindUnbrokenSymmetries(const SConfig &StandardVEVConfig, ...)               ######
######                                                                               ######
######   Version: 06.07.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Orbifold          : the orbifold of the vev-config "OriginalVEVConfig"   ######
######   2) OriginalVEVConfig : the SConfig object "OriginalVEVConfig" must contain  ######
######                          the gauge group at the orbifold point                ######
######   3) Weights           : the weights of the fields that have a vev            ######
######   4) VEVConfig         : the SConfig object "VEVConfig" contains the gauge    ######
######                          group in this vev-configuration                      ######
######   5) Print             : a CPrint object for the output                       ######
######   6) info              : print info?                                          ######
######                                                                               ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Identifies the unbroken gauge symmetries in the vev-config                  ######
######   "OriginalVEVConfig" and stores the result in "VEVConfig".                   ######
######################################################################################## */
bool CAnalyseModel::FindUnbrokenSymmetries(const COrbifold &Orbifold, const SConfig &OriginalVEVConfig, const vector<CVector> &Weights, SConfig &VEVConfig, CPrint &Print, bool info)
{
	vector<vector<double> > UnbrokenSimpleRoots;
	vector<vector<double> > UnbrokenU1Directions;

	vector<vector<double> > BrokenSimpleRoots;
	vector<string>          BrokenGaugeGroups;

	const CGaugeGroup &OriginalGaugeGroup = OriginalVEVConfig.SymmetryGroup.GaugeGroup;

	const size_t number_of_factors = OriginalGaugeGroup.factor.size();
	const size_t number_of_U1s     = OriginalGaugeGroup.u1directions.size();

	unsigned i = 0;
	unsigned j = 0;
	unsigned k = 0;
	unsigned l = 0;

	/*unsigned total_rank = number_of_U1s;
  for (i = 0; i < number_of_factors; ++i)
    total_rank += OriginalGaugeGroup.factor[i].rank;

  if (total_rank != 16)
    return false;*/

	VEVConfig = OriginalVEVConfig;

	const size_t dim = 16;

	const double prec = 0.0001;

	const size_t s1 = Weights.size();  // number of fields with vev
	bool sym_unbroken = true;
	double sp = 0.0;

	const vector<double>    WeightLine(dim, 0.0);
	vector<vector<double> > WeightMatrix(s1, WeightLine);

	// begin: collect all trivially unbroken simple roots
	{
		for (i = 0; i < number_of_factors; ++i)
		{
			const vector<vector<double> > &ggf_SimpleRoots = OriginalGaugeGroup.factor[i].simpleroots;

			const size_t q1 = ggf_SimpleRoots.size();
			for (j = 0; j < q1; ++j)
			{
				const vector<double> &double_SimpleRoot = ggf_SimpleRoots[j];

				sym_unbroken = true;
				for (k = 0; sym_unbroken && (k < s1); ++k)
				{
					sp = 0.0;
					const vector<double> &Weight = Weights[k];
					for (l = 0; l < dim; ++l)
						sp += (Weight[l] * double_SimpleRoot[l]);

					if (fabs(sp) >= prec)
						sym_unbroken = false;
				}
				if (sym_unbroken)
					UnbrokenSimpleRoots.push_back(double_SimpleRoot);
				else
				{
					VEVConfig.SymmetryGroup.GGs_AdditionalLabels[i] = "broken";
					BrokenSimpleRoots.push_back(double_SimpleRoot);
					BrokenGaugeGroups.push_back(OriginalGaugeGroup.factor[i].algebra);
				}
			}
		}
	}
	// end: collect all trivially unbroken simple roots

	bool anom_U1_unbroken = false;
	vector<string> new_U1s_AdditionalLabels;
	vector<unsigned> Select_these_U1s;

	// begin: collect all trivially unbroken U1 Directions
	{
		for (i = 0; i < number_of_U1s; ++i)
		{
			const vector<double> &double_U1Direction = OriginalGaugeGroup.u1directions[i];

			sym_unbroken = true;

			for (k = 0; sym_unbroken && (k < s1); ++k)
			{
				sp = 0.0;
				const vector<double> &Weight = Weights[k];
				for (l = 0; l < dim; ++l)
					sp += (Weight[l] * double_U1Direction[l]);

				if (fabs(sp) >= prec)
					sym_unbroken = false;
			}
			if (sym_unbroken)
			{
				if (OriginalVEVConfig.SymmetryGroup.IsFirstU1Anomalous && (i == 0))
					anom_U1_unbroken = true;

				UnbrokenU1Directions.push_back(double_U1Direction);

				if ((find(OriginalVEVConfig.SymmetryGroup.observable_sector_U1s.begin(), OriginalVEVConfig.SymmetryGroup.observable_sector_U1s.end(), i) != OriginalVEVConfig.SymmetryGroup.observable_sector_U1s.end())
						|| (OriginalVEVConfig.SymmetryGroup.U1s_AdditionalLabels[i] != ""))
					Select_these_U1s.push_back(new_U1s_AdditionalLabels.size());

				new_U1s_AdditionalLabels.push_back(OriginalVEVConfig.SymmetryGroup.U1s_AdditionalLabels[i]);
			}
		}
	}
	// end: collect all trivially unbroken U1 Directions

	const size_t t1 = UnbrokenSimpleRoots.size();
	const size_t t2 = BrokenSimpleRoots.size();
	const size_t t3 = UnbrokenU1Directions.size();

	if (info)
	{
		(*Print.out) << "\n  " << Print.cbegin << "Unbroken gauge group:" << Print.cend << "\n";
		(*Print.out) << "    " << Print.cbegin << "Simple roots of non-Abelian factors:" << Print.cend << "\n";
		if (t1 == 0)
			(*Print.out) << "      " << Print.cbegin << "none" << Print.cend << "\n";
		else
		{
			for (i = 0; i < t1; ++i)
			{
				(*Print.out) << "      ";
				Print.PrintRational(UnbrokenSimpleRoots[i], E8xE8);
				(*Print.out) << "\n";
			}
		}

		(*Print.out) << "    " << Print.cbegin << "U(1) generators:" << Print.cend << "\n";
		if (t3 == 0)
			(*Print.out) << "      " << Print.cbegin << "none" << Print.cend << "\n";
		else
		{
			for (i = 0; i < t3; ++i)
			{
				(*Print.out) << "      ";
				Print.PrintRational(UnbrokenU1Directions[i], E8xE8);
				(*Print.out) << "\n";
			}
		}

		(*Print.out) << "\n  " << Print.cbegin << "Other non-Abelian factors:" << Print.cend << "\n";
		(*Print.out) << "    " << Print.cbegin << "Simple roots of broken non-Abelian factors:" << Print.cend << "\n";
		if (t2 == 0)
			(*Print.out) << "      " << Print.cbegin << "none" << Print.cend << "\n";
		else
		{
			for (i = 0; i < t2; ++i)
			{
				(*Print.out) << "      ";
				Print.PrintRational(BrokenSimpleRoots[i], E8xE8);
				(*Print.out) << "\n";
			}
		}
	}

	vector<vector<double> > Orthogonal2SimpleRoots;
	Find_Basis_Of_Orthogonal_Space(BrokenSimpleRoots, UNSPECIFIED_LATTICE, 16, Orthogonal2SimpleRoots);

	for (i = 0; i < s1; ++i)
		Orthogonal2SimpleRoots.push_back(Weights[i]);
	Orthogonal2SimpleRoots = findBasis<double>(Orthogonal2SimpleRoots);
	vector<vector<double> > NewSimpleRoots;
	Find_Basis_Of_Orthogonal_Space(Orthogonal2SimpleRoots, UNSPECIFIED_LATTICE, 16, NewSimpleRoots);

	const size_t t4 = NewSimpleRoots.size();

	if (info)
	{
		(*Print.out) << "    " << Print.cbegin << "Leading to simple roots of unbroken subgroups:" << Print.cend << "\n";
		if (t4 == 0)
			(*Print.out) << "      " << Print.cbegin << "none" << Print.cend << "\n";
		else
		{
			for (i = 0; i < t4; ++i)
			{
				(*Print.out) << "      ";
				Print.PrintRational(NewSimpleRoots[i], E8xE8);
				(*Print.out) << "\n";
			}
		}
	}

	vector<vector<double> > OrthogonalBasis;
	for (i = 0; i < s1; ++i)
		OrthogonalBasis.push_back(Weights[i]);

	OrthogonalBasis.insert(OrthogonalBasis.end(), UnbrokenSimpleRoots.begin(), UnbrokenSimpleRoots.end());
	OrthogonalBasis.insert(OrthogonalBasis.end(), NewSimpleRoots.begin(), NewSimpleRoots.end());
	OrthogonalBasis.insert(OrthogonalBasis.end(), UnbrokenU1Directions.begin(), UnbrokenU1Directions.end());

	OrthogonalBasis = findBasis<double>(OrthogonalBasis);

	vector<vector<double> > UnbrokenU1s;
	Find_Basis_Of_Orthogonal_Space(OrthogonalBasis, UNSPECIFIED_LATTICE, 16, UnbrokenU1s);

	const size_t t6 = UnbrokenU1s.size();

	if (info)
	{
		(*Print.out) << "\n  " << Print.cbegin << "Further unbroken U(1) generators:" << Print.cend << "\n";
		if (t6 == 0)
			(*Print.out) << "      " << Print.cbegin << "none" << Print.cend << "\n";
	}

	for (i = 0; i < t6; ++i)
	{
		const vector<double> &UnbrokenU1 = UnbrokenU1s[i];

		if (info)
		{
			(*Print.out) << "      ";
			Print.PrintRational(UnbrokenU1, E8xE8);
		}

		bool not_from_sr = true;
		double sp = 0.0;

		for (j = 0; not_from_sr && (j < t2); ++j)
		{
			const vector<double> &BrokenSimpleRoot = BrokenSimpleRoots[j];

			sp = 0.0;
			for (k = 0; k < dim; ++k)
				sp += UnbrokenU1[k] * BrokenSimpleRoot[k];

			if (fabs(sp) >= prec)
			{
				if (info)
					(*Print.out) << "  " << Print.cbegin << "Unbroken U(1) subgroup of non-Abelian factor." << Print.cend;

				new_U1s_AdditionalLabels.push_back("from" + BrokenGaugeGroups[j]);
				not_from_sr = false;
			}
		}
		if (not_from_sr)
			new_U1s_AdditionalLabels.push_back("");

		if (info)
			(*Print.out) << endl;
	}

	vector<vector<double> > &NewU1Directions = VEVConfig.SymmetryGroup.GaugeGroup.u1directions;
	NewU1Directions.clear();
	NewU1Directions.insert(NewU1Directions.end(), UnbrokenU1Directions.begin(), UnbrokenU1Directions.end());
	NewU1Directions.insert(NewU1Directions.end(), UnbrokenU1s.begin(), UnbrokenU1s.end());

	// begin: find broken U1s
	vector<vector<double> > FindBrokenU1s = NewU1Directions;
	FindBrokenU1s.insert(FindBrokenU1s.end(), NewSimpleRoots.begin(), NewSimpleRoots.end());
	FindBrokenU1s.insert(FindBrokenU1s.end(), UnbrokenSimpleRoots.begin(), UnbrokenSimpleRoots.end());

	vector<vector<double> > BrokenU1s;
	Find_Basis_Of_Orthogonal_Space(FindBrokenU1s, UNSPECIFIED_LATTICE, 16, BrokenU1s);

	const size_t t7 = BrokenU1s.size();

	for (i = 0; i < t7; ++i)
		new_U1s_AdditionalLabels.push_back("broken");

	if (info)
	{
		(*Print.out) << "\n  " << Print.cbegin << "Broken U(1) generators:" << Print.cend << "\n";
		if (t7 == 0)
			(*Print.out) << "      " << Print.cbegin << "none" << Print.cend << "\n";
		else
		{
			for (i = 0; i < t7; ++i)
			{
				(*Print.out) << "      ";
				Print.PrintRational(BrokenU1s[i], E8xE8);
				(*Print.out) << "\n";
			}
		}
	}
	// end: find broken U1s

	NewU1Directions.insert(NewU1Directions.end(), BrokenU1s.begin(), BrokenU1s.end());

	if (!anom_U1_unbroken)
	{
		VEVConfig.SymmetryGroup.IsFirstU1Anomalous = false;
		VEVConfig.SymmetryGroup.D0_FI_term = 0.0;
	}

	VEVConfig.SymmetryGroup.observable_sector_U1s.clear();
	if (NewU1Directions.size() == 1)
		VEVConfig.SymmetryGroup.observable_sector_U1s.push_back(0);
	else
		VEVConfig.SymmetryGroup.observable_sector_U1s = Select_these_U1s;

	VEVConfig.SymmetryGroup.U1s_AdditionalLabels = new_U1s_AdditionalLabels;
	Orbifold.Config_ComputeNewU1Charges(VEVConfig);
	return true;
}



/* ########################################################################################
######   FindInvariantConfig(...)                                                    ######
######                                                                               ######
######   Version: 28.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Orbifold          : the orbifold of the vev-config "OriginalVEVConfig"   ######
######   2) Background        : a U(1) gauge background                              ######
######   3) OriginalVEVConfig : the SConfig object "OriginalVEVConfig" must contain  ######
######                          the gauge group at the orbifold point                ######
######   4) NewVEVConfig      : the resulting vev-config                             ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Computes the spectrum when a U(1) gauge background has been turned on.      ######
######################################################################################## */
bool CAnalyseModel::FindInvariantConfig(const COrbifold &Orbifold, const CVector &Background, const SConfig &OriginalVEVConfig, SConfig &NewVEVConfig)
{
	if (Background.GetSize() != 16)
	{
		cout << "\n  Warning in bool CAnalyseModel::FindInvariantConfig(...) : \"Background\" ill-defined. Return false." << endl;
		return false;
	}
	if (Background.IsZero())
		return true;

	const double prec = 0.0001;

	unsigned i = 0;
	unsigned j = 0;
	unsigned k = 0;
	unsigned l = 0;
	unsigned m = 0;

	size_t s1 = 0;
	size_t s2 = 0;
	size_t s3 = 0;

	bool not_charged1 = true;
	bool not_charged2 = true;

	// check that "Background" is not charged with respect to two non-Ableian gauge group factors
	const CGaugeGroup &OrigGaugeGroup = OriginalVEVConfig.SymmetryGroup.GaugeGroup;
	s1 = OrigGaugeGroup.factor.size();
	for (i = 0; (not_charged1 || not_charged2) && (i < s1); ++i)
	{
		const vector<vector<double> > &ggf_SimpleRoots = OrigGaugeGroup.factor[i].simpleroots;

		s2 = ggf_SimpleRoots.size();
		for (j = 0; j < s2; ++j)
		{
			if (fabs(ggf_SimpleRoots[j] * Background) > prec)
			{
				if (not_charged1)
				{
					not_charged1 = false;
					break;
				}
				else
				{
					not_charged2 = false;
					break;
				}
			}
		}
	}
	if (!not_charged1 && !not_charged2)
		return false;

	const vector<CField>  &Fields  = OriginalVEVConfig.Fields;
	const vector<CSector> &Sectors = Orbifold.GetSectors();

	SSymmetryGroup &NewSymmetryGroup = NewVEVConfig.SymmetryGroup;
	CGaugeGroup    &GaugeGroup       = NewSymmetryGroup.GaugeGroup;

	vector<vector<double> > Roots;

	// begin: find unbroken gauge group
	s1 = Fields.size();
	for (i = 0; i < s1; ++i)
	{
		const CField &Field = Fields[i];
		if (Field.Multiplet == Vector)
		{
			s2 = Field.GetNumberOfLMWeights();
			for (j = 0; j < s2; ++j)
			{
				const CVector &Weight = Field.GetLMWeight(j, Sectors);
				if (fabs(Weight * Background) < prec)
					Roots.push_back(Weight);
			}
		}
	}
	if (Roots.size() != 0)
		GaugeGroup = determineAlgebra(Roots);
	else
	{
		GaugeGroup.algebra = "";
		GaugeGroup.factor.clear();
	}
	GaugeGroup.u1directions.clear();
	// end: find unbroken gauge group

	// begin: create U(1) directions
	vector<vector<double> > SimpleRoots;

	const size_t number_of_factors = GaugeGroup.factor.size();
	for (i = 0; i < number_of_factors; ++i)
	{
		const vector<vector<double> > &ggf_SimpleRoots = GaugeGroup.factor[i].simpleroots;
		SimpleRoots.insert(SimpleRoots.end(), ggf_SimpleRoots.begin(), ggf_SimpleRoots.end());
	}
	GaugeGroup.u1directions.push_back(Background);

	if (SimpleRoots.size() != 15)
	{
		SimpleRoots.push_back(Background);

		vector<vector<double> > D_U1Generators;
		const SelfDualLattice GaugeLattice = Orbifold.OrbifoldGroup.GetLattice();

		if (!Find_Basis_Of_Orthogonal_Space(SimpleRoots, GaugeLattice, 16, D_U1Generators))
		{
			cout << "\n  Warning in bool CAnalyseModel::FindInvariantConfig(...) : Return false." << endl;
			return false;
		}

		GaugeGroup.u1directions.insert(GaugeGroup.u1directions.end(), D_U1Generators.begin(), D_U1Generators.end());
	}
	const size_t number_of_U1s = GaugeGroup.u1directions.size();
	// end: create U(1) directions

	NewSymmetryGroup.observable_sector_GGs.clear();
	NewSymmetryGroup.observable_sector_U1s.clear();
	NewSymmetryGroup.GGs_AdditionalLabels.clear();
	NewSymmetryGroup.U1s_AdditionalLabels.clear();

	for (i = 0; i < number_of_factors; ++i)
	{
		NewSymmetryGroup.observable_sector_GGs.push_back(i);
		NewSymmetryGroup.GGs_AdditionalLabels.push_back("");
	}

	for (i = 0; i < number_of_U1s; ++i)
	{
		NewSymmetryGroup.observable_sector_U1s.push_back(i);
		NewSymmetryGroup.U1s_AdditionalLabels.push_back("");
	}
	NewSymmetryGroup.U1s_AdditionalLabels[0] = "flux";

	NewSymmetryGroup.Position_of_and_in_GaugeGroup = 0;
	NewSymmetryGroup.IsFirstU1Anomalous = true;

	vector<CVector>                  List_U1Charges;
	vector<vector<bool> >            List_NonAbelianRep;
	vector<int>                      List_Multiplicities;

	vector<unsigned>                 tmp_Indices;
	vector<vector<unsigned> >        List_Indices;

	long double  Charge = 0.0;
	CVector      U1Charges(number_of_U1s);
	vector<bool> NonAbelianRep(number_of_factors, false);
	string       new_Label = "";

	NewVEVConfig.Fields.clear();

	// begin: decompose the fields
	s1 = Fields.size();
	for (i = 0; i < s1; ++i)
	{
		const CField &Field = Fields[i];

		List_Indices.clear();
		List_U1Charges.clear();
		List_NonAbelianRep.clear();

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// begin: sort all weights according to their "U1Charges" and non-Abelian charges "NonAbelianRep"
		s2 = Field.GetNumberOfLMWeights();
		for (j = 0; j < s2; ++j)
		{
			const CVector &weight = Field.GetLMWeight(j, Sectors);

			// begin: compute the U(1) charges
			for (k = 0; k < number_of_U1s; ++k)
			{
				const doubleVector &U1Direction = GaugeGroup.u1directions[k];

				// compute the U(1) charge
				Charge = 0.0;
				for (l = 0; l < 16; ++l)
					Charge += U1Direction[l] * weight[l];

				RoundCharge(Charge);
				U1Charges[k] = (double)Charge;
			}
			// end: compute the U(1) charges

			// begin: is the weight charged under the k-th gauge group factor
			NonAbelianRep.assign(number_of_factors, false);
			for (k = 0; k < number_of_factors; ++k)
			{
				const vector<vector<double> > &SimpleRoots = GaugeGroup.factor[k].simpleroots;
				s3 = SimpleRoots.size();
				for (l = 0; !NonAbelianRep[k] && (l < s3); ++l)
				{
					const vector<double> &SimpleRoot = SimpleRoots[l];

					Charge = 0.0;
					for (m = 0; m < 16; ++m)
						Charge += SimpleRoot[m] * weight[m];

					if (fabs(Charge) > prec)
						NonAbelianRep[k] = true;
				}
			}
			// end: is the weight charged under the k-th gauge group factor

			// begin: sort
			bool unknown = true;
			s3 = List_U1Charges.size();
			for (k = 0; unknown && (k < s3); ++k)
			{
				if ((List_U1Charges[k] == U1Charges) && (List_NonAbelianRep[k] == NonAbelianRep))
				{
					List_Indices[k].push_back(Field.WeightIndices[j]);
					unknown = false;
				}
			}
			if (unknown)
			{
				List_U1Charges.push_back(U1Charges);
				List_NonAbelianRep.push_back(NonAbelianRep);

				tmp_Indices.clear();
				tmp_Indices.push_back(Field.WeightIndices[j]);
				List_Indices.push_back(tmp_Indices);
			}
			// end: sort
		}
		// end: sort all weights according to their "U1Charges" and non-Abelian charges "NonAbelianRep"
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////

		ostringstream os;
		os << Field.Numbers[NewVEVConfig.use_Labels];
		new_Label = Field.Labels[NewVEVConfig.use_Labels] + os.str();

		s2 = List_U1Charges.size();
		for (j = 0; j < s2; ++j)
		{
			CField NewField = Field;

			NewField.WeightIndices = List_Indices[j];
			NewField.U1Charges     = List_U1Charges[j];
			NewField.VEVs.Assign(NewField.WeightIndices.size(), 0.0);

			if (s2 != 1)
				NewField.Labels[NewVEVConfig.use_Labels] = new_Label;

			NewVEVConfig.Fields.push_back(NewField);
		}
	}
	// end: decompose the fields

	vector<vector<double> > Weights;
	cirrep<double>          Irrep;
	vector<int>             HW;
	SDimension              DimOfRep;
	CState                  tmpState;

	unsigned         Index = 0;
	bool             LabelNotKnown = true;
	vector<string>   KnownLabels;
	vector<unsigned> Counter;

	// begin: compute new non-Abelian charges
	s1 = NewVEVConfig.Fields.size();
	for (i = 0; i < s1; ++i)
	{
		CField &Field = NewVEVConfig.Fields[i];

		Field.Dimensions.clear();
		Field.HighestWeights_DL.clear();

		Weights.clear();
		s2 = Field.GetNumberOfLMWeights();
		for (j = 0; j < s2; ++j)
			Weights.push_back(Field.GetLMWeight(j, Sectors));

		Irrep = findHighestWeight<double>(GaugeGroup, Weights);

		for (j = 0; j < number_of_factors; ++j)
		{
			const gaugeGroupFactor<double> &ggf = GaugeGroup.factor[j];
			HW = findDynkinLabels(ggf, Irrep.highestweight);
			if (!tmpState.DetermineDimension(ggf, HW, DimOfRep))
				return false;

			Field.Dimensions.push_back(DimOfRep);
			Field.HighestWeights_DL.push_back(HW);
		}

		const string &CurrentLabel = Field.Labels[NewVEVConfig.use_Labels];

		LabelNotKnown = true;
		s2 = KnownLabels.size();
		for (j = 0; LabelNotKnown && (j < s2); ++j)
		{
			if (KnownLabels[j] == CurrentLabel)
			{
				Index = j;
				LabelNotKnown = false;
			}
		}
		if (LabelNotKnown)
		{
			Index = Counter.size();
			Counter.push_back(1);
			KnownLabels.push_back(CurrentLabel);
		}

		Field.Numbers[NewVEVConfig.use_Labels] = Counter[Index];
		++Counter[Index];
	}
	// end: compute new non-Abelian charges

	return true;
}

