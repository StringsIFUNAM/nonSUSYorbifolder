#include "crandommodel.h"
#include "globalfunctions.h"
#include "cprint.h"

#define CHECKERROR true


// used by:
// bool CRandomModel::Initiate(...)
// translate the integer "number" to a binary number
void binary(int number, vector<unsigned> &bin)
{
	int remainder;

	if(number <= 1)
	{
		bin.push_back(number);
		return;
	}

	remainder = number%2;
	binary(number >> 1, bin);
	bin.push_back(remainder);
}



/* ########################################################################################
######   CRandomModel(const SelfDualLattice &Lattice)                                ######
######                                                                               ######
######   Version: 16.09.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Lattice : E8xE8 or Spin32                                                ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard constructor of a CRandomModel object. Creates the E8xE8 or Spin32  ###### 
######   roots needed to create random shifts and Wilson lines. To initialize the    ######
######   random model creation call "Initiate(...)".                                 ######
######################################################################################## */
CRandomModel::CRandomModel(const SelfDualLattice &Lattice)
{
	// begin: create the E8 x E8' or SO(32) gauge group
	const CVector Null_Shift(16);
	S_OscillatorExcitation Excitation;
	Excitation.NumberOperator = 0.0;
	Excitation.ZeroPointEnergy = 2.0;
	CMasslessHalfState Roots10D(LeftMover, Excitation);
	Roots10D.SolveMassEquation(Null_Shift, Lattice);

	if (Roots10D.Weights.size() != 480)
	{
		cout << "\n  Warning in CRandomModel::CRandomModel(...) : check the 10D gauge group. Return." << endl;
		return;
	}
	this->LatticeVectors = Roots10D.Weights;

	// end: create the E8 x E8' or SO(32) gauge group

	this->MAX_Emergency_Exit = 250000;
}



/* ########################################################################################
######   ~CRandomModel()                                                             ######
######                                                                               ######
######   Version: 17.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard destructor of a CWilsonLine object.                                ######
######################################################################################## */
CRandomModel::~CRandomModel()
{
}

/* ########################################################################################
######   Obtain_Symmetry_Factors(const vector<unsigned> &oldfactors, ...) const      ######
######                                                                               ######
######   Version: 17.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) oldfactors  : the permutation symmetry factors without "Shift_or_WL"     ######
######   2) Shift_or_WL : a CVector object                                           ######
######   output:                                                                     ######
######   return value   : the new permutation symmetry factors respecting            ######
######                    "oldfactors" and the permutation symmetry of "Shift_or_WL" ######
###########################################################################################
######   description:                                                                ######
######   Analyzes the permutation symmetry of the vector "Shift_or_WL" starting from ######
######   the permutation symmetry stored in "oldfactors". Returns the result.        ######
######################################################################################## */
vector<unsigned> CRandomModel::Obtain_Symmetry_Factors(const vector<unsigned> &oldfactors, const CVector &Shift_or_WL) const
{
	vector<unsigned> factors;
	unsigned sum = 0, id_elem = 0, j, k;
	double current = 0.0, tmp;

	const size_t s2 = oldfactors.size();
	const double prec = 0.0001;

	for (j = 0; j < s2; ++j)
	{
		current = 0.0;
		id_elem = 0;
		for (k = sum; k < sum + oldfactors[j]; ++k)
		{
			tmp = Shift_or_WL[k];
			if (k == sum) current = tmp;
			if (fabs(current - tmp) < prec)
				id_elem++;
			else
			{
				factors.push_back(id_elem);
				id_elem = 1;
				current = tmp;
			}
			if (k == sum + oldfactors[j] - 1)
				factors.push_back(id_elem);
		}
		sum += oldfactors[j];
	}
	return factors;
}



/* ########################################################################################
######   Initiate(COrbifoldGroup &OrbifoldGroup, ...)                                ######
######                                                                               ######
######   Version: 29.09.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) OrbifoldGroup               : need for the relations between the Wilson  ######
######                                    lines, e.g. W_1 = W_2                      ######
######   2) UseOrigShiftsAndWilsonLines : a vector of boolean in the order V_1, V_2, ######
######                                    W_1,...,W_6; true if the original V or W   ######
######                                    shall be taken for the new model, false    ######
######                                    otherwise                                  ######
######   3) UnbrokenRoots               : the new model shall leave these roots      ######
######                                    unbroken                                   ######
######   output:                                                                     ######
######   return value                   : finished successfully?                     ######
###########################################################################################
######   description:                                                                ######
######   Initialize the random model creation. Now one can call "CreateRandom(...)"  ######
######   in the class "COrbifoldGroup" to create new models.                         ######
######################################################################################## */
bool CRandomModel::Initiate(const COrbifoldGroup &OrbifoldGroup, const vector<bool> &UseOrigShiftsAndWilsonLines, const vector<CVector> &UnbrokenRoots)
{
	if ((OrbifoldGroup.GetModelIndependent_CheckStatus() != CheckedAndGood) || (OrbifoldGroup.GetModularInvariance_CheckStatus() != CheckedAndGood))
	{
		cout << "\n  Warning in bool CRandomModel::Initiate(...) : Original orbifold group ill-defined. Return false." << endl;
		return false;
	}

	srand ( time(NULL) );

	const CSpaceGroup          &SpaceGroup = OrbifoldGroup.GetSpaceGroup();
	const vector<CTwistVector> &Twists     = SpaceGroup.GetTwists();
	const SelfDualLattice      Lattice     = OrbifoldGroup.GetLattice();

	this->v1v2 = Twists[1] * Twists[2];				//hacking here!!!
	this->v0v1 = Twists[0] * Twists[1];
	this->v0v2 = Twists[0] * Twists[2];

	unsigned max_length = 8;
	unsigned max_binary = 256; // = 2^(8-2)			//why not 2^8=256 ??
	if (Lattice == Spin32)
	{
		max_length = 16;
		max_binary = 16384; // = 2^(16-2)
	}
	unsigned max_length_minus_one = max_length - 1;
	unsigned max_length_plus_one  = max_length + 1;

	unsigned i = 0;
	unsigned j = 0;

	size_t s1 = 0;
	size_t s2 = 0;

	this->VectorialBlocks.clear();
	this->SpinorialBlocks.clear();

	vector<vector<vector<double> > > VectorialBlocks;
	vector<vector<vector<double> > > SpinorialBlocks;

	vector<vector<double> > ListOfPossibilities;
	vector<vector<double> > VectorialBlocksOfLength_a1;
	vector<vector<double> > SpinorialBlocksOfLength_a1;
	vector<double> Possibility(max_length,0);
	vector<double> VectorialPossibility;
	vector<double> SpinorialPossibility;
	vector<unsigned> bin;

	int a1 = 0;
	int a2 = 0;
	int a3 = 0;
	int a4 = 0;

	// Old routines
	// create vectors
	ListOfPossibilities.clear();
	  for (a1 = -4; a1 <= 4; ++a1)
	  {
	    Possibility[0] = a1;
	    for (a2 = -3; a2 <= 3; ++a2)
	    {
	      Possibility[1] = a2;
	      for (a3 = 0; a3 < max_binary; ++a3)
	      {
	        bin.clear();
	        binary(a3, bin);
	        s1 = bin.size();
	        for (a4 = 2; a4 < max_length-s1; ++a4)
	          Possibility[a4] = 0;
	        for (a4 = 0; a4 < s1; ++a4)
	          Possibility[max_length_minus_one-a4] = bin[s1-a4-1];

	        ListOfPossibilities.push_back(Possibility);
	      }
	    }
	  }
	  s1 = ListOfPossibilities.size();


	/*ListOfPossibilities.clear();
	for (a1 = -6; a1 <= 12; ++a1) // original upper bound: 4
	{
		Possibility[0] = a1;
		for (a2 = -5; a2 <= 11; ++a2) // original upper bound: 3
		{
			Possibility[1] = a2;

			for (int a22 = -4; a22 <= 10; ++a22) // originally did not exist
			{                                    // originally did not exist
				Possibility[2] = a22;                // originally did not exist

				for (int a23 = -3; a23 <= 10; ++a23) // originally did not exist
				{                                    // originally did not exist
					Possibility[3] = a23;                // originally did not exist

					int start_a4 = 4;
					if (a22 == 0)
					{
						if (a23 == 0)
							start_a4 = 2;
						else
							start_a4 = 3;
					}
					
					for (a3 = 0; a3 < max_binary; ++a3)
					{
						bin.clear();
						binary(a3, bin);
						s1 = bin.size();
						for (a4 = start_a4; a4 < max_length-s1; ++a4)  // originally start from 2
							Possibility[a4] = 0;
						for (a4 = 0; a4 < s1; ++a4)
							Possibility[max_length_minus_one-a4] = bin[s1-a4-1];

						ListOfPossibilities.push_back(Possibility);
					}
				}
			}
		}
	}
	s1 = ListOfPossibilities.size();*/

	vector<double> PossibleOrders;
	PossibleOrders.push_back(-1.0);
	PossibleOrders.push_back(-1.0);
	PossibleOrders.push_back(2.0);
	PossibleOrders.push_back(3.0);
	PossibleOrders.push_back(4.0);
	PossibleOrders.push_back(-1.0);
	PossibleOrders.push_back(6.0);
	PossibleOrders.push_back(7.0);
	PossibleOrders.push_back(8.0);
	PossibleOrders.push_back(-1.0);
	PossibleOrders.push_back(-1.0);
	PossibleOrders.push_back(-1.0);
	PossibleOrders.push_back(12.0);
	const size_t p1 = PossibleOrders.size();

	for (i = 0; i < p1; ++i)
	{
		const double Order = PossibleOrders[i];

		VectorialBlocks.clear();
		SpinorialBlocks.clear();

		if (Order != -1)
		{
			for (a1 = 1; a1 < max_length_plus_one; ++a1)
			{
				VectorialBlocksOfLength_a1.clear();
				SpinorialBlocksOfLength_a1.clear();
				VectorialPossibility.assign(a1,0);
				SpinorialPossibility.assign(a1,0);

				for (a2 = 0; a2 < s1; ++a2)
				{
					const vector<double> &Possibility = ListOfPossibilities[a2];
					for (a3 = 0; a3 < a1; ++a3)
					{
						VectorialPossibility[a3] = Possibility[a3]/Order;
						SpinorialPossibility[a3] = Possibility[a3];

						if (a3 == 0)
							SpinorialPossibility[a3] += 0.5;
						else
							SpinorialPossibility[a3] -= 0.5;

						SpinorialPossibility[a3] /= Order;
					}

					stable_sort(VectorialPossibility.begin(), VectorialPossibility.end());
					stable_sort(SpinorialPossibility.begin(), SpinorialPossibility.end());

					if (find(VectorialBlocksOfLength_a1.begin(), VectorialBlocksOfLength_a1.end(), VectorialPossibility) == VectorialBlocksOfLength_a1.end())
						VectorialBlocksOfLength_a1.push_back(VectorialPossibility);

					if (find(SpinorialBlocksOfLength_a1.begin(), SpinorialBlocksOfLength_a1.end(), SpinorialPossibility) == SpinorialBlocksOfLength_a1.end())
						SpinorialBlocksOfLength_a1.push_back(SpinorialPossibility);
				}
				VectorialBlocks.push_back(VectorialBlocksOfLength_a1);
				SpinorialBlocks.push_back(SpinorialBlocksOfLength_a1);
			}
		}
		this->VectorialBlocks.push_back(VectorialBlocks);
		this->SpinorialBlocks.push_back(SpinorialBlocks);
	}

	this->UseOriginalVectors = UseOrigShiftsAndWilsonLines;
	this->UnbrokenRoots      = UnbrokenRoots;

	if (this->UseOriginalVectors.size() != 9)				//hacking here!!!
			{
		cout << "\n  Warning in bool CRandomModel::Initiate(...) : UseOriginalVectors does not contain 9 entries. Return false." << endl;
		return false;
			}

	vector<CVector> OrigVector;							//hacking here!!! random generation only for V_1, V_2, Witten shift V_0 the same
	if (this->UseOriginalVectors[1])
		OrigVector.push_back(OrbifoldGroup.GetShift(1));
	if (this->UseOriginalVectors[2])
		OrigVector.push_back(OrbifoldGroup.GetShift(2));
	for (i = 0; i < 6; ++i)
	{
		if (this->UseOriginalVectors[i+3])
			OrigVector.push_back(OrbifoldGroup.GetWilsonLines().GetWilsonLine(i));
	}

	const vector<vector<unsigned> > &WL_Relations     = SpaceGroup.GetWL_Relations();
	const vector<unsigned>          &WL_AllowedOrders = SpaceGroup.GetWL_AllowedOrders();

	// values of CreateOrIdentifyWL:
	// -2      : create WL
	// -1      : do not create WL
	// i, else : identify WL with i-th WL
	this->CreateOrIdentifyWL.assign(6, -2);
	for (i = 0; i < 6; ++i)
	{
		if ((WL_AllowedOrders[i] == 1) || this->UseOriginalVectors[i+3])
			this->CreateOrIdentifyWL[i] = -1;
	}

	s1 = WL_Relations.size();
	for (i = 0; i < s1; ++i)
	{
		const vector<unsigned> &WL_Relation = WL_Relations[i];

		s2 = WL_Relation.size();
		for (j = 1; j < s2; ++j)
		{
			if (this->CreateOrIdentifyWL[WL_Relation[j]] == -2)
				this->CreateOrIdentifyWL[WL_Relation[j]] = WL_Relation[0];
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// begin: compute the original symmetry factors
	this->OrigFactors.clear();
	if (OrbifoldGroup.GetLattice() == E8xE8)
	{
		this->OrigFactors.push_back(8);
		this->OrigFactors.push_back(8);
	}
	else
		this->OrigFactors.push_back(16);

	unsigned NumberOfOrigShifts = OrigVector.size();
	for (i = 0; i < NumberOfOrigShifts; ++i)
		this->OrigFactors = Obtain_Symmetry_Factors(this->OrigFactors, OrigVector[i]);

	const size_t NumberOfUnbrokenRoots = this->UnbrokenRoots.size();
	for (i = 0; i < NumberOfUnbrokenRoots; ++i)
		this->OrigFactors = Obtain_Symmetry_Factors(this->OrigFactors, this->UnbrokenRoots[i]);
	// end: compute the original symmetry factors
	////////////////////////////////////////////////////////////////////////////////////////////////////////


	if (SpaceGroup.ShiftsWL_ScalarProductFactors.size() == 0)
	{
		cout << "\n  Warning in bool CRandomModel::Initiate(...) : The prefactors for the modular invariance conditions were not computed in the space group. Return false." << endl;
		return false;
	}

	return true;
}



/* ########################################################################################
######   CheckUnbrokenRoots(const CVector &Vector) const                             ######
######                                                                               ######
######   Version: 16.05.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Vector    : a shift or a Wilson line                                     ######
######   output:                                                                     ######
######   return value : are the roots of "UnbrokenRoots" unbroken?                   ######
###########################################################################################
######   description:                                                                ######
######   Checks whether the CVector object "Vector" leaves the roots (stored in the  ######
######   private member variable "UnbrokenRoots") unbroken, i.e. if the scalar       ######
######   product is integer.                                                         ######
######################################################################################## */
bool CRandomModel::CheckUnbrokenRoots(const CVector &Vector) const
{
	const size_t s1 = this->UnbrokenRoots.size();
	for (unsigned i = 0; i < s1; ++i)
	{
		if (!is_integer(Vector * this->UnbrokenRoots[i]))
			return false;
	}
	return true;
}



/* ########################################################################################
######   &GetVectorialBlocksOfOrder(const unsigned &Order) const                     ######
######                                                                               ######
######   Version: 16.05.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   Order        : order of the shift or Wilson line                            ######
######   output:                                                                     ######
######   return value : building blocks of various lengths for shifts/ Wilson lines  ######
###########################################################################################
######   description:                                                                ######
######   Gives access to the building blocks of various lengths.                     ######
######################################################################################## */
const vector<vector<vector<double> > > &CRandomModel::GetVectorialBlocksOfOrder(const unsigned &Order) const
{
#ifdef CHECKERROR
	if (Order > 12)
	{
		cout << "\n  Warning in const vector<vector<vector<double> > > &CRandomModel::GetVectorialBlocksOfOrder(...) const : Index i out of range. Set i = 0." << endl;
		return this->VectorialBlocks[0];
	}
#endif

	return this->VectorialBlocks[Order];
}



/* ########################################################################################
######   &GetSpinorialBlocksOfOrder(const unsigned &Order) const                     ######
######                                                                               ######
######   Version: 16.05.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   Order        : order of the shift or Wilson line                            ######
######   output:                                                                     ######
######   return value : building blocks of various lengths for shifts/ Wilson lines  ######
###########################################################################################
######   description:                                                                ######
######   Gives access to the building blocks of various lengths.                     ######
######################################################################################## */
const vector<vector<vector<double> > > &CRandomModel::GetSpinorialBlocksOfOrder(const unsigned &Order) const
{
#ifdef CHECKERROR
	if (Order > 12)
	{
		cout << "\n  Warning in const vector<vector<vector<double> > > &CRandomModel::GetSpinorialBlocksOfOrder(...) const : Index i out of range. Set i = 0." << endl;
		return this->SpinorialBlocks[0];
	}
#endif

	return this->SpinorialBlocks[Order];
}
