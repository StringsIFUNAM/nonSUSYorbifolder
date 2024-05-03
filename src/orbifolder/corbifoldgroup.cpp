#include <iostream> 
#include <cstdlib>

#include "corbifoldgroup.h"
#include "globalfunctions.h"
#include "cprint.h"
//#include "lattice.hpp"
#include "crandommodel.h"

using std::cout;
using std::endl;
using std::exit;
using std::flush;
using std::vector;


/* ########################################################################################
######   COrbifoldGroup()                                                            ######
######                                                                               ######
######   Version: 10.08.2011                                                         ######
######   Check-Level: 2                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard constructor of a COrbifoldGroup object. No content is specified.   ######
######################################################################################## */
COrbifoldGroup::COrbifoldGroup()
: UseFreelyActingWL(false),
  ModelIndependent_CheckStatus(NotChecked),
  ModularInvariance_CheckStatus(NotChecked),
  OrbifoldGroup_CheckStatus(NotChecked)
{
	CShiftVector ZeroShift;
	this->Shifts.assign(3, ZeroShift);			//hacking here!!!

	this->NumberOfSupersymmetry = -1;
	this->SetDiscreteTorsionToZero();
}



/* ########################################################################################
######   COrbifoldGroup(const CSpaceGroup &SpaceGroup)                               ######
######                                                                               ######
######   Version: 10.08.2011                                                         ######
######   Check-Level: 2                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) SpaceGroup : a CSpaceGroup object                                        ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Constructor of a COrbifoldGroup object. Creates an orbifold group based on  ######
######   the space group "SpaceGroup".                                               ######
######################################################################################## */
COrbifoldGroup::COrbifoldGroup(const CSpaceGroup &SpaceGroup)
: UseFreelyActingWL(false),
  SpaceGroup(SpaceGroup),
  ModelIndependent_CheckStatus(NotChecked),
  ModularInvariance_CheckStatus(NotChecked),
  OrbifoldGroup_CheckStatus(NotChecked)
{
	CShiftVector ZeroShift;
	this->Shifts.assign(3, ZeroShift);			//hacking here!!!

	this->SetDiscreteTorsionToZero();

	this->CreateModelIndependentPart();
	if (this->ModelIndependent_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in COrbifoldGroup::COrbifoldGroup(...): Model independent part of the orbifold group can not be created." << endl;
}



/* ########################################################################################
######   COrbifoldGroup(...)                                                         ######
######                                                                               ######
######   Version: 10.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) SpaceGroup  : a CSpaceGroup object                                       ######
######   1) Shifts      : two CShiftVector objects                                   ######
######   1) Wilsonlines : a CWilsonLines object                                      ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Constructor of a COrbifoldGroup object. Creates an orbifold group based on  ######
######   the space group "SpaceGroup" using the gauge embedding specified by "Shifts"######
######   and "Wilsonlines".                                                          ######
######################################################################################## */
COrbifoldGroup::COrbifoldGroup(const CSpaceGroup &SpaceGroup, const vector<CShiftVector> &Shifts, const CWilsonLines &Wilsonlines)
: UseFreelyActingWL(false),
  Shifts(Shifts),
  WilsonLines(Wilsonlines),
  SpaceGroup(SpaceGroup),
  ModelIndependent_CheckStatus(NotChecked),
  ModularInvariance_CheckStatus(NotChecked),
  OrbifoldGroup_CheckStatus(NotChecked)
{
	if (this->Shifts.size() != 3)		//hacking here!!!
	{
		cout << "\n  Warning in COrbifoldGroup::COrbifoldGroup(...) : three shifts need to be defined. Use V_1 = 0, V_2 = 0  for Z_M orbifolds and V_2 = 0 for Z_MxZ_N orbifolds." << endl;
		return;
	}

	this->FreelyActingWilsonLine.SetLattice(this->Shifts[0].Lattice);
	this->SetDiscreteTorsionToZero();

	this->CreateModelIndependentPart();
	if (this->ModelIndependent_CheckStatus != CheckedAndGood)
	{
		cout << "\n  Warning in COrbifoldGroup::COrbifoldGroup(...) : Model independent part of the orbifold group can not be created." << endl;
		return;
	}

	CPrint Print(Tstandard, &cout);
	this->CreateModelDependentPart(Print, true);
	if (this->OrbifoldGroup_CheckStatus != CheckedAndGood)
	{
		cout << "\n  Warning in COrbifoldGroup::COrbifoldGroup(...) : Orbifold group ill-defined." << endl;
		return;
	}
}



/* ########################################################################################
######   ~COrbifoldGroup()                                                           ######
######                                                                               ######
######   Version: 10.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard destructor of a COrbifoldGroup object.                             ######
######################################################################################## */
COrbifoldGroup::~COrbifoldGroup()
{
}



/* ########################################################################################
######   CreateRandom(const CRandomModel &RandomModel, bool CreateBrotherModel)      ######
######                                                                               ######
######   Version: 31.01.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) RandomModel        : a CRandomModel object                               ######
######   2) CreateBrotherModel : create a brother model (add lattice vectors to      ######
######                           shifts and Wilson lines)                            ######
######   output:                                                                     ######
######   return value          : random orbifold group created succesfully?          ######
###########################################################################################
######   description:                                                                ######
######   Randomly creates a modular invariant orbifold group.                        ######
###########################################################################################
######   update:                                                                     ######
######   The shift V_1 (i.e. "this->Shifts[0]") was not created.                     ######
######################################################################################## */
bool COrbifoldGroup::CreateRandom(const CRandomModel &RandomModel, bool CreateBrotherModel)
{
	bool verbose = false;
	unsigned Emergency_Exit = 0;
	const unsigned &MAX_Emergency_Exit = RandomModel.GetMAX_Emergency_Exit();

	bool modular_invariant = true;
	unsigned i = 0;
	unsigned j = 0;

	const vector<vector<double> > &ScalarProductFactors = this->SpaceGroup.ShiftsWL_ScalarProductFactors;
	/*for (int i=0; i<ScalarProductFactors.size(); i++) {
		for (int j=0; j<ScalarProductFactors.size(); j++) {
			cout<<ScalarProductFactors[i][j]<<" ";
		}
		cout<<endl;
	}*/

	vector<unsigned>        factors            = RandomModel.GetOrigFactors();
	const vector<bool>     &UseOriginalVectors = RandomModel.GetUseOriginalVectors();

	const bool ZMxZN = this->SpaceGroup.IsZMxZN();
	const bool ZMxZNxZK = this->SpaceGroup.IsZMxZNxZK();

	// values of CreateOrIdentifyWL:
	// -2      : create WL
	// -1      : do not create WL
	// i, else : identify WL with i-th WL

	vector<CWilsonLine> new_WLs;
	for (i = 0; i < LatticeDim; ++i)
		new_WLs.push_back(this->WilsonLines.GetWilsonLine(i));

	// create shift V_1
	/*if (!UseOriginalVectors[0])				//hacking here!!!		V_0 will be the Witten shift, no random generation
	{
		const double v1_sqr = this->SpaceGroup.GetTwist(0).GetSqrTo(4);
		const unsigned M = (unsigned)ScalarProductFactors[0][0];

		Emergency_Exit = 0;
		while (Emergency_Exit < MAX_Emergency_Exit)
		{
			++Emergency_Exit;

			// assume that the new shift is ok/modular invariant
			modular_invariant = true;
			if (CreateBrotherModel)
			{
				if (!this->Shifts[0].CreateRandomBrother(RandomModel.GetLatticeVectors(), v1_sqr, M) || !RandomModel.CheckUnbrokenRoots(this->Shifts[0]))
					modular_invariant = false;
			}
			else
			{
				if (!this->Shifts[0].CreateRandom(M, v1_sqr, factors, RandomModel.GetVectorialBlocksOfOrder(M), RandomModel.GetSpinorialBlocksOfOrder(M)) || !RandomModel.CheckUnbrokenRoots(this->Shifts[0]))
					modular_invariant = false;
			}

			if (modular_invariant)
			{
				// begin: check modular invariance
				// 1) with second shift
				if (ZMxZN && UseOriginalVectors[1])
				{
					// x (V_1 V_2 - v_1 v_2) = 0 mod 2
					if (!is_even(ScalarProductFactors[1][0] * ((this->Shifts[0] * this->Shifts[1]) - RandomModel.Get_v1v2())))
						modular_invariant = false;
				}
				// 2) with Wilson lines that need not be generated
				for (j = 0; modular_invariant && (j < LatticeDim); ++j)
				{
					if (RandomModel.GetCreateOrIdentifyWL(j) == -1)
					{
						// x V_1 W_j = 0 mod 2
						if (!is_even(ScalarProductFactors[0][j+2] * (this->Shifts[0] * new_WLs[j])))
							modular_invariant = false;
					}
				}
				// end: check modular invariance

				if (modular_invariant)
					break;
			}
		}
		if (Emergency_Exit == MAX_Emergency_Exit)
			return false;
		else
			factors = RandomModel.Obtain_Symmetry_Factors(factors, this->Shifts[0]);
	}*/

	// create shift V_1
	if ( (ZMxZN or ZMxZNxZK) && !UseOriginalVectors[1])
	{
		const double v1_sqr = this->SpaceGroup.GetTwist(1).GetSqrTo(4);
		const unsigned N = (unsigned)ScalarProductFactors[1][1];

		Emergency_Exit = 0;
		while (Emergency_Exit < MAX_Emergency_Exit)
		{
			++Emergency_Exit;

			// assume that the new shift is ok/modular invariant
			modular_invariant = true;
			if (CreateBrotherModel)
			{
				if (!this->Shifts[1].CreateRandomBrother(RandomModel.GetLatticeVectors(), v1_sqr, N) || !RandomModel.CheckUnbrokenRoots(this->Shifts[1]))
					modular_invariant = false;
			}
			else
			{
				if (!this->Shifts[1].CreateRandom(N, v1_sqr, factors, RandomModel.GetVectorialBlocksOfOrder(N), RandomModel.GetSpinorialBlocksOfOrder(N)) || !RandomModel.CheckUnbrokenRoots(this->Shifts[1]))
					modular_invariant = false;
			}

			if (modular_invariant)
			{
				// begin: check modular invariance
				// 1) with first shift
				// x (V_0 V_1 - v_0 v_1) = 0 mod 2
				if (!is_even(ScalarProductFactors[1][0] * ((this->Shifts[0] * this->Shifts[1]) - RandomModel.Get_v0v1()))) {
					modular_invariant = false;
				}

				// 2) additional constraints from modular invariance/orbifold periodicity, maybe too strong
				// e_8^T V_a = 0 mod 4		//hacking here!!!

				/*double sum = 0.0;
				CVector V_1=this->Shifts[1];
				for (int l=0; l<8; l++) {
					sum += V_1[l];
				}
				if ( !is_even(sum*0.5) ) {
					modular_invariant = false;
				}
				sum = 0.0;
				for (int l=8; l<16; l++) {
					sum += V_1[l];
				}
				if ( !is_even(sum*0.5) ) {
					modular_invariant = false;
				}*/

				// 3) with Wilson lines that need not be generated
				for (j = 0; modular_invariant && (j < LatticeDim); ++j)
				{
					if (RandomModel.GetCreateOrIdentifyWL(j) == -1)
					{
						// x V_1 W_j = 0 mod 2
						if (!is_even(ScalarProductFactors[1][j+3] * (this->Shifts[1] * new_WLs[j])))
							modular_invariant = false;
					}
				}
				// end: check modular invariance

				if (modular_invariant)
					break;
			}
		}
		if (Emergency_Exit == MAX_Emergency_Exit)
			return false;
		else
			factors = RandomModel.Obtain_Symmetry_Factors(factors, this->Shifts[1]);
	}

	// create shift V_2
	if (ZMxZNxZK && !UseOriginalVectors[2])
	{
		const double v2_sqr = this->SpaceGroup.GetTwist(2).GetSqrTo(4);
		const unsigned K = (unsigned)ScalarProductFactors[2][2];

		Emergency_Exit = 0;
		while (Emergency_Exit < MAX_Emergency_Exit)
		{
			++Emergency_Exit;

			// assume that the new shift is ok/modular invariant
			modular_invariant = true;
			if (CreateBrotherModel)
			{
				if (!this->Shifts[2].CreateRandomBrother(RandomModel.GetLatticeVectors(), v2_sqr, K) || !RandomModel.CheckUnbrokenRoots(this->Shifts[2]))
					modular_invariant = false;
			}
			else
			{
				if (!this->Shifts[2].CreateRandom(K, v2_sqr, factors, RandomModel.GetVectorialBlocksOfOrder(K), RandomModel.GetSpinorialBlocksOfOrder(K)) || !RandomModel.CheckUnbrokenRoots(this->Shifts[2]))
					modular_invariant = false;
			}

			if (modular_invariant)
			{
				// begin: check modular invariance
				// 1a) with first shift
				// x (V_0 V_2 - v_0 v_2) = 0 mod 2
				if (!is_even(ScalarProductFactors[2][0] * ((this->Shifts[0] * this->Shifts[2]) - RandomModel.Get_v0v2()))) {
					modular_invariant = false;
				}
				// 1b) with second shift
				// x (V_1 V_2 - v_1 v_2) = 0 mod 2
				if (!is_even(ScalarProductFactors[2][1] * ((this->Shifts[1] * this->Shifts[2]) - RandomModel.Get_v1v2()))) {
					modular_invariant = false;
				}

				// 2) additional constraints from modular invariance/orbifold periodicity, maybe too strong
				// e_8^T V_a = 0 mod 4		//hacking here!!!

				/*double sum = 0.0;
					CVector V_2=this->Shifts[1];
					for (int l=0; l<8; l++) {
						sum += V_2[l];
					}
					if ( !is_even(sum*0.5) ) {
						modular_invariant = false;
					}
					sum = 0.0;
					for (int l=8; l<16; l++) {
						sum += V_2[l];
					}
					if ( !is_even(sum*0.5) ) {
						modular_invariant = false;
					}*/

				// 3) with Wilson lines that need not be generated
				for (j = 0; modular_invariant && (j < LatticeDim); ++j)
				{
					if (RandomModel.GetCreateOrIdentifyWL(j) == -1)
					{
						// x V_2 W_j = 0 mod 2
						if (!is_even(ScalarProductFactors[2][j+3] * (this->Shifts[2] * new_WLs[j])))
							modular_invariant = false;
					}
				}
				// end: check modular invariance

				if (modular_invariant)
					break;
			}
		}
		if (Emergency_Exit == MAX_Emergency_Exit)
			return false;
		else
			factors = RandomModel.Obtain_Symmetry_Factors(factors, this->Shifts[2]);
	}

	const vector<unsigned>          &WL_AllowedOrders = this->SpaceGroup.GetWL_AllowedOrders();
	const vector<vector<unsigned> > &WL_Relations     = this->SpaceGroup.GetWL_Relations();

	// run through the six Wilson lines
	for (i = 0; i < LatticeDim; ++i)
	{
		// if the Wilson line needs to be generated
		if (RandomModel.GetCreateOrIdentifyWL(i) == -2)
		{
			CWilsonLine &WL_i = new_WLs[i];
			const unsigned O = WL_AllowedOrders[i];

			Emergency_Exit = 0;
			while (Emergency_Exit < MAX_Emergency_Exit)
			{
				++Emergency_Exit;

				// assume that the new Wilson line is ok/modular invariant
				modular_invariant = true;
				if (CreateBrotherModel)
				{
					if (!WL_i.CreateRandomBrother(RandomModel.GetLatticeVectors()) || !RandomModel.CheckUnbrokenRoots(WL_i))
						modular_invariant = false;
				}
				else
				{
					if (!WL_i.CreateRandom(O, factors, RandomModel.GetVectorialBlocksOfOrder(O), RandomModel.GetSpinorialBlocksOfOrder(O)) || !RandomModel.CheckUnbrokenRoots(WL_i))
						modular_invariant = false;
				}

				if (modular_invariant)
				{
					// begin: check modular invariance
					// 1) with shifts
					// x V_0 W_i = 0 mod 2
					if (!is_even(ScalarProductFactors[i+3][0] * (WL_i * this->Shifts[0])))
						modular_invariant = false;

					if ( (ZMxZN or ZMxZNxZK) && modular_invariant)
					{
						// x V_1 W_i = 0 mod 2
						if (is_even(ScalarProductFactors[i+3][1] * (WL_i * this->Shifts[1]))) {
							if (ZMxZNxZK) {
								// x V_2 W_i = 0 mod 2
								if (!is_even(ScalarProductFactors[i+3][2] * (WL_i * this->Shifts[2])))
									modular_invariant = false;
							}
						}
						else {
							modular_invariant = false;
						}
					}
					// 2) additional constraints from modular invariance/orbifold periodicity, maybe too strong
					// e_8^T W_a = 0 mod 4			//hacking here!!!

					/*double sum = 0.0;
					CVector WL=WL_i;
					for (int l=0; l<8; l++) {
						sum += WL[l];
					}
					if ( !is_even(sum*0.5) ) {
						modular_invariant = false;
					}
					sum = 0.0;
					for (int l=8; l<16; l++) {
						sum += WL[l];
					}
					if ( !is_even(sum*0.5) ) {
						modular_invariant = false;
					}*/

					// 3) with previously generated Wilson lines
					for (j = 0; modular_invariant && (j <= i); ++j)
					{
						if (RandomModel.GetCreateOrIdentifyWL(j) == -2)
						{
							// x W_i W_j = 0 mod 2
							if (!is_even(ScalarProductFactors[i+3][j+3] * (WL_i * new_WLs[j])))
								modular_invariant = false;
						}
					}
					// 4) with Wilson lines that need not be generated
					for (j = 0; modular_invariant && (j < LatticeDim); ++j)
					{
						if (RandomModel.GetCreateOrIdentifyWL(j) == -1)
						{
							// x W_i W_j = 0 mod 2
							if (!is_even(ScalarProductFactors[i+3][j+3] * (WL_i * new_WLs[j])))
								modular_invariant = false;
						}
					}
					// end: check modular invariance

					if (modular_invariant)
						break;
				}
			}
			if (Emergency_Exit == MAX_Emergency_Exit)
			{
				if (CreateBrotherModel)
					return false;
				else
					WL_i.SetToZero();
			}
			else
				factors = RandomModel.Obtain_Symmetry_Factors(factors, WL_i);
		}
	}

	for (i = 0; i < LatticeDim; ++i)
	{
		if (RandomModel.GetCreateOrIdentifyWL(i) >= 0)
			new_WLs[i] = new_WLs[RandomModel.GetCreateOrIdentifyWL(i)];
	}

	for (i = 0; i < LatticeDim; ++i)
		this->WilsonLines.SetWilsonLine(i, new_WLs[i]);

	this->WilsonLines.Check(WL_Relations, WL_AllowedOrders);

	if (this->UseFreelyActingWL)
	{
		if ((this->WilsonLines.GetWilsonLine(1) != this->WilsonLines.GetWilsonLine(3)) || (this->WilsonLines.GetWilsonLine(1) != this->WilsonLines.GetWilsonLine(5)))
		{
			cout << "Warning in bool CRandomModel::Create() : freely acting Wilson line failed!" << endl;
			return false;
		}
		this->FreelyActingWilsonLine = this->WilsonLines.GetWilsonLine(1) * 1.5;
	}

	this->ModularInvariance_CheckStatus = NotChecked;
	this->OrbifoldGroup_CheckStatus     = NotChecked;

	CPrint Print(Tstandard, &cout);
	this->CreateModelDependentPart(Print, true);
	if (this->OrbifoldGroup_CheckStatus != CheckedAndGood)
	{
	   if(verbose)
	   {
		cout << "\n  Warning in bool CRandomModel::Create() : Orbifold group ill-defined." << endl;
		return false;
	   }
	}
	return true;
}



/* ########################################################################################
######   CreateRandom_shifts(const CRandomModel &RandomModel, bool CreateBrotherModel)      ######
######                                                                               ######
######   Version: 31.01.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) RandomModel        : a CRandomModel object                               ######
######   2) CreateBrotherModel : create a brother model (add lattice vectors to      ######
######                           shifts and Wilson lines)                            ######
######   output:                                                                     ######
######   return value          : random orbifold group created succesfully?          ######
###########################################################################################
######   description:                                                                ######
######   Randomly creates a modular invariant orbifold group.                        ######
###########################################################################################
######   update:                                                                     ######
######   The shift V_1 (i.e. "this->Shifts[0]") was not created.                     ######
######################################################################################## */
bool COrbifoldGroup::CreateRandom_shifts(const CRandomModel &RandomModel, bool CreateBrotherModel)
{
	unsigned Emergency_Exit = 0;
	const unsigned &MAX_Emergency_Exit = RandomModel.GetMAX_Emergency_Exit();

	bool modular_invariant = true;
	unsigned i = 0;
	unsigned j = 0;

	const vector<vector<double> > &ScalarProductFactors = this->SpaceGroup.ShiftsWL_ScalarProductFactors;

	vector<unsigned>        factors            = RandomModel.GetOrigFactors();
	const vector<bool>     &UseOriginalVectors = RandomModel.GetUseOriginalVectors();

	const bool ZMxZN = this->SpaceGroup.IsZMxZN();
	const bool ZMxZNxZK = this->SpaceGroup.IsZMxZNxZK();

	// values of CreateOrIdentifyWL:
	// -2      : create WL
	// -1      : do not create WL
	// i, else : identify WL with i-th WL

	this->WilsonLines.SetToZero(E8xE8);					//Set Wilson lines to zero

	// create shift V_1
	if ( (ZMxZN or ZMxZNxZK) && !UseOriginalVectors[1])
	{
		const double v1_sqr = this->SpaceGroup.GetTwist(1).GetSqrTo(4);
		const unsigned N = (unsigned)ScalarProductFactors[1][1];

		Emergency_Exit = 0;
		while (Emergency_Exit < MAX_Emergency_Exit)
		{
			++Emergency_Exit;

			// assume that the new shift is ok/modular invariant
			modular_invariant = true;
			if (CreateBrotherModel)
			{
				if (!this->Shifts[1].CreateRandomBrother(RandomModel.GetLatticeVectors(), v1_sqr, N) || !RandomModel.CheckUnbrokenRoots(this->Shifts[1]))
					modular_invariant = false;
			}
			else
			{
				if (!this->Shifts[1].CreateRandom(N, v1_sqr, factors, RandomModel.GetVectorialBlocksOfOrder(N), RandomModel.GetSpinorialBlocksOfOrder(N)) || !RandomModel.CheckUnbrokenRoots(this->Shifts[1]))
					modular_invariant = false;
			}

			if (modular_invariant)
			{
				// begin: check modular invariance
				// 1) with first shift
				// x (V_0 V_1 - v_0 v_1) = 0 mod 2
				if (!is_even(ScalarProductFactors[1][0] * ((this->Shifts[0] * this->Shifts[1]) - RandomModel.Get_v0v1()))) {
					modular_invariant = false;
				}

				// 2) additional constraints from modular invariance/orbifold periodicity, maybe too strong
				// e_8^T V_a = 0 mod 4		//hacking here!!!

				/*double sum = 0.0;
				CVector V_1=this->Shifts[1];
				for (int l=0; l<8; l++) {
					sum += V_1[l];
				}
				if ( !is_even(sum*0.5) ) {
					modular_invariant = false;
				}
				sum = 0.0;
				for (int l=8; l<16; l++) {
					sum += V_1[l];
				}
				if ( !is_even(sum*0.5) ) {
					modular_invariant = false;
				}*/

				// end: check modular invariance

				if (modular_invariant)
					break;
			}
		}
		if (Emergency_Exit == MAX_Emergency_Exit)
			return false;
		else
			factors = RandomModel.Obtain_Symmetry_Factors(factors, this->Shifts[1]);
	}

	// create shift V_2
	if (ZMxZNxZK && !UseOriginalVectors[2])
	{
		const double v2_sqr = this->SpaceGroup.GetTwist(2).GetSqrTo(4);
		const unsigned K = (unsigned)ScalarProductFactors[2][2];

		Emergency_Exit = 0;
		while (Emergency_Exit < MAX_Emergency_Exit)
		{
			++Emergency_Exit;

			// assume that the new shift is ok/modular invariant
			modular_invariant = true;
			if (CreateBrotherModel)
			{
				if (!this->Shifts[2].CreateRandomBrother(RandomModel.GetLatticeVectors(), v2_sqr, K) || !RandomModel.CheckUnbrokenRoots(this->Shifts[2]))
					modular_invariant = false;
			}
			else
			{
				if (!this->Shifts[2].CreateRandom(K, v2_sqr, factors, RandomModel.GetVectorialBlocksOfOrder(K), RandomModel.GetSpinorialBlocksOfOrder(K)) || !RandomModel.CheckUnbrokenRoots(this->Shifts[2]))
					modular_invariant = false;
			}

			if (modular_invariant)
			{
				// begin: check modular invariance
				// 1a) with first shift
				// x (V_0 V_2 - v_0 v_2) = 0 mod 2
				if (!is_even(ScalarProductFactors[2][0] * ((this->Shifts[0] * this->Shifts[2]) - RandomModel.Get_v0v2()))) {
					modular_invariant = false;
				}
				// 1b) with second shift
				// x (V_1 V_2 - v_1 v_2) = 0 mod 2
				if (!is_even(ScalarProductFactors[2][1] * ((this->Shifts[1] * this->Shifts[2]) - RandomModel.Get_v1v2()))) {
					modular_invariant = false;
				}

				// 2) additional constraints from modular invariance/orbifold periodicity, maybe too strong
				// e_8^T V_a = 0 mod 4		//hacking here!!!

				/*double sum = 0.0;
					CVector V_2=this->Shifts[1];
					for (int l=0; l<8; l++) {
						sum += V_2[l];
					}
					if ( !is_even(sum*0.5) ) {
						modular_invariant = false;
					}
					sum = 0.0;
					for (int l=8; l<16; l++) {
						sum += V_2[l];
					}
					if ( !is_even(sum*0.5) ) {
						modular_invariant = false;
					}*/

				// end: check modular invariance

				if (modular_invariant)
					break;
			}
		}
		if (Emergency_Exit == MAX_Emergency_Exit)
			return false;
		else
			factors = RandomModel.Obtain_Symmetry_Factors(factors, this->Shifts[2]);
	}


	this->ModularInvariance_CheckStatus = NotChecked;
	this->OrbifoldGroup_CheckStatus     = NotChecked;

	CPrint Print(Tstandard, &cout);
	this->CreateModelDependentPart(Print, true);
	if (this->OrbifoldGroup_CheckStatus != CheckedAndGood)
	{
		cout << "\n  Warning in bool CRandomModel::Create() : Orbifold group ill-defined." << endl;
		return false;
	}
	return true;
}



/* ########################################################################################
######   CreateRandom_internal(const CRandomModel &RandomModel, bool CreateBrotherModel)      ######
######                                                                               ######
######   Version: 31.01.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) RandomModel        : a CRandomModel object                               ######
######   2) CreateBrotherModel : create a brother model (add lattice vectors to      ######
######                           shifts and Wilson lines)                            ######
######   output:                                                                     ######
######   return value          : random orbifold group created succesfully?          ######
###########################################################################################
######   description:                                                                ######
######   Randomly creates a modular invariant orbifold group.                        ######
###########################################################################################
######   update:                                                                     ######
######   The shift V_1 (i.e. "this->Shifts[0]") was not created.                     ######
######################################################################################## */
bool COrbifoldGroup::CreateRandom_internal(const CRandomModel &RandomModel, bool CreateBrotherModel)
{
	unsigned Emergency_Exit = 0;
	const unsigned &MAX_Emergency_Exit = RandomModel.GetMAX_Emergency_Exit();

	bool modular_invariant = true;
	unsigned i = 0;
	unsigned j = 0;

	const vector<vector<double> > &ScalarProductFactors = this->SpaceGroup.ShiftsWL_ScalarProductFactors;
	/*for (int i=0; i<ScalarProductFactors.size(); i++) {
		for (int j=0; j<ScalarProductFactors.size(); j++) {
			cout<<ScalarProductFactors[i][j]<<" ";
		}
		cout<<endl;
	}*/

	vector<unsigned>        factors            = RandomModel.GetOrigFactors();
	const vector<bool>     &UseOriginalVectors = RandomModel.GetUseOriginalVectors();

	const bool ZMxZN = this->SpaceGroup.IsZMxZN();
	const bool ZMxZNxZK = this->SpaceGroup.IsZMxZNxZK();

	// values of CreateOrIdentifyWL:
	// -2      : create WL
	// -1      : do not create WL
	// i, else : identify WL with i-th WL

	vector<CWilsonLine> new_WLs;
	for (i = 0; i < LatticeDim; ++i)
		new_WLs.push_back(this->WilsonLines.GetWilsonLine(i));

	//hacking here!!! No shift generation

	const vector<unsigned>          &WL_AllowedOrders = this->SpaceGroup.GetWL_AllowedOrders();
	const vector<vector<unsigned> > &WL_Relations     = this->SpaceGroup.GetWL_Relations();

	// run through the six Wilson lines
	for (i = 0; i < LatticeDim; ++i)
	{
		// if the Wilson line needs to be generated
		if (RandomModel.GetCreateOrIdentifyWL(i) == -2)
		{
			CWilsonLine &WL_i = new_WLs[i];
			const unsigned O = WL_AllowedOrders[i];

			Emergency_Exit = 0;
			while (Emergency_Exit < MAX_Emergency_Exit)
			{
				++Emergency_Exit;

				// assume that the new Wilson line is ok/modular invariant
				modular_invariant = true;
				if (CreateBrotherModel)
				{
					if (!WL_i.CreateRandomBrother(RandomModel.GetLatticeVectors()) || !RandomModel.CheckUnbrokenRoots(WL_i))
						modular_invariant = false;
				}
				else
				{
					if (!WL_i.CreateRandom(O, factors, RandomModel.GetVectorialBlocksOfOrder(O), RandomModel.GetSpinorialBlocksOfOrder(O)) || !RandomModel.CheckUnbrokenRoots(WL_i))
						modular_invariant = false;
				}

				if (modular_invariant)
				{
					// begin: check modular invariance
					// 1) with shifts
					// x V_0 W_i = 0 mod 2
					if (!is_even(ScalarProductFactors[i+3][0] * (WL_i * this->Shifts[0])))
						modular_invariant = false;

					if ( (ZMxZN or ZMxZNxZK) && modular_invariant)
					{
						// x V_1 W_i = 0 mod 2
						if (is_even(ScalarProductFactors[i+3][1] * (WL_i * this->Shifts[1]))) {
							if (ZMxZNxZK) {
								// x V_2 W_i = 0 mod 2
								if (!is_even(ScalarProductFactors[i+3][2] * (WL_i * this->Shifts[2])))
									modular_invariant = false;
							}
						}
						else {
							modular_invariant = false;
						}
					}
					// 2) additional constraints from modular invariance/orbifold periodicity, maybe too strong
					// e_8^T W_a = 0 mod 4			//hacking here!!!

					/*double sum = 0.0;
					CVector WL=WL_i;
					for (int l=0; l<8; l++) {
						sum += WL[l];
					}
					if ( !is_even(sum*0.5) ) {
						modular_invariant = false;
					}
					sum = 0.0;
					for (int l=8; l<16; l++) {
						sum += WL[l];
					}
					if ( !is_even(sum*0.5) ) {
						modular_invariant = false;
					}*/

					// 3) with previously generated Wilson lines
					for (j = 0; modular_invariant && (j <= i); ++j)
					{
						if (RandomModel.GetCreateOrIdentifyWL(j) == -2)
						{
							// x W_i W_j = 0 mod 2
							if (!is_even(ScalarProductFactors[i+3][j+3] * (WL_i * new_WLs[j])))
								modular_invariant = false;
						}
					}
					// 4) with Wilson lines that need not be generated
					for (j = 0; modular_invariant && (j < LatticeDim); ++j)
					{
						if (RandomModel.GetCreateOrIdentifyWL(j) == -1)
						{
							// x W_i W_j = 0 mod 2
							if (!is_even(ScalarProductFactors[i+3][j+3] * (WL_i * new_WLs[j])))
								modular_invariant = false;
						}
					}
					// end: check modular invariance

					if (modular_invariant)
						break;
				}
			}
			if (Emergency_Exit == MAX_Emergency_Exit)
			{
				if (CreateBrotherModel)
					return false;
				else
					WL_i.SetToZero();
			}
			else
				factors = RandomModel.Obtain_Symmetry_Factors(factors, WL_i);
		}
	}

	for (i = 0; i < LatticeDim; ++i)
	{
		if (RandomModel.GetCreateOrIdentifyWL(i) >= 0)
			new_WLs[i] = new_WLs[RandomModel.GetCreateOrIdentifyWL(i)];
	}

	for (i = 0; i < LatticeDim; ++i)
		this->WilsonLines.SetWilsonLine(i, new_WLs[i]);

	this->WilsonLines.Check(WL_Relations, WL_AllowedOrders);

	if (this->UseFreelyActingWL)
	{
		if ((this->WilsonLines.GetWilsonLine(1) != this->WilsonLines.GetWilsonLine(3)) || (this->WilsonLines.GetWilsonLine(1) != this->WilsonLines.GetWilsonLine(5)))
		{
			cout << "Warning in bool CRandomModel::Create() : freely acting Wilson line failed!" << endl;
			return false;
		}
		this->FreelyActingWilsonLine = this->WilsonLines.GetWilsonLine(1) * 1.5;
	}

	this->ModularInvariance_CheckStatus = NotChecked;
	this->OrbifoldGroup_CheckStatus     = NotChecked;

	CPrint Print(Tstandard, &cout);
	this->CreateModelDependentPart(Print, true);
	if (this->OrbifoldGroup_CheckStatus != CheckedAndGood)
	{
		cout << "\n  Warning in bool CRandomModel::Create() : Orbifold group ill-defined." << endl;
		return false;
	}
	return true;
}




/* ########################################################################################
######   GetShiftVector(...) const                                                   ######
######                                                                               ######
######   Version: 10.08.2011                                                         ######
######   Check-Level: 2                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) constructing_Element : a CSpaceGroupElement object                       ######
######   output:                                                                     ######
######   1) result               : the local shift corresponding to the constructing ######
######                             element                                           ######
######   return value            : finished successfully?                            ######
###########################################################################################
######   description:                                                                ######
######   Returns the local shift (k V_1 + l V_2 + n_alpha W_alpha).                  ######
######################################################################################## */
bool COrbifoldGroup::GetShiftVector(const CSpaceGroupElement &constructing_Element, CShiftVector &result) const
{
	if (this->ModularInvariance_CheckStatus != CheckedAndGood)
	{
		cout << "\n  Waring in bool COrbifoldGroup::GetShiftVector(...) : Shift and Wilson lines not modular invariant. Return false." << endl;
		CShiftVector NullShift;
		result = NullShift;
		return false;
	}

	if (this->SpaceGroup.IsZMxZN())
	{
		result =   (this->Shifts[0] * constructing_Element.Get_m())
            																 + (this->Shifts[1] * constructing_Element.Get_n())
            																 + (this->WilsonLines * constructing_Element.GetLatticeElement());
		return true;
	}
	if (this->SpaceGroup.IsZMxZNxZK())
	{
		result =   (this->Shifts[0] * constructing_Element.Get_m())
	            																 + (this->Shifts[1] * constructing_Element.Get_n())
	            																 + (this->Shifts[2] * constructing_Element.Get_k())
	            																 + (this->WilsonLines * constructing_Element.GetLatticeElement());
		return true;
	}

	result = (this->Shifts[0] * constructing_Element.Get_k()) + (this->WilsonLines * constructing_Element.GetLatticeElement());
	return true;
}



/* ########################################################################################
######   GetTwistVector(...) const                                                   ######
######                                                                               ######
######   Version: 04.08.2011                                                         ######
######   Check-Level: 2                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) constructing_Element : a CSpaceGroupElement object                       ######
######   output:                                                                     ######
######   1) result               : the local twist corresponding to the constructing ######
######                             element                                           ######
######   return value            : finished successfully?                            ######
###########################################################################################
######   description:                                                                ######
######   Returns the local twist (k v_1 + l v_2).                                    ######
######################################################################################## */
bool COrbifoldGroup::GetTwistVector(const CSpaceGroupElement &constructing_Element, CTwistVector &result) const
{
	const vector<CTwistVector> &Twists = this->SpaceGroup.GetTwists();

	if (this->SpaceGroup.IsZMxZN())
	{
		result = (Twists[0] * constructing_Element.Get_m()) + (Twists[1] * constructing_Element.Get_n());
		return true;
	}
	if (this->SpaceGroup.IsZMxZNxZK())
	{
		result = (Twists[0] * constructing_Element.Get_m()) + (Twists[1] * constructing_Element.Get_n()) + (Twists[2] * constructing_Element.Get_k());
		return true;
	}

	result = Twists[0] * constructing_Element.Get_k();
	return true;
}



/* ########################################################################################
######   CheckModularInvariance(CPrint &Print, bool info)                            ######
######                                                                               ######
######   Version: 17.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Print     : if "info" is true print the output to this CPrint object     ######
######   2) info      : print info?                                                  ######
######   output:                                                                     ######
######   return value : is this orbifold group modular invariant?                    ######
###########################################################################################
######   description:                                                                ######
######   Checks whether this orbifold group is modular invariant and sets the        ######
######   private member variable "ModularInvariance_CheckStatus" accordingly.        ######
######################################################################################## */
bool COrbifoldGroup::CheckModularInvariance(CPrint &Print, bool info)
{
	bool verbose = false;
	const vector<vector<double> > &SP_factors = this->SpaceGroup.ShiftsWL_ScalarProductFactors;
	if (SP_factors.size() == 0)
	{
		cout << "\n  Warning in bool COrbifoldGroup::CheckModularInvariance(...) : the prefactors for the modular invariance conditions were not computed in the space group. Return false." << endl;
		return false;
	}

	bool IsModularInvariant = true;

	const unsigned M = this->SpaceGroup.GetM();
	const unsigned N = this->SpaceGroup.GetN();
	const unsigned K = this->SpaceGroup.GetK();

	const CShiftVector &ZM_Shift = this->Shifts[0];
	const CShiftVector &ZN_Shift = this->Shifts[1];
	const CShiftVector &ZK_Shift = this->Shifts[2];

	const vector<CTwistVector> &Twists = this->SpaceGroup.GetTwists();

	const CTwistVector &ZM_Twist = Twists[0];
	const CTwistVector &ZN_Twist = Twists[1];
	const CTwistVector &ZK_Twist = Twists[2];

	const int Order_ZM_Twist = ZM_Twist.OrderOfTwist();
	const int Order_ZM_Shift = ZM_Shift.OrderOfShift();

	double sp  = 0.0;

	unsigned i = 0;
	unsigned j = 0;

	//CTwistVector PG_Element;
	//vector<bool> Rotation_In_ComplexPlane(3, false);

	//  ZM  //////////////////////////////////////////////////////////////////////////////////////////////////////
	//check the order of the twist
	if ((Order_ZM_Twist != 1) && (Order_ZM_Twist != M))
	{
		if (info)
		{
			(*Print.out) << "  Warning in bool COrbifoldGroup::CheckModularInvariance(...) : Check order of the Z_M twist.\n";
			(*Print.out) << "  Order_ZM_Twist = " << Order_ZM_Twist << endl;
		}
		IsModularInvariant = false;
	}

	//check the order of the shift
	if ((M % Order_ZM_Shift) != 0)
	{
		if (info)
		{
			(*Print.out) << "  Warning in bool COrbifoldGroup::CheckModularInvariance(...) : Check order of the Z_M shift.\n";
			(*Print.out) << "  Order_ZM_Shift = " << Order_ZM_Shift << endl;
		}
		IsModularInvariant = false;
	}

	// M(V_0) (V_0^2 - v_0^2) = 0 mod 2
	sp = M * ((ZM_Shift * ZM_Shift) - (ZM_Twist * ZM_Twist));
	if (!is_even(sp))
	{
		if (info)
			(*Print.out) << "\n  " << Order_ZM_Shift << " (V_M^2 - v_M^2) = " << sp << " = 0 mod 2 failed!" << flush;
		IsModularInvariant = false;
	}

	const bool ZMxZN = this->SpaceGroup.IsZMxZN();
	const bool ZMxZNxZK = this->SpaceGroup.IsZMxZNxZK();
	//  ZM x ZN  or  ZM x ZN x ZK  /////////////////////////////////////////////////////////////////////////////////////////////////
	if (ZMxZN || ZMxZNxZK)		//hacking here!!!
	{
		const int Order_ZN_Twist = ZN_Twist.OrderOfTwist();
		const int Order_ZN_Shift = ZN_Shift.OrderOfShift();

		// check the order of the twist
		if ((Order_ZN_Twist != 1) && (Order_ZN_Twist != N))
		{
			if (info)
			{
				(*Print.out) << "  Warning in bool COrbifoldGroup::CheckModularInvariance(...) : Check order of the Z_N twist.\n";
				(*Print.out) << "  Order_ZN_Twist = " << Order_ZN_Twist << endl;
			}
			IsModularInvariant = false;
		}

		//check the order of the shift
		if ((N % Order_ZN_Shift) != 0)
		{
			if (info)
			{
				(*Print.out) << "  Warning in bool COrbifoldGroup::CheckModularInvariance(...) : Check order of the Z_N shift.\n";
				(*Print.out) << "  Order_ZN_Shift = " << Order_ZN_Shift << endl;
			}
			IsModularInvariant = false;
		}

		// N(V_1) (V_1^2 - v_1^2) = 0 mod 2
		sp = N * ((ZN_Shift * ZN_Shift) - (ZN_Twist * ZN_Twist));
		if (!is_even(sp))
		{
			if (info)
				(*Print.out) << "\n  " << Order_ZN_Shift << " (V_N^2 - v_N^2) = " << sp << " = 0 mod 2 failed!" << flush;
			IsModularInvariant = false;
		}

		// GGT(V_0, V_1) (V_0 x V_1 - v_0 x v_1) = 0 mod 2
		sp = SP_factors[0][1] * ((ZN_Shift * ZM_Shift) - (ZN_Twist * ZM_Twist));
		if (!is_even(sp))
		{
			if (info)
				(*Print.out) << "\n  " << SP_factors[0][1] << " (V_N x V_M - v_N x v_M) = " << sp << " = 0 mod 2 failed!" << flush;
			IsModularInvariant = false;
		}
		//  ZM x ZN x ZK /////////////////////////////////////////////////////////////////////////////////////////////////
		if (ZMxZNxZK) {
			const int Order_ZK_Twist = ZK_Twist.OrderOfTwist();
			const int Order_ZK_Shift = ZK_Shift.OrderOfShift();

			// check the order of the twist
			if ((Order_ZK_Twist != 1) && (Order_ZK_Twist != K))
			{
				if (info)
				{
					(*Print.out) << "  Warning in bool COrbifoldGroup::CheckModularInvariance(...) : Check order of the Z_K twist.\n";
					(*Print.out) << "  Order_ZK_Twist = " << Order_ZK_Twist << endl;
				}
				IsModularInvariant = false;
			}

			//check the order of the shift
			if ((K % Order_ZK_Shift) != 0)
			{
				if (info)
				{
					(*Print.out) << "  Warning in bool COrbifoldGroup::CheckModularInvariance(...) : Check order of the Z_K shift.\n";
					(*Print.out) << "  Order_ZK_Shift = " << Order_ZK_Shift << endl;
				}
				IsModularInvariant = false;
			}

			// K(V_2) (V_2^2 - v_2^2) = 0 mod 2
			sp = K * ((ZK_Shift * ZK_Shift) - (ZK_Twist * ZK_Twist));
			if (!is_even(sp))
			{
				if (info)
					(*Print.out) << "\n  " << Order_ZK_Shift << " (V_N^2 - v_N^2) = " << sp << " = 0 mod 2 failed!" << flush;
				IsModularInvariant = false;
			}

			// GGT(V_1, V_2) (V_1 x V_2 - v_1 x v_2) = 0 mod 2
			sp = SP_factors[1][2] * ((ZK_Shift * ZN_Shift) - (ZK_Twist * ZN_Twist));
			if (!is_even(sp))
			{
				if (info)
					(*Print.out) << "\n  " << SP_factors[0][1] << " (V_N x V_M - v_N x v_M) = " << sp << " = 0 mod 2 failed!" << flush;
				IsModularInvariant = false;
			}
			// GGT(V_0, V_2) (V_0 x V_2 - v_0 x v_2) = 0 mod 2
			sp = SP_factors[0][2] * ((ZK_Shift * ZM_Shift) - (ZK_Twist * ZM_Twist));
			if (!is_even(sp))
			{
				if (info)
					(*Print.out) << "\n  " << SP_factors[0][1] << " (V_N x V_M - v_N x v_M) = " << sp << " = 0 mod 2 failed!" << flush;
				IsModularInvariant = false;
			}

		}
	}

	//  Wilsonlines  ///////////////////////////////////////////////////////////////////////////////////////////////
	// run through all six Wilson lines
	for(i = 0; i < LatticeDim; ++i)
	{
		const CWilsonLine &W_i = this->WilsonLines.GetWilsonLine(i);

		// first check the order of the Wilson line
		if (!W_i.GetIs_Zero())
		{
			// then check the scalar product with the shifts
			// GGT(V_0, W_i) V_0 x W_i = 0 mod 2
			sp = SP_factors[i+3][0] * (W_i * ZM_Shift);

			if (!is_even(sp))
			{
				if (info)
					(*Print.out) << "\n  " << SP_factors[i+3][0] << " V_M" << " x W_" << i+1 << " = " << sp << " = 0 mod 2 failed!" << flush;
				IsModularInvariant = false;
			}


			if (ZMxZN || ZMxZNxZK)				//hacking here!!!
			{
				// GGT(V_1, W_i) V_1 x W_i = 0 mod 2
				sp = SP_factors[i+3][1] * (W_i * ZN_Shift);

				if (!is_even(sp))
				{
					if (info)
						(*Print.out) << "\n  " << SP_factors[i+3][1] << " V_N" << " x W_" << i+1 << " = " << sp << " = 0 mod 2 failed!" << flush;
					IsModularInvariant = false;
				}
				if (ZMxZNxZK) {
					sp = SP_factors[i+3][2] * (W_i * ZK_Shift);

					if (!is_even(sp))
					{
						if (info)
							(*Print.out) << "\n  " << SP_factors[i+3][2] << " V_K" << " x W_" << i+1 << " = " << sp << " = 0 mod 2 failed!" << flush;
						IsModularInvariant = false;
					}
				}
			}

			// and the scalar product with the other wilson lines
			for(j = i; j < LatticeDim; ++j)
			{
				const CWilsonLine &W_j = this->WilsonLines.GetWilsonLine(j);

				if (!W_j.GetIs_Zero())
				{
					sp = SP_factors[i+3][j+3] * W_i * W_j;		//hacking here!!! first three entries in SP_factors are for ZMxZNxZK shfits

					// begin: get the factor "ggt"
					/*if (i == j)
            ggt = SP_factors[i+2][i+2];
          else
          {
            // ??????????????????????????????????????????????????????????????????????????
            // Is there a point group element which rotates one WL, but not the other?
            // (for both point groups ZM and ZM x ZN)
            bool stop = false;

            for (k = 0; !stop && (k < M); ++k)
            {
              for (l = 0; !stop && (l < N); ++l)
              {
                PG_Element = (ZM_Twist * k) + (ZN_Twist * l);

                Rotation_In_ComplexPlane.assign(3, false);

                // run through the twist
                for (a = 1; a < 4; ++a)
                {
                  tmp = fabs(PG_Element[a]);

                  if (fabs(round2(tmp) - tmp) > prec)
                    Rotation_In_ComplexPlane[a-1] = true;
                }
                // if such an element exists:
                if (( Rotation_In_ComplexPlane[(unsigned)floor(((double)i)/2.0)] && !Rotation_In_ComplexPlane[(unsigned)floor(((double)j)/2.0)])
                 || (!Rotation_In_ComplexPlane[(unsigned)floor(((double)i)/2.0)] &&  Rotation_In_ComplexPlane[(unsigned)floor(((double)j)/2.0)]))
                {
                  stop = true;
                }
              }
            }
            // if such an element does not exist:
            if (!stop && (Allowed_Order_W_i == 2))
              ggt = Allowed_Order_W_i * Allowed_Order_W_j;
            else
              ggt = SP_factors[i+2][j+2];
          }*/
					// end: get the factor "ggt"

					//if ((i == 0) && (j == 0))
					//  (*Print.out) << "\n  " << SP_factors[i+2][j+2] << " W_" << i+1 << " x W_" << j+1 << " = " << sp << " = 0 mod 2" << flush;

					if (!is_even(sp))
					{
						if (info)
						if(verbose)
						{
							(*Print.out) << "\n  " << SP_factors[i+3][j+3] << " W_" << i+1 << " x W_" << j+1 << " = " << sp << " = 0 mod 2 failed!" << flush;
						IsModularInvariant = false;
					    } 
					}
				}
			}
		}
	}

	if (this->UseFreelyActingWL)
	{
		//out << "modular invariance with freely acting WL:" << endl;
		//out << "V_1 * W_tau = " << this->ZM_constructing_Element.Shift * this->FreelyActingWilsonLine << endl;
		//out << "V_2 * W_tau = " << this->ZN_constructing_Element.Shift * this->FreelyActingWilsonLine << endl;
		for (i = 0; i < LatticeDim; ++i)
		{
			sp = 2.0 * (this->WilsonLines.GetWilsonLine(i) * this->FreelyActingWilsonLine);
			if (!is_even(sp))
			{
				if (info)
					(*Print.out) << "\n  2 W_" << i+1 << " * W_tau = " << sp << " = 0 mod 2 failed!" << flush;
				IsModularInvariant = false;
			}
		}

		sp = 4.0 * (this->FreelyActingWilsonLine * this->FreelyActingWilsonLine);
		if (!is_even(sp))
		{
			if (info)
				(*Print.out) << "\n  4 W_tau * W_tau = " << sp << " = 0 mod 2 failed!" << flush;
			IsModularInvariant = false;
		}
	}
	if (!IsModularInvariant && info)
		(*Print.out) << endl;

	if (IsModularInvariant)
		this->ModularInvariance_CheckStatus = CheckedAndGood;
	else
		this->ModularInvariance_CheckStatus = CheckedAndFailed;

	return IsModularInvariant;
}



/* ########################################################################################
######   CreateModelIndependentPart()                                                ######
######                                                                               ######
######   Version: 10.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Sets the twist parts of the private member variables "Elements" and         ######
######   "Centralizer" and computes the orbifold invariant supercharges.             ######
######################################################################################## */
bool COrbifoldGroup::CreateModelIndependentPart()
{
	this->ModelIndependent_CheckStatus  = NotChecked;
	this->ModularInvariance_CheckStatus = NotChecked;
	this->OrbifoldGroup_CheckStatus     = NotChecked;

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// begin: create the constructing elements
	this->Elements.clear();

	const SelfDualLattice Lattice = this->Shifts[0].Lattice;
	COrbifoldGroupElement Element(Lattice);

	unsigned i = 0;
	unsigned j = 0;
	size_t s2 = 0;

	// run through all sectors
	size_t s1 = this->SpaceGroup.GetNumberOfSectors();
	for (i = 0; i < s1; ++i)
	{
		// run through all elements of the current sector
		s2 = this->SpaceGroup.GetNumberOfElementsOfSector(i);
		for (j = 0; j < s2; ++j)
		{
			Element.SGElement = this->SpaceGroup.GetElement(i,j);

			this->GetTwistVector(Element.SGElement, Element.Twist);
			this->Elements.push_back(Element);
		}
	}
	if (this->Elements.size() == 0)
	{
		cout << "\n  Warning in bool COrbifoldGroup::CreateModelIndependentPart(): Constructing elements can not be created. Return false." << endl;
		return false;
	}
	// end: create the constructing elements
	////////////////////////////////////////////////////////////////////////////////////////////////////////

	const bool ZMxZN = this->SpaceGroup.IsZMxZN();
	const bool ZMxZNxZK = this->SpaceGroup.IsZMxZNxZK();

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// begin: find number of supersymmetry
	this->InvariantSupercharges.clear();
	//This needs to be changed to include the twisted supersymmetric E8xE8 !!!!!!
	CTwistVector Supercharge1(-0.5,-0.5,-0.5, 0.5);
	CTwistVector Supercharge2(-0.5,-0.5, 0.5,-0.5);
	CTwistVector Supercharge3(-0.5, 0.5,-0.5,-0.5);
	CTwistVector Supercharge4( 0.5,-0.5,-0.5,-0.5);

	const CTwistVector &ZM_Twist = this->SpaceGroup.GetTwist(0);
	const CTwistVector &ZN_Twist = this->SpaceGroup.GetTwist(1);
	const CTwistVector &ZK_Twist = this->SpaceGroup.GetTwist(2);
	if (ZMxZN)
	{
		if (Supercharge1.IsScalarProductInteger(ZM_Twist) && Supercharge1.IsScalarProductInteger(ZN_Twist))
			this->InvariantSupercharges.push_back(Supercharge1);
		if (Supercharge2.IsScalarProductInteger(ZM_Twist) && Supercharge2.IsScalarProductInteger(ZN_Twist))
			this->InvariantSupercharges.push_back(Supercharge2);
		if (Supercharge3.IsScalarProductInteger(ZM_Twist) && Supercharge3.IsScalarProductInteger(ZN_Twist))
			this->InvariantSupercharges.push_back(Supercharge3);
		if (Supercharge4.IsScalarProductInteger(ZM_Twist) && Supercharge4.IsScalarProductInteger(ZN_Twist))
			this->InvariantSupercharges.push_back(Supercharge4);
	}
	else if (ZMxZNxZK)
	{
		if (Supercharge1.IsScalarProductInteger(ZM_Twist) && Supercharge1.IsScalarProductInteger(ZN_Twist) && Supercharge1.IsScalarProductInteger(ZK_Twist))
			this->InvariantSupercharges.push_back(Supercharge1);
		if (Supercharge2.IsScalarProductInteger(ZM_Twist) && Supercharge2.IsScalarProductInteger(ZN_Twist) && Supercharge2.IsScalarProductInteger(ZK_Twist))
			this->InvariantSupercharges.push_back(Supercharge2);
		if (Supercharge3.IsScalarProductInteger(ZM_Twist) && Supercharge3.IsScalarProductInteger(ZN_Twist) && Supercharge3.IsScalarProductInteger(ZK_Twist))
			this->InvariantSupercharges.push_back(Supercharge3);
		if (Supercharge4.IsScalarProductInteger(ZM_Twist) && Supercharge4.IsScalarProductInteger(ZN_Twist) && Supercharge4.IsScalarProductInteger(ZK_Twist))
			this->InvariantSupercharges.push_back(Supercharge4);
	}
	else
	{
		if (Supercharge1.IsScalarProductInteger(ZM_Twist))
			this->InvariantSupercharges.push_back(Supercharge1);
		if (Supercharge2.IsScalarProductInteger(ZM_Twist))
			this->InvariantSupercharges.push_back(Supercharge2);
		if (Supercharge3.IsScalarProductInteger(ZM_Twist))
			this->InvariantSupercharges.push_back(Supercharge3);
		if (Supercharge4.IsScalarProductInteger(ZM_Twist))
			this->InvariantSupercharges.push_back(Supercharge4);
	}
	this->NumberOfSupersymmetry = this->InvariantSupercharges.size();
	// end: find number of supersymmetry
	////////////////////////////////////////////////////////////////////////////////////////////////////////


	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// begin: create centralizer
	this->Centralizer.clear();
	vector<COrbifoldGroupElement> tmp_Centralizer;

	COrbifoldGroupElement Centralizer_Element(Lattice);

	const vector<vector<CSpaceGroupElement> > &SG_Centralizer = this->SpaceGroup.GetCentralizer();

	s1 = SG_Centralizer.size();
	if (this->Elements.size() != s1)
	{
		cout << "\n  Warning in bool COrbifoldGroup::CreateModelIndependentPart(): Centralizer can not be created. Return false." << endl;
		return false;
	}

	for (i = 0; i < s1; ++i)
	{
		const vector<CSpaceGroupElement> &Current_Centralizer =  SG_Centralizer[i];

		s2 = Current_Centralizer.size();
		if (s2 == 0)
		{
			cout << "\n  Warning in bool COrbifoldGroup::CreateModelIndependentPart(): Centralizer can not be created. Return false." << endl;
			return false;
		}

		tmp_Centralizer.clear();

		for (j = 0; j < s2; ++j)
		{
			Centralizer_Element.SGElement = Current_Centralizer[j];

			this->GetTwistVector(Centralizer_Element.SGElement, Centralizer_Element.Twist);
			tmp_Centralizer.push_back(Centralizer_Element);
		}
		this->Centralizer.push_back(tmp_Centralizer);
	}
	// end: create centralizer

	this->ModelIndependent_CheckStatus = CheckedAndGood;
	return true;
}



/* ########################################################################################
######   CreateModelDependentPart(CPrint &Print, bool PrintWarningModularInvariance) ######
######                                                                               ######
######   Version: 10.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Print                         : a CPrint object                          ######
######   2) PrintWarningModularInvariance : needed for "CheckModularInvariance(...)" ######
######   output:                                                                     ######
######   return value                     : finished successfully?                   ######
###########################################################################################
######   description:                                                                ######
######   Creates the gauge embeddding of the private member variables "Elements" and ######
######   "Centralizer" and checks modular invariance.                                ######
######################################################################################## */
bool COrbifoldGroup::CreateModelDependentPart(CPrint &Print, bool PrintWarningModularInvariance)
{
	if (this->OrbifoldGroup_CheckStatus == CheckedAndGood)
	{
		(*Print.out) << "\n  Warning in bool COrbifoldGroup::CreateModelDependentPart(...) : Orbifold group does not need to be created. Return true." << endl;
		return true;
	}
	this->ModularInvariance_CheckStatus = NotChecked;
	this->OrbifoldGroup_CheckStatus     = NotChecked;

	if (this->ModelIndependent_CheckStatus != CheckedAndGood)
	{
		this->Reset(true, false, false, false);
		if (this->ModelIndependent_CheckStatus != CheckedAndGood)
		{
			(*Print.out) << "\n  Warning in bool COrbifoldGroup::CreateModelDependentPart(...) : Model independent part of the orbifold group ill-defined. Return false." << endl;
			return false;
		}
	}


	// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// begin: check modular invariance
	this->CheckModularInvariance(Print, PrintWarningModularInvariance);
	if (this->ModularInvariance_CheckStatus != CheckedAndGood)
		return false;
	// end: check modular invariance
	// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	size_t s1 = 0;
	size_t s2 = 0;
	unsigned i = 0;
	unsigned j = 0;

	const SelfDualLattice Lattice = this->Shifts[0].Lattice;

	// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// begin: create local shifts of constructing elements
	s1 = this->Elements.size();
	if (s1 == 0)
	{
		(*Print.out) << "\n  Warning in bool COrbifoldGroup::CreateModelDependentPart(...) : Local shifts of the constructing elements can not be created." << endl;
		return false;
	}

	// run through all constructing elements
	for (i = 0; i < s1; ++i)
	{
		COrbifoldGroupElement &Element = this->Elements[i];
		this->GetShiftVector(Element.SGElement, Element.Shift);
		Element.Shift.Lattice = Lattice;
	}
	// end: create local shifts of constructing elements
	// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// begin: create local shifts of centralizer elements
	if (this->Centralizer.size() != s1)
	{
		(*Print.out) << "\n  Warning in bool COrbifoldGroup::CreateModelDependentPart(...) : There is not one centralizer for each constructing element. Return false." << endl;
		return false;
	}

	// run through all centralizer
	for (i = 0; i < s1; ++i)
	{
		vector<COrbifoldGroupElement> &tmp_Centralizer = this->Centralizer[i];
		s2 = tmp_Centralizer.size();
		if (s2 == 0)
		{
			(*Print.out) << "\n  Warning in bool COrbifoldGroup::CreateModelDependentPart(...) : Current centralizer is empty. Return false." << endl;
			return false;
		}

		for (j = 0; j < s2; ++j)
		{
			COrbifoldGroupElement &Centralizer_Element = tmp_Centralizer[j];
			this->GetShiftVector(Centralizer_Element.SGElement, Centralizer_Element.Shift);
			Centralizer_Element.Shift.Lattice = Lattice;
		}
	}
	// end: create local shifts of centralizer elements
	// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	this->OrbifoldGroup_CheckStatus = CheckedAndGood;
	return true;
}



/* ########################################################################################
######   CreateSubGroup(...) const                                                   ######
######                                                                               ######
######   Version: 27.01.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Print               : a CPrint object where possible warnings shall be   ######
######                            printed                                            ######
######   2) TwistFactors1       : sets new order of the Z_N or ZN x ZM sub orbifold  ######
######   3) TwistFactors2       : sets new order of the Z_N or ZN x ZM sub orbifold  ######
######   4) use_WL_in_FixedTori : allow for Wilson lines in fixed tori               ######
######   output:                                                                     ######
######   1) SubOrbifoldGroup    : the sub group of this orbifold group as specified  ######
######                            by "TwistFactors1" and "TwistFactors2"             ######
######   return value           : finished successfully?                             ######
###########################################################################################
######   description:                                                                ######
######   If this orbifold group is Z_N: "TwistFactors1" contains one number k and    ######
######   the sub group has a twist k v_1. "TwistFactors2" contains one number l and  ######
######   the sub group has a second twist l v_2.                                     ######
######   If this orbifold group is Z_N x Z_M: "TwistFactors1" contains two numbers   ######
######   (k_1, k_2) and the sub group has a twist k_1 v_1 + k_2 v_2. "TwistFactors2" ######
######   contains two numbers (l_1, l_2) and the sub group has a second twist        ######
######   l_1 v_1 + l_2 v_2.                                                          ######
######################################################################################## */
bool COrbifoldGroup::CreateSubGroup(CPrint &Print, vector<unsigned> &TwistFactors1, vector<unsigned> &TwistFactors2, bool use_WL_in_FixedTori, COrbifoldGroup &SubOrbifoldGroup) const
{
	if (this->OrbifoldGroup_CheckStatus != CheckedAndGood)
	{
		(*Print.out) << "\n  Warning in bool COrbifoldGroup::CreateSubGroup(...) : Orbifold group ill-defined. Return false." << endl;
		return false;
	}

	unsigned i = 0;

	const SelfDualLattice Lattice = this->Shifts[0].Lattice;

	// begin: create the space group of the sub-orbifold
	CSpaceGroup SubGroup;
	if (!this->SpaceGroup.CreateSubGroup(Print, TwistFactors1, TwistFactors2, use_WL_in_FixedTori, SubGroup))
		return false;

	if (!SubGroup.Check())
	{
		(*Print.out) << "\n  Warning in bool COrbifoldGroup::CreateSubGroup(...) : Space group ill-defined. Return false." << endl;
		return false;
	}

	CShiftVector Shift1(Lattice);
	CShiftVector Shift2(Lattice);

	if (SubGroup.IsZMxZN())
	{
		if (this->SpaceGroup.IsZMxZN())
		{
			Shift1 = (this->Shifts[0] * TwistFactors1[0]) + (this->Shifts[1] * TwistFactors1[1]);
			Shift2 = (this->Shifts[0] * TwistFactors2[0]) + (this->Shifts[1] * TwistFactors2[1]);
		}
		else
		{
			Shift1 = this->Shifts[0] * TwistFactors1[0];
			Shift2 = this->Shifts[0] * TwistFactors2[0];
		}
	}
	else
	{
		if (this->SpaceGroup.IsZMxZN())
			Shift1 = (this->Shifts[0] * TwistFactors1[0]) + (this->Shifts[1] * TwistFactors1[1]);
		else
			Shift1 = this->Shifts[0] * TwistFactors1[0];
	}

	// begin: create the orbifold group SubOrbifoldGroup of the sub-orbifold
	// begin: use only those Wilson lines that are not invariant under the sub orbfifold twists
	CWilsonLines SubOrb_WLs(Lattice);
	const vector<unsigned> &S_WL_AllowedOrders = SubGroup.GetWL_AllowedOrders();
	for (i = 0; i < LatticeDim; ++i)
	{
		if (S_WL_AllowedOrders[i] != 1)
			SubOrb_WLs.SetWilsonLine(i, this->WilsonLines.GetWilsonLine(i));
	}
	SubOrb_WLs.Check(SubGroup.GetWL_Relations(), SubGroup.GetWL_AllowedOrders());
	// end: use only those Wilson lines that are not invariant under the sub orbfifold twists

	vector<CShiftVector> NewShifts;
	NewShifts.push_back(Shift1);
	NewShifts.push_back(Shift2);

	COrbifoldGroup OrbifoldGroup(SubGroup, NewShifts, SubOrb_WLs);
	SubOrbifoldGroup = OrbifoldGroup;
	// end: create the orbifold group SubOrbifoldGroup of the sub-orbifold
	SubOrbifoldGroup.Label = this->Label + "_Sub";

	{
		std::ostringstream os;
		os << TwistFactors1[0];

		if (this->SpaceGroup.IsZMxZN())
			os << TwistFactors1[1];

		if (SubOrbifoldGroup.GetSpaceGroup().IsZMxZN())
		{
			os << "_";
			os << TwistFactors2[0];

			if (this->SpaceGroup.IsZMxZN())
				os << TwistFactors2[1];
		}
		SubOrbifoldGroup.Label += os.str();
	}

	return true;
}



/* ########################################################################################
######   LoadOrbifoldGroup(std::ifstream &in, string &ProgramFilename)               ######
######                                                                               ######
######   Version: 28.02.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) in              : ifstream object containing the model file              ######
######   2) ProgramFilename : optionally, the filename of the program file to be     ######
######                        executed after the model has been loaded               ######
######   output:                                                                     ######
######   return value       : finished successfully?                                 ######
###########################################################################################
######   description:                                                                ######
######   Loads an orbifold group from an ifstream object.                            ######
######################################################################################## */
bool COrbifoldGroup::LoadOrbifoldGroup(std::ifstream &in, string &ProgramFilename)
{
	this->ModelIndependent_CheckStatus  = NotChecked;
	this->ModularInvariance_CheckStatus = NotChecked;
	this->OrbifoldGroup_CheckStatus     = NotChecked;

	this->UseFreelyActingWL = false;
	this->SetDiscreteTorsionToZero();
	this->InvariantSupercharges.clear();
	this->NumberOfSupersymmetry = -1;
	this->LoadedU1Generators.clear();

	this->Label     = "";
	ProgramFilename = "";

	string currentline = "";

	SelfDualLattice Lattice = UNSPECIFIED_LATTICE;

	// begin: find next model
	bool beginmodel_found = false;
	while (GetSaveLine(in, currentline))
	{
		if (currentline == "begin model")
		{
			beginmodel_found = true;
			break;
		}
	}
	if (!beginmodel_found)
		return false;
	// end: find next model

	bool SyntaxError = false;
	vector<bool> NecessaryDataRead(3, false);

	while (GetSaveLine(in, currentline) && (currentline != "end model"))
	{
		if (currentline.substr(0,2) == "//")
		{
			//ignore comments
		}
		else
			if (currentline.substr(0,6) == "Label:")
				this->Label = currentline.substr(6,string::npos);
			else
				if (currentline.substr(0,11) == "SpaceGroup:")
				{
					if (this->SpaceGroup.LoadSpaceGroup(currentline.substr(11,string::npos)))
						NecessaryDataRead[0] = true;
					else
					{
						cout << "\n  Warning in bool COrbifoldGroup::LoadOrbifoldGroup(...) : geometry-file \"" << currentline.substr(11,string::npos) << "\" is corrupt. Return false." << endl;
						return false;
					}
				}
				else
					if (currentline.substr(0,8) == "Lattice:")
					{
						if (currentline.substr(8,string::npos) == "E8xE8")
							Lattice = E8xE8;
						else
							if (currentline.substr(8,string::npos) == "Spin32")
								Lattice = Spin32;
							else
							{
								cout << "\n  Warning in bool COrbifoldGroup::LoadOrbifoldGroup(...) : cannot load the lattice. Return false." << endl;
								return false;
							}
						NecessaryDataRead[1] = true;
						this->WilsonLines.SetLattice(Lattice);
						this->FreelyActingWilsonLine.SetLattice(Lattice);
					}
					else
						if (currentline.substr(0,23) == "Shifts and Wilsonlines:")
						{
							if (Lattice == UNSPECIFIED_LATTICE)
							{
								cout << "\n  Warning in bool COrbifoldGroup::LoadOrbifoldGroup(...) : lattice not defined. Set to E8xE8 and continue." << endl;
								Lattice = E8xE8;
							}

							if (!this->Shifts[0].LoadShiftVector(Lattice, in)
									|| !this->Shifts[1].LoadShiftVector(Lattice, in)
									|| !this->Shifts[2].LoadShiftVector(Lattice, in)		//hacking here!!!
									|| !this->WilsonLines.LoadWilsonLines(Lattice, in))
							{
								cout << "\n  Warning in bool COrbifoldGroup::LoadOrbifoldGroup(...) : cannot load shifts or Wilsonlines. Return false." << endl;
								return false;
							}

							NecessaryDataRead[2] = true;
						}
						else
							if (currentline.substr(0,17) == "Program filename:")
								ProgramFilename = currentline.substr(17,string::npos);
							else
								if (currentline.substr(0,21) == "use freely acting WL:")
								{
									if (Lattice == UNSPECIFIED_LATTICE)
									{
										cout << "\n  Warning in bool COrbifoldGroup::LoadOrbifoldGroup(...) : lattice not defined. Set to E8xE8 and continue." << endl;
										Lattice = E8xE8;
									}

									CShiftVector tmp(Lattice);
									if (!tmp.LoadShiftVector(Lattice, in))
									{
										cout << "\n  Warning in bool COrbifoldGroup::LoadOrbifoldGroup(...) : cannot load freely acting Wilsonline. Return false." << endl;
										return false;
									}
									this->FreelyActingWilsonLine.SetLattice(Lattice);
									this->FreelyActingWilsonLine = tmp;
									this->UseFreelyActingWL = true;
								}
								else
									if (currentline.substr(0,19) == "begin U1 generators")
									{
										if (this->LoadedU1Generators.size() != 0)
										{
											cout << "\n  Warning in bool COrbifoldGroup::LoadOrbifoldGroup(...) : U(1) generators not empty - now cleared." << endl;
											this->LoadedU1Generators.clear();
										}
										if (Lattice == UNSPECIFIED_LATTICE)
										{
											cout << "\n  Warning in bool COrbifoldGroup::LoadOrbifoldGroup(...) : lattice not defined. Set to E8xE8 and continue." << endl;
											Lattice = E8xE8;
										}

										GetSaveLine(in, currentline);

										bool loadU1s = (currentline.substr(0,14) == "number of U1s:");
										if (loadU1s)
										{
											currentline = currentline.substr(14,string::npos);
											loadU1s = (currentline.find_first_not_of("0123456789") == string::npos);
										}

										if (!loadU1s)
										{
											cout << "\n  Warning in bool COrbifoldGroup::LoadOrbifoldGroup(...) : number of U(1) generators not specified. Return false." << endl;
											return false;
										}
										const unsigned number_of_U1s = (unsigned)atoi(currentline.c_str());

										CShiftVector tmp(Lattice);
										for (unsigned i = 0; i < number_of_U1s; ++i)
										{
											if (!tmp.LoadShiftVector(Lattice, in))
											{
												cout << "\n  Warning in bool COrbifoldGroup::LoadOrbifoldGroup(...) : cannot load U(1) generator. Return false." << endl;
												return false;
											}
											this->LoadedU1Generators.push_back(tmp);
										}
									}
									else
										if (currentline.substr(0,16) == "discrete torsion")
										{
											SyntaxError = false;
											vector<rational<int> > gdt;

											GDT_Parameters DiscreteTorsionFromFile;

											// load a
											GetSaveLine(in, currentline);
											convert_string_to_vector_of_rational(currentline, gdt);
											if (gdt.size() != 1)
												SyntaxError = true;
											else
												DiscreteTorsionFromFile.a = gdt[0];

											// load b
											GetSaveLine(in, currentline);
											convert_string_to_vector_of_rational(currentline, gdt);
											if (gdt.size() != LatticeDim)
												SyntaxError = true;
											else
												DiscreteTorsionFromFile.b = gdt;

											// load c
											GetSaveLine(in, currentline);
											convert_string_to_vector_of_rational(currentline, gdt);
											if (gdt.size() != LatticeDim)
												SyntaxError = true;
											else
												DiscreteTorsionFromFile.c = gdt;

											// load d
											unsigned i = 0;
											for (i = 0; !SyntaxError && (i < LatticeDim); ++i)
											{
												GetSaveLine(in, currentline);
												convert_string_to_vector_of_rational(currentline, gdt);
												if (gdt.size() != LatticeDim)
													SyntaxError = true;
												else
													DiscreteTorsionFromFile.d[i] = gdt;
											}

											if (SyntaxError || !this->SetDiscreteTorsion(DiscreteTorsionFromFile))
											{
												cout << "\n  Warning in bool COrbifoldGroup::LoadOrbifoldGroup(...) : cannot load generalized discrete torsion. Set it to zero and continue." << endl;
												this->SetDiscreteTorsionToZero();
											}
											else
												this->DiscreteTorsion.UseGDT = true;
										}
		//else
	}

	if (find(NecessaryDataRead.begin(), NecessaryDataRead.end(), false) != NecessaryDataRead.end())
	{
		vector<string> ErrorMessages(3,"");
		ErrorMessages[0] = "space group is";
		ErrorMessages[1] = "the 16-dim. lattice is";
		ErrorMessages[2] = "shifts and Wilsonlines are";

		for (unsigned i = 0; i < 3; ++i)
		{
			if (!NecessaryDataRead[i])
				cout << "\n  Cannot load model from file:  " << ErrorMessages[i] << " not specified." << endl;
		}
		return false;
	}

	if (!this->SpaceGroup.Check())
	{
		cout << "\n  Cannot load model from file:  Space group ill-defined." << endl;
		return false;
	}

	if (!this->WilsonLines.Check(this->SpaceGroup.GetWL_Relations(), this->SpaceGroup.GetWL_AllowedOrders()))
	{
		cout << "\n  Cannot load model from file: Relations between Wilson lines not fulfilled." << endl;
		return false;
	}


	this->CreateModelIndependentPart();
	if (this->ModelIndependent_CheckStatus != CheckedAndGood)
	{
		cout << "\n  Warning in bool COrbifoldGroup::LoadOrbifoldGroup(...): Model independent part of the orbifold group can not be created." << endl;
		return false;
	}

	CPrint Print(Tstandard, &cout);
	this->CreateModelDependentPart(Print, true);
	if (this->OrbifoldGroup_CheckStatus != CheckedAndGood)
	{
		cout << "\n  Warning in bool COrbifoldGroup::LoadOrbifoldGroup(...): Orbifold group ill-defined." << endl;
		return false;
	}

	this->OrbifoldGroup_CheckStatus = CheckedAndGood;

	return true;
}



/* ########################################################################################
######   Reset(bool ResetSpaceGroup)                                                 ######
######                                                                               ######
######   Version: 10.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) CreateModelIndependentPart : use "true" if the space group has changed   ######
######   2) ResetShifts                : set the shifts to zero                      ######
######   3) ResetWilsonLines           : set the Wilson lines to zero                ######
######   4) ResetTorsion               : set discrete torsion to zero                ######
######   output:                                                                     ######
######   return value                  : finished successfully?                      ######
###########################################################################################
######   description:                                                                ######
######   Resets the orbifold group.                                                  ######
######################################################################################## */
bool COrbifoldGroup::Reset(bool CreateModelIndependentPart, bool ResetShifts, bool ResetWilsonLines, bool ResetTorsion)
{
	// keep the Label and the SpaceGroup

	// use old lattice
	SelfDualLattice Lattice = E8xE8;
	if (this->Shifts[0].Lattice != UNSPECIFIED_LATTICE)
		Lattice = this->Shifts[0].Lattice;

	if (ResetShifts)
	{
		CShiftVector NullVector(Lattice);
		this->Shifts[0] = NullVector;
		this->Shifts[1] = NullVector;
	}

	if (ResetWilsonLines)
	{
		CWilsonLine null(Lattice);
		this->UseFreelyActingWL = false;
		this->FreelyActingWilsonLine = null;
		this->WilsonLines.SetToZero(Lattice);
	}

	if (ResetTorsion)
		this->SetDiscreteTorsionToZero();

	this->LoadedU1Generators.clear();

	if (ResetShifts || ResetWilsonLines)
		this->ModularInvariance_CheckStatus = NotChecked;

	this->OrbifoldGroup_CheckStatus = NotChecked;

	if (CreateModelIndependentPart)
	{
		this->ModelIndependent_CheckStatus = NotChecked;

		this->Elements.clear();
		this->Centralizer.clear();

		this->NumberOfSupersymmetry = -1;
		this->InvariantSupercharges.clear();

		this->CreateModelIndependentPart();
		if (this->ModelIndependent_CheckStatus != CheckedAndGood)
			return false;
	}

	return true;
}



/* ########################################################################################
######   LocalGaugeGroupOnFixedTorus(std::ostream &out)                              ######
######                                                                               ######
######   Version: 10.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) out : print the output to "out"                                          ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Identifies fixed tori in the various twisted sectors and computes the gauge ######
######   group in the 6D orbifold GUT limit.                                         ######
######################################################################################## */
void COrbifoldGroup::LocalGaugeGroupOnFixedTorus(std::ostream &out) const
{
	// "Is_Ti_Fixed" contains the information whether the i-th torus is fixed or not
	vector<bool>   Is_Ti_Fixed(3, false);

	// "Ti_GaugeGroup" contains the gauge group of the i-th torus as a string
	vector<string> Ti_GaugeGroup;

	// Get the Wilsonlines
	const CWilsonLine  &W1 = this->WilsonLines.GetWilsonLine(0);
	const CWilsonLine  &W2 = this->WilsonLines.GetWilsonLine(1);
	const CWilsonLine  &W3 = this->WilsonLines.GetWilsonLine(2);
	const CWilsonLine  &W4 = this->WilsonLines.GetWilsonLine(3);
	const CWilsonLine  &W5 = this->WilsonLines.GetWilsonLine(4);
	const CWilsonLine  &W6 = this->WilsonLines.GetWilsonLine(5);

	const SelfDualLattice Lattice = this->Shifts[0].Lattice;

	// begin: create the E8 x E8' or SO(32) gauge group
	const CVector Null_Shift(16);
	S_OscillatorExcitation Excitation;
	Excitation.NumberOperator = 0.0;
	Excitation.ZeroPointEnergy = 2.0;
	CMasslessHalfState GaugeGroup_10D(LeftMover, Excitation);
	GaugeGroup_10D.SolveMassEquation(Null_Shift, Lattice);

	if (GaugeGroup_10D.Weights.size() != 480)
	{
		cout << "\n  Warning in void COrbifoldGroup::LocalGaugeGroupOnFixedTorus(...) const : check the 10D gauge group. Return." << endl;
		return;
	}
	// end: create the E8 x E8' or SO(32) gauge group

	const size_t s1 = this->Elements.size();

	// i = 0 corresponds to the untwisted sector
	for (unsigned i = 1; i < s1; ++i)
	{
		const COrbifoldGroupElement &Element = this->Elements[i];
		const CTwistVector          &Twist   = Element.Twist;

		// find the fixed branes
		for (unsigned j = 0; j < 3; ++j)
		{
			if ((!Is_Ti_Fixed.at(j)) && is_integer(Twist.at(j+1)))
			{
				Is_Ti_Fixed.at(j) = true;

				// begin: create the shift vector of the (m,n,k) sector
				CSpaceGroupElement Label(Element.SGElement.Get_m(), Element.SGElement.Get_n(), Element.SGElement.Get_k());

				CShiftVector Shift(Lattice);
				this->GetShiftVector(Label, Shift);
				// end: create the shift vector of the (k,l) sector

				// calculate the local gauge group on the fixed brane
				vector<CVector> Six_D_Bulk_Gauge_Group;
				for (unsigned k = 0; k < 480; ++k)
				{
					const CVector &tmp = GaugeGroup_10D.Weights[k];

					double sp1 = tmp * Shift;
					double sp2 = 0;
					double sp3 = 0;
					double sp4 = 0;
					double sp5 = 0;

					switch (j)
					{
					case 0:
					{
						if (!W3.GetIs_Zero()) sp2 = tmp * W3;
						if (!W4.GetIs_Zero()) sp3 = tmp * W4;
						if (!W5.GetIs_Zero()) sp4 = tmp * W5;
						if (!W6.GetIs_Zero()) sp5 = tmp * W6;
						break;
					}
					case 1:
					{
						if (!W1.GetIs_Zero()) sp2 = tmp * W1;
						if (!W2.GetIs_Zero()) sp3 = tmp * W2;
						if (!W5.GetIs_Zero()) sp4 = tmp * W5;
						if (!W6.GetIs_Zero()) sp5 = tmp * W6;
						break;
					}
					case 2:
					{
						if (!W1.GetIs_Zero()) sp2 = tmp * W1;
						if (!W2.GetIs_Zero()) sp3 = tmp * W2;
						if (!W3.GetIs_Zero()) sp4 = tmp * W3;
						if (!W4.GetIs_Zero()) sp5 = tmp * W4;
						break;
					}
					}

					if (is_integer(sp1) && is_integer(sp2) && is_integer(sp3) && is_integer(sp4) && is_integer(sp5))
						Six_D_Bulk_Gauge_Group.push_back(tmp);
				}

				vector<vector<double> > GaugeGroup2_p1;
				vector<vector<double> >  GaugeGroup2_p2;

				const size_t s3 = Six_D_Bulk_Gauge_Group.size();

				// E8 x E8'
				if (Lattice == E8xE8)
				{
					for (unsigned k = 0; k < s3; ++k)
					{
						vector<double> Root_p1;
						vector<double> Root_p2;

						const CVector &tmp = Six_D_Bulk_Gauge_Group[k];
						for (unsigned l = 0; l < 8; ++l)
						{
							Root_p1.push_back(tmp[l]);
							Root_p2.push_back(tmp[l + 8]);
						}
						GaugeGroup2_p1.push_back(Root_p1);
						GaugeGroup2_p2.push_back(Root_p2);
					}
				}
				// Spin32
				else
				{
					for (unsigned k = 0; k < s3; ++k)
					{
						vector<double> Root;

						const CVector &tmp = Six_D_Bulk_Gauge_Group[k];
						for (unsigned l = 0; l < 16; ++l)
							Root.push_back(tmp[l]);

						GaugeGroup2_p1.push_back(Root);
					}
				}

				const gaugeGroup<double> Group2_p1 = determineAlgebra(GaugeGroup2_p1);

				std::ostringstream tmp;
				tmp << " " << j+1 << "-th torus: " << Group2_p1.algebra;

				if (Lattice == E8xE8)
				{
					const gaugeGroup<double> Group2_p2 = determineAlgebra(GaugeGroup2_p2);
					tmp << " and " << Group2_p2.algebra;
				}

				Ti_GaugeGroup.push_back(tmp.str());
			}
		}
	}

	// print
	const size_t s3 = Ti_GaugeGroup.size();
	if (s3 == 0)
		out << "This Orbifold has no fixed tori.\n";
	else
	{
		out << "Local gauge group on i-th fixed torus:\n";
		for (unsigned i = 0; i < s3; ++i)
			out << Ti_GaugeGroup[i] << "\n";
	}
	out << endl;
}



/* ########################################################################################
######   void COrbifoldGroup::PrintToFile(...) const                                 ######
######                                                                               ######
######   Version: 04.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) out : save the model file of this orbifold group to "out"                ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Saves this orbifold group in a model file specified by "out".               ######
######################################################################################## */
void COrbifoldGroup::PrintToFile(std::ostream &out) const
{
	const SelfDualLattice Lattice = this->Shifts[0].Lattice;

	CPrint Print(Tstandard, &out);
	out << "begin model\n";
	out << "Label:" << this->Label << "\n"; 
	out << "SpaceGroup:" << this->SpaceGroup.GeometryFilename << "\n";

	if (Lattice == E8xE8)
		out << "Lattice:E8xE8\n";
	else
		out << "Lattice:Spin32\n";

	if (this->UseFreelyActingWL)
	{
		out << "use freely acting WL:\n";
		Print.PrintRational(this->FreelyActingWilsonLine);
		out << "\n";
	}

	out << "Shifts and Wilsonlines:\n";
	unsigned i = 0;
	for (i = 0; i < 3; ++i)						//hacking here!!!
	{
		Print.PrintRational(this->Shifts[i]);
		out << "\n";
	}

	for (i = 0; i < LatticeDim; ++i)
	{
		Print.PrintRational(this->WilsonLines.GetWilsonLine(i));
		out << "\n";
	}

	if (this->DiscreteTorsion.UseGDT)
	{
		bool print = true;
		if ((this->DiscreteTorsion.b.size() != LatticeDim) || (this->DiscreteTorsion.c.size() != LatticeDim) || (this->DiscreteTorsion.d.size() != LatticeDim))
			print = false;

		for (i = 0; print && (i < LatticeDim); ++i)
		{
			if (this->DiscreteTorsion.d[i].size() != LatticeDim)
				print = false;
		}

		if (!print)
		{
			cout << "Warning in void COrbifoldGroup::PrintToFile(...) const : cannot print discrete torsion to file." << endl;
			out << "end model" << endl;
			return;
		}

		out << "discrete torsion\n";
		Print.PrintRational(this->DiscreteTorsion.a);
		out << "\n";

		for (i = 0; i < LatticeDim; ++i)
		{
			Print.PrintRational(this->DiscreteTorsion.b[i]);
			out << " ";
		}
		out << "\n";

		for (i = 0; i < LatticeDim; ++i)
		{
			Print.PrintRational(this->DiscreteTorsion.c[i]);
			out << " ";
		}
		out << "\n";

		unsigned j = 0;
		for (i = 0; i < LatticeDim; ++i)
		{
			const rationalVector &d_i = this->DiscreteTorsion.d[i];
			for (j = 0; j < LatticeDim; ++j)
			{
				Print.PrintRational(d_i[j]);
				out << " ";
			}
			out << "\n";
		}
	}
	out << "end model" << endl;
	out << endl;
}



/* ########################################################################################
######   SetDiscreteTorsion(const GDT_Parameters &DiscreteTorsion)                   ######
######                                                                               ######
######   Version: 04.09.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) DiscreteTorsion : a GDT_Parameters object                                ######
######   output:                                                                     ######
######   return value       : finished successfully?                                 ######
###########################################################################################
######   description:                                                                ######
######   Sets the private member variable "DiscreteTorsion".                         ######
######################################################################################## */
bool COrbifoldGroup::SetDiscreteTorsion(const GDT_Parameters &DiscreteTorsion)
{
	if ((DiscreteTorsion.b.size() != LatticeDim) || (DiscreteTorsion.c.size() != LatticeDim) || (DiscreteTorsion.d.size() != LatticeDim))
	{
		cout << "\n  Warning in bool COrbifoldGroup::SetDiscreteTorsion: Discrete torsion (b,c or d) ill-defined. Return false." << endl;
		return false;
	}

	unsigned i = 0;
	unsigned j = 0;
	for (i = 0; i < LatticeDim; ++i)
	{
		if (DiscreteTorsion.d[i].size() != LatticeDim)
		{
			cout << "\n  Warning in bool COrbifoldGroup::SetDiscreteTorsion: Discrete torsion (d) ill-defined. Return false." << endl;
			return false;
		}

		for (j = 0; j < LatticeDim; ++j)
		{
			if (DiscreteTorsion.d[i][j] != -DiscreteTorsion.d[j][i])
			{
				cout << "\n  Warning in bool COrbifoldGroup::SetDiscreteTorsion: Discrete torsion phase d_ab not antisymmetric. Return false." << endl;
				return false;
			}
		}
	}
	this->DiscreteTorsion = DiscreteTorsion;
	this->DiscreteTorsion.UseGDT = false;

	if (DiscreteTorsion.a != 0)
	{
		this->DiscreteTorsion.UseGDT = true;
		return true;
	}

	for (i = 0; i < LatticeDim; ++i)
	{
		if (DiscreteTorsion.b[i] != 0)
		{
			this->DiscreteTorsion.UseGDT = true;
			return true;
		}
	}

	for (i = 0; i < LatticeDim; ++i)
	{
		if (DiscreteTorsion.c[i] != 0)
		{
			this->DiscreteTorsion.UseGDT = true;
			return true;
		}
	}

	for (i = 0; i < LatticeDim; ++i)
	{
		for (j = 0; j < LatticeDim; ++j)
		{
			if (DiscreteTorsion.d[i][j] != 0)
			{
				this->DiscreteTorsion.UseGDT = true;
				return true;
			}
		}
	}

	return true;
}



/* ########################################################################################
######   SetDiscreteTorsionToZero()                                                  ######
######                                                                               ######
######   Version: 04.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Sets discrete torsion phases to zero.                                       ######
######################################################################################## */
void COrbifoldGroup::SetDiscreteTorsionToZero()
{
	this->DiscreteTorsion.UseGDT = false;
	this->DiscreteTorsion.a = 0;
	this->DiscreteTorsion.b.assign(LatticeDim,0);
	this->DiscreteTorsion.c.assign(LatticeDim,0);

	rationalVector d_line(LatticeDim,0);
	this->DiscreteTorsion.d.assign(LatticeDim, d_line);
}



/* ########################################################################################
######   AccessShift(const unsigned &i)                                              ######
######                                                                               ######
######   Version: 10.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) i         : index i = 0,1                                                ######
######   output:                                                                     ######
######   return value : either the shift V_1 or V_2                                  ######
###########################################################################################
######   description:                                                                ######
######   Get access to the shifts: V_1 for i = 0 and V_2 for i = 1.                  ######
######################################################################################## */
CShiftVector &COrbifoldGroup::AccessShift(const unsigned &i)
{
	this->ModularInvariance_CheckStatus = NotChecked;
	this->OrbifoldGroup_CheckStatus     = NotChecked;

	if (i > 2)  
	{
		cout << "\n  Warning in CShiftVector &COrbifoldGroup::AccessShift(...) : Index i out of range. Set i = 0." << endl;
		return this->Shifts[0];
	}
	return this->Shifts[i];
}



/* ########################################################################################
######   AccessSpaceGroup()                                                          ######
######                                                                               ######
######   Version: 10.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : the CSpaceGroup object of this orbifold group                ######
###########################################################################################
######   description:                                                                ######
######   Get access to the space group.                                              ######
######################################################################################## */
CSpaceGroup &COrbifoldGroup::AccessSpaceGroup()
{
	this->SpaceGroup.SG_NotChecked();

	this->ModelIndependent_CheckStatus  = NotChecked;
	this->ModularInvariance_CheckStatus = NotChecked;
	this->OrbifoldGroup_CheckStatus     = NotChecked;

	return this->SpaceGroup;
}



/* ########################################################################################
######   AccessWilsonLines()                                                         ######
######                                                                               ######
######   Version: 10.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : the CWilsonLines object of this orbifold group               ######
###########################################################################################
######   description:                                                                ######
######   Get access to the Wilson lines.                                             ######
######################################################################################## */
CWilsonLines &COrbifoldGroup::AccessWilsonLines()
{
	this->WilsonLines.WL_NotChecked();

	this->ModularInvariance_CheckStatus = NotChecked;
	this->OrbifoldGroup_CheckStatus     = NotChecked;

	return this->WilsonLines;
}



/* ########################################################################################
######   GetElement(const unsigned &i) const                                         ######
######                                                                               ######
######   Version: 10.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) i         : index i from 0,..., number of elements                       ######
######   output:                                                                     ######
######   return value : the i-th COrbifoldGroupElement object                        ######
###########################################################################################
######   description:                                                                ######
######   Get (constant) access to the "i"-th constructing element.                   ######
######################################################################################## */
const COrbifoldGroupElement &COrbifoldGroup::GetElement(const unsigned &i) const
{
	if (this->OrbifoldGroup_CheckStatus != CheckedAndGood)
	{
		cout << "\n  Warning in const COrbifoldGroupElement &COrbifoldGroup::GetElement(...) const : Orbifold group ill-defined." << endl;
		return this->Elements[0];
	}
	if (i >= this->Elements.size())
	{
		cout << "\n  Warning in const COrbifoldGroupElement &COrbifoldGroup::GetElement(...) const : Index i out of range. Set i = 0." << endl;
		return this->Elements[0];
	}
	return this->Elements[i];
}



/* ########################################################################################
######   GetElements() const                                                         ######
######                                                                               ######
######   Version: 10.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : all constructing elements                                    ######
###########################################################################################
######   description:                                                                ######
######   Get (constant) access to all constructing element.                          ######
######################################################################################## */
const vector<COrbifoldGroupElement> &COrbifoldGroup::GetElements() const
{
	if (this->OrbifoldGroup_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const vector<COrbifoldGroupElement> &COrbifoldGroup::GetElements() const : Orbifold group ill-defined." << endl;
	return this->Elements;
}



/* ########################################################################################
######   GetCentralizer(const unsigned &i) const                                     ######
######                                                                               ######
######   Version: 12.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : all centralizer elements                                     ######
###########################################################################################
######   description:                                                                ######
######   Get (constant) access to all centralizer element.                           ######
######################################################################################## */
const vector<COrbifoldGroupElement> &COrbifoldGroup::GetCentralizer(const unsigned &i) const
{
	if (this->OrbifoldGroup_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const vector<COrbifoldGroupElement> &COrbifoldGroup::GetCentralizer(...) const : Orbifold group ill-defined." << endl;

	if (i >= this->Centralizer.size())
	{
		cout << "\n  Warning in const vector<COrbifoldGroupElement> &COrbifoldGroup::GetCentralizer(...) const : Index i out of range. Set i = 0." << endl;
		return this->Centralizer[0];
	}
	return this->Centralizer[i];
}



/* ########################################################################################
######   GetShift(const unsigned &i) const                                           ######
######                                                                               ######
######   Version: 12.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) i         : index i = 0,1                                                ######
######   output:                                                                     ######
######   return value : either the shift V_1 or V_2                                  ######
###########################################################################################
######   description:                                                                ######
######   Get (constant) access to the shifts: V_0 for i = 0, V_1 for i = 1 and V_2 for i = 2      ######
######################################################################################## */
const CShiftVector &COrbifoldGroup::GetShift(const unsigned &i) const
{
	if (this->OrbifoldGroup_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const CShiftVector &COrbifoldGroup::GetShift(...) const : Orbifold group ill-defined." << endl;

	if (i > 2)
	{
		cout << "\n  Warning in const CShiftVector &COrbifoldGroup::GetShift(...) const : Index i out of range. Set i = 0." << endl;
		return this->Shifts[0];
	}
	return this->Shifts[i];
}



/* ########################################################################################
######   GetShifts() const                                                           ######
######                                                                               ######
######   Version: 10.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : the shifts V_1 and V_2                                       ######
###########################################################################################
######   description:                                                                ######
######   Get (constant) access to the shifts.                                        ######
######################################################################################## */
const vector<CShiftVector> &COrbifoldGroup::GetShifts() const
{
	if (this->OrbifoldGroup_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const vector<CShiftVector> &COrbifoldGroup::GetShifts() const : Orbifold group ill-defined." << endl;
	return this->Shifts;
}



/* ########################################################################################
######   GetWilsonLines() const                                                      ######
######                                                                               ######
######   Version: 10.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : the set of six Wilson lines stored in a CWilsonLines object  ######
###########################################################################################
######   description:                                                                ######
######   Get (constant) access to the Wilson lines.                                  ######
######################################################################################## */
const CWilsonLines &COrbifoldGroup::GetWilsonLines() const
{
	if (this->OrbifoldGroup_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const CWilsonLines &COrbifoldGroup::GetWilsonLines() const : Orbifold group ill-defined." << endl;
	return this->WilsonLines;
}



/* ########################################################################################
######   GetDiscreteTorsion() const                                                  ######
######                                                                               ######
######   Version: 10.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : the discrete torsion phases                                  ######
###########################################################################################
######   description:                                                                ######
######   Get (constant) access to the discrete torsion phases.                       ######
######################################################################################## */
const GDT_Parameters &COrbifoldGroup::GetDiscreteTorsion() const
{
	if (this->OrbifoldGroup_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const GDT_Parameters &COrbifoldGroup::GetDiscreteTorsion() const : Orbifold group ill-defined." << endl;
	return this->DiscreteTorsion;
}



/* ########################################################################################
######   GetSpaceGroup() const                                                       ######
######                                                                               ######
######   Version: 10.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : the space group                                              ######
###########################################################################################
######   description:                                                                ######
######   Get (constant) access to the space group.                                   ######
######################################################################################## */
const CSpaceGroup &COrbifoldGroup::GetSpaceGroup() const
{
	if (this->SpaceGroup.GetCheckStatus() != CheckedAndGood)
		cout << "\n  Warning in const CSpaceGroup &COrbifoldGroup::GetSpaceGroup() const : Space group ill-defined." << endl;
	return this->SpaceGroup;
}



/* ########################################################################################
######   GetInvariantSupercharges() const                                            ######
######                                                                               ######
######   Version: 04.09.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : the set of orbifold invariant supercharges                   ######
###########################################################################################
######   description:                                                                ######
######   Get (constant) access to the set of orbifold invariant supercharges.        ######
######################################################################################## */
const vector<CVector> &COrbifoldGroup::GetInvariantSupercharges() const
{
	if (this->ModelIndependent_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const vector<CVector> &COrbifoldGroup::GetInvariantSupercharges() const : Orbifold group ill-defined." << endl;
	return this->InvariantSupercharges;
}



/* ########################################################################################
######   GetNumberOfSupersymmetry() const                                            ######
######                                                                               ######
######   Version: 10.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : the number of supersymmetries                                ######
###########################################################################################
######   description:                                                                ######
######   Get (constant) access to the number of supersymmetries.                     ######
######################################################################################## */
const int &COrbifoldGroup::GetNumberOfSupersymmetry() const
{
	if (this->OrbifoldGroup_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const int &COrbifoldGroup::GetNumberOfSupersymmetry() const : Orbifold group ill-defined." << endl;
	return this->NumberOfSupersymmetry;
}

