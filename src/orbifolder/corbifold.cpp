#include "corbifold.h"
#include "globalfunctions.h"
#include "clinalg.h"
#include "cprint.h"
#include "chugeint.h"
#include "clatticevector.h"
#include "cprompt.h"

#include <iostream>
#include <sstream>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <sys/stat.h>

#define CHECKERROR true

extern unsigned SELFDUALLATTICE;

using namespace std;



/* ########################################################################################
######   COrbifold()                                                                 ######
######                                                                               ######
######   Version: 05.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard constructor of a COrbifold object. No content is specified.        ######
######################################################################################## */
COrbifold::COrbifold()
{
	this->Orbifold_CheckStatus = NotChecked;

	string tmp = "StandardConfig";
	this->Config_Clear(this->StandardConfig, tmp);
}



/* ########################################################################################
######   COrbifold(const COrbifoldGroup &OrbifoldGroup)                              ######
######                                                                               ######
######   Version: 19.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) OrbifoldGroup : a COrbifoldGroup object specifying the space group and   ######
######                      its gauge embedding                                      ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Constructor of a COrbifold object. The complete massless orbifold spectrum  ######
######   is computed and stored in the standard vev-config "StandardConfig".         ######
######################################################################################## */
COrbifold::COrbifold(const COrbifoldGroup &OrbifoldGroup)
: OrbifoldGroup(OrbifoldGroup)
{
	this->Orbifold_CheckStatus = NotChecked;

	string tmp = "StandardConfig";
	this->Config_Clear(this->StandardConfig, tmp);
	this->StandardConfig.InvariantSupercharges = this->OrbifoldGroup.GetInvariantSupercharges();

	COrbifoldCore OrbifoldCore(this->OrbifoldGroup, this->Sectors);
	this->Create();
}



/* ########################################################################################
######   COrbifold(const COrbifoldGroup &OrbifoldGroup, ...)                         ######
######                                                                               ######
######   Version: 19.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) OrbifoldGroup : a COrbifoldGroup object specifying the space group and   ######
######                      its gauge embedding                                      ######
######   2) CoreSpectrum  : created using the class "COrbifoldCore"                  ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Constructor of a COrbifold object. "CoreSpectrum" contains the part of the  ######
######   spectrum that does not depend on the gauge embedding, i.e. the right-moving ######
######   part. Using this, the complete massless orbifold spectrum is computed and   ######
######   stored in the standard vev-config "StandardConfig".                         ######
######################################################################################## */
COrbifold::COrbifold(const COrbifoldGroup &OrbifoldGroup, const vector<CSector> &CoreSpectrum)
: OrbifoldGroup(OrbifoldGroup)
{
	this->Orbifold_CheckStatus = NotChecked;

	string tmp = "StandardConfig";
	this->Config_Clear(this->StandardConfig, tmp);
	this->StandardConfig.InvariantSupercharges = this->OrbifoldGroup.GetInvariantSupercharges();

	this->Sectors = CoreSpectrum;
	this->Create();
}



/* ########################################################################################
######   ~COrbifold()                                                                ######
######                                                                               ######
######   Version: 31.03.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard destructor of a COrbifold object.                                  ######
######################################################################################## */
COrbifold::~COrbifold()
{
}



/* ########################################################################################
######   Reset(bool CreateModelIndependentPart, ...)                                 ######
######                                                                               ######
######   Version: 08.09.2011                                                         ######
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
######   Resets the orbifold.                                                        ######
######################################################################################## */
bool COrbifold::Reset(bool CreateModelIndependentPart, bool ResetShifts, bool ResetWilsonLines, bool ResetTorsion)
{
	string tmp = "StandardConfig";
	this->Config_Clear(this->StandardConfig, tmp);
	if (CreateModelIndependentPart)
	{
		this->OrbifoldGroup.Reset(true, ResetShifts, ResetWilsonLines, ResetTorsion);

		if (this->OrbifoldGroup.GetModelIndependent_CheckStatus() == CheckedAndGood)
			this->StandardConfig.InvariantSupercharges = this->OrbifoldGroup.GetInvariantSupercharges();
	}
	else
	{
		if (this->OrbifoldGroup.GetModelIndependent_CheckStatus() == CheckedAndGood)
			this->StandardConfig.InvariantSupercharges = this->OrbifoldGroup.GetInvariantSupercharges();

		this->OrbifoldGroup.Reset(false, ResetShifts, ResetWilsonLines, ResetTorsion);
	}

	size_t s2 = 0;
	unsigned j = 0;

	const size_t s1 = this->Sectors.size();
	for (unsigned i = 0; i < s1; ++i)
	{
		CSector &Sector = this->Sectors[i];

		s2 = Sector.GetNumberOfFixedBranes();
		for (j = 0; j < s2; ++j)
			Sector.AccessFixedBrane(j).Reset();
	}

	this->Orbifold_CheckStatus = NotChecked;

	return true;
}



/* ########################################################################################
######   Create()                                                                    ######
######                                                                               ######
######   Version: 28.09.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Creates the model-dependent part of the orbifold, i.e. the left-movers, and ######
######   the complete spectrum of massless strings.                                  ######
######################################################################################## */
bool COrbifold::Create()
{
	if (this->OrbifoldGroup.GetModularInvariance_CheckStatus() != CheckedAndGood)
	{
		cout << "\n  Warning in bool COrbifold::Create() : Modular invariance failed. Return false." << endl;
		return false;
	}

	size_t number_of_sectors = this->Sectors.size();
	if (number_of_sectors == 0)
	{
		cout << "\n  Warning in bool COrbifold::Create() : \"Sectors\" is empty. Calling COrbifoldCore." << endl;

		COrbifoldCore OrbifoldCore(this->OrbifoldGroup, this->Sectors);
		number_of_sectors = this->Sectors.size();
		if (number_of_sectors == 0)
		{
			cout << "\n  Warning in bool COrbifold::Create() : \"Sectors\" is still empty. Return false." << endl;
			return false;
		}
	}

	unsigned i = 0;
	unsigned j = 0;

	size_t s1 = 0;
	size_t s2 = 0;

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// begin: define the gamma elements
	vector<COrbifoldGroupElement> Gamma_Centralizer;

	const vector<CSpaceGroupElement> &SG_Generators_Twist = this->OrbifoldGroup.GetSpaceGroup().GetSG_Generators_Twist();
	const vector<CSpaceGroupElement> &SG_Generators_Shift = this->OrbifoldGroup.GetSpaceGroup().GetSG_Generators_Shift();

	const SelfDualLattice Lattice = this->OrbifoldGroup.GetLattice();
	CShiftVector Shift(Lattice);
	CTwistVector Twist;

	s1 = SG_Generators_Twist.size();
	for (i = 0; i < s1; ++i)
	{
		const CSpaceGroupElement &SGElement = SG_Generators_Twist[i];

		this->OrbifoldGroup.GetShiftVector(SGElement, Shift);
		this->OrbifoldGroup.GetTwistVector(SGElement, Twist);

		COrbifoldGroupElement new_Element(SGElement, Shift, Twist);
		Gamma_Centralizer.push_back(new_Element);

	}
	s1 = SG_Generators_Shift.size();
	for (i = 0; i < s1; ++i)
	{
		const CSpaceGroupElement &SGElement = SG_Generators_Shift[i];

		this->OrbifoldGroup.GetShiftVector(SGElement, Shift);
		this->OrbifoldGroup.GetTwistVector(SGElement, Twist);

		COrbifoldGroupElement new_Element(SGElement, Shift, Twist);
		Gamma_Centralizer.push_back(new_Element);
	}
	// end: define the gamma elements
	////////////////////////////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// begin: create the left-movers and the states on the fixed branes
	for (i = 0; i < number_of_sectors; ++i)
	{
		CSector &Sector = this->Sectors[i];

		if (Sector.GetRightMovers().size() == 0)
		{
			cout << "\n  Warning in bool COrbifold::Create() : Right-movers of the " << i << "-th sector have not been created. Return false." << endl;
			return false;
		}

		s2 = Sector.GetNumberOfFixedBranes();
		if (s2 == 0)
		{
			cout << "\n  Warning in bool COrbifold::Create() : " << i << "-th sector has no fixed branes. Return false." << endl;
			return false;
		}

		for (j = 0; j < s2; ++j) {
			if(Sector.AccessFixedBrane(j).Create(this->OrbifoldGroup, Sector, Gamma_Centralizer) != true) {						//hacking here!!!
				//cout << "\n  Warning in bool COrbifold::Create() : " << i << "-th sector has problem in the " << j <<  "-th fixed brane. Return false." << endl;
				this->Orbifold_CheckStatus = CheckedAndFailed;
				return false;
			}
		}
	}
	// end: create the left-movers and the states on the fixed branes
	////////////////////////////////////////////////////////////////////////////////////////////////////////

	// number_of_sectors != 0 check before
	if (this->Sectors[0].GetNumberOfFixedBranes() == 0)
	{
		cout << "\n  Warning in bool COrbifold::Create() : Untwisted sector is empty. Return false." << endl;
		return false;
	}

	if (Lattice == E8xE8)
		SELFDUALLATTICE = 1;
	else
		if (Lattice == Spin32)
			SELFDUALLATTICE = 2;
		else
		{
			SELFDUALLATTICE = 1;
			cout << "\n  Self-dual lattice not defined. Set to E8xE8.\n" << endl;
		}

	// find the gauge group and create the anomalous U(1) generator
	if (!this->FindGaugeGroup(this->Sectors[0].GetFixedBrane(0), this->StandardConfig, true))
	{
		cout << "\n  Warning in bool COrbifold::Create() : Gauge group not found. Return false." << endl;
		return false;
	}

	this->StandardConfig.SymmetryGroup.observable_sector_GGs.clear();
	this->StandardConfig.SymmetryGroup.observable_sector_U1s.clear();
	this->StandardConfig.SymmetryGroup.GGs_AdditionalLabels.clear();
	this->StandardConfig.SymmetryGroup.U1s_AdditionalLabels.clear();

	CGaugeGroup &GaugeGroup = this->StandardConfig.SymmetryGroup.GaugeGroup;
	s1 = GaugeGroup.factor.size();

	for (i = 0; i < s1; ++i)
	{
		this->StandardConfig.SymmetryGroup.observable_sector_GGs.push_back(i);
		this->StandardConfig.SymmetryGroup.GGs_AdditionalLabels.push_back("");
	}

	s1 = GaugeGroup.u1directions.size();
	for (i = 0; i < s1; ++i)
	{
		this->StandardConfig.SymmetryGroup.observable_sector_U1s.push_back(i);
		this->StandardConfig.SymmetryGroup.U1s_AdditionalLabels.push_back("");
	}

	//create tachyonic representations, disable if tachyon-free model construction
	this->TachyonicStandardConfig = this->StandardConfig;
	this->TachyonicCreateRepresentations();

	if (this->CreateRepresentations() && this->CheckCPPartner())
		this->Orbifold_CheckStatus = CheckedAndGood;

	// Load tachyonic reps into StandardConfig for Spectra comparison to work
	this->StandardConfig.Fields.insert(this->StandardConfig.Fields.end(), this->TachyonicStandardConfig.Fields.begin(), this->TachyonicStandardConfig.Fields.end());

	this->YukawaCouplings.Initiate(this->OrbifoldGroup.GetSpaceGroup(), this->StandardConfig);

	return true;
}



/* ########################################################################################
######   FindGaugeGroup(const CFixedBrane &UntwistedSector, SConfig &VEVConfig, ...) ######
######                                                                               ######
######   Version: 24.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) UntwistedSector            : the CFixedBrane object corresponding to the ######
######                                   untwisted sector where the gauge bosons     ######
######                                   come from                                   ######
######   2) VEVConfig                  : the gauge group is stored here              ######
######   3) CreateAnomalousU1Generator : shall the anomalous U(1) generator be       ######
######                                   created?                                    ######
######   output:                                                                     ######
######   return value                  : finished successfully?                      ######
###########################################################################################
######   description:                                                                ######
######   Find the gauge bosons in "UntwistedSector" and determine the gauge group.   ######
######   The result is stored in "VEVConfig".                                        ######
######################################################################################## */
bool COrbifold::FindGaugeGroup(const CFixedBrane &UntwistedSector, SConfig &VEVConfig, bool CreateAnomalousU1Generator)
{
	const double prec = 0.001;

	CGaugeGroup &GaugeGroup = VEVConfig.SymmetryGroup.GaugeGroup;
	unsigned    &PosAnd     = VEVConfig.SymmetryGroup.Position_of_and_in_GaugeGroup;

	const SelfDualLattice GaugeLattice = this->OrbifoldGroup.GetLattice();

	unsigned i = 0, j = 0, k = 0, l = 0, m = 0;
	bool evs_ok = true;

	size_t s0 = 0;
	size_t s1 = 0;
	size_t s2 = 0;
	size_t s3 = 0;
	size_t s4 = 0;

	vector<vector<double> > SimpleRoots;

	// search the gauge bosons in the untwisted sector
	s0 = UntwistedSector.GetNumberOfInvariantStates();
	for (m = 0; m < s0; ++m)
	{
		const CHalfState     &LeftMover   = UntwistedSector.GetInvariantState(m).GetLeftMover();
		const vector<double> &Eigenvalues = LeftMover.Eigenvalues;

		s1 = this->OrbifoldGroup.GetCentralizer(UntwistedSector.Getconstructing_Element()).size();
		if (s1 == 0)
		{
			cout << "\n  Warning in bool COrbifold::FindGaugeGroup(...) const: Eigenvalues are not computed yet. Return false." << endl;
			return false;
		}

		evs_ok = true;
		for (i = 0; evs_ok && (i < s1); ++i)
		{
			if (fabs(Eigenvalues[i] - round(Eigenvalues[i])) > prec)
				evs_ok = false;
		}
		if (evs_ok)
		{
			// get the set of weights for this state
			const vector<CVector>  &all_Weights      = UntwistedSector.GetMasslessLeftMover(LeftMover.GetIndex()).Weights;
			const vector<unsigned> &unsigned_Weights = LeftMover.Weights;

			s2 = unsigned_Weights.size();
			if (s2 > 1)
			{
				vector<vector<double> > Roots;
				for (i = 0; i < s2; ++i)
					Roots.push_back(all_Weights.at(unsigned_Weights[i]));

				GaugeGroup = determineAlgebra(Roots);

				const vector<gaugeGroupFactor<double> > copy_factors = GaugeGroup.factor;

				// collect all simple roots
				s1 = copy_factors.size();
				for (i = 0; i < s1; ++i)
				{
					const vector<vector<double> > &ggf_SimpleRoots = copy_factors[i].simpleroots;
					SimpleRoots.insert(SimpleRoots.end(), ggf_SimpleRoots.begin(), ggf_SimpleRoots.end());
				}

				// reorder the gauge group factors depending on their origin from the first or second E_8
				if (GaugeLattice == Spin32)
					PosAnd = 16;
				else
				{
					GaugeGroup.factor.clear();

					PosAnd = 0;

					bool only_from_first_E8  = true;
					bool only_from_second_E8 = true;

					// from first E_8
					for (i = 0; i < s1; ++i)
					{
						const gaugeGroupFactor<double> &ggf            = copy_factors[i];
						const vector<double>           &ggf_SimpleRoot = ggf.simpleroots[0];

						only_from_first_E8 = true;
						for (j = 8; only_from_first_E8 && (j < 16); ++j)
						{
							if (fabs(ggf_SimpleRoot[j]) > prec)
								only_from_first_E8 = false;
						}
						if (only_from_first_E8)
						{
							GaugeGroup.factor.push_back(ggf);
							++PosAnd;
						}
					}
					// from second E_8
					for (i = 0; i < s1; ++i)
					{
						const gaugeGroupFactor<double> &ggf            = copy_factors[i];
						const vector<double>           &ggf_SimpleRoot = ggf.simpleroots[0];

						only_from_second_E8 = true;
						for (j = 0; only_from_second_E8 && (j < 8); ++j)
						{
							if (fabs(ggf_SimpleRoot[j]) > prec)
								only_from_second_E8 = false;
						}
						if (only_from_second_E8)
							GaugeGroup.factor.push_back(ggf);
					}
					if (s1 != GaugeGroup.factor.size())
					{
						cout << "\n  Warning in bool COrbifold::FindGaugeGroup(...) const: Cannot sort the gauge group factors. Number of factors before/after sorting: " << s1 << "/" << GaugeGroup.factor.size() << ". Return false." << endl;
						return false;
					}
					// reorder the string "GaugeGroup.algebra" and insert the word "and"
					if ((PosAnd != 0)  && (PosAnd != s1))
					{
						GaugeGroup.algebra = "";
						for (i = 0; i < PosAnd; ++i)
						{
							GaugeGroup.algebra += GaugeGroup.factor[i].algebra;
							if (i + 1 < PosAnd)
								GaugeGroup.algebra += " + ";
						}
						GaugeGroup.algebra += " and ";
						for (i = PosAnd; i < s1; ++i)
						{
							GaugeGroup.algebra += GaugeGroup.factor[i].algebra;
							if (i + 1 < s1)
								GaugeGroup.algebra += " + ";
						}
					}
				}
				break;
			}
		}
	}
	const unsigned number_of_srs = SimpleRoots.size();
	const unsigned number_of_U1s = 16 - number_of_srs;

	size_t number_of_LoadedU1s = this->OrbifoldGroup.LoadedU1Generators.size();
	if ((number_of_LoadedU1s != 0) && (number_of_LoadedU1s > number_of_U1s))
	{
		cout << "\n  Warning in bool COrbifold::FindGaugeGroup(...) const: Two many U(1) generators loaded. Create new ones." << endl;
		number_of_LoadedU1s = 0;
	}

	VEVConfig.SymmetryGroup.IsFirstU1Anomalous = false;
	VEVConfig.SymmetryGroup.D0_FI_term = 0.0;

	if (number_of_U1s == 0)
		return true;

	// create the U1 directions
	GaugeGroup.u1directions.clear();

	vector<vector<double> > D_U1Generators;

	bool createU1generators = true;
	if (number_of_LoadedU1s != 0)
	{
		createU1generators = false;

		double sp = 0.0;

		for (i = 0; !createU1generators && (i < number_of_LoadedU1s); ++i)
		{
			const CVector &U1Generator = this->OrbifoldGroup.LoadedU1Generators[i];
			D_U1Generators.push_back(U1Generator);

			// begin: check the scalar products of each U(1) generator with all simple roots
			for (j = 0; !createU1generators && (j < number_of_srs); ++j)
			{
				const vector<double> &SimpleRoot = SimpleRoots[j];

				sp = 0.0;
				for (k = 0; k < 16; ++k)
					sp += SimpleRoot[k] * U1Generator[k];

				if (fabs(sp) > prec)
					createU1generators = true;
			}
			// end: check the scalar products of each U(1) generator with all simple roots
		}
		if (!createU1generators && (findBasis<double>(D_U1Generators).size() != number_of_LoadedU1s))
			createU1generators = true;

		if (createU1generators)
		{
			cout << "\n  Warning in bool COrbifold::FindGaugeGroup(...) const: Cannot use loaded U(1) generators. Create new ones." << endl;
			this->OrbifoldGroup.LoadedU1Generators.clear();
		}
		else
		{
			// use the loaded U(1) generators
			CreateAnomalousU1Generator = false;
			GaugeGroup.u1directions = D_U1Generators;

			// if not all U(1) generators have been loaded the missing ones must be created
			if (number_of_LoadedU1s != number_of_U1s)
			{
				SimpleRoots.insert(SimpleRoots.end(), D_U1Generators.begin(), D_U1Generators.end());
				createU1generators = true;
			}
		}
	}

	// begin: create the anomalous U(1) generator
	//        using \sum p_sh = const. t_anom
	if (CreateAnomalousU1Generator && (this->OrbifoldGroup.GetNumberOfSupersymmetry() == 0))
	{
		CVector Anomalous_U1_Direction(16);
		unsigned number_of_left_chiral = 0;
		unsigned m = 0;

		s1 = this->Sectors.size();
		for (i = 0; i < s1; ++i)
		{
			const CSector                    &Sector              = this->Sectors[i];
			const vector<CMasslessHalfState> &MasslessRightMovers = Sector.GetMasslessRightMovers();

			s2 = Sector.GetNumberOfFixedBranes();
			for (j = 0; j < s2; ++j)
			{
				const CFixedBrane &FixedBrane = Sector.GetFixedBrane(j);

				s3 = FixedBrane.GetNumberOfInvariantStates();
				for (k = 0; k < s3; ++k)
				{
					const CState &State = FixedBrane.GetInvariantState(k);

					const vector<CVector>  &RM_Weights       = MasslessRightMovers[State.GetRightMover().GetIndex()].Weights;
					const vector<unsigned> &RM_WeightIndices = State.GetRightMover().Weights;

					// is the state left-chiral
					number_of_left_chiral = 0;
					s4 = RM_WeightIndices.size();
					for (l = 0; l < s4; ++l)
					{
						if (fabs(RM_Weights.at(RM_WeightIndices[l])[0] + 0.5) < prec)		//hacking here!!!
							++number_of_left_chiral;
					}
					if (number_of_left_chiral != 0)
					{
						const vector<CVector>  &LM_Weights       = FixedBrane.GetMasslessLeftMover(State.GetLeftMover().GetIndex()).Weights;
						const vector<unsigned> &LM_WeightIndices = State.GetLeftMover().Weights;

						s4 = LM_WeightIndices.size();
						for (l = 0; l < s4; ++l)
						{
							for (m = 0; m < number_of_left_chiral; ++m)
								Anomalous_U1_Direction += LM_Weights.at(LM_WeightIndices[l]);
						}
					}
				}
			}
		}

		for (i = 0; !VEVConfig.SymmetryGroup.IsFirstU1Anomalous && (i < 16); ++i)
		{
			if (fabs(Anomalous_U1_Direction[i]) > prec)
				VEVConfig.SymmetryGroup.IsFirstU1Anomalous = true;
		}

		if (VEVConfig.SymmetryGroup.IsFirstU1Anomalous)
		{
			CVector tmp = Anomalous_U1_Direction * (1.0/12.0);
			SimpleRoots.push_back(tmp);
			GaugeGroup.u1directions.push_back(tmp);
			VEVConfig.SymmetryGroup.IsFirstU1Anomalous = true;
			VEVConfig.SymmetryGroup.D0_FI_term = 12.0 * tmp.GetSqrTo(16);

			/*cout << "Anomalous U(1) generator found:" << endl;
      CPrint Print(Tstandard, &cout);
      Print.PrintVector(tmp);
      cout << endl;*/
		}
	}
	// end: create the anomalous U(1) generator

	if (createU1generators)
	{
		if ((GaugeGroup.u1directions.size() == 0) && (number_of_U1s == 16))
		{
			// gauge group is U(1)^16
			vector<double> u1direction(16,0);
			for (i = 0; i < 16; ++i)
			{
				u1direction.assign(16,0);
				u1direction[i] = 2;
				GaugeGroup.u1directions.push_back(u1direction);
			}
			PosAnd = 0;
			return true;
		}

		// find the directions orthogonal to the simple roots and the loaded U(1) directions / anomalous U(1) direction
		D_U1Generators.clear();
		if (!Find_Basis_Of_Orthogonal_Space(SimpleRoots, GaugeLattice, 16, D_U1Generators))
		{
			cout << "\n  Warning in bool COrbifold::FindGaugeGroup(...) const: Could not find new basis for other U1 directions. Return false." << endl;
			return false;
		}

		GaugeGroup.u1directions.insert(GaugeGroup.u1directions.end(), D_U1Generators.begin(), D_U1Generators.end());
	}
	return true;
}



/* ########################################################################################
######   CreateRepresentations()                                                     ######
######                                                                               ######
######   Version: 29.02.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Computes the U(1) charges and the non-Abelian representations for all       ######
######   massless strings. Creates the fields of the vev-config "StandardConfig".    ######
######################################################################################## */
bool COrbifold::CreateRepresentations()
{
	unsigned i = 0;
	unsigned j = 0;
	unsigned k = 0;
	unsigned l = 0;
	size_t s1 = 0;
	size_t s2 = 0;
	size_t s3 = 0;
	size_t s4 = 0;

	vector<CField> &Fields = this->StandardConfig.Fields;
	Fields.clear();

	unsigned Counter = 0;
	vector<unsigned> FieldCounters(15,1);

	s1 = this->Sectors.size();
	for(i = 0; i < s1; ++i)
	{
		CSector &Sector = this->Sectors[i];

		s2 = Sector.GetNumberOfFixedBranes();
		for (j = 0; j < s2; ++j)
		{
			CFixedBrane &FixedBrane = Sector.AccessFixedBrane(j);

			s3 = FixedBrane.GetNumberOfInvariantStates();
			for (k = 0; k < s3; ++k)
			{
				CState &State = FixedBrane.AccessInvariantState(k);
				State.SetInternalIndex(i,j,k);

				if (!State.CreateRepresentations(FixedBrane, this->StandardConfig, FieldCounters))
				{
					cout << "\n  Warning in bool COrbifold::CreateRepresentations(): Return false." << endl;
					return false;
				}
				s4 = State.GetNumberOfFieldIndices();
				for (l = 0; l < s4; ++l)
				{
					CField &Field = Fields[State.GetFieldIndex(l)];
					Field.SetInternalIndex(i,j,k);
					Field.OriginStandardConfig = Counter;
					++Counter;
				}
			}
		}
	}

	return true;
}



/* ########################################################################################
######   TachyonicCreateRepresentations()                                                     ######
######                                                                               ######
######   Version: 29.02.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Computes the U(1) charges and the non-Abelian representations for all       ######
######   massless strings. Creates the fields of the vev-config "StandardConfig".    ######
######################################################################################## */
bool COrbifold::TachyonicCreateRepresentations()
{
	unsigned i = 0;
	unsigned j = 0;
	unsigned k = 0;
	unsigned l = 0;
	size_t s1 = 0;
	size_t s2 = 0;
	size_t s3 = 0;
	size_t s4 = 0;

	vector<CField> &Fields = this->TachyonicStandardConfig.Fields;
	Fields.clear();

	unsigned Counter = 0;
	vector<unsigned> FieldCounters(15,1);

	s1 = this->Sectors.size();
	for(i = 0; i < s1; ++i)
	{
		CSector &Sector = this->Sectors[i];

		s2 = Sector.GetNumberOfFixedBranes();
		for (j = 0; j < s2; ++j)
		{
			CFixedBrane &FixedBrane = Sector.AccessFixedBrane(j);

			s3 = FixedBrane.TachyonicGetNumberOfInvariantStates();
			for (k = 0; k < s3; ++k)
			{
				CState &State = FixedBrane.TachyonicAccessInvariantState(k);
				State.SetInternalIndex(i,j,k);

				if (!State.TachyonicCreateRepresentations(FixedBrane, this->TachyonicStandardConfig, FieldCounters))
				{
					cout << "\n  Warning in bool COrbifold::TachyonicCreateRepresentations(): Return false." << endl;
					return false;
				}
				s4 = State.GetNumberOfFieldIndices();
				for (l = 0; l < s4; ++l)
				{
					CField &Field = Fields[State.GetFieldIndex(l)];
					Field.SetInternalIndex(i,j,k);
					Field.OriginStandardConfig = Counter;
					++Counter;
				}
			}

			//Create tachyons from excited right-mover level
			s4 = FixedBrane.ExcitedTachyonicGetNumberOfInvariantStates();
			for (k = 0; k < s4; ++k)
			{
				CState &State = FixedBrane.ExcitedTachyonicAccessInvariantState(k);
				State.SetInternalIndex(i,j,(k+s3));									//Continue counting from s3 for k-counter!

				if (!State.ExcitedTachyonicCreateRepresentations(FixedBrane, this->TachyonicStandardConfig, FieldCounters))
				{
					cout << "\n  Warning in bool COrbifold::TachyonicCreateRepresentations(): Return false." << endl;
					return false;
				}
				s4 = State.GetNumberOfFieldIndices();
				for (l = 0; l < s4; ++l)
				{
					CField &Field = Fields[State.GetFieldIndex(l)];
					Field.SetInternalIndex(i,j,(k+s3));								//Continue counting from s3 for k-counter!
					Field.OriginStandardConfig = Counter;							//Increase counter with R-excited tachyonic fields from the same FixedBrane
					++Counter;
				}
			}
		}
	}

	return true;
}



/* ########################################################################################
######   &AccessFixedBrane(const CSpaceGroupElement &Element, bool &FixedBraneFound) ######
######                                                                               ######
######   Version: 19.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Element         : a constructing space group element                     ######
######   output:                                                                     ######
######   2) FixedBraneFound : fixed point/brane with constructing element "Element"  ######
######                        found?                                                 ######
######   return value       : fixed point/brane with constructing element "Element"  ######
###########################################################################################
######   description:                                                                ######
######   Access the fixed point/brane with constructing element "Element".           ######
######################################################################################## */
CFixedBrane &COrbifold::AccessFixedBrane(const CSpaceGroupElement &Element, bool &FixedBraneFound)
{
	size_t s2 = 0;
	unsigned j = 0;

	const size_t s1 = this->Sectors.size();
	for (unsigned i = 0; i < s1; ++i)
	{
		CSector &Sector = this->Sectors[i];

		if ((Element.Get_m() == Sector.Get_m()) && (Element.Get_n() == Sector.Get_n()) && (Element.Get_k() == Sector.Get_k()))
		{
			s2 = Sector.GetNumberOfFixedBranes();
			for (j = 0; j < s2; ++j)
			{
				if (Element == Sector.GetFixedBrane(j).GetSGElement())
				{
					FixedBraneFound = true;
					return Sector.AccessFixedBrane(j);
				}
			}
		}
	}

	FixedBraneFound = false;
	return this->Sectors[0].AccessFixedBrane(0);
}



/* ########################################################################################
######   &AccessSector(const unsigned &i)                                            ######
######                                                                               ######
######   Version: 19.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) i         : index i = 0,..., number of sectors-1                         ######
######   output:                                                                     ######
######   return value : the "i"-th (un)twisted sector                                ######
###########################################################################################
######   description:                                                                ######
######   Access "i"-th (un)twisted sector.                                           ######
######################################################################################## */
CSector &COrbifold::AccessSector(const unsigned &i)
{
#ifdef CHECKERROR
	if (i >= this->Sectors.size())
	{
		cout << "\n  Warning in CSector &COrbifold::AccessSector(...) : Index i out of range. Set i = 0." << endl;
		return this->Sectors[0];
	}
#endif

	return this->Sectors[i];
}



/* ########################################################################################
######   &GetSectors() const                                                         ######
######                                                                               ######
######   Version: 19.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : all (un)twisted sectors                                      ######
###########################################################################################
######   description:                                                                ######
######   Constant access to all (un)twisted sectors.                                ######
######################################################################################## */
const vector<CSector> &COrbifold::GetSectors() const
{
#ifdef CHECKERROR
	if (this->Orbifold_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const vector<CSector> &COrbifold::GetSectors() const : Orbifold ill-defined." << endl;
#endif

	return this->Sectors;
}



/* ########################################################################################
######   &GetSector(const unsigned &i) const                                         ######
######                                                                               ######
######   Version: 19.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) i         : index i = 0,..., number of sectors-1                         ######
######   output:                                                                     ######
######   return value : the "i"-th (un)twisted sector                                ######
###########################################################################################
######   description:                                                                ######
######   Constant access to the "i"-th (un)twisted sector.                           ######
######################################################################################## */
const CSector &COrbifold::GetSector(const unsigned &i) const
{
#ifdef CHECKERROR
	if (this->Orbifold_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const CSector &COrbifold::GetSector(...) const : Orbifold ill-defined." << endl;

	if (i >= this->Sectors.size())
	{
		cout << "\n  Warning in const CSector &COrbifold::GetSector(...) const : Index i out of range. Set i = 0." << endl;
		return this->Sectors[0];
	}
#endif

	return this->Sectors[i];
}



/* ########################################################################################
######   CheckCPPartner() const                                                      ######
######                                                                               ######
######   Version: 31.01.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : do all fields have a CP partner?                             ######
###########################################################################################
######   description:                                                                ######
######   For all left-chiral massless strings (with constructing element g) find the ######
######   right-chiral CP partners (with constructing element g^(-1)).                ######
######################################################################################## */
bool COrbifold::CheckCPPartner() const
{
	bool CheckOk = true;

	const bool SaveProblematicModels = true;

	const vector<CField> &Fields = this->StandardConfig.Fields;
	const size_t f1 = Fields.size();

	const vector<gaugeGroupFactor<double> > &factors = this->StandardConfig.SymmetryGroup.GaugeGroup.factor;

	const size_t number_of_factors = factors.size();

	const CSpaceGroup &SpaceGroup = this->OrbifoldGroup.GetSpaceGroup();

	const bool N1SUSY = (this->OrbifoldGroup.GetNumberOfSupersymmetry() == 1);

	const unsigned M = SpaceGroup.GetM();
	const unsigned N = SpaceGroup.GetN();
	const unsigned K = SpaceGroup.GetK();

	unsigned i = 0;
	unsigned j = 0;
	unsigned k = 0;

	bool realrep = true;

	SUSYMultiplet PartnerSUSYMultiplet = NOT_DEF_SUSY;
	vector<bool> CCPartnerFound(f1, false);
	bool CC = true;
	FPCoordinates fp_i;
	FPCoordinates fp_j;
	CVector BasisVector(6);
	CSpaceGroupElement SGElement_j_inverse;

	if (OrbifoldGroup.GetNumberOfSupersymmetry() == 0) {					//hacking here!!!

		for (i = 0; i < f1; ++i) {

			if (!CCPartnerFound[i])
			{
				const CField             &Field_i     = Fields[i];
				const CSpaceGroupElement &SGElement_i = Field_i.SGElement;
				CC = false;

				switch (Field_i.Multiplet) {
				case Gauge:
				{
					PartnerSUSYMultiplet = bGauge;
					break;
				}
				case bGauge:
				{
					PartnerSUSYMultiplet = Gauge;
					break;
				}
				case Scalar:
				{
					PartnerSUSYMultiplet = bScalar;
					break;
				}
				case bScalar:
				{
					PartnerSUSYMultiplet = Scalar;
					break;
				}
				case LeftFermi:
				{
					PartnerSUSYMultiplet = RightFermi;
					break;
				}
				case RightFermi:
				{
					PartnerSUSYMultiplet = LeftFermi;
					break;
				}
				case moduli:
				{
					PartnerSUSYMultiplet = bmoduli;
					break;
				}
				case bmoduli:
				{
					PartnerSUSYMultiplet = moduli;
					break;
				}
				case Gravity:
				{
					PartnerSUSYMultiplet = GravityCC;
					break;
				}
				case GravityCC:
				{
					PartnerSUSYMultiplet = Gravity;
					break;
				}
				default:
				{
					cout << "\n  Warning in bool COrbifold::CheckCPPartner() const: SUSYMultiplet not defined. Return false." << endl;
					return false;
				}
				}
				for (j = i+1; !CC && (j < f1); ++j)
				{
					if (!CCPartnerFound[j])
					{
						const CField             &Field_j     = Fields[j];
						const CSpaceGroupElement &SGElement_j = Field_j.SGElement;

						// check multiplet matching and twisted sector
						if ((Field_j.Multiplet == PartnerSUSYMultiplet) && ((SGElement_i.Get_m() + SGElement_j.Get_m()) % M == 0)
								&& ((SGElement_i.Get_n() + SGElement_j.Get_n()) % N == 0) && ((SGElement_i.Get_k() + SGElement_j.Get_k()) % K == 0) && (Field_i.U1Charges + Field_j.U1Charges).IsZero())
						{
							// assume that the field "Field_j" is the charge conjugate partner of "Field_i"
							CC = true;

							// check U(1) charges
							if (CC) {
								for (k = 0; CC && (k < number_of_factors); ++k)
								{
									if (!AreHighestWeights_CC(Field_i.HighestWeights_DL[k], Field_j.HighestWeights_DL[k], factors[k])) {
										CC = false;
									}
								}
								if (CC)
								{
									// check the constructing elements
									SGElement_j_inverse = SpaceGroup.SG_Inverse(SGElement_j);
									if (SGElement_i != SGElement_j_inverse)
										CC = SpaceGroup.SG_FromSameConjugationClass(SGElement_i, SGElement_j_inverse);
									if (CC)
									{
										CCPartnerFound[i] = true;
										CCPartnerFound[j] = true;
									}
								}
							}
						}
					}
				}
				if (!CC)
				{
					CheckOk = false;
					cout << "\n  Warning: left and right chiral parts of the spectrum differ:" << endl;
					cout << "  No partner for: ";
					CPrint Print(Tstandard, &cout);
					Print.PrintLabel(Field_i);
					cout << " - ";
					Print.PrintRep(Field_i, this->StandardConfig.SymmetryGroup, true);
					cout << " localized at ";
					Print.PrintSGElement(SGElement_i);
					cout << endl;
				}
			}
		}
	}
	else {											//SuSy case, not optimized for ZMxZNxZK...

		for (i = 0; i < f1; ++i)
		{
			if (!CCPartnerFound[i])
			{
				const CField             &Field_i     = Fields[i];
				const CSpaceGroupElement &SGElement_i = Field_i.SGElement;

				CC = false;
				// begin: find half-hyper for N=2
				if (!N1SUSY && Field_i.U1Charges.IsZero() && ((2 * SGElement_i.Get_m()) % M == 0) && ((2 * SGElement_i.Get_n()) % N == 0))
				{
					realrep = true;
					for (j = 0; realrep && (j < number_of_factors); ++j)
					{
						if (!AreHighestWeights_CC(Field_i.HighestWeights_DL[j], Field_i.HighestWeights_DL[j], factors[j]))
							realrep = false;
					}
					if (realrep)
						CC = true;
				}
				// end: find half-hyper for N=2

				switch (Field_i.Multiplet)
				{
				case LeftChiral:
				{
					PartnerSUSYMultiplet = RightChiral;
					break;
				}
				case RightChiral:
				{
					PartnerSUSYMultiplet = LeftChiral;
					break;
				}
				case Vector:
				{
					PartnerSUSYMultiplet = VectorCC;
					break;
				}
				case VectorCC:
				{
					PartnerSUSYMultiplet = Vector;
					break;
				}
				case Hyper:
				{
					PartnerSUSYMultiplet = Hyper;
					break;
				}
				case Halfhyper:
				{
					PartnerSUSYMultiplet = Halfhyper;
					break;
				}
				case Gravity:
				{
					PartnerSUSYMultiplet = GravityCC;
					break;
				}
				case GravityCC:
				{
					PartnerSUSYMultiplet = Gravity;
					break;
				}
				case LCModulus:
				{
					PartnerSUSYMultiplet = RCModulus;
					break;
				}
				case RCModulus:
				{
					PartnerSUSYMultiplet = LCModulus;
					break;
				}
				default:
				{
					cout << "\n  Warning in bool COrbifold::CheckCPPartner() const: SUSYMultiplet not defined. Return false." << endl;
					return false;
				}
				}

				for (j = i+1; !CC && (j < f1); ++j)
				{
					if (!CCPartnerFound[j])
					{
						const CField             &Field_j     = Fields[j];
						const CSpaceGroupElement &SGElement_j = Field_j.SGElement;

						CC = false;

						// check SUSY type and twisted sector
						if ((Field_j.Multiplet == PartnerSUSYMultiplet)
								&& ((SGElement_i.Get_m() + SGElement_j.Get_k()) % M == 0)
								&& ((SGElement_i.Get_n() + SGElement_j.Get_n()) % N == 0)
								&& (Field_i.U1Charges + Field_j.U1Charges).IsZero())
						{
							// assume that the field "Field_j" is the charge conjugate partner of "Field_i"
							CC = true;

							// check R charges for N=1 SUSY
							if (N1SUSY)
							{
								if ((!(Field_i.q_sh + Field_j.q_sh).IsZero()) || (!(Field_i.GetState(this->Sectors).GetOsciContribution() + Field_j.GetState(this->Sectors).GetOsciContribution()).IsZero()))
									CC = false;
							}

							// check U(1) charges
							if (CC)
							{
								for (k = 0; CC && (k < number_of_factors); ++k)
								{
									if (!AreHighestWeights_CC(Field_i.HighestWeights_DL[k], Field_j.HighestWeights_DL[k], factors[k]))
										CC = false;
								}

								if (CC)
								{
									// check the constructing elements
									if (N1SUSY)
									{
										SGElement_j_inverse = SpaceGroup.SG_Inverse(SGElement_j);
										if (SGElement_i != SGElement_j_inverse)
											CC = SpaceGroup.SG_FromSameConjugationClass(SGElement_i, SGElement_j_inverse);
									}

									if (CC)
									{
										CCPartnerFound[i] = true;
										CCPartnerFound[j] = true;
									}
								}
							}
						}
					}
				}
				if (!CC)
				{
					CheckOk = false;
					cout << "\n  Warning: left and right chiral parts of the spectrum differ:" << endl;
					cout << "  No partner for: ";
					CPrint Print(Tstandard, &cout);
					Print.PrintLabel(Field_i);
					cout << " - ";
					Print.PrintRep(Field_i, this->StandardConfig.SymmetryGroup, true);
					cout << " localized at ";
					Print.PrintSGElement(SGElement_i);
					cout << endl;

					if (SaveProblematicModels)
					{
						string Filename = "ProblematicModels_NoCPPartner.txt";
						ofstream tmp_out;
						tmp_out.open(Filename.data(), ofstream::app | ios::ate);
						this->OrbifoldGroup.PrintToFile(tmp_out);
						tmp_out.close();
					}
				}
			}
		}
	}

	return CheckOk;
}



/* ########################################################################################
######   CheckAnomaly(SConfig &VEVConfig, ...) const                                 ######
######                                                                               ######
######   Version: 28.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) VEVConfig    : contains all fields                                       ######
######   2) GaugeIndices : contains the quadratic and cubic indices of non-Abelian   ######
######                     representations                                           ######
######   3) Print        : if "info" is true print the output to this CPrint object  ######
######   4) info         : print info?                                               ######
######   output:                                                                     ######
######   return value    : finished succesfully and spectrum is free of anomalies?   ######
###########################################################################################
######   description:                                                                ######
######   Checks that all gauge anomalies of the spectrum of "VEVConfig" vanish or    ######
######   can be canceled by the Green-Schwarz mechanism.                             ######
######################################################################################## */
bool COrbifold::CheckAnomaly(SConfig &VEVConfig, const CGaugeIndices &GaugeIndices, CPrint &Print, bool info, double AddFactor) const
{
	if (fabs(AddFactor) < 0.0001)
	{
		(*Print.out) << "\n  Warning in bool COrbifold::CheckAnomaly(...) const: \"AddFactor\" ill-defined. Return false." << endl;
		return false;
	}

	bool AllAnomaliesCancel = true;

	const vector<CField> &Fields = VEVConfig.Fields;
	const size_t f1 = Fields.size();

	if (f1 == 0)
	{
		(*Print.out) << "\n  Warning in bool COrbifold::CheckAnomaly(...) const: Spectrum is empty. Return true." << endl;
		return true;
	}

	unsigned i = 0;
	unsigned j = 0;
	unsigned k = 0;

	const CGaugeGroup          &GaugeGroup    = VEVConfig.SymmetryGroup.GaugeGroup;
	const vector<doubleVector> &U1_Directions = GaugeGroup.u1directions;
	const size_t number_of_factors = GaugeGroup.factor.size();
	const size_t number_of_U1s     = U1_Directions.size();

	// for N=2 SUSY check n(hyper) - n(vector) = 244
	if (this->OrbifoldGroup.GetNumberOfSupersymmetry() == 2)
	{
		unsigned mult = 1;
		unsigned number_of_hhypers = 0;
		unsigned number_of_vectors = number_of_U1s;
		// {NOT_DEF_SUSY = 0, AnyKind, LeftChiral, RightChiral, Vector, VectorCC, Hyper, Halfhyper, Gravity, GravityCC};
		for (i = 0; i < f1; ++i)
		{
			const CField &Field = Fields[i];
			mult = 1;
			for (j = 0; j < number_of_factors; ++j)
				mult *= abs(Field.Dimensions[j].Dimension);

			switch (Field.Multiplet)
			{
			case Gravity:
				break;
			case GravityCC:
				break;
			case VectorCC:
				break;
			case Vector:
				number_of_vectors += mult;
				break;
			case Hyper:
				number_of_hhypers += 2 * mult;
				break;
			case Halfhyper:
				number_of_hhypers += mult;
				break;
			default:
			{
				(*Print.out) << "\n  Warning in bool COrbifold::CheckAnomaly(...) const: N=2 SUSY multiplet not defined. Return false." << endl;
				return false;
			}
			}
		}
		if (info)
		{
			(*Print.out) << "\n";
			(*Print.out) << "  " << Print.cbegin << "check n(hyper) - n(vector) = 244:" << Print.cend << "\n";
			(*Print.out) << "    " << Print.cbegin << "n(halfhyper) = " << number_of_hhypers << Print.cend << "\n";
			(*Print.out) << "    " << Print.cbegin << "n(vector)    = " << number_of_vectors << Print.cend << endl;
		}
		if ((number_of_hhypers/2) - number_of_vectors != 244)
		{
			(*Print.out) << "\n  " << Print.cbegin << "4d N=2 gauge theory has anomalies." << Print.cend << endl;
			AllAnomaliesCancel = false;
		}
	}

	// N=1 case
	const vector<vector<int> > &CubicIndices = GaugeIndices.GetCubicIndices();

	vector<CHugeInt> total_Dimension(f1, 0);
	// begin: compute the total dimensions of all fields representations
	// e.g.
	// gauge group SO(10) x SU(3) x U(1) x E_8'
	// 2 (16,3,1)_+1
	// 1 (16,3,1)_-1  -> 3x16x3 = 144 particles

	vector<rational<CHugeInt> >          tmp_RationalU1Charges(number_of_U1s, CHugeInt(0));
	vector<vector<rational<CHugeInt> > > RationalU1Charges;

	CHugeInt tmp = 1;
	for (i = 0; i < f1; ++i)
	{
		tmp_RationalU1Charges.assign(number_of_U1s, CHugeInt(0));

		const CField &Field = Fields[i];
		if (Field.Multiplet == LeftFermi)								//hacking here and below!!!!
		{
			const CVector &U1Charges = Field.U1Charges;
			for (j = 0; j < number_of_U1s; ++j)
				tmp_RationalU1Charges[j] = D2RatHugeInt(U1Charges[j]);

			const RepVector &Dimensions = Field.Dimensions;
			tmp = 1;
			for (j = 0; j < number_of_factors; ++j)
				tmp *= (CHugeInt)abs(Dimensions[j].Dimension);

			total_Dimension[i] = tmp;
		}
		RationalU1Charges.push_back(tmp_RationalU1Charges);
	}
	// end: compute the total dimension of all fields representations

	rational<CHugeInt> sum_anomaly = CHugeInt(0);

	bool AnomalyInfoPrinted = false;
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// begin: SU(N) - anomaly and global SU(2) anomaly
	if (info)
		(*Print.out) << "\n  " << Print.cbegin << "SU(N)^3 anomaly:" << Print.cend << "\n";

	for (i = 0; i < number_of_factors; ++i)
	{
		const gaugeGroupFactor<double> &ggf = GaugeGroup.factor[i];
		const string   &algebra = ggf.algebra;
		const unsigned Rank     = ggf.rank;

		// SU(N) algebra
		if (algebra[0] == 'A')
		{
			// check: global SU(2) - anomaly
			if (Rank == 1)
			{
				CHugeInt number_of_su2 = CHugeInt(0);

				// begin: run through all fields
				for (k = 0; k < f1; ++k)
				{
					const CField &Field = Fields[k];
					if (Field.Multiplet == LeftFermi)
					{
						if (abs(Field.Dimensions[i].Dimension) == 2)
							number_of_su2 += total_Dimension[k] / CHugeInt(2);
					}
				}
				// end: run through all fields

				if (info)
				{
					AnomalyInfoPrinted = true;
					(*Print.out) << "    " << Print.cbegin << i+1 << "-th gauge group factor A1: #(doublets) = " << number_of_su2 << Print.cend << "\n";
				}

				// if SU(2)_i is anomalous
				if ((number_of_su2 % 2) != CHugeInt(0))
				{
					(*Print.out) << "\n  " << Print.cbegin << "SU(2)_" << i+1 << " gauge theory has anomalies." << Print.cend << endl;
					AllAnomaliesCancel = false;
				}
			}
			// check: SU(N) - anomaly
			else
			{
				const vector<int> &CubicIndices_of_AN = CubicIndices[Rank];

				sum_anomaly = CHugeInt(0);

				// begin: run through all fields
				for (k = 0; k < f1; ++k)
				{
					const CField &Field = Fields[k];
					if (Field.Multiplet == LeftFermi)
					{
						const int Dimension   = abs(Field.Dimensions[i].Dimension);
						const int Cubic_Index = CubicIndices_of_AN[(unsigned)Dimension];

						if (Cubic_Index == -1)
						{
							(*Print.out) << "\n  Warning in bool COrbifold::CheckAnomaly(...) const.\nCubic index of rep " << Dimension << " of " << algebra[0] << Rank << " is not known. Return false." << endl;
							return false;
						}

						if (Cubic_Index != 0)
						{
							if (Fields[k].Dimensions[i].Dimension < 0)
								sum_anomaly -= (total_Dimension[k] / Dimension) * Cubic_Index;
							else
								sum_anomaly += (total_Dimension[k] / Dimension) * Cubic_Index;
						}
					}
				}
				// end: run through all fields

				if (info)
				{
					AnomalyInfoPrinted = true;
					(*Print.out) << "    " << Print.cbegin << i+1 << "-th gauge group factor " << algebra << ": tr cubic    = " << sum_anomaly << Print.cend << "\n";
				}

				// if SU(N)_i is anomalous
				if (sum_anomaly != CHugeInt(0))
				{
					(*Print.out) << "\n  " << Print.cbegin << "SU(N)_" << i+1 << " gauge theory has anomalies." << Print.cend << endl;
					if (!info)
						(*Print.out) << "    " << Print.cbegin << i+1 << "-th gauge group factor is " << algebra << ": tr cubic    = " << sum_anomaly << Print.cend << "\n";
					AllAnomaliesCancel = false;
				}
			}
		}
	}
	// begin: SU(N) - anomaly and global SU(2) anomaly
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if (info && !AnomalyInfoPrinted)
		(*Print.out) << "    " << Print.cbegin << "no SU(N) factors" << Print.cend << "\n";

	if (number_of_U1s == 0)
		return true;

	rational<CHugeInt> X = CHugeInt(1);
	if (AddFactor != 1.0)
		X = D2RatHugeInt(AddFactor);

	rational<CHugeInt> tmp_i = CHugeInt(0);
	vector<rational<CHugeInt> > U1LengthSquare(number_of_U1s, CHugeInt(0));
	for (i = 0; i < number_of_U1s; ++i)
	{
		const doubleVector &U1_direction_i = U1_Directions[i];
		tmp_i = CHugeInt(0); //tmp_i = 0;
		for (j = 0; j < 16; ++j)
			tmp_i += D2RatHugeInt(U1_direction_i[j]) * D2RatHugeInt(U1_direction_i[j]);
		U1LengthSquare[i] = tmp_i;
	}

	bool MultipleAnomalousU1s = false;
	bool NoAnomaly = true;

	rational<CHugeInt> sum_Qi       = CHugeInt(0);
	rational<CHugeInt> sum_Q0       = CHugeInt(0);
	rational<CHugeInt> sum_QiQjQj   = CHugeInt(0);
	rational<CHugeInt> sum_Qi_three = CHugeInt(0);

	vector<rational<CHugeInt> > sum_Qis;

	if (info)
		(*Print.out) << "\n  " << Print.cbegin << "U(1)_i - U(1)_i - U(1)_i, U(1)_i - U(1)_j - U(1)_j and U(1)_i - grav. - grav. anomalies:" << Print.cend << "\n";

	// begin: pure U(1) and mixed grav. anomalies
	for (i = 0; i < number_of_U1s; ++i)
	{
		sum_Qi       = CHugeInt(0);
		sum_Qi_three = CHugeInt(0);

		// begin: run through all fields
		for (j = 0; j < f1; ++j)
		{
			const CField &Field = Fields[j];
			if (Field.Multiplet == LeftFermi)
			{
				const rational<CHugeInt> &U1_Charge = RationalU1Charges[j][i];

				sum_Qi       += U1_Charge * total_Dimension[j];
				sum_Qi_three += U1_Charge * U1_Charge * U1_Charge * total_Dimension[j];
			}
		}
		// end: run through all fields

		if (info)
		{
			(*Print.out) << "    " << Print.cbegin << "tr Q_" << i+1 << "       = " << sum_Qi << Print.cend << "\n";
			(*Print.out) << "    " << Print.cbegin << "tr Q_" << i+1 << "^3     = " << sum_Qi_three << Print.cend << "\n";
		}
		sum_Qis.push_back(sum_Qi);

		// if U(1)_i is anomalous
		if ((sum_Qi != CHugeInt(0)) || (sum_Qi_three != CHugeInt(0)))
		{
			if (i == 0)
				sum_Q0 = sum_Qi;
			else
				MultipleAnomalousU1s = true;

			NoAnomaly = false;
			if (VEVConfig.SymmetryGroup.IsFirstU1Anomalous && (D2RatHugeInt(VEVConfig.SymmetryGroup.D0_FI_term) != sum_Qi))
			{
				(*Print.out) << "\n  Warning in bool COrbifold::CheckAnomaly(...) const: FI term not correct. Return false." << endl;
				return false;
			}

			// begin check: tr Q_i^3 = 1/4 * |t_i|^2 tr Q_i
			if (sum_Qi_three != (U1LengthSquare[i] * X * sum_Qi)/CHugeInt(4))
			{
				(*Print.out) << "\n  " << Print.cbegin << "Anomalies are not universal: check tr Q_" << i+1 << "^3 = 1/4 * |t_" << i+1 << "|^2 tr Q_" << i+1 << "." << Print.cend << "\n";
				(*Print.out) << "    " << Print.cbegin << "tr Q_" << i+1 << "   = " << sum_Qi << Print.cend << "\n";
				(*Print.out) << "    " << Print.cbegin << "tr Q_" << i+1 << "^3 = " << sum_Qi_three << Print.cend << "\n";
				(*Print.out) << "    " << Print.cbegin << "|t_" << i+1 << "|^2  = " << U1LengthSquare[i] << Print.cend << endl;
				AllAnomaliesCancel = false;
			}
			// end check: tr Q_i^3 = 1/4 * |t_i|^2 tr Q_i
		}

		for (j = 0; j < number_of_U1s; ++j)
		{
			if (i != j)
			{
				sum_QiQjQj = CHugeInt(0);

				// begin: run through all fields
				for (k = 0; k < f1; ++k)
				{
					const CField &Field = Fields[k];
					if (Field.Multiplet == LeftFermi)
					{
						const vector<rational<CHugeInt> > &U1Charges = RationalU1Charges[k];

						sum_QiQjQj += U1Charges[i] * U1Charges[j] * U1Charges[j] * total_Dimension[k];
					}
				}
				// end: run through all fields

				if (info)
					cout << "    " << Print.cbegin << "tr Q_" << i+1 << " Q_" << j+1 << "^2 = " << sum_QiQjQj << Print.cend << "\n";

				if (sum_QiQjQj == CHugeInt(0))
				{
					if (!NoAnomaly && (i == 0) && (j != 0))
					{
						(*Print.out) << "\n  " << Print.cbegin << "Anomalies are not universal: first U(1) is anomalous, but tr Q_" << j+1 << "^2 Q_anom = 0." << Print.cend << endl;
						AllAnomaliesCancel = false;
					}
				}
				else
				{
					if (NoAnomaly)
					{
						(*Print.out) << "\n  " << Print.cbegin << "Anomalies are not universal: tr Q_" << j+1 << "^2 Q_" << i+1 << " = " << sum_QiQjQj << " != 0, but no U(1) is anomalous." << Print.cend << endl;
						AllAnomaliesCancel = false;
					}

					// begin check: 1/(2 |t_j|^2) tr Q_j^2 Q_i = 1/24 * tr Q_i
					if (VEVConfig.SymmetryGroup.IsFirstU1Anomalous && (i == 0) && (sum_QiQjQj/(CHugeInt(2) * U1LengthSquare[j]) != ((X * sum_Qi) / CHugeInt(24))))
					{
						(*Print.out) << "\n  " << Print.cbegin << "Anomalies are not universal: check 1/(2 |t_" << j+1<< "|^2) tr Q_" << j+1 << "^2 Q_" << i+1 << " = 1/24 * tr Q_" << i+1 << "." << Print.cend << "\n";
						(*Print.out) << "    " << Print.cbegin << "tr Q_" << i+1 << "       = " << sum_Qi << Print.cend << "\n";
						(*Print.out) << "    " << Print.cbegin << "tr Q_" << i+1 << " Q_" << j+1 << "^2 = " << sum_QiQjQj << Print.cend << "\n";
						(*Print.out) << "    " << Print.cbegin << "|t_" << i+1 << "|^2 = " << U1LengthSquare[i] << Print.cend << "\n";
						(*Print.out) << "    " << Print.cbegin << "|t_" << j+1 << "|^2 = " << U1LengthSquare[j] << Print.cend << endl;
						AllAnomaliesCancel = false;
					}
					// end check: 1/(2 |t_j|^2) tr Q_j^2 Q_i = 1/24 * tr Q_i
				}
			}
		}
	}
	// end: pure U(1) and mixed grav. anomalies

	if (!MultipleAnomalousU1s && (sum_Q0 != CHugeInt(0)))
	{
		VEVConfig.SymmetryGroup.D0_FI_term = RatHugeInt2D(sum_Q0);
		VEVConfig.SymmetryGroup.IsFirstU1Anomalous = true;

		if ((this->OrbifoldGroup.LoadedU1Generators.size() == 0) && ((X * sum_Q0) != (CHugeInt(12) * U1LengthSquare[0])))
		{
			(*Print.out) << "  Warning in bool COrbifold::CheckAnomaly(...) const: FI term not given by 12|t_anom|^2. Return false." << endl;
			AllAnomaliesCancel = false;
		}
	}

	unsigned Index = 0;

	// SU(N) - SU(N) - U(1)_j  diagram
	if (info)
		(*Print.out) << "\n  " << Print.cbegin << "G - G - U(1)_j anomalies for G = SU(N), SO(N) or E_N:" << Print.cend << "\n";

	// loop: all non-abelian gauge group factors
	for (i = 0; i < number_of_factors; ++i)
	{
		const gaugeGroupFactor<double> &ggf = GaugeGroup.factor[i];

		// loop: all U(1) charges
		for (j = 0; j < number_of_U1s; ++j)
		{
			sum_anomaly = CHugeInt(0);

			// begin: run through all fields
			for (k = 0; k < f1; ++k)
			{
				const CField &Field = Fields[k];
				if (Field.Multiplet == LeftFermi)
				{
					const unsigned &Dimension = abs(Field.Dimensions[i].Dimension);
					if (!GaugeIndices.GetQuadraticIndex(ggf, Dimension, Index))
						return false;

					sum_anomaly += total_Dimension[k] * RationalU1Charges[k][j] * CHugeInt(Index) / CHugeInt(Dimension); //sum_anomaly += total_Dimension[k] * RationalU1Charges[k][j] * Index / Dimension;
				}
			}
			// end: run through all fields

			if (info)
				(*Print.out) << "    " << Print.cbegin << "tr l(rep. of " << ggf.algebra << ") Q_" << j+1 << " = " << sum_anomaly << Print.cend << "\n";

			if (sum_anomaly == CHugeInt(0))
			{
				if (!NoAnomaly && (j == 0))
				{
					(*Print.out) << "\n  " << Print.cbegin << "Anomalies are not universal: first U(1) is anomalous, but tr Q_" << j+1 << " l(rep. of " << ggf.algebra << ") = 0." << endl;
					AllAnomaliesCancel = false;
				}
			}
			else
			{
				if (NoAnomaly)
				{
					(*Print.out) << "\n  " << Print.cbegin << "Anomalies are not universal: first U(1) is not anomalous, but tr Q_" << j+1 << " l(rep. of " << ggf.algebra << ") = " << sum_anomaly << "." << endl;
					AllAnomaliesCancel = false;
				}

				// begin check: 12 tr index Q = tr Q
				if (CHugeInt(12) * sum_anomaly != X * sum_Qis[j])
				{
					(*Print.out) << "\n  " << Print.cbegin << "Anomalies are not universal: check 12 tr l(rep. of " << ggf.algebra << ") Q_" << j+1 << " = tr Q_" << j+1 << ".\n";
					(*Print.out) << "    " << Print.cbegin << "tr Q_" << j+1 << " = " << sum_Qis[j] << Print.cend << "\n";
					(*Print.out) << "    " << Print.cbegin << "tr l(rep. of " << ggf.algebra << ") Q_" << j+1 << " = " << sum_anomaly << Print.cend << endl;
					AllAnomaliesCancel = false;
				}
				// end check: 12 tr index Q = tr Q
			}
		}
	}

	if (info)
	{
		if (AllAnomaliesCancel)
			(*Print.out) << "\n  " << Print.cbegin << "All anomalies are universal, i.e. their ratios are OK." << Print.cend << "\n" << endl;
		else
			(*Print.out) << "\n  " << Print.cbegin << "Anomalies are not universal, i.e. their ratios are not OK." << Print.cend << "\n" << endl;
	}

	return AllAnomaliesCancel;
}



/* ########################################################################################
######   CheckDiscreteAnomaly(...)                                                   ######
######                                                                               ######
######   Version: 29.02.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) VEVConfig         : contains all fields                                  ######
######   2) GaugeIndices      : contains the quadratic and cubic indices of non-     ######
######                          Abelian representations                              ######
######   4) Print             : if "info" is true print the output to this CPrint    ######
######                          object                                               ######
######   5) info              : print info?                                          ######
######   output:                                                                     ######
######   3) anomalous_element : the space group element whose gauge embedding is     ######
######                          t_anom, i.e. the generator of the anomalous U(1)     ######
######   return value         : finished succesfully and spectrum is free of         ######
######                          anomalies?                                           ######
###########################################################################################
######   description:                                                                ######
######   Checks that all gauge anomalies of the spectrum of "VEVConfig" vanish or    ######
######   can be canceled by the Green-Schwarz mechanism.                             ######
######################################################################################## */
bool COrbifold::CheckDiscreteAnomaly(const SConfig &VEVConfig, const CGaugeIndices &GaugeIndices, CSpaceGroupElement &anomalous_element, CPrint &Print, bool info) const
{
	bool AnomaliesNonUniversal = false;

	const vector<CField>          &Fields       = VEVConfig.Fields;
	const CGaugeGroup             &GaugeGroup   = VEVConfig.SymmetryGroup.GaugeGroup;
	const vector<vector<double> > &U1Directions = GaugeGroup.u1directions;

	// Set some general variables
	CVector AnomU1(16);
	if (VEVConfig.SymmetryGroup.IsFirstU1Anomalous)
		AnomU1 = GaugeGroup.u1directions[0];

	const size_t number_of_factors = GaugeGroup.factor.size();
	const size_t number_of_U1s     = U1Directions.size();

	if (number_of_factors == 0)
		return true;

	const SelfDualLattice Lattice = this->OrbifoldGroup.GetShift(0).GetLattice();

	unsigned i = 0;
	unsigned j = 0;
	unsigned k = 0;

	unsigned anomalous_element_k = 0;
	unsigned anomalous_element_l = 0;
	rationalVector anomalous_element_n(6,0);

	rational<int> DiscreteCharge = 0;
	rational<int> AnomalyPolynomial;
	rational<int> GravAnomalyPolynomial;
	rational<int> first_AnomalyPolynomial;

	rational<int> Normalization  = 1;

	unsigned Index = 0;
	rational<int> Factor = 1;
	rational<int> tmp    = 1;

	rational<int> tmp1 = 1;
	rational<int> tmp2 = 1;


	const CSpaceGroup          &SpaceGroup = this->OrbifoldGroup.GetSpaceGroup();
	const vector<CTwistVector> &Twists     = SpaceGroup.GetTwists();
	const bool                 &ZMxZN      = SpaceGroup.IsZMxZN();

	const vector<SDiscreteSymmetry> &DiscreteNonRSymmetries = SpaceGroup.DiscreteNonRSymmetries;
	const vector<SDiscreteSymmetry> &DiscreteRSymmetries    = SpaceGroup.DiscreteRSymmetries;
	const vector<SModularSymmetry>  &ModularSymmetries      = SpaceGroup.ModularSymmetries;

	const size_t number_of_R_symmetries    = DiscreteRSymmetries.size();
	const size_t number_of_NonR_symmetries = DiscreteNonRSymmetries.size();
	const size_t number_of_mod_symmetries  = ModularSymmetries.size();

	const size_t number_of_symmetries = number_of_NonR_symmetries + number_of_R_symmetries;

	if ((number_of_symmetries == 0) && (number_of_mod_symmetries == 0))
	{
		if (info)
			(*Print.out) << "\n  " << Print.cbegin << "No discrete symmetry defined in the geometry file. Hence, cannot check discrete anomaly." << Print.cend << "\n";
		return true;
	}

	const vector<SDiscreteSymmetry> *DiscreteSymmetries = NULL;
	unsigned shift = 0;

	bool IsR = true;
	string Rstring = "";
	rational<int> Order = 0;
	bool SmallestIndex = true;

	rational<int> Test;

	// begin: compute the N=2 contributions
	vector<vector<rational<int> > > BetaFunctionCoefficients;
	vector<bool> N2_Sector_Exists(3, false);

	if ((number_of_R_symmetries != 0) || (number_of_mod_symmetries != 0))
	{
		CAnalyseModel Analyse;
		Analyse.ComputeN2BetaFunctionCoefficient(this, VEVConfig, GaugeIndices, N2_Sector_Exists, BetaFunctionCoefficients);
	}
	// end: compute the N=2 contributions


	//for (i = 0; i < number_of_symmetries; ++i)
	for (i = 0; i < number_of_NonR_symmetries; ++i)
	{
		if (i < number_of_NonR_symmetries)
		{
			DiscreteSymmetries = &DiscreteNonRSymmetries;
			shift = 0;
		}
		else
		{
			DiscreteSymmetries = &DiscreteRSymmetries;
			shift = number_of_NonR_symmetries;
		}

		const SDiscreteSymmetry &DiscreteSymmetry = DiscreteSymmetries->at(i - shift);

		Order = DiscreteSymmetry.Order;
		IsR = (DiscreteSymmetry.SuperpotentialCharge != 0);

		if (IsR)
			Rstring = "^R";
		else
			Rstring = "";

		SmallestIndex = false;

		if (info)
		{
			(*Print.out) << "\n  ";
			Print.PrintDiscreteSymmetry(DiscreteSymmetry);
			(*Print.out) << "\n    -----------------------------" << endl;
		}

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// begin: Z_N (R) - grav - grav
		if (IsR) // -21 from gravitinos
			GravAnomalyPolynomial = int(-21 + number_of_U1s) * (DiscreteSymmetry.SuperpotentialCharge/2); //GravAnomalyPolynomial = (-21 + number_of_U1s) * (DiscreteSymmetry.SuperpotentialCharge/2);
		else
			GravAnomalyPolynomial = 0;

		// begin: run through all (moduli, vector and) left-chiral fields
		for( vector<CField>::const_iterator it_field = Fields.begin(); it_field != Fields.end(); ++it_field)
		{
			// modulus and vector multiplet
			if (IsR)
			{
				const RepVector &Dimensions = it_field->Dimensions;

				Factor = 1;
				for (k = 0; k < number_of_factors; ++k)
					Factor *= abs(Dimensions[k].Dimension);

				if ((it_field->Multiplet == Vector) || (it_field->Multiplet == LCModulus))
					GravAnomalyPolynomial += Factor;
			}

			// left-chiral multiplet
			if (it_field->Multiplet == LeftChiral)
			{
				const RepVector &Dimensions = it_field->Dimensions;

				Factor = 1;
				for (k = 0; k < number_of_factors; ++k)
					Factor *= abs(Dimensions[k].Dimension);

				GravAnomalyPolynomial += Factor * it_field->GetDiscreteCharge(DiscreteSymmetry);
			}
		}
		// end: run through all left-chiral fields

		if (info)
			(*Print.out) << "    A_{Z_" << Order.numerator() << Rstring << " - grav - grav} = " << GravAnomalyPolynomial << "\n";
		// end: Z_N (R) - grav - grav
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		/*if (!IsR && (number_of_U1s != 0))
    {
    unsigned start = 0;
    CVector U1Gen(16);
    int den;
    int norm = 1;
    rational<int> charge;
    rational<int> U1Squared;

    if (VEVConfig.SymmetryGroup.IsFirstU1Anomalous)
      start = 1;
      AnomU1 = GaugeGroup.u1directions[0];

    if (start < number_of_U1s)
      (*Print.out) << "    -----------------------------" << endl;

    for (j = start; j < number_of_U1s; ++j)
    {
      norm = 1;
      for( vector<CField>::const_iterator it_field = Fields.begin(); it_field != Fields.end(); ++it_field)
      {
        den = (D2Rat(it_field->U1Charges[j]) * rational<int>(norm,1)).denominator();
        if (den != 1)
          norm *= den; 
      }
      //cout << "\n\n  norm = " << norm << endl << endl;

      U1Gen = GaugeGroup.u1directions[j];
      U1Squared = D2Rat(U1Gen.GetSqrTo(16));

      cout << "norm U(1) gen. = " << U1Squared << endl;

      AnomalyPolynomial = 0;
      for( vector<CField>::const_iterator it_field = Fields.begin(); it_field != Fields.end(); ++it_field)
      {
        // left-chiral multiplet
        if (it_field->Multiplet == LeftChiral)
        {
          const RepVector &Dimensions = it_field->Dimensions;

          //cout << "rep = ";
          //Print.PrintRep(Dimensions);
          //cout << " and U(1) charge = " << it_field->U1Charges[j] << " = " << D2Rat(it_field->U1Charges[j]) << " contributes with charge ";

          charge = D2Rat(it_field->U1Charges[j]);
          //cout << charge << "";
          if (charge != 0)
          {
            Factor = 1;

            for (k = 0; k < number_of_factors; ++k)
              Factor *= abs(Dimensions[k].Dimension);

            //cout << "hier" << endl;
            AnomalyPolynomial += Factor * charge * it_field->GetDiscreteCharge(DiscreteSymmetry) * it_field->GetDiscreteCharge(DiscreteSymmetry);
            //cout << " Z_N = " << it_field->GetDiscreteCharge(DiscreteSymmetry);
            //cout << " Factor = " << Factor;
          }
          //cout << "\n  -> " << AnomalyPolynomial << endl;
          // end: compute the charges under the discrete ZN symmetries
        }
      }
      AnomalyPolynomial = AnomalyPolynomial/U1Squared;
      // end: run through all left-chiral fields

      if (info)
        (*Print.out) << "    A_{Z_" << Order.numerator() << Rstring << " - U(1)_" << j+1 << " - U(1)_" << j+1 << "} = " << AnomalyPolynomial << "  ";

      if (j == start)
        first_AnomalyPolynomial = AnomalyPolynomial;

      tmp = (12 * AnomalyPolynomial - GravAnomalyPolynomial)/Order;
      if (tmp.denominator() != 1)
      {
        AnomaliesNonUniversal = true;
        if (info)
          (*Print.out) << " Z_N - grav - grav anomaly not universal";
      }

      tmp = (AnomalyPolynomial - first_AnomalyPolynomial)/Order;
      if (tmp.denominator() != 1)
      {
        AnomaliesNonUniversal = true;
        if (info)
          (*Print.out) << " Z_N - U(1) - U(1) anomaly not universal";
      }
      if (info)
        (*Print.out) << endl;
    }
  }*/

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// begin: Z_N (R) - G - G
		// run through all gauge group factors
		for (j = 0; j < number_of_factors; ++j)
		{
			if (info && ((j == 0) || (j == VEVConfig.SymmetryGroup.Position_of_and_in_GaugeGroup)))
				(*Print.out) << "    -----------------------------" << endl;

			const gaugeGroupFactor<double> &ggf = GaugeGroup.factor[j];

			if (IsR)
			{
				if (!GaugeIndices.GetQuadraticIndexAdj(ggf, Index))
					return false;

				AnomalyPolynomial = int(Index) * (DiscreteSymmetry.SuperpotentialCharge/2); //AnomalyPolynomial = Index * (DiscreteSymmetry.SuperpotentialCharge/2);
			}
			else
				AnomalyPolynomial = 0;

			Test = 0;

			// begin: run through all left-chiral fields
			for( vector<CField>::const_iterator it_field = Fields.begin(); it_field != Fields.end(); ++it_field)
			{
				// left-chiral multiplet
				if (it_field->Multiplet == LeftChiral)
				{
					const RepVector &Dimensions = it_field->Dimensions;

					// begin: compute the charges under the discrete ZN symmetries
					if (!GaugeIndices.GetQuadraticIndex(ggf, abs(Dimensions[j].Dimension), Index))
						return false;

					if (Index != 0)
					{
						Factor = Index;

						for (k = 0; k < number_of_factors; ++k)
						{
							if (k != j)
								Factor *= abs(Dimensions[k].Dimension);
						}

						if (IsR)
						{
							if (it_field->GetDiscreteCharge(DiscreteSymmetry).denominator() != 1)
								Test += Factor * (it_field->GetDiscreteCharge(DiscreteSymmetry) - DiscreteSymmetry.SuperpotentialCharge/2);
						}

						if (IsR)
							AnomalyPolynomial += Factor * (it_field->GetDiscreteCharge(DiscreteSymmetry) - DiscreteSymmetry.SuperpotentialCharge/2);
						else
							AnomalyPolynomial += Factor * it_field->GetDiscreteCharge(DiscreteSymmetry);
					}
					// end: compute the charges under the discrete ZN symmetries
				}
			}
			// end: run through all left-chiral fields

			unsigned RIndex = 0;
			if (IsR)
			{
				bool first_entry_not_found = true;
				for (k = 1; k < 4; ++k)
				{
					if (DiscreteSymmetry.ChargeOperator[k].numerator() != 0)
					{
						if (first_entry_not_found)
						{
							RIndex = k;
							first_entry_not_found = false;
						}
						else
						{
							if (info)
								(*Print.out) << "\n  " << Print.cbegin << "R symmetry ill-defined." << Print.cend << "\n";
							return false;
						}
					}
				}
				if (first_entry_not_found)
				{
					if (info)
						(*Print.out) << "\n  " << Print.cbegin << "R symmetry ill-defined." << Print.cend << "\n";
					return false;
				}
			}

			if (IsR)
				AnomalyPolynomial -= BetaFunctionCoefficients[RIndex-1][j];

			if (info)
				(*Print.out) << "    A_{Z_" << Order.numerator() << Rstring << " - " << ggf.algebra << " - " << ggf.algebra << "} = " << AnomalyPolynomial << "  ";

			// begin: create the anomalous space-group element
			if (j == 0)
			{
				if (!IsR && (AnomalyPolynomial != 0))
				{
					tmp = AnomalyPolynomial;
					while (tmp < 0)
						tmp += Order;
					while (tmp >= Order)
						tmp -= Order;

					const CSpaceGroupElement &SymmetryGenerator = DiscreteSymmetry.SG_SymmetryGenerator;

					anomalous_element_k += (tmp * rational<int>(SymmetryGenerator.Get_m(),1)).numerator();		//hacking here!!! not optimized for ZMxZNxZK case...
					anomalous_element_l += (tmp * rational<int>(SymmetryGenerator.Get_n(),1)).numerator();

					for (k = 0; k < 6; ++k)
						anomalous_element_n[k] += tmp * SymmetryGenerator.Get_n(k);
				}

				first_AnomalyPolynomial = AnomalyPolynomial;
			}
			// end: create the anomalous space-group element

			if (IsR && (j == 0) && (AnomalyPolynomial != 0))
			{
				if ((ZMxZN && is_integer(Twists[0][RIndex]+0.5) && is_integer(Twists[1][RIndex]+0.5)) || (!ZMxZN && is_integer(Twists[0][RIndex]+0.5)))
				{
					AnomaliesNonUniversal = true;
					if (info)
						(*Print.out) << " R-anomaly not zero in Z2 plane";
				}
			}

			// begin: check universality
			if (!IsR)
			{
				tmp = (12 * AnomalyPolynomial - GravAnomalyPolynomial)/Order;
				if (tmp.denominator() != 1)
				{
					AnomaliesNonUniversal = true;
					if (info)
						(*Print.out) << " Z_N - grav - grav anomaly not universal";
				}

				tmp = (AnomalyPolynomial - first_AnomalyPolynomial)/Order;
				if (tmp.denominator() != 1)
				{
					AnomaliesNonUniversal = true;
					if (info)
						(*Print.out) << " Z_N - G - G anomaly not universal";
				}
			}
			else
			{
				/*tmp = (12 * AnomalyPolynomial - GravAnomalyPolynomial)/Order;
        if (tmp.denominator() != 1)
        {
          cout << tmp << endl;
          cout << "\n  Warning in bool COrbifold::CheckDiscreteAnomaly(...) const: Mixed discrete R - grav. - grav. anomaly is not universal. Return false." << endl;
          //return false;
      }*/

				tmp = (AnomalyPolynomial - first_AnomalyPolynomial)/Order;
				if (tmp.denominator() != 1)
				{
					AnomaliesNonUniversal = true;
					if (info)
						(*Print.out) << " Z_N^R - G - G anomaly not universal";
				}
			}
			// end: check universality
			if (info)
				(*Print.out) << endl;
		}
		if (info)
			(*Print.out) << "    -----------------------------" << endl;
		// end: Z_N (R) - G - G
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	}

	if (number_of_mod_symmetries != 0)
	{
		rational<int> first_delta_GS = 0;
		rational<int> delta_GS = 0;

		for (i = 0; i < number_of_mod_symmetries; ++i)
		{
			const SModularSymmetry &ModularSymmetry = ModularSymmetries[i];

			if (info)
			{
				(*Print.out) << "\n  ";
				Print.PrintModularSymmetry(ModularSymmetry);
				(*Print.out) << endl;
			}

			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// begin: modular - grav - grav
			GravAnomalyPolynomial = 21 + 1 - number_of_U1s;

			// begin: run through all left-chiral, vector and modulus fields
			for( vector<CField>::const_iterator it_field = Fields.begin(); it_field != Fields.end(); ++it_field)
			{
				// vector multiplet
				if (it_field->Multiplet == Vector)
				{
					const RepVector &Dimensions = it_field->Dimensions;

					Factor = 1;
					for (k = 0; k < number_of_factors; ++k)
						Factor *= abs(Dimensions[k].Dimension);

					GravAnomalyPolynomial -= Factor;
				}
				else
					// moduli
					if (it_field->Multiplet == LCModulus)
						GravAnomalyPolynomial += it_field->GetModularWeight(ModularSymmetry, this->Sectors);
					else
						// left-chiral multiplet
						if (it_field->Multiplet == LeftChiral)
						{
							const RepVector &Dimensions = it_field->Dimensions;

							Factor = 1;
							for (k = 0; k < number_of_factors; ++k)
								Factor *= abs(Dimensions[k].Dimension);

							GravAnomalyPolynomial += Factor * (1 + (2 * it_field->GetModularWeight(ModularSymmetry, this->Sectors)));
						}
			}
			GravAnomalyPolynomial /= 24;
			// end: run through all left-chiral, vector and modulus fields
			// end: modular - grav - grav
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			if (info)
			{
				(*Print.out) << "    -----------------------------" << endl;
				(*Print.out) << "    A_{" << ModularSymmetry.Label << " - grav - grav}/24 = " << GravAnomalyPolynomial;

				if (N2_Sector_Exists[ModularSymmetry.Index-1])
					(*Print.out) << " (can not be checked)";
				(*Print.out) << "\n    -----------------------------" << endl;
			}

			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// begin: modular - G - G
			// run through all gauge group factors
			for (j = 0; j < number_of_factors; ++j)
			{
				if (info && (j == VEVConfig.SymmetryGroup.Position_of_and_in_GaugeGroup))
					(*Print.out) << "    -----------------------------" << endl;

				const gaugeGroupFactor<double> &ggf = GaugeGroup.factor[j];

				if (!GaugeIndices.GetQuadraticIndexAdj(ggf, Index))
					return false;

				AnomalyPolynomial = -rational<int>(Index,2);

				// begin: run through all left-chiral fields
				for( vector<CField>::const_iterator it_field = Fields.begin(); it_field != Fields.end(); ++it_field)
				{
					// left-chiral multiplet
					if (it_field->Multiplet == LeftChiral)
					{
						const RepVector &Dimensions = it_field->Dimensions;

						// begin: compute the charges under the modular symmetries
						if (!GaugeIndices.GetQuadraticIndex(ggf, abs(Dimensions[j].Dimension), Index))
							return false;

						if (Index != 0)
						{
							Factor = Index;

							for (k = 0; k < number_of_factors; ++k)
							{
								if (k != j)
									Factor *= abs(Dimensions[k].Dimension);
							}
							Factor /= 2;

							AnomalyPolynomial += Factor * (1 + (2 * it_field->GetModularWeight(ModularSymmetry, this->Sectors)));
						}
						// end: compute the charges under the modular symmetries
					}
				}
				// end: run through all left-chiral fields

				delta_GS = AnomalyPolynomial - BetaFunctionCoefficients[ModularSymmetry.Index-1][j];

				if (j == 0)
					first_delta_GS = delta_GS;

				if (info)
				{
					(*Print.out) << "      A_{" << ModularSymmetry.Label << " - " << ggf.algebra << " - " << ggf.algebra << "} = " << AnomalyPolynomial << "\n";
					(*Print.out) << "      |P_" << ModularSymmetry.Index << "|/|P| b^" << ModularSymmetry.Index << "_{N=2} = " << BetaFunctionCoefficients[ModularSymmetry.Index-1][j] << "\n";
					(*Print.out) << "    delta_GS^" << ModularSymmetry.Index << "(" << ModularSymmetry.Label << ", " << ggf.algebra << ") = " << delta_GS << "  ";
				}

				// begin: check universality
				if (delta_GS != first_delta_GS)
				{
					AnomaliesNonUniversal = true;
					(*Print.out) << " delta_GS not universal";
				}
				if (!N2_Sector_Exists[ModularSymmetry.Index-1] && (delta_GS != GravAnomalyPolynomial))
				{
					AnomaliesNonUniversal = true;
					(*Print.out) << " grav. anomaly not universal";
				}
				if ((j == 0) && (delta_GS != 0))
				{
					if ((ZMxZN && is_integer(Twists[0][ModularSymmetry.Index]+0.5) && is_integer(Twists[1][ModularSymmetry.Index]+0.5)) || (!ZMxZN && is_integer(Twists[0][ModularSymmetry.Index]+0.5)))
					{
						AnomaliesNonUniversal = true;
						(*Print.out) << " delta_GS not zero in Z2 plane";
					}
				}
				// end: check universality
				(*Print.out) << endl;
			}
			if (info)
				(*Print.out) << "    -----------------------------" << endl;
			// end: modular - G - G
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		}
	}

	anomalous_element.Set_m(anomalous_element_k);
	anomalous_element.Set_n(anomalous_element_l);
	anomalous_element.Set_n_alpha(anomalous_element_n);

	if (!VEVConfig.SymmetryGroup.IsFirstU1Anomalous || (number_of_U1s == 0))
	{
		if (AnomaliesNonUniversal)
		{
			if (info)
				(*Print.out) << "\n  Warning: Anomalies are not universal!" << endl;
			return false;
		}

		return true;
	}

	CShiftVector AnomalousShift(Lattice);
	this->OrbifoldGroup.GetShiftVector(anomalous_element, AnomalousShift);
	const CLatticeVector diff = AnomU1 - AnomalousShift;

	bool Error = false;
	if (Lattice == E8xE8)
	{
		if (!diff.From_E8_Lattice(1) || !diff.From_E8_Lattice(2))
			Error = true;
	}
	else
	{
		if (!diff.From_Spin32_Lattice())
			Error = true;
	}

	if (info)
	{
		(*Print.out) << "\n  Anomalous space-group element:\n    ";
		Print.PrintSGElement(anomalous_element);
		(*Print.out) << "\n  The corresponding shift:\n    ";
		Print.PrintRational(AnomalousShift, Lattice);
		(*Print.out) << "\n  Difference to the anomalous U(1) generator:\n    ";
		Print.PrintRational(diff, Lattice);
		(*Print.out) << endl;

		if (Error)
		{
			AnomaliesNonUniversal = true;
			(*Print.out) << "\n  Warning: Embedding of anomalous space-group element is not the anomalous U(1) generator!" << endl;
		}
	}

	if (AnomaliesNonUniversal)
	{
		if (info)
			(*Print.out) << "\n  Warning: Anomalies are not universal!" << endl;
		return false;
	}

	return true;
}



/* ########################################################################################
######   Config_Clear(SConfig &VEVConfig, string &Label)                             ######
######                                                                               ######
######   Version: 19.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) VEVConfig : reset the vev-config "VEVConfig"                             ######
######   2) Label     : new label of the vev-config "VEVConfig"                      ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Reset the vev-config "VEVConfig".                                           ######
######################################################################################## */
bool COrbifold::Config_Clear(SConfig &VEVConfig, string &Label)
{
	VEVConfig.ConfigLabel  = Label;
	VEVConfig.ConfigNumber = 1;

	VEVConfig.InvariantSupercharges.clear();

	VEVConfig.use_Labels = 0;
	VEVConfig.Fields.clear();
	VEVConfig.NamesOfSetsOfFields.clear();
	VEVConfig.SetsOfFields.clear();

	VEVConfig.NamesOfMonomials.clear();
	VEVConfig.SetOfMonomials.clear();

	VEVConfig.SymmetryGroup.D0_FI_term = 0.0;
	VEVConfig.SymmetryGroup.Position_of_and_in_GaugeGroup = 0;

	VEVConfig.SymmetryGroup.observable_sector_GGs.clear();
	VEVConfig.SymmetryGroup.observable_sector_U1s.clear();
	VEVConfig.SymmetryGroup.GGs_AdditionalLabels.assign(VEVConfig.SymmetryGroup.GaugeGroup.factor.size(), "");
	VEVConfig.SymmetryGroup.U1s_AdditionalLabels.assign(VEVConfig.SymmetryGroup.GaugeGroup.u1directions.size(), "");

	VEVConfig.SymmetryGroup.IsFirstU1Anomalous = false;

	return true;
}



/* ########################################################################################
######   Config_SetU1Direction(SConfig &NewU1VEVConfig,  ...) const                  ######
######                                                                               ######
######   Version: 06.09.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) NewU1VEVConfig : vev-config where the basis of U(1) generators will be   ######
######                       changed                                                 ######
######   2) U1Direction    : new U(1) generator                                      ######
######   3) pos_of_U1      : index of new U(1) generator in the list of all U(1)     ######
######                       generators of "NewU1VEVConfig"                          ######
######   output:                                                                     ######
######   return value      : finished successfully?                                  ######
###########################################################################################
######   description:                                                                ######
######   Set the "pos_of_U1"-th U(1) generator of the vev-config "NewU1VEVConfig" to ######
######   "U1Direction" and recalculate the U(1) charges.                             ######
######################################################################################## */
bool COrbifold::Config_SetU1Direction(SConfig &NewU1VEVConfig, const CVector &U1Direction, unsigned pos_of_U1) const
{
	const double prec = 0.0001;

	if (U1Direction.GetSize() != 16)
	{
		cout << "\n  Warning in bool COrbifold::Config_SetU1Direction(...): Generator must have 16 entries. Return false." << endl;
		return false;
	}

	CGaugeGroup &GaugeGroup = NewU1VEVConfig.SymmetryGroup.GaugeGroup;

	const unsigned number_of_U1s = GaugeGroup.u1directions.size();

	if (number_of_U1s == 0)
	{
		cout << "\n  Warning in bool COrbifold::Config_SetU1Direction(...): Gauge group does not contain U(1) factors. Return false." << endl;
		return false;
	}

	unsigned i = 0;
	unsigned j = 0;

	// begin: collect all simple roots
	vector<vector<double> > SimpleRootsOrig;

	const size_t number_of_factors = GaugeGroup.factor.size();
	for (i = 0; i < number_of_factors; ++i)
	{
		const vector<vector<double> > &ggf_SimpleRoots = GaugeGroup.factor[i].simpleroots;
		SimpleRootsOrig.insert(SimpleRootsOrig.end(), ggf_SimpleRoots.begin(), ggf_SimpleRoots.end());
	}
	const size_t s1 = SimpleRootsOrig.size();
	// end: collect all simple roots

	const bool AnomalousU1 = NewU1VEVConfig.SymmetryGroup.IsFirstU1Anomalous;

	if (AnomalousU1)
	{
		if (number_of_U1s == 1)
		{
			cout << "\n  Warning in bool COrbifold::Config_SetU1Direction(...): Gauge group only contains the anomalous U1. Return false." << endl;
			return false;
		}

		if (pos_of_U1 == 0)
		{
			cout << "\n  Warning in bool COrbifold::Config_SetU1Direction(...): first U(1) is anomalous and hence can not be changed. Increasing \"pos_of_U1\" by one." << endl;
			++pos_of_U1;
		}
	}

	if (pos_of_U1 >= number_of_U1s)
	{
		cout << "\n  Warning in bool COrbifold::Config_SetU1Direction(...): index of U(1) out of range. Return false." << endl;
		return false;
	}

	// begin: check if new U1 direction is orthogonal to the simple roots and to the unchanged U1 directions
	CVector unchanged_U1Direction(16);
	bool orthogonal = true;
	for (i = 0; i < pos_of_U1; ++i)
	{
		unchanged_U1Direction = GaugeGroup.u1directions[i];

		if (fabs(unchanged_U1Direction * U1Direction) > prec)
			orthogonal = false;
	}

	if (!orthogonal)
	{
		cout << "\n  Warning in bool COrbifold::Config_SetU1Direction(...): i-th U(1) generator is not orthogonal to the j-th ones, for j < i. Return false." << endl;
		return false;
	}
	double sp = 0.0;
	for (i = 0; orthogonal && (i < s1); ++i)
	{
		vector<double> &SimpleRoot = SimpleRootsOrig[i];

		sp = 0.0;
		for (j = 0; j < 16; ++j)
			sp += SimpleRoot[j] * U1Direction[j];

		if (fabs(sp) > prec)
			orthogonal = false;
	}

	if (!orthogonal)
	{
		cout << "\n  Warning in bool COrbifold::Config_SetU1Direction(...): New U(1) direction is not orthogonal to the simple roots. Return false." << endl;
		return false;
	}
	// end: check if new U1 direction is orthogonal to the simple roots and to the unchanged U1 directions

	vector<doubleVector> &U1_Directions = GaugeGroup.u1directions;
	U1_Directions[pos_of_U1] = U1Direction;

	// the non-anomalous U(1) generators are orthogonal to the simple roots and to the old U(1) generator
	vector<vector<double> > SimpleRoots = SimpleRootsOrig;

	for (i = 0; i < pos_of_U1; ++i)
		SimpleRoots.push_back(GaugeGroup.u1directions[i]);

	SimpleRoots.push_back(U1Direction);

	// create the non-anomalous U1 directions
	if (SimpleRoots.size() == 16)
	{
		// compute the new U(1) charges of all states and create new summaries
		this->Config_ComputeNewU1Charges(NewU1VEVConfig);

		return true;
	}

	vector<vector<double> > new_U1_Directions;
	if (!Find_Basis_Of_Orthogonal_Space(SimpleRoots, this->OrbifoldGroup.GetLattice(), 16, new_U1_Directions))
	{
		cout << "\n  Warning in bool COrbifold::Config_SetU1Direction(...). Could not find new basis for other U1 directions. Return false." << endl;
		return false;
	}

	const size_t s2 = new_U1_Directions.size();

	if (s2 + pos_of_U1 + 1 != number_of_U1s)
	{
		cout << "\n  Warning in bool COrbifold::Config_SetU1Direction(...). Could not find new basis for other U1 directions. Return false." << endl;
		return false;
	}

	for (i = 0; i < s2; ++i)
	{
		if (i+pos_of_U1+1 >= U1_Directions.size())
		{
			cout << "\n  Warning in bool COrbifold::Config_SetU1Direction(...). Index out of range. Return false." << endl;
			return false;
		}

		U1_Directions[i+pos_of_U1+1] = new_U1_Directions[i];
	}

	// compute the new U(1) charges of all states and create new summaries
	this->Config_ComputeNewU1Charges(NewU1VEVConfig);
	return true;
}



/* ########################################################################################
######   Config_ComputeNewU1Charges(SConfig &NewU1VEVConfig) const                   ######
######                                                                               ######
######   Version: 25.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) NewU1VEVConfig : a vev-config with a new basis of U(1) generators        ######
######   output:                                                                     ######
######   return value      : finished successfully?                                  ######
###########################################################################################
######   description:                                                                ######
######   Recalculate the U(1) charges using the U(1) generators of "NewU1VEVConfig"  ######
######   and store the new charges there.                                            ######
######################################################################################## */
bool COrbifold::Config_ComputeNewU1Charges(SConfig &NewU1VEVConfig) const
{
	const size_t s1 = this->Sectors.size();
	size_t s2 = 0;
	size_t s3 = 0;
	unsigned j = 0;
	unsigned k = 0;

	for (unsigned i = 0; i < s1; ++i)
	{
		const CSector &Sector = this->Sectors[i];

		s2 = Sector.GetNumberOfFixedBranes();
		for (j = 0; j < s2; ++j)
		{
			const CFixedBrane &FixedBrane = Sector.GetFixedBrane(j);

			s3 = FixedBrane.GetNumberOfInvariantStates();
			for (k = 0; k < s3; ++k)
			{
				const CState          &State       = FixedBrane.GetInvariantState(k);
				const vector<CVector> &all_Weights = FixedBrane.GetMasslessLeftMover(State.GetLeftMover().GetIndex()).Weights;

				if (!State.RecalculateU1Charges(all_Weights, NewU1VEVConfig))
				{
					cout << "\n  Warning in bool COrbifold::Config_ComputeNewU1Charges(...): Return false." << endl;
					return false;
				}
			}
		}
	}

	return true;
}

