#include <iostream>
#include <cstdlib>

#include "linalg.hpp"
#include "io.hpp"

#include "cfixedpoint.h"
#include "corbifold.h"
#include "globalfunctions.h"

#include "cprint.h"

#define CHECKERROR true

using std::cout;
using std::endl;
using std::exit;
using std::vector;



/* ########################################################################################
######   CFixedBrane()                                                               ######
######                                                                               ######
######   Version: 16.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard constructor of a CFixedBrane object. No content is specified.      ######
######################################################################################## */
CFixedBrane::CFixedBrane()
: Index_SGElement(0)
{
	this->FixedBraneLabel = "";
	this->FixedBrane_CheckStatus = NotChecked;
}



/* ########################################################################################
######   CFixedBrane()                                                               ######
######                                                                               ######
######   Version: 16.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) SGElement       : the constructing element of this fixed point           ######
######   2) Index_SGElement : the index of the constructing element in the list      ######
######                        "Elements" in the class "COrbifoldGroup"               ######
######   3) FixedBraneLabel : the label of this fixed point                          ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Constructor of a CFixedBrane object.                                        ######
######################################################################################## */
CFixedBrane::CFixedBrane(const CSpaceGroupElement &SGElement, const unsigned &Index_SGElement, const string &FixedBraneLabel)
{
	this->Index_SGElement        = Index_SGElement;
	this->SGElement              = SGElement;
	this->FixedBraneLabel        = FixedBraneLabel;
	this->FixedBrane_CheckStatus = NotChecked;
}



/* ########################################################################################
######   ~CFixedBrane()                                                              ######
######                                                                               ######
######   Version: 04.05.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard destructor of a CFixedBrane object.                                ######
######################################################################################## */
CFixedBrane::~CFixedBrane()
{
}



/* ########################################################################################
######   Create(const COrbifoldGroup &OrbifoldGroup, const CSector &Sector, ...)     ######
######                                                                               ######
######   Version: 11.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) OrbifoldGroup     : containes for example the centralizer                ######
######   2) Sector            : contains right-movers and the oscillator excitations ######
######   3) Gamma_Centralizer : a set of orbifold group elements for the gamma       ######
######                          phases                                               ######
######   output:                                                                     ######
######   return value         : finished successfully?                               ######
###########################################################################################
######   description:                                                                ######
######   Creates the left-movers and combines them with the right-movers to create   ######
######   all orbifold invariant twisted states localized at this fixed point.        ######
######################################################################################## */
bool CFixedBrane::Create(const COrbifoldGroup &OrbifoldGroup, const CSector &Sector, const vector<COrbifoldGroupElement> &Gamma_Centralizer)
{

	if  (this->CreateMasslessLeftMover(OrbifoldGroup.GetElement(this->Index_SGElement), Sector.GetLM_Excitations())
			&& this->SortByEigenvalue(OrbifoldGroup, Sector.GetLM_all_Oscillators(), Gamma_Centralizer)
			&& this->CreateStates(Sector.GetRightMovers(), OrbifoldGroup.GetCentralizer(this->Index_SGElement), Gamma_Centralizer)
			&& this->FindSUSYMultiplets(Sector, OrbifoldGroup.GetInvariantSupercharges()))
	{
		if (Sector.TachyonicGetRightMovers().size() != 0) {
			if ( this->TachyonicCreateLeftMover(OrbifoldGroup.GetElement(this->Index_SGElement), Sector.TachyonicGetLM_Excitations())
					&& this->TachyonicSortByEigenvalue(OrbifoldGroup, Sector.GetLM_all_Oscillators(), Gamma_Centralizer)
					&& this->TachyonicCreateStates(Sector.TachyonicGetRightMovers(), OrbifoldGroup.GetCentralizer(this->Index_SGElement), Gamma_Centralizer)
					&& this->TachyonicFindSUSYMultiplets(Sector, OrbifoldGroup.GetInvariantSupercharges()))
			{
				this->FixedBrane_CheckStatus = CheckedAndGood;
				return true;
			}
		}
		this->FixedBrane_CheckStatus = CheckedAndGood;
		return true;
	}
	return false;
}



/* ########################################################################################
######   Reset()                                                                     ######
######                                                                               ######
######   Version: 19.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Clear the content of this fixed point.                                      ######
######################################################################################## */
bool CFixedBrane::Reset()
{
	this->MasslessLeftMovers.clear();
	this->LeftMovers.clear();
	this->InvariantStates.clear();

	this->FixedBrane_CheckStatus = NotChecked;

	return true;
}



/* ########################################################################################
######   CreateMasslessLeftMover(...)                                                ######
######                                                                               ######
######   Version: 12.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) constructing_Element : the COrbifoldGroupElement object containing the   ######
######                             local shift V_loc                                 ######
######   2) Excitations          : vector of S_OscillatorExcitation objects, one     ######
######                             element is oscillator number 0, i.e. no excitation######
######   output:                                                                     ######
######   return value            : finished successfully?                            ######
###########################################################################################
######   description:                                                                ######
######   Solves the equation for massless left-movers and saves the result in the    ######
######   private member variable "MasslessLeftMovers".                               ######
######################################################################################## */
bool CFixedBrane::CreateMasslessLeftMover(const COrbifoldGroupElement &constructing_Element, const vector<S_OscillatorExcitation> &Excitations)
{
	// begin: check the local modular invariance condition
	const double mod_inv = constructing_Element.Twist.OrderOfTwist() * (constructing_Element.Shift.GetSqrTo(16) - constructing_Element.Twist.GetSqrTo(4));
	if (!is_even(mod_inv)) {
		cout << "\n  Warning in bool CFixedBrane::CreateMasslessHalfState(...): local condition "<<constructing_Element.Twist.OrderOfTwist()<<" (V_loc^2 - v_loc^2) = " << mod_inv << " = 0 mod 2 failed. Return false." << endl;
		return false;
	}
	// end: check the local modular invariance condition

	const size_t s1 = Excitations.size();
	if (s1 == 0)
	{
		cout << "\n  Warning in bool CFixedBrane::CreateMasslessHalfState(...): Set of excitations is empty. Return false." << endl;
		return false;
	}
	if (this->MasslessLeftMovers.size() != 0)
	{
		cout << "\n  Warning in bool CFixedBrane::CreateMasslessHalfState(...): Massless left-movers have been created before - now cleared." << endl;
		this->MasslessLeftMovers.clear();
	}

	const CShiftVector    &Shift   = constructing_Element.Shift;
	const SelfDualLattice &Lattice = Shift.Lattice;

	// run through the excitations of the left moving state
	for (unsigned i = 0; i < s1; ++i)
	{
		CMasslessHalfState HalfState(LeftMover, Excitations[i]);

		if (HalfState.SolveMassEquation(Shift, Lattice))
			this->MasslessLeftMovers.push_back(HalfState);
	}
	return true;
}



/* ########################################################################################
######   TachyonicCreateLeftMover(...)                                                ######
######                                                                               ######
######   Version: 12.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) constructing_Element : the COrbifoldGroupElement object containing the   ######
######                             local shift V_loc                                 ######
######   2) Excitations          : vector of S_OscillatorExcitation objects, one     ######
######                             element is oscillator number 0, i.e. no excitation######
######   output:                                                                     ######
######   return value            : finished successfully?                            ######
###########################################################################################
######   description:                                                                ######
######   Solves the equation for tachyonic left-movers and saves the result in the    ######
######   private member variable "TachyonicMassLeftMovers" or "ExcitedTachyonicMassLeftMovers".  ######
######################################################################################## */
bool CFixedBrane::TachyonicCreateLeftMover(const COrbifoldGroupElement &constructing_Element, const vector<vector<S_OscillatorExcitation> > &Excitations)
{
	// begin: check the local modular invariance condition
	const double mod_inv = constructing_Element.Twist.OrderOfTwist() * (constructing_Element.Shift.GetSqrTo(16) - constructing_Element.Twist.GetSqrTo(4));
	if (!is_even(mod_inv)) {
		cout << "\n  Warning in bool CFixedBrane::CreateMasslessHalfState(...): local condition "<<constructing_Element.Twist.OrderOfTwist()<<" (V_loc^2 - v_loc^2) = " << mod_inv << " = 0 mod 2 failed. Return false." << endl;
		return false;
	}
	// end: check the local modular invariance condition

	unsigned i = 0;
	const CShiftVector    &Shift   = constructing_Element.Shift;
	const SelfDualLattice &Lattice = Shift.Lattice;

	// run through the excitations of the left moving state
	const size_t s1 = Excitations.size();
	if (s1==1)															//Only lowest R-moving tachyonic q_i, N_R=0
	{
		const size_t s2 = Excitations[0].size();
		for (i = 0; i < s2; ++i)
		{
			CMasslessHalfState HalfState(LeftMover, Excitations[0][i]);

			if (HalfState.SolveMassEquation(Shift, Lattice))
			{
				this->TachyonicMassLeftMovers.push_back(HalfState);
			}
		}
	}
	else if (s1>1)														//Also excited R-moving tachyonic level
	{
		const size_t s2_0 = Excitations[0].size();
		for (i = 0; i < s2_0; ++i)
		{
			CMasslessHalfState HalfState(LeftMover, Excitations[0][i]);

			if (HalfState.SolveMassEquation(Shift, Lattice))
				this->TachyonicMassLeftMovers.push_back(HalfState);
		}

		const size_t s2_1 = Excitations[1].size();
		for (i = 0; i < s2_1; ++i)
		{
			CMasslessHalfState HalfState(LeftMover, Excitations[1][i]);

			if (HalfState.SolveMassEquation(Shift, Lattice))
				this->ExcitedTachyonicMassLeftMovers.push_back(HalfState);
		}
	}

	return true;
}




/* ########################################################################################  hacking here!!!
######   CreateTachyonLeftMover(...)                                                ######
######                                                                               ######
######   Version: 12.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) constructing_Element : the COrbifoldGroupElement object containing the   ######
######                             local shift V_loc                                 ######
######   2) Excitations          : vector of S_OscillatorExcitation objects, one     ######
######                             element is oscillator number 0, i.e. no excitation######
######   output:                                                                     ######
######   return value            : finished successfully?                            ######
###########################################################################################
######   description:                                                                ######
######   Solves the equation for tachyonic left-movers and checks if the state   	 ######
######   survives the orbifold and GSO projections, if yes it returns false.         ######
######################################################################################## */
/*bool CFixedBrane::CreateTachyonicLeftMover(const unsigned &windex, const COrbifoldGroup &OrbifoldGroup, const COrbifoldGroupElement &constructing_Element, const vector<vector<S_OscillatorExcitation> > &Excitations, const vector<CModedOscillator> &all_tachyonicOscillators, const CTachyonHalfState &TachyonicRightMovers)
{
		const size_t s1 = Excitations.size();

		if (s1 == 0)							//if no level-matched tachyons are in principle possible, return true
		{
			//cout << "\n  Warning in bool CFixedBrane::CreateMasslessHalfState(...): Set of excitations is empty. Return false." << endl;
			return true;
		}
		if (this->TachyonicLeftMovers.size() != 0)
		{
			cout << "\n  Warning in bool CFixedBrane::CreateTachyonicLeftMover(...): Massless left-movers have been created before - now cleared." << endl;
			this->MasslessLeftMovers.clear();
		}

		const CShiftVector    &Shift   = constructing_Element.Shift;
		const SelfDualLattice &Lattice = Shift.Lattice;

		// create negative mass left-movers
		const CTwistVector &v_g = constructing_Element.Twist;			//take twist from the constructing element
		for(unsigned q_index = 0; q_index < s1; ++q_index) {			//index q_index refers to m_R^2(q_i)<0
			const CVector   &RMcharge = TachyonicRightMovers.Weights[q_index];		//take q_i
			const size_t s2 = Excitations[q_index].size();
			for(unsigned i = 0; i < s2; ++i) {								//run through the LM.excitations that are level-matched with right-moving tachyons, index i refers to possible (vacuum energies+N) for given q_i
				if (Excitations[q_index][i].ZeroPointEnergy>0) {			//check if in principle a tachyonic level-matching is possible
					CMasslessHalfState LM_HalfState(LeftMover, Excitations[q_index][i]);
					if (LM_HalfState.SolveMassEquation(Shift, Lattice)) {  		//if there are M_L^2=M_R^2<0
						cout<<"Level-matched tachyon!\n";

						const size_t s3 = LM_HalfState.Weights.size();						//take number of level-matched left-movers

						//compute the orbifold phase
						// Set the precision
						const double prec = 0.0001;
						// begin: define some constants and variables
						const vector<COrbifoldGroupElement> &Centralizer    = OrbifoldGroup.GetCentralizer(this->Index_SGElement);
						const size_t centralizer_size = Centralizer.size();
						unsigned pos1 = 0;
						unsigned pos2 = 0;
						bool     OK   = true;
						double   Eigenvalue = 0;
						// end: define some constants and variables

						// make a copy of the centralizer
						vector<COrbifoldGroupElement> tmp_Centralizer = Centralizer;

						const size_t   number_of_eigenvalues = tmp_Centralizer.size();
						vector<double> LM_Eigenvalues(number_of_eigenvalues, 0.0);
						vector<double> OsciEigenvalues(number_of_eigenvalues, 0.0);

						// begin: compute the transformation of the oscillators
						const vector<unsigned> &OscillatorIndices = Excitations[q_index][i].OscillatorIndices;
						if (OscillatorIndices.size() != 0) {
							for (int j = 0; j < number_of_eigenvalues; ++j) {
								const CTwistVector &SxG_Twist = tmp_Centralizer[j].Twist;
								Eigenvalue = 0.0;
								for (int k=0; k<OscillatorIndices.size(); k++)		//run through the moded LMoscillators of given sector and pic those giving the Excitations[q_index][i] zero-point energy
									Eigenvalue += all_tachyonicOscillators[OscillatorIndices[k]].GetTransformationProperty(SxG_Twist);
								 OsciEigenvalues[j]=Eigenvalue;
							}
						}
						// end: compute the transformation of the oscillators

						// run through all tachyonic p_j half-states
						for (int j= 0; j < s3; ++j)
						{
							const CVector   &LMcharge = LM_HalfState.Weights[j];				//take p_j

							// calculate the scalar-products (p_sh x V_h)
							for (int k = 0; k < number_of_eigenvalues; ++k)
							{
								const COrbifoldGroupElement &SxG_Element = tmp_Centralizer[k];

								Eigenvalue = SxG_Element.Shift * LMcharge;

								// add the contribution from the vacuum phase
								Eigenvalue -= 1.0/2.0 * ((constructing_Element.Shift * SxG_Element.Shift) - (constructing_Element.Twist * SxG_Element.Twist));

								// begin: add the transformation of the oscillators
								if (OscillatorIndices.size() != 0)
									Eigenvalue += OsciEigenvalues[k];
								// end: add the transformation of the oscillators

								RoundDouble(Eigenvalue);
								Eigenvalue -= floor(Eigenvalue);

								LM_Eigenvalues[k] = Eigenvalue;
							}

							// For the left-movers, two eigenvalues with respect to two centralizer elements
							// have to be the same if the elements belong to the same twisted sector.
							// Otherwise there cannot be a right-mover to form an invariant state.
							for (pos1 = 0; OK && (pos1 < centralizer_size); ++pos1)
							{
								const CSpaceGroupElement &Sector1 = Centralizer[pos1].SGElement;
								for (pos2 = pos1+1; OK && (pos2 < centralizer_size); ++pos2)
								{
									const CSpaceGroupElement &Sector2 = Centralizer[pos2].SGElement;

									// do both elements belong to same twisted sector?
									if ((Sector1.Get_m() == Sector2.Get_m()) && (Sector1.Get_n() == Sector2.Get_n()) && (Sector1.Get_k() == Sector2.Get_k()))
									{
										if (fabs(LM_Eigenvalues[pos1] - LM_Eigenvalues[pos2]) > prec)
											return true;													//no tachyonic invariant state possible
									}
								}
							}

							// right mover "tensor" left mover must be invariant under the centralizer
							bool invariant=true;													//assume invariant tachyonic state
							for (int k = 0; k < number_of_eigenvalues; ++k)							//run through centralizer elements
							{
								double tmp = RMcharge*Centralizer[k].Twist;							//q_sh * v_h[k]
								RoundDouble(tmp);
								tmp -= floor(tmp);
								// begin: check that the state is invariant under the constructing element
								if (fabs(LM_Eigenvalues[k] - tmp) > prec)
								{
									invariant=false;
								}
								// end: check that the state is invariant under the constructing element
							}

							if (invariant==true) {													//if orbifold invariant tachyon found
								//cout << "Warning in bool CFixedBrane::CreateTachyonicLeftMover(...): There is at least one orbifold invariant tachyonic state in spectrum. Return false.\n";
								for (int l=0; l<4; l++) {
									cout<<RMcharge[l]<<" ";
								}
								cout << endl;
								for (int l=0; l<16; l++) {
									cout<<LMcharge[l]<<" ";
								}
								cout << endl;
								return true;
							}
						}
					}
				}
			}
		}
		return true;
}*/


/* ########################################################################################
######   FindSUSYMultiplets(const CSector &Sector, ...)                              ######
######                                                                               ######
######   Version: 11.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Sector                : reference to the untwisted or twisted sector     ######
######                              this fixed brane belongs to                      ######
######   2) InvariantSupercharges : list of orbifold invariant supercharges          ######
######   output:                                                                     ######
######   return value             : finished successfully?                           ######
###########################################################################################
######   description:                                                                ######
######   For each CState object of "InvariantStates" analyzes the right-movers to    ######
######   find out which SUSY multiplets are contained for this state.                ######
######################################################################################## */
bool CFixedBrane::FindSUSYMultiplets(const CSector &Sector, const vector<CVector> &InvariantSupercharges)
{
	const size_t NumberOfSUSY = InvariantSupercharges.size();

	// begin: recursive counting
	vector<vector<unsigned> > AllNumbers;		//need initialization for N=0????
	if (NumberOfSUSY == 1)
	{
		vector<unsigned> tmp(1,0);
		AllNumbers.push_back(tmp);
		tmp[0] = 1;  // add the supercharge
		AllNumbers.push_back(tmp);
		tmp[0] = 2;  // subtract the supercharge
		AllNumbers.push_back(tmp);
	}
	else if (NumberOfSUSY > 1) 					//hacking here!
	{
		vector<unsigned> MaxDigits(NumberOfSUSY, 3);
		vector<unsigned> currentNumber(NumberOfSUSY, 0);
		RecursiveCounting(currentNumber, 0, MaxDigits[0], MaxDigits, AllNumbers);
	}
	// end: recursive counting

	unsigned i = 0;
	unsigned j = 0;
	size_t   s1 = 0;
	unsigned Osci_Index = 0;

	const size_t s2 = this->InvariantStates.size();
	for (i = 0; i < s2; ++i)
	{
		CState &State = this->InvariantStates[i];
		State.FindSUSYMultiplets(Sector, InvariantSupercharges, AllNumbers);			//AllNumber is relevant for susy cases

		CVector OsciContribution(4);
		// begin: compute the contribution of the oscillator excitations
		const vector<unsigned> &OscillatorIndices = this->MasslessLeftMovers[State.GetLeftMover().GetIndex()].Excitation.OscillatorIndices;
		s1 = OscillatorIndices.size();
		for (j = 0; j < s1; ++j)
		{
			const CModedOscillator &ModedOscillator = Sector.GetLM_Oscillator(OscillatorIndices[j]);
			Osci_Index = ModedOscillator.GetIndex();				//direction in which oscillator points

			if (Osci_Index < 4)
			{
				if (ModedOscillator.GetComplex())
					OsciContribution[Osci_Index] += 1.0;
				else
					OsciContribution[Osci_Index] -= 1.0;
			}
		}
		State.SetOsciContribution(OsciContribution);
		// end: compute the contribution of the oscillator excitations
	}
	return true;
}



/* ########################################################################################
######   TachyonicFindSUSYMultiplets(const CSector &Sector, ...)                              ######
######                                                                               ######
######   Version: 11.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Sector                : reference to the untwisted or twisted sector     ######
######                              this fixed brane belongs to                      ######
######   2) InvariantSupercharges : list of orbifold invariant supercharges          ######
######   output:                                                                     ######
######   return value             : finished successfully?                           ######
###########################################################################################
######   description:                                                                ######
######   For each CState object of "TachyonicInvariantStates" analyzes the right-movers to    ######
######   find out which SUSY multiplets are contained for this state.                ######
######################################################################################## */
bool CFixedBrane::TachyonicFindSUSYMultiplets(const CSector &Sector, const vector<CVector> &InvariantSupercharges)
{
	const size_t NumberOfSUSY = InvariantSupercharges.size();
	vector<vector<unsigned> > AllNumbers;		//need initialization for N=0????

	if (NumberOfSUSY != 0)
	{
		cout << "\n  Warning in bool CFixedBrane::TachyonicFindSUSYMultiplets(...) : There are no tachyons in SuSy theory. Return true." << endl;
		return false;
	}

	unsigned i = 0;
	unsigned j = 0;
	size_t   s1 = 0;
	unsigned Osci_Index = 0;

	size_t s2 = this->TachyonicInvariantStates.size();

	for (i = 0; i < s2; ++i)
	{
		CState &State = this->TachyonicInvariantStates[i];
		State.TachyonicFindSUSYMultiplets(Sector, InvariantSupercharges, AllNumbers);			//AllNumber is relevant for susy cases

		CVector OsciContribution(4);
		// begin: compute the contribution of the oscillator excitations //needs attention!!!!
		const vector<unsigned> &OscillatorIndices = this->TachyonicMassLeftMovers[State.GetLeftMover().GetIndex()].Excitation.OscillatorIndices;
		s1 = OscillatorIndices.size();
		for (j = 0; j < s1; ++j)
		{
			const CModedOscillator &ModedOscillator = Sector.GetLM_Oscillator(OscillatorIndices[j]);
			Osci_Index = ModedOscillator.GetIndex();				//direction in which oscillator points

			if (Osci_Index < 4)
			{
				if (ModedOscillator.GetComplex())
					OsciContribution[Osci_Index] += 1.0;
				else
					OsciContribution[Osci_Index] -= 1.0;
			}
		}
		State.SetOsciContribution(OsciContribution);
		// end: compute the contribution of the oscillator excitations
	}

	s2 = this->ExcitedTachyonicInvariantStates.size();

	for (i = 0; i < s2; ++i)
	{
		CState &State = this->ExcitedTachyonicInvariantStates[i];
		State.TachyonicFindSUSYMultiplets(Sector, InvariantSupercharges, AllNumbers);			//AllNumber is relevant for susy cases

		CVector OsciContribution(4);
		// begin: compute the contribution of the oscillator excitations //needs attention!!!!
		const vector<unsigned> &OscillatorIndices = this->ExcitedTachyonicMassLeftMovers[State.GetLeftMover().GetIndex()].Excitation.OscillatorIndices;
		s1 = OscillatorIndices.size();
		for (j = 0; j < s1; ++j)
		{
			const CModedOscillator &ModedOscillator = Sector.GetLM_Oscillator(OscillatorIndices[j]);
			Osci_Index = ModedOscillator.GetIndex();				//direction in which oscillator points

			if (Osci_Index < 4)
			{
				if (ModedOscillator.GetComplex())
					OsciContribution[Osci_Index] += 1.0;
				else
					OsciContribution[Osci_Index] -= 1.0;
			}
		}
		State.SetOsciContribution(OsciContribution);
		// end: compute the contribution of the oscillator excitations
	}
	return true;
}



/* ########################################################################################
######   CreateStates(const vector<CHalfState> &RightMovers, ...)                    ######
######                                                                               ######
######   Version: 17.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) RightMovers       : list of all right-moving states of the current       ######
######                          (un-)twisted sector                                  ######
######   2) Centralizer       : list of centralizer elements; the tensor product of  ######
######                          right- and left-movers must be invariant with        ######
######                          respect to these orbifold group elements             ######
######   3) Gamma_Centralizer : for computing the gamma phases                       ######
######   output:                                                                     ######
######   return value         : finished successfully?                               ######
###########################################################################################
######   description:                                                                ######
######   Creates the orbifold invariant states and saves them in "InvariantStates".  ######
######################################################################################## */
bool CFixedBrane::CreateStates(const vector<CHalfState> &RightMovers, const vector<COrbifoldGroupElement> &Centralizer, const vector<COrbifoldGroupElement> &Gamma_Centralizer)
{
	const size_t c1 = Centralizer.size();
	if (c1 == 0)
	{
		cout << "\n  Warning in bool CFixedBrane::CreateStates(...) : Centralizer is empty. Return false." << endl;
		return false;
	}

	// Set the precision
	const double prec = 0.0001;

	// begin: define some constants and variables
	const size_t s1 = this->LeftMovers.size();
	const size_t s2 = RightMovers.size();

	double tmp = 0.0;
	bool   invariant = true;
	bool   check_constructingElement = true;

	unsigned j = 0;
	unsigned k = 0;
	// end: define some constants and variables

	// for a left mover there may be MORE THAN ONE right mover to construct an invariant state!
	for(unsigned i = 0; i < s1; ++i)
	{
		const CHalfState   &LM             = this->LeftMovers[i];
		const doubleVector &LM_eigenvalues = LM.Eigenvalues;

		// search right mover
		for (j = 0; j < s2; ++j)
		{
			const CHalfState   &RM             = RightMovers[j];
			const doubleVector &RM_eigenvalues = RM.Eigenvalues;

			if (RM_eigenvalues.size() != 3)							//hacking here!!!
			{
				cout << "\n  Warning in bool CFixedBrane::CreateStates(...) : Eigenvalues of the right-mover are ill-defined. Return false." << endl;
				return false;
			}

			invariant = true;
			check_constructingElement = true;

			// right mover "tensor" left mover must be invariant under the centralizer
			for (k = 0; invariant && (k < c1); ++k)
			{
				const CSpaceGroupElement &CentralizerElement = Centralizer[k].SGElement;

				tmp = (CentralizerElement.Get_m() * RM_eigenvalues[0]) + (CentralizerElement.Get_n() * RM_eigenvalues[1]) + (CentralizerElement.Get_k() * RM_eigenvalues[2]);
				RoundDouble(tmp);
				tmp -= floor(tmp);

				// begin: check that the state is invariant under the constructing element
				if (check_constructingElement && (CentralizerElement == this->SGElement))
				{
					check_constructingElement = false;

					if (fabs(LM_eigenvalues[k] - tmp) > prec)
					{
						cout << "\n  Warning in bool CFixedBrane::CreateStates(...) : State is not invariant under constructing element. Return false." << endl;
						return false;
					}
				}
				// end: check that the state is invariant under the constructing element

				if (fabs(LM_eigenvalues[k] - tmp) > prec)
					invariant = false;
			}

			if (invariant)
			{
				//create new state
				this->InvariantStates.push_back(CState());
				CState &new_State = this->InvariantStates[this->InvariantStates.size()-1];
				new_State.Create(LM, RM, Gamma_Centralizer);
			}
		}
	}
	return true;
}



/* ########################################################################################
######   TachyonicCreateStates(const vector<CHalfState> &RightMovers, ...)                    ######
######                                                                               ######
######   Version: 17.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) RightMovers       : list of all right-moving states of the current       ######
######                          (un-)twisted sector                                  ######
######   2) Centralizer       : list of centralizer elements; the tensor product of  ######
######                          right- and left-movers must be invariant with        ######
######                          respect to these orbifold group elements             ######
######   3) Gamma_Centralizer : for computing the gamma phases                       ######
######   output:                                                                     ######
######   return value         : finished successfully?                               ######
###########################################################################################
######   description:                                                                ######
######   Creates the orbifold invariant states and saves them in "InvariantStates".  ######
######################################################################################## */
bool CFixedBrane::TachyonicCreateStates(const vector<CHalfState> &RightMovers, const vector<COrbifoldGroupElement> &Centralizer, const vector<COrbifoldGroupElement> &Gamma_Centralizer)
{
	const size_t c1 = Centralizer.size();
	if (c1 == 0)
	{
		cout << "\n  Warning in bool CFixedBrane::CreateStates(...) : Centralizer is empty. Return false." << endl;
		return false;
	}

	// Set the precision
	const double prec = 0.00001;

	// begin: define some constants and variables
	const size_t s1 = this->TachyonicLeftMovers.size();
	const size_t s2 = RightMovers.size();					//s2=1 means only ground tachyonic level, s2>1 means also excited tachyon

	double tmp = 0.0;
	bool   invariant = true;
	bool   check_constructingElement = true;

	unsigned j = 0;
	unsigned k = 0;
	// end: define some constants and variables

	//search for orbifold invariant tachyons from the R-moving lowest tachyonic level
	// for a left mover there may be MORE THAN ONE right mover to construct an invariant state!
	for(unsigned i = 0; i < s1; ++i)
	{
		const CHalfState   &LM             = this->TachyonicLeftMovers[i];
		const doubleVector &LM_eigenvalues = LM.Eigenvalues;

		// search right mover
		for (j = 0; j < 1; ++j)							//assuming first right-mover is tachyonic lowest level
		{
			const CHalfState   &RM             = RightMovers[j];
			const doubleVector &RM_eigenvalues = RM.Eigenvalues;

			if (RM_eigenvalues.size() != 3)							//hacking here!!!
			{
				cout << "\n  Warning in bool CFixedBrane::CreateStates(...) : Eigenvalues of the right-mover are ill-defined. Return false." << endl;
				return false;
			}

			invariant = true;
			check_constructingElement = true;

			// right mover "tensor" left mover must be invariant under the centralizer
			for (k = 0; invariant && (k < c1); ++k)
			{
				const CSpaceGroupElement &CentralizerElement = Centralizer[k].SGElement;

				tmp = (CentralizerElement.Get_m() * RM_eigenvalues[0]) + (CentralizerElement.Get_n() * RM_eigenvalues[1]) + (CentralizerElement.Get_k() * RM_eigenvalues[2]);
				RoundDouble(tmp);
				tmp -= floor(tmp);

				// begin: check that the state is invariant under the constructing element
				if (check_constructingElement && (CentralizerElement == this->SGElement))
				{
					check_constructingElement = false;

					if (fabs(LM_eigenvalues[k] - tmp) > prec)
					{
						cout << "\n  Warning in bool CFixedBrane::TachyonicCreateStates(...) : State is not invariant under constructing element. Return false." << endl;
						return false;
					}
				}
				// end: check that the state is invariant under the constructing element

				if (fabs(LM_eigenvalues[k] - tmp) > prec)
					invariant = false;
			}

			if (invariant)
			{
				//create new state
				this->TachyonicInvariantStates.push_back(CState());
				CState &new_State = this->TachyonicInvariantStates[this->TachyonicInvariantStates.size()-1];
				new_State.Create(LM, RM, Gamma_Centralizer);
			}
		}
	}

	//search for orbifold invariant tachyons from the R-moving excited tachyonic level
	// for a left mover there may be MORE THAN ONE right mover to construct an invariant state!
	if (s2>1)								//There may be two excited R-moving tachyons corresponding to the SAME negative mass level
	{
		const size_t s1 = this->ExcitedTachyonicLeftMovers.size();
		for(unsigned i = 0; i < s1; ++i)
		{
			const CHalfState   &LM             = this->ExcitedTachyonicLeftMovers[i];
			const doubleVector &LM_eigenvalues = LM.Eigenvalues;

			// search right mover
			for (j = 1; j < s2; ++j)							//assuming second right-mover is tachyonic excited level
			{
				const CHalfState   &RM             = RightMovers[j];
				const doubleVector &RM_eigenvalues = RM.Eigenvalues;

				if (RM_eigenvalues.size() != 3)							//hacking here!!!
				{
					cout << "\n  Warning in bool CFixedBrane::CreateStates(...) : Eigenvalues of the right-mover are ill-defined. Return false." << endl;
					return false;
				}

				invariant = true;
				check_constructingElement = true;

				// right mover "tensor" left mover must be invariant under the centralizer
				for (k = 0; invariant && (k < c1); ++k)
				{
					const CSpaceGroupElement &CentralizerElement = Centralizer[k].SGElement;

					tmp = (CentralizerElement.Get_m() * RM_eigenvalues[0]) + (CentralizerElement.Get_n() * RM_eigenvalues[1]) + (CentralizerElement.Get_k() * RM_eigenvalues[2]);
					RoundDouble(tmp);
					tmp -= floor(tmp);

					// begin: check that the state is invariant under the constructing element
					if (check_constructingElement && (CentralizerElement == this->SGElement))
					{
						check_constructingElement = false;

						if (fabs(LM_eigenvalues[k] - tmp) > prec)
						{
							cout << "\n  Warning in bool CFixedBrane::TachyonicCreateStates(...) : State is not invariant under constructing element. Return false." << endl;
							return false;
						}
					}
					// end: check that the state is invariant under the constructing element

					if (fabs(LM_eigenvalues[k] - tmp) > prec)
						invariant = false;
				}

				if (invariant)
				{
					//create new state
					this->ExcitedTachyonicInvariantStates.push_back(CState());
					CState &new_State = this->ExcitedTachyonicInvariantStates[this->ExcitedTachyonicInvariantStates.size()-1];
					new_State.Create(LM, RM, Gamma_Centralizer);
				}
			}
		}
	}
	return true;
}



/* ########################################################################################
######   GetFieldIndices(const vector<CField> &Fields, ...) const                    ######
######                                                                               ######
######   Version: 31.01.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Fields       : vector of SField objects where to look for "Multiplet"    ######
######   2) Multiplet    : the SUSY type (e.g. left-chiral) to look for              ######
######   3) FieldIndices : list of internal field-indices from "Fields" of SUSY type ######
######                     "Multiplet"                                               ######
######   output:                                                                     ######
######   return value    : finished successfully?                                    ######
###########################################################################################
######   description:                                                                ######
######   Finds all fields from "Fields" with SUSY type "Multiplet" and stores the    ######
######   corresponding field indices in "FieldIndices". For example, find all        ######
######   "LeftChiral" fields of this fixed point / fixed brane.                      ######
######################################################################################## */
bool CFixedBrane::GetFieldIndices(const vector<CField> &Fields, const SUSYMultiplet &Multiplet, vector<unsigned> &FieldIndices) const
{
	const size_t s1 = this->InvariantStates.size();
	for (unsigned i = 0; i < s1; ++i)
		this->InvariantStates[i].GetFieldIndices(Fields, Multiplet, FieldIndices);

	return true;
}



/* ########################################################################################
######   GetFieldIndices(const vector<CField> &Fields, ...) const                    ######
######                                                                               ######
######   Version: 31.01.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Fields       : vector of SField objects where to look for "Multiplet"    ######
######   2) Multiplet    : the SUSY type (e.g. left-chiral) to look for              ######
######   3) FieldIndices : list of internal field-indices from "Fields" of SUSY type ######
######                     "Multiplet"                                               ######
######   output:                                                                     ######
######   return value    : finished successfully?                                    ######
###########################################################################################
######   description:                                                                ######
######   Finds all fields from "Fields" with SUSY type "Multiplet" and stores the    ######
######   corresponding field indices in "FieldIndices". For example, find all        ######
######   "LeftChiral" fields of this fixed point / fixed brane.                      ######
######################################################################################## */
bool CFixedBrane::TachyonicGetFieldIndices(const vector<CField> &Fields, const SUSYMultiplet &Multiplet, vector<unsigned> &FieldIndices) const
{
	const size_t s1 = this->TachyonicInvariantStates.size();
	for (unsigned i = 0; i < s1; ++i)
		this->TachyonicInvariantStates[i].GetFieldIndices(Fields, Multiplet, FieldIndices);

	return true;
}



/* ########################################################################################
######   SortByEigenvalue(const COrbifoldGroup &OrbifoldGroup, ...)                  ######
######                                                                               ######
######   Version: 05.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) OrbifoldGroup     : containes the discrete torsion parameters and the    ######
######                          centralizer                                          ######
######   2) all_Oscillators   : the set of all possible oscillator excitations at    ######
######                          this fixed point                                     ######
######   3) Gamma_Centralizer : set of COrbifoldGroupElement objects used to         ######
######                          compute the gamma phases, i.e. the phases between    ######
######                          terms of the linear combinations of identified states######
######   output:                                                                     ######
######   return value         : finished successfully?                               ######
###########################################################################################
######   description:                                                                ######
######   Creates the private member variable "LeftMovers" by sorting the left-movers ######
######   with respect to the centralizer eigenvalues.                                ######
######################################################################################## */
bool CFixedBrane::SortByEigenvalue(const COrbifoldGroup &OrbifoldGroup, const vector<CModedOscillator> &all_Oscillators, const vector<COrbifoldGroupElement> &Gamma_Centralizer)
{
	// Set the precision
	const double prec = 0.001;

	// begin: Set discrete torsion
	const GDT_Parameters &DiscreteTorsion = OrbifoldGroup.GetDiscreteTorsion();

	const rational<int>  &a = DiscreteTorsion.a;
	const rationalVector &b = DiscreteTorsion.b;
	const rationalVector &c = DiscreteTorsion.c;
	const rationalMatrix &d = DiscreteTorsion.d;
	// end: Set discrete torsion

	// begin: define some constants and variables
	const COrbifoldGroupElement         &constr_Element = OrbifoldGroup.GetElement(this->Index_SGElement);
	const vector<COrbifoldGroupElement> &Centralizer    = OrbifoldGroup.GetCentralizer(this->Index_SGElement);
	const size_t centralizer_size = Centralizer.size();
	unsigned i = 0;
	unsigned j = 0;
	unsigned k = 0;
	unsigned l = 0;
	unsigned m = 0;

	unsigned pos1 = 0;
	unsigned pos2 = 0;
	bool     OK   = true;

	double   Eigenvalue = 0;
	// end: define some constants and variables

	const bool ZMxZN = OrbifoldGroup.GetSpaceGroup().IsZMxZN();

	// make a copy of the centralizer
	vector<COrbifoldGroupElement> tmp_Centralizer = Centralizer;

	// begin: define the gamma elements
	if (Gamma_Centralizer.size() != 0)
		tmp_Centralizer.insert(tmp_Centralizer.end(), Gamma_Centralizer.begin(), Gamma_Centralizer.end());
	// end: define the gamma elements

	// begin: mod out freely acting Z2
	if (OrbifoldGroup.UseFreelyActingWL)
	{
		CLatticeElement NullElement;

		if ((constr_Element.SGElement.Get_m() == 0) && (constr_Element.SGElement.Get_n() == 0) && (constr_Element.SGElement.Get_k() == 0) && (constr_Element.SGElement.GetLatticeElement() == NullElement))
		{
			// run through all set of weights
			for (vector<CMasslessHalfState>::iterator MHState = this->MasslessLeftMovers.begin(); MHState != this->MasslessLeftMovers.end(); ++MHState)
			{
				vector<CVector> &Weights = MHState->Weights;

				// get next weight of the current list
				for( vector<CVector>::iterator Weight = Weights.begin();
						Weight != Weights.end(); ++Weight)
				{
					double p_tau = (*Weight) * OrbifoldGroup.FreelyActingWilsonLine;

					if ( fabs(p_tau - roundf(p_tau)) > prec)
						Weights.erase(Weight--);
				}
				if (Weights.size() == 0)
					this->MasslessLeftMovers.erase(MHState--);
			}
		}
	}
	// end: mod out freely acting Z2

	size_t s1 = this->MasslessLeftMovers.size();
	size_t s2 = 0;
	size_t s3 = 0;

	bool with_excitation       = false;
	bool Are_Eigenvalues_Equal = true;
	bool Eigenvalues_NotKnown  = true;

	rational<int>  RationalEigenvalue(0);
	const size_t   number_of_eigenvalues = tmp_Centralizer.size();
	vector<double> OsciEigenvalues(number_of_eigenvalues, 0.0);
	vector<double> Eigenvalues(number_of_eigenvalues, 0.0);

	// run through all massless half-states
	for (i = 0; i < s1; ++i)
	{
		const CMasslessHalfState &MasslessHalfState = this->MasslessLeftMovers[i];
		const vector<CVector>    &Weights           = MasslessHalfState.Weights;

		// begin: compute the transformation of the oscillators
		with_excitation = (MasslessHalfState.Excitation.OscillatorIndices.size() != 0);
		if (with_excitation)
		{
			const vector<unsigned> &OscillatorIndices = MasslessHalfState.Excitation.OscillatorIndices;
			s2 = OscillatorIndices.size();

			for (j = 0; j < number_of_eigenvalues; ++j)
			{
				const CTwistVector &SxG_Twist = tmp_Centralizer[j].Twist;

				Eigenvalue = 0.0;
				for (k = 0; k < s2; ++k)
					Eigenvalue += all_Oscillators[OscillatorIndices[k]].GetTransformationProperty(SxG_Twist);

				OsciEigenvalues[j] = Eigenvalue;
			}
		}
		// end: compute the transformation of the oscillators

		s2 = Weights.size();
		// get next weight of the current half-state
		for (j = 0; j < s2; ++j)
		{
			const CVector &Weight = Weights[j];

			// calculate the scalar-products (p_sh x V_h)
			for (k = 0; k < number_of_eigenvalues; ++k)
			{
				const COrbifoldGroupElement &SxG_Element = tmp_Centralizer[k];

				Eigenvalue = SxG_Element.Shift * Weight;

				// begin: discrete torsion				//hacking here, disable discrete torsion for N0 orbifolds
				// (because of the orbifold core, discrete torsion is implemented for the left movers)
				/*if (DiscreteTorsion.UseGDT)
				{
					rational<int> tmp = 0;

					const CSpaceGroupElement &g = constr_Element.SGElement;
					const CSpaceGroupElement &h = SxG_Element.SGElement;

					if (ZMxZN)
					{
						tmp += a * (g.Get_k() * h.Get_l() - g.Get_l() * h.Get_k());

						for (l = 0; l < 6; ++l)
							tmp += c[l] * (((rational<int>(g.Get_l())) * h.Get_n(l)) - (g.Get_n(l) * (rational<int>(h.Get_l()))));
					}

					for (l = 0; l < 6; ++l)
						tmp += b[l] * (((rational<int>(g.Get_k())) * h.Get_n(l)) - (g.Get_n(l) * (rational<int>(h.Get_k()))));

					for (l = 0; l < 6; ++l)
					{
						const rationalVector &tmp_d     = d[l];
						const rational<int>  &g_n_alpha = g.Get_n(l);
						for (m = 0; m < 6; ++m)
							tmp += tmp_d[m] * g_n_alpha * h.Get_n(m);
					}
					Eigenvalue += (double)tmp.numerator()/(double)tmp.denominator();
				}*/
				// end: discrete torsion

				// add the contribution from the vacuum phase
				Eigenvalue -= 1.0/2.0 * ((constr_Element.Shift * SxG_Element.Shift) - (constr_Element.Twist * SxG_Element.Twist));

				// begin: add the transformation of the oscillators
				if (with_excitation)
					Eigenvalue += OsciEigenvalues[k];
				// end: add the transformation of the oscillators

				RoundDouble(Eigenvalue);
				Eigenvalue -= floor(Eigenvalue);

				Eigenvalues[k] = Eigenvalue;
			}

			// For the left-movers, two eigenvalues with respect to two centralizer elements
			// have to be the same if the elements belong to the same twisted sector.
			// Otherwise there cannot be a right-mover to form an invariant state.
			OK = true;
			for (pos1 = 0; OK && (pos1 < centralizer_size); ++pos1)
			{
				const CSpaceGroupElement &Sector1 = Centralizer[pos1].SGElement;
				for (pos2 = pos1+1; OK && (pos2 < centralizer_size); ++pos2)
				{
					const CSpaceGroupElement &Sector2 = Centralizer[pos2].SGElement;

					// do both elements belong to same twisted sector?
					if ((Sector1.Get_m() == Sector2.Get_m()) && (Sector1.Get_n() == Sector2.Get_n()) && (Sector1.Get_k() == Sector2.Get_k()))
					{
						if (fabs(Eigenvalues[pos1] - Eigenvalues[pos2]) > prec)
							OK = false;
					}
				}
			}

			if (OK)
			{
				// check whether this combination of eigenvalues is known or not
				Eigenvalues_NotKnown = true;

				s3 = this->LeftMovers.size();
				for (k = 0; Eigenvalues_NotKnown && (k < s3); ++k)
				{
					CHalfState &LeftMover = this->LeftMovers[k];
					if (LeftMover.GetIndex() == i)
					{
						// if the set of eigenvalues is known
						Are_Eigenvalues_Equal = true;
						for (l = 0; Are_Eigenvalues_Equal && (l < number_of_eigenvalues); ++l)
						{
							if (fabs(Eigenvalues[l] - LeftMover.Eigenvalues[l]) > prec)
								Are_Eigenvalues_Equal = false;
						}
						// add current weight to the corresponding right-mover
						if (Are_Eigenvalues_Equal)
						{
							LeftMover.Weights.push_back(j);
							Eigenvalues_NotKnown = false;
						}
					}
				}

				// if the current set of eigenvalues is not known
				// create a new LeftMover
				if (Eigenvalues_NotKnown)
				{
					CHalfState new_LeftMover(LeftMover, i);
					this->LeftMovers.push_back(new_LeftMover);
					this->LeftMovers.rbegin()->Eigenvalues = Eigenvalues;
					this->LeftMovers.rbegin()->Weights.push_back(j);
				}
			}
		}
	}
	return true;
}



/* ########################################################################################
######   TachyonicSortByEigenvalue(const COrbifoldGroup &OrbifoldGroup, ...)                  ######
######                                                                               ######
######   Version: 05.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) OrbifoldGroup     : containes the discrete torsion parameters and the    ######
######                          centralizer                                          ######
######   2) all_Oscillators   : the set of all possible oscillator excitations at    ######
######                          this fixed point                                     ######
######   3) Gamma_Centralizer : set of COrbifoldGroupElement objects used to         ######
######                          compute the gamma phases, i.e. the phases between    ######
######                          terms of the linear combinations of identified states######
######   output:                                                                     ######
######   return value         : finished successfully?                               ######
###########################################################################################
######   description:                                                                ######
######   Creates the private member variable "LeftMovers" by sorting the left-movers ######
######   with respect to the centralizer eigenvalues.                                ######
######################################################################################## */
bool CFixedBrane::TachyonicSortByEigenvalue(const COrbifoldGroup &OrbifoldGroup, const vector<CModedOscillator> &all_Oscillators, const vector<COrbifoldGroupElement> &Gamma_Centralizer)
{
	// Set the precision
	const double prec = 0.001;

	// begin: Set discrete torsion
	const GDT_Parameters &DiscreteTorsion = OrbifoldGroup.GetDiscreteTorsion();

	const rational<int>  &a = DiscreteTorsion.a;
	const rationalVector &b = DiscreteTorsion.b;
	const rationalVector &c = DiscreteTorsion.c;
	const rationalMatrix &d = DiscreteTorsion.d;
	// end: Set discrete torsion

	// begin: define some constants and variables
	const COrbifoldGroupElement         &constr_Element = OrbifoldGroup.GetElement(this->Index_SGElement);
	const vector<COrbifoldGroupElement> &Centralizer    = OrbifoldGroup.GetCentralizer(this->Index_SGElement);
	const size_t centralizer_size = Centralizer.size();
	unsigned i = 0;
	unsigned j = 0;
	unsigned k = 0;
	unsigned l = 0;
	unsigned m = 0;

	unsigned pos1 = 0;
	unsigned pos2 = 0;
	bool     OK   = true;

	double   Eigenvalue = 0;
	// end: define some constants and variables

	const bool ZMxZN = OrbifoldGroup.GetSpaceGroup().IsZMxZN();

	// make a copy of the centralizer
	vector<COrbifoldGroupElement> tmp_Centralizer = Centralizer;

	// begin: define the gamma elements
	if (Gamma_Centralizer.size() != 0)
		tmp_Centralizer.insert(tmp_Centralizer.end(), Gamma_Centralizer.begin(), Gamma_Centralizer.end());
	// end: define the gamma elements

	// begin: mod out freely acting Z2
	if (OrbifoldGroup.UseFreelyActingWL)
	{
		CLatticeElement NullElement;

		if ((constr_Element.SGElement.Get_m() == 0) && (constr_Element.SGElement.Get_n() == 0) && (constr_Element.SGElement.Get_k() == 0) && (constr_Element.SGElement.GetLatticeElement() == NullElement))
		{
			// run through all set of weights
			for (vector<CMasslessHalfState>::iterator MHState = this->TachyonicMassLeftMovers.begin(); MHState != this->TachyonicMassLeftMovers.end(); ++MHState)
			{
				vector<CVector> &Weights = MHState->Weights;

				// get next weight of the current list
				for( vector<CVector>::iterator Weight = Weights.begin();
						Weight != Weights.end(); ++Weight)
				{
					double p_tau = (*Weight) * OrbifoldGroup.FreelyActingWilsonLine;

					if ( fabs(p_tau - roundf(p_tau)) > prec)
						Weights.erase(Weight--);
				}
				if (Weights.size() == 0)
					this->TachyonicMassLeftMovers.erase(MHState--);
			}
		}
	}
	// end: mod out freely acting Z2

	size_t s1 = this->TachyonicMassLeftMovers.size();
	size_t s2 = 0;
	size_t s3 = 0;

	bool with_excitation       = false;
	bool Are_Eigenvalues_Equal = true;
	bool Eigenvalues_NotKnown  = true;

	rational<int>  RationalEigenvalue(0);
	const size_t   number_of_eigenvalues = tmp_Centralizer.size();
	vector<double> OsciEigenvalues(number_of_eigenvalues, 0.0);
	vector<double> Eigenvalues(number_of_eigenvalues, 0.0);

	// run through all tachyonic half-states
	for (i = 0; i < s1; ++i)
	{
		const CMasslessHalfState &MasslessHalfState = this->TachyonicMassLeftMovers[i];
		const vector<CVector>    &Weights           = MasslessHalfState.Weights;

		// begin: compute the transformation of the oscillators
		with_excitation = (MasslessHalfState.Excitation.OscillatorIndices.size() != 0);
		if (with_excitation)
		{
			const vector<unsigned> &OscillatorIndices = MasslessHalfState.Excitation.OscillatorIndices;
			s2 = OscillatorIndices.size();

			for (j = 0; j < number_of_eigenvalues; ++j)
			{
				const CTwistVector &SxG_Twist = tmp_Centralizer[j].Twist;

				Eigenvalue = 0.0;
				for (k = 0; k < s2; ++k)
					Eigenvalue += all_Oscillators[OscillatorIndices[k]].GetTransformationProperty(SxG_Twist);

				OsciEigenvalues[j] = Eigenvalue;
			}
		}
		// end: compute the transformation of the oscillators

		s2 = Weights.size();
		// get next weight of the current half-state
		for (j = 0; j < s2; ++j)
		{
			const CVector &Weight = Weights[j];

			// calculate the scalar-products (p_sh x V_h)
			for (k = 0; k < number_of_eigenvalues; ++k)
			{
				const COrbifoldGroupElement &SxG_Element = tmp_Centralizer[k];
				//cout<<SxG_Element.Twist<<endl;

				Eigenvalue = SxG_Element.Shift * Weight;

				// add the contribution from the vacuum phase
				Eigenvalue -= 1.0/2.0 * ((constr_Element.Shift * SxG_Element.Shift) - (constr_Element.Twist * SxG_Element.Twist));

				// begin: add the transformation of the oscillators
				if (with_excitation)
					Eigenvalue += OsciEigenvalues[k];
				// end: add the transformation of the oscillators

				RoundDouble(Eigenvalue);
				Eigenvalue -= floor(Eigenvalue);

				Eigenvalues[k] = Eigenvalue;
			}

			// For the left-movers, two eigenvalues with respect to two centralizer elements
			// have to be the same if the elements belong to the same twisted sector.
			// Otherwise there cannot be a right-mover to form an invariant state.
			OK = true;
			for (pos1 = 0; OK && (pos1 < centralizer_size); ++pos1)
			{
				const CSpaceGroupElement &Sector1 = Centralizer[pos1].SGElement;
				for (pos2 = pos1+1; OK && (pos2 < centralizer_size); ++pos2)
				{
					const CSpaceGroupElement &Sector2 = Centralizer[pos2].SGElement;

					// do both elements belong to same twisted sector?
					if ((Sector1.Get_m() == Sector2.Get_m()) && (Sector1.Get_n() == Sector2.Get_n()) && (Sector1.Get_k() == Sector2.Get_k()))
					{
						if (fabs(Eigenvalues[pos1] - Eigenvalues[pos2]) > prec)
							OK = false;
					}
				}
			}

			if (OK)
			{
				// check whether this combination of eigenvalues is known or not
				Eigenvalues_NotKnown = true;

				s3 = this->TachyonicLeftMovers.size();
				for (k = 0; Eigenvalues_NotKnown && (k < s3); ++k)
				{
					CHalfState &LeftMover = this->TachyonicLeftMovers[k];
					if (LeftMover.GetIndex() == i)
					{
						// if the set of eigenvalues is known
						Are_Eigenvalues_Equal = true;
						for (l = 0; Are_Eigenvalues_Equal && (l < number_of_eigenvalues); ++l)
						{
							if (fabs(Eigenvalues[l] - LeftMover.Eigenvalues[l]) > prec)
								Are_Eigenvalues_Equal = false;
						}
						// add current weight to the corresponding right-mover
						if (Are_Eigenvalues_Equal)
						{
							LeftMover.Weights.push_back(j);
							Eigenvalues_NotKnown = false;
						}
					}
				}

				// if the current set of eigenvalues is not known
				// create a new LeftMover
				if (Eigenvalues_NotKnown)
				{
					CHalfState new_LeftMover(LeftMover, i);
					this->TachyonicLeftMovers.push_back(new_LeftMover);
					this->TachyonicLeftMovers.rbegin()->Eigenvalues = Eigenvalues;
					this->TachyonicLeftMovers.rbegin()->Weights.push_back(j);
				}
			}
		}
	}

	s1 = this->ExcitedTachyonicMassLeftMovers.size();
	s2 = 0;
	s3 = 0;

	//rational<int>  RationalEigenvalue(0);							//why do we need that??
	fill(OsciEigenvalues.begin(), OsciEigenvalues.end(), 0.0);
	fill(Eigenvalues.begin(), Eigenvalues.end(), 0.0);

	// run through all massless half-states
	for (i = 0; i < s1; ++i)
	{
		const CMasslessHalfState &MasslessHalfState = this->ExcitedTachyonicMassLeftMovers[i];
		const vector<CVector>    &Weights           = MasslessHalfState.Weights;

		// begin: compute the transformation of the oscillators
		with_excitation = (MasslessHalfState.Excitation.OscillatorIndices.size() != 0);
		if (with_excitation)
		{
			const vector<unsigned> &OscillatorIndices = MasslessHalfState.Excitation.OscillatorIndices;
			s2 = OscillatorIndices.size();

			for (j = 0; j < number_of_eigenvalues; ++j)
			{
				const CTwistVector &SxG_Twist = tmp_Centralizer[j].Twist;

				Eigenvalue = 0.0;
				for (k = 0; k < s2; ++k)
					Eigenvalue += all_Oscillators[OscillatorIndices[k]].GetTransformationProperty(SxG_Twist);

				OsciEigenvalues[j] = Eigenvalue;
			}
		}
		// end: compute the transformation of the oscillators

		s2 = Weights.size();
		// get next weight of the current half-state
		for (j = 0; j < s2; ++j)
		{
			const CVector &Weight = Weights[j];

			// calculate the scalar-products (p_sh x V_h)
			for (k = 0; k < number_of_eigenvalues; ++k)
			{
				const COrbifoldGroupElement &SxG_Element = tmp_Centralizer[k];
				//cout<<SxG_Element.Twist<<endl;

				Eigenvalue = SxG_Element.Shift * Weight;

				// add the contribution from the vacuum phase
				Eigenvalue -= 1.0/2.0 * ((constr_Element.Shift * SxG_Element.Shift) - (constr_Element.Twist * SxG_Element.Twist));

				// begin: add the transformation of the oscillators
				if (with_excitation)
					Eigenvalue += OsciEigenvalues[k];
				// end: add the transformation of the oscillators

				RoundDouble(Eigenvalue);
				Eigenvalue -= floor(Eigenvalue);

				Eigenvalues[k] = Eigenvalue;
			}

			// For the left-movers, two eigenvalues with respect to two centralizer elements
			// have to be the same if the elements belong to the same twisted sector.
			// Otherwise there cannot be a right-mover to form an invariant state.
			OK = true;
			for (pos1 = 0; OK && (pos1 < centralizer_size); ++pos1)
			{
				const CSpaceGroupElement &Sector1 = Centralizer[pos1].SGElement;
				for (pos2 = pos1+1; OK && (pos2 < centralizer_size); ++pos2)
				{
					const CSpaceGroupElement &Sector2 = Centralizer[pos2].SGElement;

					// do both elements belong to same twisted sector?
					if ((Sector1.Get_m() == Sector2.Get_m()) && (Sector1.Get_n() == Sector2.Get_n()) && (Sector1.Get_k() == Sector2.Get_k()))
					{
						if (fabs(Eigenvalues[pos1] - Eigenvalues[pos2]) > prec)
							OK = false;
					}
				}
			}

			if (OK)
			{
				// check whether this combination of eigenvalues is known or not
				Eigenvalues_NotKnown = true;

				s3 = this->ExcitedTachyonicLeftMovers.size();
				for (k = 0; Eigenvalues_NotKnown && (k < s3); ++k)
				{
					CHalfState &LeftMover = this->ExcitedTachyonicLeftMovers[k];
					if (LeftMover.GetIndex() == i)
					{
						// if the set of eigenvalues is known
						Are_Eigenvalues_Equal = true;
						for (l = 0; Are_Eigenvalues_Equal && (l < number_of_eigenvalues); ++l)
						{
							if (fabs(Eigenvalues[l] - LeftMover.Eigenvalues[l]) > prec)
								Are_Eigenvalues_Equal = false;
						}
						// add current weight to the corresponding right-mover
						if (Are_Eigenvalues_Equal)
						{
							LeftMover.Weights.push_back(j);
							Eigenvalues_NotKnown = false;
						}
					}
				}

				// if the current set of eigenvalues is not known
				// create a new LeftMover
				if (Eigenvalues_NotKnown)
				{
					CHalfState new_LeftMover(LeftMover, i);
					this->ExcitedTachyonicLeftMovers.push_back(new_LeftMover);
					this->ExcitedTachyonicLeftMovers.rbegin()->Eigenvalues = Eigenvalues;
					this->ExcitedTachyonicLeftMovers.rbegin()->Weights.push_back(j);
				}
			}
		}
	}
	return true;
}




/* ########################################################################################
######   AccessInvariantState(const unsigned &i)                                     ######
######                                                                               ######
######   Version: 16.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) i         : index i = 0,1,...                                            ######
######   output:                                                                     ######
######   return value : the "i"-th CState object from "InvariantStates"              ######
###########################################################################################
######   description:                                                                ######
######   Get access to the i-th invariant state of this fixed point.                 ######
######################################################################################## */
CState &CFixedBrane::AccessInvariantState(const unsigned &i)
{
#ifdef CHECKERROR
	if (i >= this->InvariantStates.size())
	{
		cout << "\n  Warning in CState &CFixedBrane::AccessInvariantState(...) : Index i out of range. Set i = 0." << endl;
		return this->InvariantStates[0];
	}
#endif

	return this->InvariantStates[i];
}



/* ########################################################################################
######   TachyonicAccessInvariantState(const unsigned &i)                                     ######
######                                                                               ######
######   Version: 16.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) i         : index i = 0,1,...                                            ######
######   output:                                                                     ######
######   return value : the "i"-th CState object from "InvariantStates"              ######
###########################################################################################
######   description:                                                                ######
######   Get access to the i-th invariant state of this fixed point.                 ######
######################################################################################## */
CState &CFixedBrane::TachyonicAccessInvariantState(const unsigned &i)
{
#ifdef CHECKERROR
	if (i >= this->TachyonicInvariantStates.size())
	{
		cout << "\n  Warning in CState &CFixedBrane::TachyonicAccessInvariantState(...) : Index i out of range. Set i = 0." << endl;
		return this->TachyonicInvariantStates[0];
	}
#endif

	return this->TachyonicInvariantStates[i];
}



/* ########################################################################################
######   ExcitedTachyonicAccessInvariantState(const unsigned &i)                                     ######
######                                                                               ######
######   Version: 16.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) i         : index i = 0,1,...                                            ######
######   output:                                                                     ######
######   return value : the "i"-th CState object from "InvariantStates"              ######
###########################################################################################
######   description:                                                                ######
######   Get access to the i-th invariant state of this fixed point.                 ######
######################################################################################## */
CState &CFixedBrane::ExcitedTachyonicAccessInvariantState(const unsigned &i)
{
#ifdef CHECKERROR
	if (i >= this->ExcitedTachyonicInvariantStates.size())
	{
		cout << "\n  Warning in CState &CFixedBrane::ExcitedTachyonicAccessInvariantState(...) : Index i out of range. Set i = 0." << endl;
		return this->ExcitedTachyonicInvariantStates[0];
	}
#endif

	return this->ExcitedTachyonicInvariantStates[i];
}




/* ########################################################################################
######   AccessInvariantState(const unsigned &i) const                               ######
######                                                                               ######
######   Version: 16.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : the CSpaceGroupElement object of this fixed point            ######
###########################################################################################
######   description:                                                                ######
######   Get (constant) access to the constructing element of this fixed point.      ######
######################################################################################## */
const CSpaceGroupElement &CFixedBrane::GetSGElement() const
{
#ifdef CHECKERROR
	if (this->FixedBrane_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const CSpaceGroupElement &CFixedBrane::GetSGElement() const : Fixed brane ill-defined." << endl;
#endif

	return this->SGElement;
}



/* ########################################################################################
######   Getconstructing_Element() const                                             ######
######                                                                               ######
######   Version: 16.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : the CSpaceGroupElement object of this fixed point            ######
###########################################################################################
######   description:                                                                ######
######   Get (constant) access to the index of the constructing element.             ######
######################################################################################## */
const unsigned &CFixedBrane::Getconstructing_Element() const
{
#ifdef CHECKERROR
	if (this->FixedBrane_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const unsigned &CFixedBrane::Getconstructing_Element() const : Fixed brane ill-defined." << endl;
#endif

	return this->Index_SGElement;
}



/* ########################################################################################
######   GetFixedBraneLabel() const                                                  ######
######                                                                               ######
######   Version: 16.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : the label of this fixed point                                ######
###########################################################################################
######   description:                                                                ######
######   Get (constant) access to the label of the constructing element.             ######
######################################################################################## */
const string &CFixedBrane::GetFixedBraneLabel() const
{
#ifdef CHECKERROR
	if (this->FixedBrane_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const string &CFixedBrane::GetFixedBraneLabel() const : Fixed brane ill-defined." << endl;
#endif

	return this->FixedBraneLabel;
}



/* ########################################################################################
######   GetMasslessLeftMover(const unsigned &i) const                               ######
######                                                                               ######
######   Version: 16.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) i         : index i = 0,1,...                                            ######
######   output:                                                                     ######
######   return value : the "i"-th CMasslessHalfState object from                    ######
######                  "MasslessLeftMovers".                                        ######
###########################################################################################
######   description:                                                                ######
######   Get (constant) access to the i-th massless left-mover of this fixed point.  ######
######################################################################################## */
const CMasslessHalfState &CFixedBrane::GetMasslessLeftMover(const unsigned &i) const
{
#ifdef CHECKERROR
	if (this->FixedBrane_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const CMasslessHalfState &CFixedBrane::GetMasslessLeftMover(...) const : Fixed brane ill-defined." << endl;

	if (i >= this->MasslessLeftMovers.size())
	{
		cout << "\n  Warning in const CMasslessHalfState &CFixedBrane::GetMasslessLeftMover(...) const : Index i out of range. Set i = 0." << endl;
		return this->MasslessLeftMovers[0];
	}
#endif

	return this->MasslessLeftMovers[i];
}



/* ########################################################################################
######   TachyonicGetMassLeftMover(const unsigned &i) const                               ######
######                                                                               ######
######   Version: 16.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) i         : index i = 0,1,...                                            ######
######   output:                                                                     ######
######   return value : the "i"-th CMasslessHalfState object from                    ######
######                  "MasslessLeftMovers".                                        ######
###########################################################################################
######   description:                                                                ######
######   Get (constant) access to the i-th massless left-mover of this fixed point.  ######
######################################################################################## */
const CMasslessHalfState &CFixedBrane::TachyonicGetMassLeftMover(const unsigned &i) const
{
#ifdef CHECKERROR
	if (this->FixedBrane_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const CMasslessHalfState &CFixedBrane::TachyonicGetMassLeftMover(...) const : Fixed brane ill-defined." << endl;

	if (i >= this->TachyonicMassLeftMovers.size())
	{
		cout << "\n  Warning in const CMasslessHalfState &CFixedBrane::TachyonicGetMassLeftMover(...) const : Index i out of range. Set i = 0." << endl;
		return this->TachyonicMassLeftMovers[0];
	}
#endif

	return this->TachyonicMassLeftMovers[i];
}



/* ########################################################################################
######   ExcitedTachyonicGetMassLeftMover(const unsigned &i) const                               ######
######                                                                               ######
######   Version: 16.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) i         : index i = 0,1,...                                            ######
######   output:                                                                     ######
######   return value : the "i"-th CMasslessHalfState object from                    ######
######                  "MasslessLeftMovers".                                        ######
###########################################################################################
######   description:                                                                ######
######   Get (constant) access to the i-th massless left-mover of this fixed point.  ######
######################################################################################## */
const CMasslessHalfState &CFixedBrane::ExcitedTachyonicGetMassLeftMover(const unsigned &i) const
{
#ifdef CHECKERROR
	if (this->FixedBrane_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const CMasslessHalfState &CFixedBrane::ExcitedTachyonicGetMassLeftMover(...) const : Fixed brane ill-defined." << endl;

	if (i >= this->ExcitedTachyonicMassLeftMovers.size())
	{
		cout << "\n  Warning in const CMasslessHalfState &CFixedBrane::ExcitedTachyonicGetMassLeftMover(...) const : Index i out of range. Set i = 0." << endl;
		return this->ExcitedTachyonicMassLeftMovers[0];
	}
#endif

	return this->ExcitedTachyonicMassLeftMovers[i];
}



/* ########################################################################################
######   GetLeftMover(const unsigned &i) const                                       ######
######                                                                               ######
######   Version: 16.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) i         : index i = 0,1,...                                            ######
######   output:                                                                     ######
######   return value : the "i"-th CHalfState object from "LeftMovers"               ######
###########################################################################################
######   description:                                                                ######
######   Get (constant) access to the i-th left-mover of this fixed point.           ######
######################################################################################## */
const CHalfState &CFixedBrane::GetLeftMover(const unsigned &i) const
{
#ifdef CHECKERROR
	if (this->FixedBrane_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const CHalfState &CFixedBrane::GetLeftMover(...) const : Fixed brane ill-defined." << endl;

	if (i >= this->LeftMovers.size())
	{
		cout << "\n  Warning in const CHalfState &CFixedBrane::GetLeftMover(...) const : Index i out of range. Set i = 0." << endl;
		return this->LeftMovers[0];
	}
#endif

	return this->LeftMovers[i];
}



/* ########################################################################################
######   GetInvariantState(const unsigned &i) const                                  ######
######                                                                               ######
######   Version: 16.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) i         : index i = 0,1,...                                            ######
######   output:                                                                     ######
######   return value : the "i"-th CState object from "InvariantStates"              ######
###########################################################################################
######   description:                                                                ######
######   Get (constant) access to the i-th invariant state of this fixed point.      ######
######################################################################################## */
const CState &CFixedBrane::GetInvariantState(const unsigned &i) const
{
#ifdef CHECKERROR
	if (this->FixedBrane_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const CState &CFixedBrane::GetInvariantState(...) const : Fixed brane ill-defined." << endl;

	if (i >= this->InvariantStates.size())
	{
		cout << "\n  Warning in const CState &CFixedBrane::GetInvariantState(...) const : Index i out of range. Set i = 0." << endl;
		return this->InvariantStates[0];
	}
#endif

	return this->InvariantStates[i];
}



/* ########################################################################################
######   GetInvariantState(const unsigned &i) const                                  ######
######                                                                               ######
######   Version: 16.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) i         : index i = 0,1,...                                            ######
######   output:                                                                     ######
######   return value : the "i"-th CState object from "InvariantStates"              ######
###########################################################################################
######   description:                                                                ######
######   Get (constant) access to the i-th invariant state of this fixed point.      ######
######################################################################################## */
const CState &CFixedBrane::TachyonicGetInvariantState(const unsigned &i) const
{
#ifdef CHECKERROR
	if (this->FixedBrane_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const CState &CFixedBrane::TachyonicGetInvariantState(...) const : Fixed brane ill-defined." << endl;

	if (i >= this->TachyonicInvariantStates.size())
	{
		cout << "\n  Warning in const CState &CFixedBrane::TachyonicGetInvariantState(...) const : Index i out of range. Set i = 0." << endl;
		return this->TachyonicInvariantStates[0];
	}
#endif

	return this->TachyonicInvariantStates[i];
}



/* ########################################################################################
######   GetInvariantState(const unsigned &i) const                                  ######
######                                                                               ######
######   Version: 16.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) i         : index i = 0,1,...                                            ######
######   output:                                                                     ######
######   return value : the "i"-th CState object from "InvariantStates"              ######
###########################################################################################
######   description:                                                                ######
######   Get (constant) access to the i-th invariant state of this fixed point.      ######
######################################################################################## */
const CState &CFixedBrane::ExcitedTachyonicGetInvariantState(const unsigned &i) const
{
#ifdef CHECKERROR
	if (this->FixedBrane_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const CState &CFixedBrane::ExcitedTachyonicGetInvariantState(...) const : Fixed brane ill-defined." << endl;

	if (i >= this->ExcitedTachyonicInvariantStates.size())
	{
		cout << "\n  Warning in const CState &CFixedBrane::ExcitedTachyonicGetInvariantState(...) const : Index i out of range. Set i = 0." << endl;
		return this->ExcitedTachyonicInvariantStates[0];
	}
#endif

	return this->ExcitedTachyonicInvariantStates[i];
}

