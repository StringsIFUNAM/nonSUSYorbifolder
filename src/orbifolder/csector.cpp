#include "csector.h"
#include "globalfunctions.h"
#include "cprint.h"

#include <vector>
#include <iostream>

#define CHECKERROR true

using std::cout;
using std::endl;
using std::exit;
using std::vector;



/* ########################################################################################
######   CSector(...)                                                                ######
######                                                                               ######
######   Version: 12.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard constructor of a CSector object. No content is specified. The      ######
######   check variable "Sector_CheckStatus" is set to "NotChecked".                 ######
######################################################################################## */
CSector::CSector()
{
	this->Sector_CheckStatus = NotChecked;
}



/* ########################################################################################
######   ~CSector()                                                                  ######
######                                                                               ######
######   Version: 02.06.2006                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard destructor of a CSector object.                                    ######
######################################################################################## */
CSector::~CSector()
{
}



/* ########################################################################################
######   Create(unsigned k, unsigned l, ...)                                         ######
######                                                                               ######
######   Version: 12.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1,2) m,n,k              : specify the T_{m,n,k} twisted sector.             ######
######   3) constructing_Twist : the twist v_g = kv_1 + lv_2                         ######
######   4) ZMxZN_Twists       : v_1 and v_2 where v_2 = (0^4) for Z_M orbifolds     ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Creates the content of this CSector object by calling the functions:        ######
######     1) "CreateOscillators": creates all (left- and right-moving) oscillators  ######
######        of this sector.                                                        ######
######     2) "CreateOscillatorExcitations": find all (left- and right-moving)       ###### 
######        oscillator excitations such that the equations for massless left- or   ######
######        right-movers still might find solutions.                               ######
######     3) "CreateMasslessRightMover": find all solutions to the equation of      ######
######        massless right-movers 1/2 (q+v_g)^2 - 1/2 + delta c = 0.               ######
######     4) "SortByEigenvalue": sort the solutions (q+v_g) with respect to their   ######
######        eigenvalues under transformations with the twists v_1 and v_2.         ######
######################################################################################## */
bool CSector::Create(unsigned m, unsigned n , unsigned k, const CTwistVector &constructing_Twist, const vector<CTwistVector> &ZMxZNxZK_Twists)
{
	this->Sector.push_back(m);
	this->Sector.push_back(n);
	this->Sector.push_back(k);

	this->Twist = constructing_Twist;

	if (!this->CreateOscillators(ZMxZNxZK_Twists)
			|| !this->CreateOscillatorExcitations()
			|| !this->CreateMasslessRightMover()
			|| !this->SortByEigenvalue(ZMxZNxZK_Twists) )
	{
		cout << "\n  Warning in bool CSector::CreateSector(...) : Cannot create sector. Return false." << endl;
		this->Sector_CheckStatus = CheckedAndFailed;
		return false;
	}
	//hacking here to create tachyonic right moving states for this sector and LM_excitation for level-matching check
	if (m==1) {
		if ( !this->CreateTachyonicRightMover()
				|| !this->TachyonicSortByEigenvalue(ZMxZNxZK_Twists))
		{
			cout << "\n  Warning in bool CSector::CreateSector(...) : Cannot create sector. Return false." << endl;
			this->Sector_CheckStatus = CheckedAndFailed;
			return false;
		}
	}

	CVector Twist0(4);
	CVector Twist1(4);
	CVector Twist2(4);
	for (int i=0; i<4; i++) {
		Twist0[i]=m*ZMxZNxZK_Twists[0][i];
		Twist1[i]=n*ZMxZNxZK_Twists[1][i];
		Twist2[i]=k*ZMxZNxZK_Twists[2][i];
	}
	CTwistVector Twist;
	Twist.operator=(Twist0);
	Twists.push_back(Twist);			//hacking here!!! to seperate twists
	Twist.operator=(Twist1);
	Twists.push_back(Twist);
	Twist.operator=(Twist2);
	Twists.push_back(Twist);

	this->Sector_CheckStatus = CheckedAndGood;
	return true;
}



/* ########################################################################################
######   AddFixedBrane(const CFixedBrane &FixedBrane)                                ######
######                                                                               ######
######   Version: 16.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) SGElement       : the constructing space group element                   ######
######                        (\theta^k \omega^l, n_\alpha e_\alpha)                 ######
######   2) Index_SGElement : position of the constructing element in the member     ######
######                        variable "Elements" in class "OrbifoldGroup"           ######
######   3) FixedBraneLabel : a label for the new fixed point / fixed brane          ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Adds a new, empty fixed point / fixed brane to this twisted sector.         ######
######################################################################################## */
bool CSector::AddFixedBrane(const CSpaceGroupElement &SGElement, const unsigned &Index_SGElement, const string &FixedBraneLabel)
{
#ifdef CHECKERROR
	if ((SGElement.Get_m() != this->Sector[0]) || (SGElement.Get_n() != this->Sector[1]) || (SGElement.Get_k() != this->Sector[2]))
	{
		cout << "\n  Warning in bool CSector::AddFixedBrane(...) : Fixed brane does not belong to this sector. Return false." << endl;
		return false;
	}
#endif

	CFixedBrane FixedBrane(SGElement, Index_SGElement, FixedBraneLabel);
	this->FixedBranes.push_back(FixedBrane);
	return true;
}




/* ########################################################################################
######   CreateOscillators(const vector<CTwistVector> &ZMxZN_Twists)                 ######
######                                                                               ######
######   Version: 12.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) ZMxZN_Twists: v_1 and v_2 where v_2 = (0^4) for Z_M orbifolds            ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Creates all (left- and right-moving) oscillators of this sector and stores  ######
######   the result in the member variables "LM_all_Oscillators" and                 ###### 
######   "RM_all_Oscillators", respectively.                                         ######
######################################################################################## */
bool CSector::CreateOscillators(const vector<CTwistVector> &ZMxZNxZK_Twists)
{
	if ((this->LM_all_Oscillators.size() != 0) || (this->RM_all_Oscillators.size() != 0))
	{
		cout << "\n  Warning in bool CSector::CreateOscillators() : Sets of oscillators are not empty. Now cleared." << endl;
		this->LM_all_Oscillators.clear();
		this->RM_all_Oscillators.clear();
	}

	// untwisted sector
	if ((this->Sector[0] == 0) && (this->Sector[1] == 0) && (this->Sector[2] == 0))		//hacking here!!!
	{
		unsigned i = 0;
		unsigned j = 0;
		vector<unsigned> PossibleOscillatorIndices;

		for (i = 0; i < ZMxZNxZK_Twists.size(); ++i)
		{
			const CTwistVector &Twist = ZMxZNxZK_Twists[i];

			for (j = 0; j < 4; ++j)							//hacking here to include gravity!!!!
			{
				if (find(PossibleOscillatorIndices.begin(), PossibleOscillatorIndices.end(), j) == PossibleOscillatorIndices.end())
					PossibleOscillatorIndices.push_back(j);
			}
		}

		const size_t s1 = PossibleOscillatorIndices.size();
		for (i = 0; i < s1; ++i)
		{
			unsigned Index = PossibleOscillatorIndices[i];

			CModedOscillator LM_alpha(LeftMover, Index, false, this->Twist, -1);
			if (LM_alpha.GetLadderOperatorType() != CreationOperator)
			{
				cout << "\n  Warning in bool CSector::CreateOscillators(...) : Untwisted oscillators failed. Return false." << endl;
				return false;
			}
			this->LM_all_Oscillators.push_back(LM_alpha);

			CModedOscillator LM_alpha_bar(LeftMover, Index, true, this->Twist, -1);
			if (LM_alpha_bar.GetLadderOperatorType() != CreationOperator)
			{
				cout << "\n  Warning in bool CSector::CreateOscillators(...) : Untwisted oscillators failed. Return false." << endl;
				return false;
			}
			this->LM_all_Oscillators.push_back(LM_alpha_bar);
		}

		/*
    // gauge bosons of the 16 Cartans									hacking here to include the 16 Cartans!!!!
    for (unsigned Index = 10; Index < 26; ++Index)
    {
      CModedOscillator LM_alpha(LeftMover, Index, false, this->Twist, -1);
      if (LM_alpha.GetLadderOperatorType() != CreationOperator)
      {
        cout << "\n  Warning in bool CSector::CreateOscillators(...): Untwisted oscillators failed. Return false." << endl;
        return false;
      }
      this->LM_all_Oscillators.push_back(LM_alpha);
    }*/
		return true;
	}

	// twisted sector
	int m = -1;
	// create all moded oscillators
	for (unsigned Index = 1; Index < 4; ++Index)
	{
		// if the twist is non-trivial in the current plane
		if (!is_integer(this->Twist[Index]))
		{
			for (m = -1; m <= 0; ++m)
			{
				// moded oscillator: /tilde/alpha_{m - /eta^Index}^Index
				CModedOscillator LM_alpha(LeftMover, Index, false, this->Twist, m);
				if (LM_alpha.GetLadderOperatorType() == CreationOperator)
					this->LM_all_Oscillators.push_back(LM_alpha);

				CModedOscillator LM_alpha_bar(LeftMover, Index, true, this->Twist, m);
				if (LM_alpha_bar.GetLadderOperatorType() == CreationOperator)
					this->LM_all_Oscillators.push_back(LM_alpha_bar);

				CModedOscillator RM_alpha(RightMover, Index, false, this->Twist, m);
				if (RM_alpha.GetLadderOperatorType() == CreationOperator)
					this->RM_all_Oscillators.push_back(RM_alpha);

				CModedOscillator RM_alpha_bar(RightMover, Index, true, this->Twist, m);
				if (RM_alpha_bar.GetLadderOperatorType() == CreationOperator)
					this->RM_all_Oscillators.push_back(RM_alpha_bar);
			}
		}
	}
	return true;
}



/* ########################################################################################
######   CreateOscillatorExcitations()                                               ######
######                                                                               ######
######   Version: 12.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Find all (left- and right-moving) oscillator excitations such that the      ######
######   equations for massless left- and right-movers still might find solutions.   ######
######   The result is stored in the member variables "LM_Excitations" and           ######
######   "RM_Excitations".                                                           ###### 
######################################################################################## */
bool CSector::CreateOscillatorExcitations()
{
	if ((this->LM_Excitations.size() != 0) || (this->RM_Excitations.size() != 0))
	{
		cout << "\n  Warning in bool CSector::CreateOscillatorExcitations(...) : Sets of oscillators are not empty. Now cleared." << endl;
		this->LM_Excitations.clear();
		this->RM_Excitations.clear();
	}

	const size_t s1 = this->LM_all_Oscillators.size();
	const size_t s2 = this->RM_all_Oscillators.size();

	if (((s1 == 0) || (s2 == 0)) && (this->Twist.OrderOfTwist() != 1) && (this->Twist.OrderOfTwist() != 2)) //hacking here!!
	{
		cout << "\n  Warning in bool CSector::CreateOscillatorExcitations(...) : Sets of oscillators are empty. Return false.." << endl;
		return false;
	}

	// The mass equation for massless left movers reads:
	// (p + V_g)^2 = -2(N + a_L) .
	//
	// Surely, to solve the equation the following must hold:
	// (p + V_g)^2 >= 0 .
	//
	// It follows for the Number Operator:
	// N <= -a_L .
	//
	// The right movers case is analogous.
	//
	// create all subsets of moded oscillators with Number Operator N <= -a_L > 0

	S_OscillatorExcitation NoExcitation;
	NoExcitation.NumberOperator  = 0.0;

	NoExcitation.ZeroPointEnergy = -2.0 * this->Twist.Get_a_L();
	this->RecursiveCreate_LM_OscillatorExcitations(NoExcitation, 0);

	NoExcitation.ZeroPointEnergy = -2.0 * this->Twist.Get_a_R();
	this->RecursiveCreate_RM_OscillatorExcitations(NoExcitation, 0);

	/*cout << "sector (" <<Sector[0] <<", " << Sector[1]<<")\n";
	cout<<Twist.Get_a_L()<<endl;
	cout<<Twist.Get_a_R()<<endl;*/

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
######   "LeftChiral" fields of this sector.                                         ######
######################################################################################## */
bool CSector::GetFieldIndices(const vector<CField> &Fields, const SUSYMultiplet &Multiplet, vector<unsigned> &FieldIndices) const
{
	const size_t s1 = this->FixedBranes.size();
	for (unsigned i = 0; i < s1; ++i)
		this->FixedBranes[i].GetFieldIndices(Fields, Multiplet, FieldIndices);

	return true;
}



/* ########################################################################################
######   RecursiveCreate_LM_OscillatorExcitations(...)                               ######
######                                                                               ######
######   Version: 12.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) CurrentOscillatorExcitation : will be added to "LM_Excitations" if       ######
######                                    massless solutions might still exist       ######
######   2) Index                       : index in "LM_all_Oscillators"              ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Adds recursively an aditional oscillator excitation as long as there still  ######
######   might be solutions to the equation for massless left-movers.                ######
######################################################################################## */
void CSector::RecursiveCreate_LM_OscillatorExcitations(const S_OscillatorExcitation &CurrentOscillatorExcitation, unsigned Index)
{
	if (CurrentOscillatorExcitation.NumberOperator + Twist.Get_a_L() < 0.00001)
	{
		this->LM_Excitations.push_back(CurrentOscillatorExcitation);

		for (unsigned i = Index; i < this->LM_all_Oscillators.size(); ++i)
		{
			// add one moded oscillator
			S_OscillatorExcitation NewOscillatorExcitation = CurrentOscillatorExcitation;
			NewOscillatorExcitation.OscillatorIndices.push_back(i);
			NewOscillatorExcitation.NumberOperator += this->LM_all_Oscillators[i].GetNumberOperator();
			NewOscillatorExcitation.ZeroPointEnergy = -2.0 * (Twist.Get_a_L() + NewOscillatorExcitation.NumberOperator);

			this->RecursiveCreate_LM_OscillatorExcitations(NewOscillatorExcitation, i);
		}
	}
}



/* ########################################################################################
######   RecursiveCreate_tachyonicLM_OscillatorExcitations(...)                               ######
######                                                                               ######
######   Version: 12.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) CurrentOscillatorExcitation : will be added to "LM_Excitations" if       ######
######                                    massless solutions might still exist       ######
######   2) Index                       : index in "LM_all_Oscillators"              ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Adds recursively an additional oscillator excitation as long as there still  ######
######   might be solutions to the equation for massless left-movers.                ######
######################################################################################## */
void CSector::RecursiveCreate_tachyonicLM_OscillatorExcitations(const S_OscillatorExcitation &CurrentOscillatorExcitation, unsigned Index, const double &M_R)
{
	//a_L+N < 0 must be more negative than M_R^2<0 to level-match (p+V_g)^2/2 + a_L +N =M_L^2/8=M_R^2/8 since (p+V_g)^2>=0

	if (CurrentOscillatorExcitation.NumberOperator + Twist.Get_a_L() - M_R < 0.00001)
	{
		this->TachyonicLM_Excitations_tmp.push_back(CurrentOscillatorExcitation);

		for (unsigned i = Index; i < this->LM_all_Oscillators.size(); ++i)
		{
			// add one moded oscillator
			S_OscillatorExcitation NewOscillatorExcitation = CurrentOscillatorExcitation;
			NewOscillatorExcitation.OscillatorIndices.push_back(i);
			NewOscillatorExcitation.NumberOperator += this->LM_all_Oscillators[i].GetNumberOperator();
			NewOscillatorExcitation.ZeroPointEnergy = -2.0 * (Twist.Get_a_L() + NewOscillatorExcitation.NumberOperator - M_R);

			this->RecursiveCreate_tachyonicLM_OscillatorExcitations(NewOscillatorExcitation, i, M_R);
		}
	}
}



/* ########################################################################################
######   RecursiveCreate_RM_OscillatorExcitations(...)                               ######
######                                                                               ######
######   Version: 12.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) CurrentOscillatorExcitation : will be added to "RM_Excitations" if       ######
######                                    massless solutions might still exist       ######
######   2) Index                       : index in "RM_all_Oscillators"              ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Adds recursively an aditional oscillator excitation as long as there still  ######
######   might be solutions to the equation for massless right-movers.               ######
######################################################################################## */
void CSector::RecursiveCreate_RM_OscillatorExcitations(const S_OscillatorExcitation &CurrentOscillatorExcitation, unsigned Index)
{
	if (CurrentOscillatorExcitation.NumberOperator + Twist.Get_a_R() < 0.00001)
	{
		this->RM_Excitations.push_back(CurrentOscillatorExcitation);

		for (unsigned i = Index; i < this->RM_all_Oscillators.size(); ++i)
		{
			// add one moded oscillator
			S_OscillatorExcitation NewOscillatorExcitation = CurrentOscillatorExcitation;
			NewOscillatorExcitation.OscillatorIndices.push_back(i);
			NewOscillatorExcitation.NumberOperator += this->RM_all_Oscillators[i].GetNumberOperator();
			NewOscillatorExcitation.ZeroPointEnergy = -2.0 * (Twist.Get_a_R() + NewOscillatorExcitation.NumberOperator);

			this->RecursiveCreate_RM_OscillatorExcitations(NewOscillatorExcitation, i);
		}
	}
}



/* ########################################################################################
######   CreateMasslessRightMover()                                                  ######
######                                                                               ######
######   Version: 02.02.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Find all solutions to the equation of massless right-movers:                ###### 
######     1/2 (q+v_g)^2 - 1/2 + delta c = 0                                         ######
######   The result is stored in the member variable "MasslessRightMovers".          ###### 
######################################################################################## */
bool CSector::CreateMasslessRightMover()
{
	const size_t s1 = this->RM_Excitations.size();
	if (s1 == 0)
	{
		cout << "\n  Warning in bool CSector::CreateMasslessRightMover(...) : Set of excitations is empty. Return false.";
		return false;
	}
	if (this->MasslessRightMovers.size() != 0)
	{
		cout << "\n  Warning in bool CSector::CreateMasslessRightMover(...) : Massless right-movers have been created before - now cleared." << endl;
		this->MasslessRightMovers.clear();
	}
	for (unsigned i = 0; i < s1; ++i)
	{
		CMasslessHalfState HalfState(RightMover, this->RM_Excitations[i]);

		if (HalfState.SolveMassEquation(this->Twist, SO8))
			this->MasslessRightMovers.push_back(HalfState);
	}
	return true;
}



/* ########################################################################################	    //hacking here	to compute right moving tachyons
######   CreateTachyonicRightMover()                                                 ######
######                                                                               ######
######   Version: 02.02.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Find all solutions to the equation of massless right-movers:                ######
######     1/2 (q+v_g)^2 - 1/2 + delta c = 0                                         ######
######   The result is stored in the member variable "TachyonicRightMovers".         ######
######################################################################################## */
bool CSector::CreateTachyonicRightMover()
{
	const size_t s1 = this->RM_Excitations.size();

	for (unsigned i = 0; i < s1; ++i)
	{
		CTachyonHalfState TachyonHalfState(RightMover, this->RM_Excitations[i]);

		if ( TachyonHalfState.TachyonSolver(Sector[1], Sector[2], this->Twist, SO8) )
			TachyonMassRightMovers.push_back(TachyonHalfState);
	}

	S_OscillatorExcitation NoExcitation;
	NoExcitation.NumberOperator  = 0.0;
	const size_t s2 = this->TachyonMassRightMovers.size();

	//Run through right moving tachyons in this (1,n,k) sector
	for (unsigned i=0; i<s2; ++i) {

		//define left-moving vacuum shift to level-match right-moving tachyon mass
		NoExcitation.ZeroPointEnergy = -2.0 * (this->Twist.Get_a_L()-this->TachyonMassRightMovers[i].tachyonmass[0]);
		this->RecursiveCreate_tachyonicLM_OscillatorExcitations(NoExcitation, 0, this->TachyonMassRightMovers[i].tachyonmass[0]);

		this->TachyonicLM_Excitations.push_back(this->TachyonicLM_Excitations_tmp);
		this->TachyonicLM_Excitations_tmp.clear();
	}

	return true;
}



/* ########################################################################################
######   SectorHasLeftchiralRightmover() const                                       ######
######                                                                               ######
######   Version: 25.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : does this sector have a massless right-mover which in        ######
######                  principle can give left-chiral fields?                       ######    
###########################################################################################
######   description:                                                                ######
######   Goes through the list of massless right-movers (stored in the member        ######
######   variable "MasslessRightMovers") and checks whether at least one of them can ######
######   give left-chiral fields in principle.                                       ######    
######################################################################################## */
bool CSector::SectorHasLeftchiralRightmover() const
{
	unsigned j = 0;
	size_t s2 = 0;

	const size_t s1 = this->MasslessRightMovers.size();
	for (unsigned i = 0; i < s1; ++i)
	{
		const vector<CVector> &Weights = this->MasslessRightMovers[i].Weights;
		s2 = Weights.size();
		for (j = 0; j < s2; ++j)
		{
			if (fabs(Weights[j][0] + 0.5) < 0.001)
				return true;
		}
	}
	return false;
}



/* ########################################################################################
######   SortByEigenvalue(const vector<CTwistVector> &ZMxZNxZK_Twists)               ######
######                                                                               ######
######   Version: 12.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   4) ZMxZN_Twists : v_1 and v_2 where v_2 = (0^4) for Z_M orbifolds           ######
######   output:                                                                     ######
######   return value    : finished successfully?                                    ######
###########################################################################################
######   description:                                                                ######
######   Sort the solutions (q+v_g) from "MasslessRightMovers" with respect to their ######
######   eigenvalues under transformations with the twist(s) v_1 (and v_2, for       ######
######   Z_M x Z_N orbifolds) and stores the sorted result in "RightMovers".         ######
######################################################################################## */
bool CSector::SortByEigenvalue(const vector<CTwistVector> &ZMxZNxZK_Twists)
{
	// Set the precision
	const double prec = 0.001;

	const size_t number_of_eigenvalues = ZMxZNxZK_Twists.size();

	double Eigenvalue = 0;
	vector<double> OsciEigenvalues(number_of_eigenvalues, 0.0);
	vector<double> Eigenvalues(number_of_eigenvalues, 0.0);

	unsigned i = 0;
	unsigned j = 0;
	unsigned k = 0;
	unsigned l = 0;

	size_t s2 = 0;
	size_t s3 = 0;

	bool with_excitation       = false;
	bool Are_Eigenvalues_Equal = true;
	bool Eigenvalues_NotKnown  = true;

	// run through all massless right-movers
	const size_t s1 = this->MasslessRightMovers.size();
	for (i = 0; i < s1; ++i)
	{
		const CMasslessHalfState &MasslessHalfState = this->MasslessRightMovers[i];
		const vector<CVector>    &Weights           = MasslessHalfState.Weights;

		// begin: compute the transformation of the oscillators
		with_excitation = (MasslessHalfState.Excitation.OscillatorIndices.size() != 0);
		if (with_excitation)
		{
			const vector<unsigned> &OscillatorIndices = MasslessHalfState.Excitation.OscillatorIndices;
			s2 = OscillatorIndices.size();

			for (j = 0; j < number_of_eigenvalues; ++j)
			{
				const CTwistVector &SxG_Twist = ZMxZNxZK_Twists[j];

				Eigenvalue = 0.0;
				for (k = 0; k < s2; ++k)
					Eigenvalue += this->RM_all_Oscillators[OscillatorIndices[k]].GetTransformationProperty(SxG_Twist);

				OsciEigenvalues[j] = Eigenvalue;
			}
		}
		// end: compute the transformation of the oscillators

		s2 = Weights.size();
		// get next weight of the current right-mover
		for (j = 0; j < s2; ++j)
		{
			const CVector &Weight = Weights[j];

			// calculate the scalar-products (q x v)
			for (k = 0; k < number_of_eigenvalues; ++k)
			{
				Eigenvalue = ZMxZNxZK_Twists[k] * Weight;

				// begin: add the transformation of the oscillators
				if (with_excitation)
					Eigenvalue -= OsciEigenvalues[k];							//hacking here!!! wrong trafo sign for R-moving oscillator excitations in original Orbifolder
				// end: add the transformation of the oscillators

				RoundDouble(Eigenvalue);
				Eigenvalue -= floor(Eigenvalue);

				Eigenvalues[k] = Eigenvalue;
			}

			// check whether this combination of eigenvalues is known or not
			Eigenvalues_NotKnown = true;

			s3 = this->RightMovers.size();
			for (k = 0; Eigenvalues_NotKnown && (k < s3); ++k)
			{
				CHalfState &RightMover = this->RightMovers[k];
				if (RightMover.GetIndex() == i)
				{
					// if the set of eigenvalues is known
					Are_Eigenvalues_Equal = true;
					for (l = 0; Are_Eigenvalues_Equal && (l < number_of_eigenvalues); ++l)
					{
						if (fabs(Eigenvalues[l] - RightMover.Eigenvalues[l]) > prec)
							Are_Eigenvalues_Equal = false;
					}
					// add current weight to the corresponding right-mover
					if (Are_Eigenvalues_Equal)
					{
						RightMover.Weights.push_back(j);
						Eigenvalues_NotKnown = false;
					}
				}
			}

			// if the current set of eigenvalues is not known
			// create a new RightMover
			if (Eigenvalues_NotKnown)
			{
				CHalfState new_RightMover(RightMover, i);
				this->RightMovers.push_back(new_RightMover);
				this->RightMovers.rbegin()->Eigenvalues = Eigenvalues;
				this->RightMovers.rbegin()->Weights.push_back(j);
			}
		}
	}
	return true;
}



/* ########################################################################################
######   TachyonicSortByEigenvalue(const vector<CTwistVector> &ZMxZNxZK_Twists)               ######
######                                                                               ######
######   Version: 12.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   4) ZMxZN_Twists : v_1 and v_2 where v_2 = (0^4) for Z_M orbifolds           ######
######   output:                                                                     ######
######   return value    : finished successfully?                                    ######
###########################################################################################
######   description:                                                                ######
######   Sort the solutions (q+v_g) from "MasslessRightMovers" with respect to their ######
######   eigenvalues under transformations with the twist(s) v_1 (and v_2, for       ######
######   Z_M x Z_N orbifolds) and stores the sorted result in "RightMovers".         ######
######################################################################################## */
bool CSector::TachyonicSortByEigenvalue(const vector<CTwistVector> &ZMxZNxZK_Twists)
{
	// Set the precision
	const double prec = 0.001;

	const size_t number_of_eigenvalues = ZMxZNxZK_Twists.size();

	double Eigenvalue = 0;
	vector<double> OsciEigenvalues(number_of_eigenvalues, 0.0);
	vector<double> Eigenvalues(number_of_eigenvalues, 0.0);

	unsigned i = 0;
	unsigned j = 0;
	unsigned k = 0;
	unsigned l = 0;

	size_t s2 = 0;
	size_t s3 = 0;

	bool with_excitation       = false;
	bool Are_Eigenvalues_Equal = true;
	bool Eigenvalues_NotKnown  = true;

	// run through all massless right-movers
	const size_t s1 = this->TachyonMassRightMovers.size();
	for (i = 0; i < s1; ++i)
	{
		const CTachyonHalfState &TachyonHalfState = this->TachyonMassRightMovers[i];
		const vector<CVector>   &Weights          = TachyonHalfState.Weights;

		// begin: compute the transformation of the oscillators
		with_excitation = (TachyonHalfState.Excitation.OscillatorIndices.size() != 0);
		if (with_excitation)
		{
			const vector<unsigned> &OscillatorIndices = TachyonHalfState.Excitation.OscillatorIndices;
			s2 = OscillatorIndices.size();

			for (j = 0; j < number_of_eigenvalues; ++j)
			{
				const CTwistVector &SxG_Twist = ZMxZNxZK_Twists[j];

				Eigenvalue = 0.0;
				for (k = 0; k < s2; ++k)
					Eigenvalue += this->RM_all_Oscillators[OscillatorIndices[k]].GetTransformationProperty(SxG_Twist);

				OsciEigenvalues[j] = Eigenvalue;
			}
		}
		// end: compute the transformation of the oscillators

		s2 = Weights.size();
		// get next weight of the current right-mover
		for (j = 0; j < s2; ++j)
		{
			const CVector &Weight = Weights[j];

			// calculate the scalar-products (q x v)
			for (k = 0; k < number_of_eigenvalues; ++k)
			{
				Eigenvalue = ZMxZNxZK_Twists[k] * Weight;

				// begin: add the transformation of the oscillators
				if (with_excitation)
					Eigenvalue -= OsciEigenvalues[k];							//hacking here!!! wrong trafo sign for R-moving oscillator excitations in original Orbifolder
				// end: add the transformation of the oscillators

				RoundDouble(Eigenvalue);
				Eigenvalue -= floor(Eigenvalue);

				Eigenvalues[k] = Eigenvalue;
			}

			// check whether this combination of eigenvalues is known or not
			Eigenvalues_NotKnown = true;

			s3 = this->TachyonicRightMovers.size();
			for (k = 0; Eigenvalues_NotKnown && (k < s3); ++k)
			{
				CHalfState &RightMover = this->TachyonicRightMovers[k];
				if (RightMover.GetIndex() == i)
				{
					// if the set of eigenvalues is known
					Are_Eigenvalues_Equal = true;
					for (l = 0; Are_Eigenvalues_Equal && (l < number_of_eigenvalues); ++l)
					{
						if (fabs(Eigenvalues[l] - RightMover.Eigenvalues[l]) > prec)
							Are_Eigenvalues_Equal = false;
					}
					// add current weight to the corresponding right-mover
					if (Are_Eigenvalues_Equal)
					{
						RightMover.Weights.push_back(j);
						Eigenvalues_NotKnown = false;
					}
				}
			}

			// if the current set of eigenvalues is not known
			// create a new RightMover
			if (Eigenvalues_NotKnown)
			{
				CHalfState new_RightMover(RightMover, i);						//refers to the corresponding CTachyonHalfState object
				if (with_excitation)
					new_RightMover.Excited = 1;

				this->TachyonicRightMovers.push_back(new_RightMover);
				this->TachyonicRightMovers.rbegin()->Eigenvalues = Eigenvalues;
				this->TachyonicRightMovers.rbegin()->Weights.push_back(j);
			}
		}
	}
	return true;
}



/* ########################################################################################
######   AccessFixedBrane(const unsigned &i)                                         ######
######                                                                               ######
######   Version: 12.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) i         : index in the private member variable "FixedBranes"           ######
######   output:                                                                     ######
######   return value : reference to the CFixedBrane object in the i-th position     ######
######                  of "FixedBranes"                                             ######
###########################################################################################
######   description:                                                                ######
######   Allows access to the content of the private member variable "FixedBranes".  ######
######################################################################################## */
CFixedBrane &CSector::AccessFixedBrane(const unsigned &i)
{
#ifdef CHECKERROR
	if (i >= this->FixedBranes.size())
	{
		cout << "\n  Warning in CFixedBrane &CSector::AccessFixedBrane(...) : Index i out of range. Set i = 0." << endl;
		return this->FixedBranes[0];
	}
#endif

	return this->FixedBranes[i];
}



/* ########################################################################################
######   &GetFixedBrane(const unsigned &i) const                                     ######
######                                                                               ######
######   Version: 12.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) i         : index in the private member variable "FixedBranes"           ######
######   output:                                                                     ######
######   return value : constant reference to the CFixedBrane object in the i-th     ######
######                  position of "FixedBranes"                                    ######
###########################################################################################
######   description:                                                                ######
######   Allows constant access to the content of the private member variable        ######
######   "FixedBranes".                                                              ######
######################################################################################## */
const CFixedBrane &CSector::GetFixedBrane(const unsigned &i) const
{
#ifdef CHECKERROR
	if (this->Sector_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const CFixedBrane &CSector::GetFixedBrane(...) const : Sector ill-defined." << endl;

	if (i >= this->FixedBranes.size())
	{
		cout << "\n  Warning in const CFixedBrane &CSector::GetFixedBrane(...) const : Index i out of range. Set i = 0." << endl;
		return this->FixedBranes[0];
	}
#endif

	return this->FixedBranes[i];
}



/* ########################################################################################
######   &GetLM_Excitations() const                                                  ######
######                                                                               ######
######   Version: 12.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : constant reference to "LM_Excitations"                       ######
###########################################################################################
######   description:                                                                ######
######   Allows constant access to the content of the private member variable        ######
######   "LM_Excitations".                                                           ######
######################################################################################## */
const vector<S_OscillatorExcitation> &CSector::GetLM_Excitations() const
{
#ifdef CHECKERROR
	if (this->Sector_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const vector<S_OscillatorExcitation> &CSector::GetLM_Excitations() const : Sector ill-defined." << endl;
#endif

	return this->LM_Excitations;
}



/* ########################################################################################
######   &GetLM_Excitations() const                                                  ######
######                                                                               ######
######   Version: 12.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : constant reference to "LM_Excitations"                       ######
###########################################################################################
######   description:                                                                ######
######   Allows constant access to the content of the private member variable        ######
######   "LM_Excitations".                                                           ######
######################################################################################## */
const vector<vector<S_OscillatorExcitation> > &CSector::TachyonicGetLM_Excitations() const
{
#ifdef CHECKERROR
	if (this->Sector_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const vector<vector<S_OscillatorExcitation> > &CSector::TachyonicGetLM_Excitations() const : Sector ill-defined." << endl;
#endif

	return this->TachyonicLM_Excitations;

}



/* ########################################################################################
######   &GetLM_all_Oscillators() const                                              ######
######                                                                               ######
######   Version: 12.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : constant reference to "LM_all_Oscillators"                   ######
###########################################################################################
######   description:                                                                ######
######   Allows constant access to the content of the private member variable        ######
######   "LM_all_Oscillators".                                                       ######
######################################################################################## */
const vector<CModedOscillator> &CSector::GetLM_all_Oscillators() const
{
#ifdef CHECKERROR
	if (this->Sector_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const vector<CModedOscillator> &CSector::GetLM_all_Oscillators() const : Sector ill-defined." << endl;
#endif

	return this->LM_all_Oscillators;
}



/* ########################################################################################
######   &GetLM_Oscillator(const unsigned &i) const                                  ######
######                                                                               ######
######   Version: 12.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) i         : specify the i-th oscillator in "LM_all_Oscillators"          ######
######   output:                                                                     ######
######   return value : constant reference to "LM_all_Oscillators[i]"                ######
###########################################################################################
######   description:                                                                ######
######   Allows constant access to the content of the private member variable        ######
######   "LM_all_Oscillators" at the "i"-th position.                                ######
######################################################################################## */
const CModedOscillator &CSector::GetLM_Oscillator(const unsigned &i) const
{
#ifdef CHECKERROR
	if (this->Sector_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const CModedOscillator &CSector::GetLM_Oscillator(...) const : Sector ill-defined." << endl;

	if (i >= this->LM_all_Oscillators.size())
	{
		cout << "\n  Warning in const CModedOscillator &CSector::GetLM_Oscillator(...) const : Index i out of range. Set i = 0." << endl;
		return this->LM_all_Oscillators[0];
	}

	if (this->LM_all_Oscillators[i].GetOscillator_CheckStatus() != CheckedAndGood)
		cout << "\n  Warning in const CModedOscillator &CSector::GetLM_Oscillator(...) const : Oscillator ill-defined." << endl;
#endif

	return this->LM_all_Oscillators[i];
}



/* ########################################################################################
######   &GetRM_Excitations() const                                                  ######
######                                                                               ######
######   Version: 12.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : constant reference to "RM_Excitations"                       ######
###########################################################################################
######   description:                                                                ######
######   Allows constant access to the content of the private member variable        ######
######   "RM_Excitations".                                                           ######
######################################################################################## */
const vector<S_OscillatorExcitation> &CSector::GetRM_Excitations() const
{
#ifdef CHECKERROR
	if (this->Sector_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const vector<S_OscillatorExcitation> &CSector::GetRM_Excitations() const : Sector ill-defined." << endl;
#endif

	return this->RM_Excitations;
}



/* ########################################################################################
######   &GetRM_all_Oscillators() const                                              ######
######                                                                               ######
######   Version: 12.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : constant reference to "RM_all_Oscillators"                   ######
###########################################################################################
######   description:                                                                ######
######   Allows constant access to the content of the private member variable        ######
######   "RM_all_Oscillators".                                                       ######
######################################################################################## */
const vector<CModedOscillator> &CSector::GetRM_all_Oscillators() const
{
#ifdef CHECKERROR
	if (this->Sector_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const vector<CModedOscillator> &CSector::GetRM_all_Oscillators() const : Sector ill-defined." << endl;
#endif

	return this->RM_all_Oscillators;
}



/* ########################################################################################
######   &GetRM_Oscillator(const unsigned &i) const                                  ######
######                                                                               ######
######   Version: 12.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) i         : specify the i-th oscillator in "RM_all_Oscillators"          ######
######   output:                                                                     ######
######   return value : constant reference to "RM_all_Oscillators[i]"                ######
###########################################################################################
######   description:                                                                ######
######   Allows constant access to the content of the private member variable        ######
######   "RM_all_Oscillators" at the "i"-th position.                                ######
######################################################################################## */
const CModedOscillator &CSector::GetRM_Oscillator(const unsigned &i) const
{
#ifdef CHECKERROR
	if (this->Sector_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const CModedOscillator &CSector::GetRM_Oscillator(...) const : Sector ill-defined." << endl;

	if (i >= this->RM_all_Oscillators.size())
	{
		cout << "\n  Warning in const CModedOscillator &CSector::GetRM_Oscillator(...) const : Index i out of range. Set i = 0." << endl;
		return this->RM_all_Oscillators[0];
	}

	if (this->RM_all_Oscillators[i].GetOscillator_CheckStatus() != CheckedAndGood)
		cout << "\n  Warning in const CModedOscillator &CSector::GetRM_Oscillator(...) const : Oscillator ill-defined." << endl;
#endif

	return this->RM_all_Oscillators[i];
}



/* ########################################################################################
######   &GetMasslessRightMover(const unsigned &i) const                             ######
######                                                                               ######
######   Version: 12.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) i         : specify the i-th solution in "MasslessRightMovers"           ######
######   output:                                                                     ######
######   return value : constant reference to "MasslessRightMovers[i]"               ######
###########################################################################################
######   description:                                                                ######
######   Allows constant access to the content of the private member variable        ######
######   "MasslessRightMovers" at the "i"-th position.                               ######
######################################################################################## */
const CMasslessHalfState &CSector::GetMasslessRightMover(const unsigned &i) const
{
#ifdef CHECKERROR
	if (this->Sector_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const CMasslessHalfState &CSector::GetMasslessRightMover(...) const : Sector ill-defined." << endl;

	if (i >= this->MasslessRightMovers.size())
	{
		cout << "\n  Warning in const CMasslessHalfState &CSector::GetMasslessRightMover(...) const : Index i out of range. Set i = 0." << endl;
		return this->MasslessRightMovers[0];
	}
#endif

	return this->MasslessRightMovers[i];
}



/* ########################################################################################
######   &GetMasslessRightMovers() const                                             ######
######                                                                               ######
######   Version: 12.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : constant reference to "MasslessRightMovers"                  ######
###########################################################################################
######   description:                                                                ######
######   Allows constant access to the content of the private member variable        ######
######   "MasslessRightMovers".                                                      ######
######################################################################################## */
const vector<CMasslessHalfState> &CSector::GetMasslessRightMovers() const
{
#ifdef CHECKERROR
	if (this->Sector_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const vector<CMasslessHalfState> &CSector::GetMasslessRightMovers() const : Sector ill-defined." << endl;
#endif

	return this->MasslessRightMovers;
}



/* ########################################################################################
######   &GetRightMover(const unsigned &i) const                                     ######
######                                                                               ######
######   Version: 12.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) i         : specify the i-th solution in "RightMovers"                   ######
######   output:                                                                     ######
######   return value : constant reference to "RightMovers[i]"                       ######
###########################################################################################
######   description:                                                                ######
######   Allows constant access to the content of the private member variable        ######
######   "RightMovers" at the "i"-th position.                                       ######
######################################################################################## */
const CHalfState &CSector::GetRightMover(const unsigned &i) const
{
#ifdef CHECKERROR
	if (this->Sector_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const CHalfState &CSector::GetRightMover(...) const : Sector ill-defined." << endl;

	if (i >= this->RightMovers.size())
	{
		cout << "\n  Warning in const CHalfState &CSector::GetRightMover(...) const : Index i out of range. Set i = 0." << endl;
		return this->RightMovers[0];
	}
#endif

	return this->RightMovers[i];
}



/* ########################################################################################
######   &GetRightMovers() const                                                     ######
######                                                                               ######
######   Version: 12.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : constant reference to "RightMovers"                          ######
###########################################################################################
######   description:                                                                ######
######   Allows constant access to the content of the private member variable        ######
######   "RightMovers".                                                              ######
######################################################################################## */
const vector<CHalfState> &CSector::GetRightMovers() const
{
#ifdef CHECKERROR
	if (this->Sector_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const vector<CHalfState> &CSector::GetRightMovers() const : Sector ill-defined." << endl;
#endif

	return this->RightMovers;
}



/* ########################################################################################
######   &GetRightMovers() const                                                     ######
######                                                                               ######
######   Version: 12.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : constant reference to "RightMovers"                          ######
###########################################################################################
######   description:                                                                ######
######   Allows constant access to the content of the private member variable        ######
######   "RightMovers".                                                              ######
######################################################################################## */
const vector<CHalfState> &CSector::TachyonicGetRightMovers() const
{
#ifdef CHECKERROR
	if (this->Sector_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const vector<CHalfState> &CSector::GetRightMovers() const : Sector ill-defined." << endl;
#endif

	return this->TachyonicRightMovers;
}



/* ########################################################################################			//hacking here!!!
######   &GetTachyonicRightMovers() const                                            ######
######                                                                               ######
######   Version: 12.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : constant reference to "RightMovers"                          ######
###########################################################################################
######   description:                                                                ######
######   Allows constant access to the content of the private member variable        ######
######   "TachyonicRightMovers".                                                     ######
######################################################################################## */
const vector<CTachyonHalfState> &CSector::GetRTachyons() const
{
#ifdef CHECKERROR
	if (this->Sector_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const vector<CHalfState> &CSector::GetRightMovers() const : Sector ill-defined." << endl;
#endif

	return this->TachyonMassRightMovers;
}

