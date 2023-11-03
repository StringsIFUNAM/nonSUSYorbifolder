#include <iostream>

#include "cspacegroup.h"
#include "cspacegroupelement.h"
#include "globalfunctions.h"
#include "cprint.h"
#include "corbifoldgroupelement.h"
#include "clinalg.h"

using std::cout;
using std::endl;
using std::polar;



/* ########################################################################################
######   CSpaceGroup()                                                               ######
######                                                                               ######
######   Version: 09.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard constructor of a CSpaceGroup object. No content is specified.      ######
######################################################################################## */
CSpaceGroup::CSpaceGroup()
: WL_AllowedOrders(LatticeDim, 0)
{
	CTwistVector ZeroTwist;
	this->Twists.assign(3, ZeroTwist);				//hacking here!!!

	this->ZMxZN = false;
	this->ZMxZNxZK = false;
	this->M = 1;
	this->N = 1;
	this->K = 1;
	this->additional_label = "";
	this->GeometryFilename = "";
	this->lattice_label = "";

	this->SpaceGroup_CheckStatus = NotChecked;
}



/* ########################################################################################
######   CSpaceGroup(const unsigned &M, const unsigned &N)                           ######
######                                                                               ######
######   Version: 09.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) M : order of Z_M                                                         ######
######   2) N : order of Z_N; set N=1 for Z_M orbifolds                              ######
######   2) K : order of Z_K; set K=1 for Z_MxZ_N orbifolds                          ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Constructor of a CSpaceGroup object. The point group is Z_M x Z_N.          ######
######################################################################################## */
CSpaceGroup::CSpaceGroup(const unsigned &M, const unsigned &N, const unsigned &K)
: WL_AllowedOrders(LatticeDim, 0)
{
	CTwistVector ZeroTwist;
	this->Twists.assign(3, ZeroTwist);		//hacking here!!!

	this->additional_label = "";
	this->GeometryFilename = "";
	this->lattice_label    = "";

	this->SetOrder(M, N, K);		//hacking here!!!

	this->SpaceGroup_CheckStatus = NotChecked;
}



/* ########################################################################################
######   Check()                                                                     ######
######                                                                               ######
######   Version: 05.09.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : space group checked succesfully?                             ######
###########################################################################################
######   description:                                                                ######
######   Checks the consistency of this space group.                                 ######
######################################################################################## */
bool CSpaceGroup::Check()
{
	if (this->SpaceGroup_CheckStatus == CheckedAndGood)
		return true;

	if ((this->T6_Lattice.size() != LatticeDim) || (this->DualT6_Lattice.size() != LatticeDim)
			|| (this->Twists.size() != 3)              || (this->TwistMatrices.size() == 0)					//hacking here!!!
			|| (this->ShiftsWL_ScalarProductFactors.size() == 0)
			|| (this->SG_Generators_Twist.size() == 0) || (this->SG_Generators_Shift.size() == 0)
			|| (this->Sectors.size() == 0)             || (this->Centralizer.size() == 0)
			|| (this->SG_AllRotations.size() == 0)     || (this->WL_AllowedOrders.size() != LatticeDim))
	{
		cout << "\n  Warning in bool CSpaceGroup::Check() : problems with:" << endl;

		if (this->T6_Lattice.size() != LatticeDim)
			cout << "    this->T6_Lattice.size()" << endl;
		if (this->Twists.size() != 3)																		//hacking here!!!
			cout << "    this->Twists.size()" << endl;
		if (this->ShiftsWL_ScalarProductFactors.size() == 0)
			cout << "    this->ShiftsWL_ScalarProductFactors.size()" << endl;
		if (this->TwistMatrices.size() == 0)
			cout << "    this->TwistMatrices.size()" << endl;
		if (this->DualT6_Lattice.size() != LatticeDim)
			cout << "    this->DualT6_Lattice.size()" << endl;
		if (this->SG_Generators_Twist.size() == 0)
			cout << "    this->SG_Generators_Twist.size()" << endl;
		if (this->SG_Generators_Shift.size() == 0)
			cout << "    this->SG_Generators_Shift.size()" << endl;
		if (this->Sectors.size() == 0)
			cout << "    this->Sectors.size()" << endl;
		if (this->Centralizer.size() == 0)
			cout << "    this->Centralizer.size()" << endl;
		if (this->SG_AllRotations.size() == 0)
			cout << "    this->SG_AllRotations.size()" << endl;
		if (this->WL_AllowedOrders.size() != LatticeDim)
			cout << "    this->WL_AllowedOrders.size()" << endl;

		this->SpaceGroup_CheckStatus = CheckedAndFailed;
		return false;
	}

	this->SpaceGroup_CheckStatus = CheckedAndGood;
	return true;
}



/* ########################################################################################
######   Clear()                                                                     ######
######                                                                               ######
######   Version: 16.11.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Clear the content of this space group.                                      ######
######################################################################################## */
void CSpaceGroup::Clear()
{
	CTwistVector ZeroTwist;
	this->Twists.assign(3, ZeroTwist);

	this->ZMxZN = false;
	this->ZMxZNxZK = false;
	this->M = 1;
	this->N = 1;
	this->K = 1;

	this->additional_label = "";
	this->GeometryFilename = "";
	this->lattice_label = "";

	this->T6_Lattice.clear();

	this->ShiftsWL_ScalarProductFactors.clear();

	this->DiscreteNonRSymmetries.clear();
	this->DiscreteRSymmetries.clear();
	this->ModularSymmetries.clear();

	this->TwistMatrices.clear();
	this->DualT6_Lattice.clear();

	this->SG_Generators_Twist.clear();
	this->SG_Generators_Shift.clear();

	this->Sectors.clear();
	this->Centralizer.clear();

	this->SG_SetOfElements.clear();

	this->SG_AllRotations.clear();
	this->SG_AllNonStandardShifts.clear();

	this->WL_Relations.clear();
	this->WL_AllowedOrders.assign(LatticeDim, 0);

	this->SpaceGroup_CheckStatus = NotChecked;
}



/* ########################################################################################
######   SetOrder(const unsigned &M, const unsigned &N, const unsigned &K)           ######
######                                                                               ######
######   Version: 09.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) M         : the order of Z_M                                             ######
######   2) N         : the order of Z_N; set N=1 for Z_M orbifolds                  ######
######   output:                                                                     ######
######   return value : is the point group correct?                                  ######
###########################################################################################
######   description:                                                                ######
######   Set the order(s) of the point group.                                        ######
######################################################################################## */
bool CSpaceGroup::SetOrder(const unsigned &M, const unsigned &N, const unsigned &K)
{
	if (this->SpaceGroup_CheckStatus == CheckedAndGood)
		this->SpaceGroup_CheckStatus = NotChecked;

	this->M = M;
	this->N = N;
	this->K = K;
    
    

	if (M == 1)
	{
		cout << "\n  Warning in bool CSpaceGroup::SetOrder(...) : M == 1." << endl;
		if (N != 1 && K != 1)
		{
			cout << "  Warning in bool CSpaceGroup::SetOrder(...) : Set M=N and N=K." << endl;
			this->M = N;
			this->N = K;
		}
		else if (N != 1 && K == 1)
		{
			cout << "  Warning in bool CSpaceGroup::SetOrder(...) : exchange M and N." << endl;
			this->M = N;
			this->N = M;
		}
		else if (N == 1 && K != 1)
		{
			cout << "  Warning in bool CSpaceGroup::SetOrder(...) : exchange M and K." << endl;
			this->M = K;
			this->K = M;
		}
		return false;
	}

	if (this->N == 1) {
		this->ZMxZN = false;
		this->ZMxZNxZK = false;
	}
	else {
		this->ZMxZN = true;
		if (this->K == 1)
			this->ZMxZNxZK = false;
		else
			this->ZMxZNxZK = true;
	}

	return true;
}

/* ########################################################################################
######   CreateConstructingElements()                                                ######
######                                                                               ######
######   Version: 11.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Find all inequivalent constructing elements from the set "SG_SetOfElements".######
######################################################################################## */
bool CSpaceGroup::CreateConstructingElements()
{
	const double prec = 0.0001;

	if (this->Sectors.size() != 0)
	{
		cout << "\n  Warning in bool CSpaceGroup::CreateConstructingElements(): Sectors are not empty! Now cleared." << endl;
		this->Sectors.clear();
	}

	const size_t s1 = this->SG_SetOfElements.size();
	if (s1 == 0)
	{
		cout << "\n  Warning in bool CSpaceGroup::CreateConstructingElements(): Elements of the space group not defined. Return false." << endl;
		return false;
	}

	unsigned i = 0;
	unsigned j = 0;
	unsigned k = 0;

	vector<unsigned>          Sector_mnk(3, 0);
	vector<vector<unsigned> > Known_Sectors_kl;
	vector<vector<unsigned> > ::iterator pos;

	FPCoordinates FixedPoint1;
	FPCoordinates FixedPoint2;
	vector<CSpaceGroupElement> Sectors_mnk;
	bool new_element = true;
	bool positive1   = true;
	bool positive2   = true;
	bool closer      = true;
	size_t s2 = 0;

	for (i = 0; i < s1; ++i)
	{
		const CSpaceGroupElement &SG_Element = this->SG_SetOfElements[i];

		// does the new element have an associated fixed point?
		new_element = this->SG_FixedPoint(SG_Element, FixedPoint1);

		if (new_element)
		{
			Sector_mnk[0] = SG_Element.Get_m();
			Sector_mnk[1] = SG_Element.Get_n();
			Sector_mnk[2] = SG_Element.Get_k();


			pos = find(Known_Sectors_kl.begin(), Known_Sectors_kl.end(), Sector_mnk);
			if (pos == Known_Sectors_kl.end())
			{
				Known_Sectors_kl.push_back(Sector_mnk);

				Sectors_mnk.clear();
				Sectors_mnk.push_back(SG_Element);
				this->Sectors.push_back(Sectors_mnk);
				//cout << "  new constr. element: " << SG_Element.k << " " << SG_Element.l << " " << SG_Element.n_alpha[0] << " " << SG_Element.n_alpha[1] << " " << SG_Element.n_alpha[2] << " " << SG_Element.n_alpha[3] << " " << SG_Element.n_alpha[4] << " " << SG_Element.n_alpha[5] << endl;
			}
			else
			{
				vector<CSpaceGroupElement> &knownElementsOfCurrentSector = this->Sectors[distance(Known_Sectors_kl.begin(), pos)];

				// is the new element in the same conjugacy class of a known one?
				s2 = knownElementsOfCurrentSector.size();
				for (j = 0; j < s2; ++j)
				{
					CSpaceGroupElement &knownElement = knownElementsOfCurrentSector[j];

					// if the new element in the same conjugacy class
					if (this->SG_FromSameConjugationClass(knownElement, SG_Element))
					{
						//cout << "from same cc" << endl;
						new_element = false;

						// are the coordinates of the new element all positive?
						positive1 = true;
						for (k = 0; positive1 && (k < ComplexLatticeDim); ++k)
						{
							if (!FixedPoint1.FixedTorus[k] && ((FixedPoint1.Coordinates[k].real() < -prec) || (FixedPoint1.Coordinates[k].imag() < -prec)))
								positive1 = false;
						}

						if (positive1)
						{
							// are the coordinates of the known element all positive?
							this->SG_FixedPoint(knownElement, FixedPoint2);
							positive2 = true;
							for (k = 0; positive2 && (k < ComplexLatticeDim); ++k)
							{
								if (!FixedPoint2.FixedTorus[k] && ((FixedPoint2.Coordinates[k].real() < -prec) || (FixedPoint2.Coordinates[k].imag() < -prec)))
									positive2 = false;
							}

							// is the new element closer to zero than the known one?
							closer = true;
							for (k = 0; closer && (k < ComplexLatticeDim); ++k)
							{
								if (!FixedPoint1.FixedTorus[k]
								                            && ((fabs(FixedPoint1.Coordinates[k].real()) > fabs(FixedPoint2.Coordinates[k].real()) + prec)
								                            		|| (fabs(FixedPoint1.Coordinates[k].imag()) > fabs(FixedPoint2.Coordinates[k].imag()) + prec)))
									closer = false;
							}

							if (!positive2 || closer)
								knownElement = SG_Element;
						}
					}
				}
				if (new_element)
				{
					//cout << "new" << endl;
					//cout << "  new constr. element: " << SG_Element.k << " " << SG_Element.l << " " << SG_Element.n_alpha[0] << " " << SG_Element.n_alpha[1] << " " << SG_Element.n_alpha[2] << " " << SG_Element.n_alpha[3] << " " << SG_Element.n_alpha[4] << " " << SG_Element.n_alpha[5] << endl;
					knownElementsOfCurrentSector.push_back(SG_Element);
				}
			}
		}
	}
	return true;
}



/* ########################################################################################
######   CreateCentralizerElements()                                                 ######
######                                                                               ######
######   Version: 13.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   For each constructing element find all centralizer elements from the set    ######
######   "SG_SetOfElements".                                                         ######
######################################################################################## */
bool CSpaceGroup::CreateCentralizerElements()
{
	const CSpaceGroupElement NullElement;
	vector<CSpaceGroupElement> Current_Centralizer;
	const CShiftVector NullShift;
	const CTwistVector NullTwist;

	size_t s1 = this->Sectors.size();
	size_t s2 = 0;
	unsigned i = 0;
	unsigned j = 0;
	unsigned k = 0;

	size_t s3 = this->SG_SetOfElements.size();
	if (s3 == 0)
	{
		cout << "\n  Warning in bool CSpaceGroup::CreateCentralizerElements(): Elements of the space group not defined. Return false." << endl;
		return false;
	}

	for (i = 0; i < s1; ++i)
	{
		const vector<CSpaceGroupElement> &Sector = Sectors[i];

		s2 = Sector.size();
		for (j = 0; j < s2; ++j)
		{
			Current_Centralizer.clear();

			const CSpaceGroupElement &constr_Element = Sector[j];

			// if untwisted sector
			if (constr_Element.NoTwist())
			{
				Current_Centralizer.insert(Current_Centralizer.end(), this->SG_Generators_Twist.begin(), this->SG_Generators_Twist.end());
				Current_Centralizer.insert(Current_Centralizer.end(), this->SG_Generators_Shift.begin(), this->SG_Generators_Shift.end());
			}
			// if twisted sector
			else
			{
				for (k = 0; k < s3; ++k)
				{
					const CSpaceGroupElement &SG_Element = this->SG_SetOfElements[k];
					if ((SG_Element != NullElement) && this->SG_Commute(constr_Element, SG_Element))
						Current_Centralizer.push_back(SG_Element);
				}
				// !!! find basis of centralizer....
			}

			this->Centralizer.push_back(Current_Centralizer);
		}
	}
	return true;
}



/* ########################################################################################
######   CreateDualTorus()                                                           ######
######                                                                               ######
######   Version: 16.05.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Computes six vectors e_j^* such that e_j^* e_i = \delta_{ij}.               ######
######################################################################################## */
bool CSpaceGroup::CreateDualTorus()
{
	if (this->T6_Lattice.size() != LatticeDim)
	{
		cout << "\n  Warning in bool CSpaceGroup::CreateDualTorus(): " << LatticeDim << "-torus \"T6_Lattice\" not defined. Return false." << endl;
		return false;
	}

	// Set the precision
	const double prec = 0.0001;

	unsigned i = 0;
	unsigned j = 0;

	// begin: create real basis vectors of the six-torus
	CVector         BasisVector(LatticeDim);
	vector<CVector> RealT6_Lattice(LatticeDim, BasisVector);
	for (i = 0; i < LatticeDim; ++i)
	{
		const complexVector &complexLatticeVector = this->T6_Lattice[i];

		for (j = 0; j < ComplexLatticeDim; ++j)
		{
			const complex<double> &c = complexLatticeVector[j];
			RealT6_Lattice[i][2 * j]     = c.real();
			RealT6_Lattice[i][2 * j + 1] = c.imag();
		}
	}
	// end: create real basis vectors of the six-torus

	// begin: compute the metric of the six-torus
	// Metric[i][j] = Lattice_Vector[i] * Lattice_Vector[j]
	intVector row(LatticeDim, 0);
	intMatrix Metric(LatticeDim, row);

	double tmp1 = 0.0;
	double tmp2 = 0.0;
	for (i = 0; i < LatticeDim; ++i)
	{
		const CVector &e_i = RealT6_Lattice[i];

		for (j = 0; j < LatticeDim; ++j)
		{
			tmp1 = e_i * RealT6_Lattice[j];
			tmp2 = round_double_to_int(tmp1);

			if (fabs(tmp1 - tmp2) > prec)
			{
				cout << "\n  Warning bool CSpaceGroup::CreateDualTorus(): Can not compute the metric of the T" << LatticeDim << " lattice. Return false." << endl;
				return false;
			}
			Metric[i][j] = (int)tmp2;
		}
	}
	// end: compute the metric of the six-torus

	rationalMatrix Inverse_Metric = inverse(Metric);

	// begin: compute the dual basis of the six-torus
	this->DualT6_Lattice.assign(LatticeDim, BasisVector);

	for (i = 0; i < LatticeDim; ++i)
	{
		const rationalVector &rat_row = Inverse_Metric[i];

		BasisVector.Assign(LatticeDim, 0.0);
		for (j = 0; j < LatticeDim; ++j)
		{
			const rational<int> &rat_g_ij = rat_row[j];
			const double        g_ij      = ((double)rat_g_ij.numerator())/((double)rat_g_ij.denominator());
			const CVector       &e_j      = RealT6_Lattice[j];
			BasisVector = BasisVector + (g_ij * e_j);
		}
		this->DualT6_Lattice[i] = BasisVector;
	}
	// end: compute the dual basis on the six-torus

	return true;
}



/* ########################################################################################
######   CreateSetOfElements()                                                       ######
######                                                                               ######
######   Version: 03.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Starting from the space group generators "SG_Generators_Twist" and          ######
######   "SG_Generators_Shift" create a large set of space group elements.           ######
######################################################################################## */
bool CSpaceGroup::CreateSetOfElements()
{
	if (this->SG_SetOfElements.size() != 0)
	{
		cout << "\n  Warning in bool CSpaceGroup::CreateSetOfElements(): \"SG_SetOfElements\" is not empty. Now cleared." << endl;
		this->SG_SetOfElements.clear();
	}
	const size_t s1 = this->SG_Generators_Twist.size();
	const size_t s2 = this->SG_Generators_Shift.size();
	const size_t s3 = s1 + s2;

	unsigned i = 0;
	unsigned j = 0;
	unsigned k = 0;

	const vector<rational<int> > Zero_n_alpha(LatticeDim,0);
	const CSpaceGroupElement SGIdentity;
	CSpaceGroupElement SG_Element;
	CSpaceGroupElement SG_tmpResult;

	CPrint Print(Tstandard, &cout);
	vector<int> order_of_lattice_vectors1(s2, 1);
	vector<int> order_of_lattice_vectors2(s2, 1);

	for (i = 0; i < s1; ++i)
	{
		CSpaceGroupElement SG_Generator_PurTwist = this->SG_Generators_Twist[i];
		SG_Generator_PurTwist.Set_n_alpha(Zero_n_alpha);

		for (j = 0; j < s2; ++j)
		{
			const CSpaceGroupElement &SG_Generator_Shift = this->SG_Generators_Shift[j];

			SG_Element = SG_Generator_Shift;

			// k countes how often the vector has been rotated
			k = 1;
			this->SG_Multiply(SG_Generator_PurTwist, SG_Element, SG_tmpResult);
			SG_Element = SG_tmpResult;

			while ((SG_Element.GetLatticeElement() != SG_Generator_Shift.GetLatticeElement()) && (k < 24))
			{
				++k;
				this->SG_Multiply(SG_Generator_PurTwist, SG_Element, SG_tmpResult);
				SG_Element = SG_tmpResult;
			}
			if (k == 24)
			{
				if (i == 0)
					order_of_lattice_vectors1[j] = 1;
				else
					order_of_lattice_vectors2[j] = 1;
			}
			else
			{
				if (i == 0)
					order_of_lattice_vectors1[j] = k;
				else
					order_of_lattice_vectors2[j] = k;
			}
		}
	}

	vector<int> order_of_lattice_vectors(s2, 1);
	for (i = 0; i < s2; ++i)
		order_of_lattice_vectors[i] = kgV(order_of_lattice_vectors1[i], order_of_lattice_vectors2[i]);

	// begin: recursive counting
	vector<unsigned> MaxDigits(s3, 0);
	vector<int> ShiftPureTranslations(s2, 0);

	MaxDigits[0] = (unsigned)this->M;
	if (this->ZMxZN)
		MaxDigits[1] = (unsigned)this->N;

	for (i = 0; i < LatticeDim; ++i)
	{
		if (this->WL_AllowedOrders[i] == 1)
		{
			ShiftPureTranslations[i] = -order_of_lattice_vectors[i] + 1;
			MaxDigits[i + s1]        = 2 * order_of_lattice_vectors[i] - 1;
		}
		else
		{
			ShiftPureTranslations[i] = -this->WL_AllowedOrders[i] + 1;
			MaxDigits[i + s1]        = 2 * this->WL_AllowedOrders[i] - 1;
		}
	}
	for (i = LatticeDim; i < s2; ++i)
	{
		ShiftPureTranslations[i] = -order_of_lattice_vectors[i] + 1;
		MaxDigits[i + s1]        = 2 * order_of_lattice_vectors[i] - 1;
	}
	// end: recursive counting

	CSpaceGroupElement SG_tmpElement(0,0,0);		//hacking here!!!
	CSpaceGroupElement SG_Element1;
	CSpaceGroupElement SG_Element2;

	int c = 0;
	CLatticeElement LatticeVector;

	this->SG_SetOfElements.push_back(SGIdentity);
	this->SG_SetOfElements.insert(this->SG_SetOfElements.end(), this->SG_Generators_Shift.begin(), this->SG_Generators_Shift.end());

	vector<unsigned> currentNumber(s3, 0);
	do
	{
		SG_Element1 = SGIdentity;
		for (j = 0; j < s2; ++j)
		{
			c = currentNumber[j + s1] + ShiftPureTranslations[j];

			if (c != 0)
			{
				SG_tmpElement.SetLatticeElement(this->SG_Generators_Shift[j].GetLatticeElement() * c);

				this->SG_Multiply(SG_tmpElement, SG_Element1, SG_tmpResult);
				SG_Element1 = SG_tmpResult;
			}
		}

		SG_Element2 = SGIdentity;
		for (j = 0; j < s1; ++j)
		{
			c = currentNumber[j];
			if (c != 0)
			{
				const CSpaceGroupElement &SG_Generator = this->SG_Generators_Twist[j];
				for (k = 0; k < (unsigned)abs(c); ++k)
				{
					this->SG_Multiply(SG_Generator, SG_Element2, SG_tmpResult);
					SG_Element2 = SG_tmpResult;
				}
			}
		}

		this->SG_Multiply(SG_Element1, SG_Element2, SG_Element);
		if (!SG_Element.NoTwist())
		{
			if (find(this->SG_SetOfElements.begin(), this->SG_SetOfElements.end(), SG_Element) == this->SG_SetOfElements.end())
				this->SG_SetOfElements.push_back(SG_Element);

			this->SG_Multiply(SG_Element2, SG_Element1, SG_Element);
			if (find(this->SG_SetOfElements.begin(), this->SG_SetOfElements.end(), SG_Element) == this->SG_SetOfElements.end())
				this->SG_SetOfElements.push_back(SG_Element);
		}
	}
	while (NextNumber(currentNumber, MaxDigits, s3));

	if (this->SG_SetOfElements.size() == 0)
	{
		cout << "\n  Warning in bool CSpaceGroup::CreateSetOfElements(): \"SG_SetOfElements\" is empty. Return false." << endl;
		return false;
	}
	return true;
}



/* ########################################################################################
######   CreateTwistMatrices()                                                       ######
######                                                                               ######
######   Version: 16.05.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   For all twisted sectors computes the matrices of the twists in the lattice  ######
######   basis, e.g. computes M_ji with \theta e_i = M_ji e_j.                       ######
######################################################################################## */
bool CSpaceGroup::CreateTwistMatrices()
{
	if (this->T6_Lattice.size() != LatticeDim)
	{
		cout << "\n  Warning in bool CSpaceGroup::CreateTwistMatrices(): Six-torus \"T6_Lattice\" not defined. Return false." << endl;
		return false;
	}
	if (this->DualT6_Lattice.size() != 6)
	{
		cout << "\n  Warning in bool CSpaceGroup::CreateTwistMatrices(): Dual six-torus \"DualT6_Lattice\" not defined. Return false." << endl;
		return false;
	}
	if (this->TwistMatrices.size() != 0)
	{
		cout << "\n  Warning in bool CSpaceGroup::CreateTwistMatrices(): \"TwistMatrices\" is not empty. Now cleared." << endl;
		this->TwistMatrices.clear();
	}

	unsigned i = 0;
	const rationalVector tmp_line(LatticeDim, 0);
	rationalMatrix IdentityMatrix(LatticeDim, tmp_line);
	for (i = 0; i < LatticeDim; ++i)
		IdentityMatrix[i][i] = 1;

	if ( (this->ZMxZN) or (this->ZMxZNxZK) )			//the Z2W need not be considered here, since trivial in geometry
	{
		// Set the precision
		const double prec = 0.0001;
		unsigned j = 0;

		rationalMatrix ZNTwistMatrix(LatticeDim, tmp_line);
		rationalMatrix ZKTwistMatrix(LatticeDim, tmp_line);

		complex<double> e(0);
		CVector         BasisVector(LatticeDim);
		double tmp1 = 0.0;
		double tmp2 = 0.0;

		const double TwoPi = 2.0 * M_PI;

		for (i = 0; i < LatticeDim; ++i)
		{
			const complexVector &e_i = this->T6_Lattice[i];

			// compute the theta rotated basis vector of the six-torus
			for (j = 0; j < ComplexLatticeDim; ++j)
			{
				e = polar(1.0, TwoPi * this->Twists[1][j+1]) * e_i[j];			//hacking here!!!

				BasisVector[2 * j]     = e.real();
				BasisVector[2 * j + 1] = e.imag();
			}

			for (j = 0; j < LatticeDim; ++j)
			{
				tmp1 = BasisVector * this->DualT6_Lattice[j];
				tmp2 = round_double_to_int(tmp1);

				if (fabs(tmp1 - tmp2) > prec)
				{
					cout << "\n  Warning in bool CSpaceGroup::CreateTwistMatrices(): Rotated basis vector of the six-lattice is not in the lattice. Return false." << endl;
					return false;
				}
				ZNTwistMatrix[j][i] = (int)tmp2;									//hacking here!!!
			}
		}

		if (this->ZMxZNxZK)
		{
			for (i = 0; i < LatticeDim; ++i)
			{
				const complexVector &e_i = this->T6_Lattice[i];

				// compute the omega rotated basis vector of the six-torus
				for (j = 0; j < ComplexLatticeDim; ++j)
				{
					e = polar(1.0, TwoPi * this->Twists[2][j+1]) * e_i[j];			//hacking here!!!

					BasisVector[2 * j]     = e.real();
					BasisVector[2 * j + 1] = e.imag();
				}

				for (j = 0; j < LatticeDim; ++j)
				{
					tmp1 = BasisVector * this->DualT6_Lattice[j];
					tmp2 = round_double_to_int(tmp1);

					if (fabs(tmp1 - tmp2) > prec)
					{
						cout << "\n  Warning in bool CSpaceGroup::CreateTwistMatrices(): Rotated basis vector of the six-lattice is not in the lattice. Return false." << endl;
						return false;
					}
					ZKTwistMatrix[j][i] = (int)tmp2;								//hacking here!!!
				}
			}
		}

		vector<rationalMatrix> tmp1_TwistMatrices(this->K, IdentityMatrix);			//hacking here!!!
		this->TwistMatrices.assign(this->N, tmp1_TwistMatrices);						//hacking here!!!

		unsigned n = 0;
		unsigned k = 0;
		unsigned m = 0;

		// n = 0
		for (k = 1; k < this->K; ++k)
		{
			rationalMatrix &PreviousTwistMatrix = this->TwistMatrices[0][k-1];
			rationalMatrix &TwistMatrix         = this->TwistMatrices[0][k];

			for (i = 0; i < LatticeDim; ++i)
			{
				for (j = 0; j < LatticeDim; ++j)
				{
					TwistMatrix[i][j] = 0;
					for (m = 0; m < LatticeDim; ++m)
						TwistMatrix[i][j] += PreviousTwistMatrix[m][j] * ZKTwistMatrix[i][m];
				}
			}
		}
		for (n = 1; n < this->N; ++n)
		{
			for (k = 0; k < this->K; ++k)
			{
				rationalMatrix &PreviousTwistMatrix = this->TwistMatrices[n-1][k];
				rationalMatrix &TwistMatrix         = this->TwistMatrices[n][k];

				for (i = 0; i < LatticeDim; ++i)
				{
					for (j = 0; j < LatticeDim; ++j)
					{
						TwistMatrix[i][j] = 0;
						for (m = 0; m < LatticeDim; ++m)
							TwistMatrix[i][j] += PreviousTwistMatrix[m][j] * ZNTwistMatrix[i][m];
					}
				}
			}
		}
	}
	else {																			//case Z2W only
		vector<rationalMatrix> tmp1_TwistMatrices(1, IdentityMatrix);				//hacking here!!!
		this->TwistMatrices.assign(2, tmp1_TwistMatrices);
	}

	/*cout << "Twist Matrices\n";
  for (k = 0; k < this->M; ++k)
  {
    for (l = 0; l < this->N; ++l)
    {
      cout << "(k, l) = (" << k << ", " << l << ")" << endl;
      const rationalMatrix &TwistMatrix = TwistMatrices[k][l];
      for (i = 0; i < LatticeDim; ++i)
      {
        for (j = 0; j < LatticeDim; ++j)
          cout << TwistMatrix[i][j] << " ";
        cout << endl;
      }
    }
}*/
	return true;
}




/* ########################################################################################
######   FindRelationsOfWilsonLines()                                                ######
######                                                                               ######
######   Version: 05.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Tries to find all relations between the Wilson lines on the orbifold, e.g.  ######
######   W_1 = W_2 for the Z_3 orbifold.                                             ######
######################################################################################## */
bool CSpaceGroup::FindRelationsOfWilsonLines()
{
	if (this->TwistMatrices.size() != this->M)
	{
		cout << "\n  Warning in bool CSpaceGroup::FindRelationsOfWilsonLines(): TwistMatrices not defined. Return false." << endl;
		return false;
	}

	unsigned k = 0;
	unsigned l = 0;

	unsigned i = 0;
	unsigned j = 0;
	unsigned m = 0;

	const size_t s1 = this->SG_Generators_Twist.size();
	const size_t s2 = this->SG_Generators_Shift.size();

	const vector<rational<int> > ZeroRelation(LatticeDim, 0);
	rationalVector Relation(LatticeDim, 0);
	vector<vector<rational<int> > > AllRelations;

	for (k = 0; k < this->M; ++k)
	{
		const vector<rationalMatrix> &k_TwistMatrices = this->TwistMatrices[k];
		if (k_TwistMatrices.size() != this->N)
		{
			cout << "\n  Warning in bool CSpaceGroup::FindRelationsOfWilsonLines(): TwistMatrices not defined. Return false." << endl;
			return false;
		}

		for (l = 0; l < this->N; ++l)
		{
			// no relation between Wilson lines for identity rotation
			if ((k != 0) || (l != 0))
			{
				const rationalMatrix &kl_TwistMatrices = k_TwistMatrices[l];

				// begin: run through the generators of the space group
				for (i = 0; i < s1; ++i)
				{
					const CLatticeElement &n_alpha = this->SG_Generators_Twist[i].GetLatticeElement();

					Relation = n_alpha.Vector();

					for (j = 0; j < LatticeDim; ++j)
					{
						rational<int> &n_j = Relation[j];
						const rationalVector &TwistVector = kl_TwistMatrices[j];

						for (m = 0; m < LatticeDim; ++m)
							n_j -= TwistVector[m] * n_alpha[m];
					}
					if (Relation != ZeroRelation)
					{
						if (find(AllRelations.begin(), AllRelations.end(), Relation) == AllRelations.end())
							AllRelations.push_back(Relation);
					}
				}
				for (i = 0; i < s2; ++i)
				{
					const CLatticeElement &n_alpha = this->SG_Generators_Shift[i].GetLatticeElement();

					Relation = n_alpha.Vector();

					for (j = 0; j < LatticeDim; ++j)
					{
						rational<int> &n_j = Relation[j];
						const rationalVector &TwistVector = kl_TwistMatrices[j];

						for (m = 0; m < LatticeDim; ++m)
							n_j -= TwistVector[m] * n_alpha[m];
					}
					if (Relation != ZeroRelation)
					{
						if (find(AllRelations.begin(), AllRelations.end(), Relation) == AllRelations.end())
							AllRelations.push_back(Relation);
					}
				}
				// end: run through the generators of the space group
			}
		}
	}

	this->WL_AllowedOrders.assign(LatticeDim,0);
	vector<vector<rational<int> > > OrdersMatrix;

	CLinAlg<int> LA;
	LA.FindOrdersOfWilsonLines(AllRelations, OrdersMatrix);

	unsigned pos1 = 0;
	unsigned pos2 = 0;
	unsigned counter = 0;
	const size_t s3 = OrdersMatrix.size();
	for (i = 0; i < s3; ++i)
	{
		const vector<rational<int> > &ith_Order = OrdersMatrix[i];
		if (ith_Order != ZeroRelation)
		{
			counter = 0;
			for (j = 0; j < LatticeDim; ++j)
			{
				if (ith_Order[j] != 0)
				{
					if (ith_Order[j].denominator() != 1)
						counter += 2;
					else
					{
						pos1 = j;
						++counter;
					}
				}
			}
			if (counter == 1)
				this->WL_AllowedOrders[pos1] = (unsigned)ith_Order[pos1].numerator();
			else
			{
				cout << "\n  Warning in bool CSpaceGroup::FindRelationsOfWilsonLines(): Order of Wilson line not identified. Return false." << endl;
				return false;
			}
		}
	}
	for (i = 0; i < LatticeDim; ++i)
	{
		if (this->WL_AllowedOrders[i] == 0)
		{
			cout << "\n  Warning in bool CSpaceGroup::FindRelationsOfWilsonLines(): Order of Wilson line not identified. Return false." << endl;
			return false;
		}
	}

	this->WL_Relations.clear();

	bool AllOrdersOne = true;
	for (i = 0; AllOrdersOne && (i < LatticeDim); ++i)
	{
		if (this->WL_AllowedOrders[i] != 1)
			AllOrdersOne = false;
	}
	if (AllOrdersOne)
		return true;

	unsigned number_of_ones      = 0;
	unsigned number_of_minusones = 0;
	unsigned number_of_other     = 0;

	bool pos1_found = false;
	bool pos2_found = false;
	unsigned set_of_pos1 = 0;
	unsigned set_of_pos2 = 0;

	bool need_new_Identification = true;
	bool first_one_found = false;
	size_t s5 = 0;

	vector<vector<rational<int> > > DifficultRelations;
	vector<vector<rational<int> > > UsedRelations;

	const size_t s4 = AllRelations.size();
	for (i = 0; i < s4; ++i)
	{
		Relation = AllRelations[i];

		first_one_found = false;

		// begin: use the orders of the Wilson lines to simplify the relation
		for (j = 0; j < LatticeDim; ++j)
		{
			if (this->WL_AllowedOrders[j] == 1)
				Relation[j] = 0;

			if (Relation[j] != 0)
			{
				if (Relation[j] > 1)
				{
					while (Relation[j] > 1)
						Relation[j] -= this->WL_AllowedOrders[j];
				}
				else
				{
					while (Relation[j] < -1)
						Relation[j] += this->WL_AllowedOrders[j];
				}

				if (((Relation[j] == 1) || (Relation[j] == -1)) && (this->WL_AllowedOrders[j] == 2))
				{
					if (first_one_found && (Relation[j] == 1))
						Relation[j] -= 2;

					if (!first_one_found)
					{
						if (Relation[j] == -1)
							Relation[j] += 2;

						first_one_found = true;
					}
				}
			}
		}
		// end: use the orders of the Wilson lines to simplify the relation

		if (Relation != ZeroRelation)
		{
			// begin: check number of +1, -1
			number_of_ones      = 0;
			number_of_minusones = 0;
			number_of_other     = 0;
			for (j = 0; j < LatticeDim; ++j)
			{
				if (Relation[j] != 0)
				{
					if (Relation[j] == 1)
						++number_of_ones;
					else
					{
						if (Relation[j] == -1)
							++number_of_minusones;
						else
							++number_of_other;
					}
				}
			}
			// end: check number of +1, -1

			if ((number_of_ones == 1) && (number_of_minusones == 1) && (number_of_other == 0))
			{
				for (j = 0; j < LatticeDim; ++j)
				{
					if (Relation[j] == 1)
						pos1 = j;
					else
						if (Relation[j] == -1)
							pos2 = j;
				}

				pos1_found = false;
				pos2_found = false;
				need_new_Identification = true;

				s5 = this->WL_Relations.size();
				for (k = 0; need_new_Identification && (k < s5); ++k)
				{
					const vector<unsigned> &Identification = this->WL_Relations[k];

					if (find(Identification.begin(), Identification.end(), pos1) != Identification.end())
					{
						if (find(Identification.begin(), Identification.end(), pos2) != Identification.end())
							need_new_Identification = false;
						else
						{
							set_of_pos1 = k;
							pos1_found = true;
						}
					}

					if (need_new_Identification && find(Identification.begin(), Identification.end(), pos2) != Identification.end())
					{
						set_of_pos2 = k;
						pos2_found = true;
					}
				}
				if (need_new_Identification)
				{
					if (pos1_found)
					{
						if (pos2_found)
						{
							vector<unsigned> &Identification1 = this->WL_Relations[set_of_pos1];
							vector<unsigned> &Identification2 = this->WL_Relations[set_of_pos2];

							Identification1.insert(Identification1.end(), Identification2.begin(), Identification2.end());
							this->WL_Relations.erase (this->WL_Relations.begin() + set_of_pos2);
							UsedRelations.push_back(Relation);
						}
						else
						{
							this->WL_Relations[set_of_pos1].push_back(pos2);
							UsedRelations.push_back(Relation);
						}
					}
					else
					{
						if (pos2_found)
						{
							this->WL_Relations[set_of_pos2].push_back(pos1);
							UsedRelations.push_back(Relation);
						}
						else
						{
							vector<unsigned> new_Identification;
							new_Identification.push_back(pos1);
							new_Identification.push_back(pos2);

							this->WL_Relations.push_back(new_Identification);
							UsedRelations.push_back(Relation);
						}
					}
				}
			}
			else
				DifficultRelations.push_back(Relation);
		}
	}

	s5 = this->WL_Relations.size();
	for (i = 0; i < s5; ++i)
	{
		vector<unsigned> &Identification = this->WL_Relations[i];
		sort(Identification.begin(), Identification.end());
	}

	const size_t s6 = DifficultRelations.size();
	if (s6 == 0)
		return true;

	vector<vector<rational<int> > > StillDifficultRelations;

	size_t s7 = 0;
	for (i = 0; i < s6; ++i)
	{
		vector<rational<int> > &DifficultRelation = DifficultRelations[i];

		for (j = 0; j < s5; ++j)
		{
			vector<unsigned> &Identification = this->WL_Relations[j];

			s7 = Identification.size();
			for (k = 1; k < s7; ++k)
			{
				DifficultRelation[Identification[0]] += DifficultRelation[Identification[k]];
				DifficultRelation[Identification[k]] = 0;
			}
		}

		// begin: use the orders of the Wilson lines to simplify the relation
		for (j = 0; j < LatticeDim; ++j)
		{
			if (this->WL_AllowedOrders[j] == 1)
				DifficultRelation[j] = 0;

			if (DifficultRelation[j] != 0)
			{
				if (DifficultRelation[j] > 1)
				{
					while (DifficultRelation[j] > 1)
						DifficultRelation[j] -= this->WL_AllowedOrders[j];
				}
				else
				{
					while (DifficultRelation[j] < -1)
						DifficultRelation[j] += this->WL_AllowedOrders[j];
				}
			}
		}
		// end: use the orders of the Wilson lines to simplify the relation

		if (DifficultRelation != ZeroRelation)
			StillDifficultRelations.push_back(DifficultRelation);
	}

	if (StillDifficultRelations.size() == 0)
		return true;

	vector<vector<rational<int> > > C = UsedRelations;
	const size_t k1 = findBasis<rational<int> >(C).size();

	C.insert(C.end(), StillDifficultRelations.begin(), StillDifficultRelations.end());
	const size_t k2 = findBasis<rational<int> >(C).size();

	if (k1 == k2)
		return true;

	cout << "\n  Warning in bool CSpaceGroup::FindRelationsOfWilsonLines(): Cannot interpret all relations between Wilson lines automatically. Return false." << endl;
	return false;
}



/* ########################################################################################
######   LoadSpaceGroup(string ifilename, bool reload)                               ######
######                                                                               ######
######   Version: 12.01.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) ifilename : filename                                                     ######
######   2) reload    : if the file is loaded again only some parts (like the        ######
######                  list of constructing elements) have to be loaded.            ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Loads the space group from file named "ifilename".                          ######
######################################################################################## */
bool CSpaceGroup::LoadSpaceGroup(string ifilename, bool reload)
{
	this->SpaceGroup_CheckStatus = NotChecked;

	// Set the precision
	const double prec = 0.0001;
	const string AllPositiveDigits = "0123456789";
	const string AllDigits         = "0123456789-";
	const string AllRealNumbers    = "0123456789.+-";
	const string AllCharacters     = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890!@#$%^&*()-_=+[{]};:',<.>/?| ";

	std::ifstream in;
	in.open(ifilename.data(), ifstream::in);
	if(!in.is_open() || !in.good())
	{
		cout << "\n  Warning! Could not find the file " << ifilename << " for the elements of the space group." << endl;
		return false;
	}

	if (!reload)
	{
		this->T6_Lattice.clear();
		this->DualT6_Lattice.clear();

		this->additional_label = "";
		this->lattice_label = "";
		this->GeometryFilename = ifilename;

		this->TwistMatrices.clear();

		CTwistVector ZeroTwist;
		this->Twists.assign(3, ZeroTwist);				//hacking here!!!

		this->SG_Generators_Twist.clear();
		this->SG_Generators_Shift.clear();

		this->DiscreteRSymmetries.clear();
		this->DiscreteNonRSymmetries.clear();
		this->ModularSymmetries.clear();

		this->SG_AllNonStandardShifts.clear();
		this->SG_AllRotations.clear();
		this->SG_SetOfElements.clear();
	}

	this->WL_AllowedOrders.assign(LatticeDim,0);
	this->WL_Relations.clear();
	this->ShiftsWL_ScalarProductFactors.clear();
	this->Sectors.clear();
	this->Centralizer.clear();

	unsigned i = 0;
	unsigned j = 0;

	size_t s1 = 0;

	string currentline = "";

	vector<rationalVector> RationalTwistVectors;
	rationalVector RationalVector;
	CSpaceGroupElement SGElement;
	vector<CSpaceGroupElement> New_Sector;

	const unsigned number_of_NecessaryData = 11;

	vector<bool> NecessaryDataRead(number_of_NecessaryData, false);
	if (reload)
	{
		NecessaryDataRead[0] = true; // point group
		NecessaryDataRead[1] = true; // ZMxZN
		NecessaryDataRead[10] = true; // ZMxZNxZK
		NecessaryDataRead[2] = true; // twist
		NecessaryDataRead[3] = true; // lattice
		NecessaryDataRead[4] = true; // space group generators with twist
		NecessaryDataRead[5] = true; // space group generators without twist
	}

	bool AutoCreateConstructingElements = false;
	bool AutoCreateCentralizerElements  = false;
	bool AutoFindRelations              = false;

	while (GetSaveLine(in, currentline))
	{
		// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if (!reload && !NecessaryDataRead[0] && (currentline == "point group"))
		{
			// get next line
			GetSaveLine(in, currentline);
			if (currentline.find_first_not_of(AllPositiveDigits) != string::npos)
			{
				cout << "\n  Warning in bool CSpaceGroup::LoadSpaceGroup(...): order M of Z_M ill-defined. Return false." << endl;
				return false;
			}
			// get the order M
			std::istringstream line1(currentline);
			line1 >> this->M;

			// get next line
			GetSaveLine(in, currentline);
			if (currentline.find_first_not_of(AllPositiveDigits) != string::npos)
			{
				cout << "\n  Warning in bool CSpaceGroup::LoadSpaceGroup(...): order N of Z_N ill-defined. Return false." << endl;
				return false;
			}
			// get the order N
			std::istringstream line2(currentline);
			line2 >> this->N;

			// get next line
			GetSaveLine(in, currentline);
			if (currentline.find_first_not_of(AllPositiveDigits) != string::npos)
			{
				cout << "\n  Warning in bool CSpaceGroup::LoadSpaceGroup(...): order K of Z_K ill-defined. Return false." << endl;
				return false;
			}
			// get the order K
			std::istringstream line3(currentline);
			line3 >> this->K;

			NecessaryDataRead[0] = true;
		}
		else
			// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			if (!reload && !NecessaryDataRead[1] && (currentline == "ZMxZN"))
			{
				// get next line
				GetSaveLine(in, currentline);
				if (currentline == "true")
				{
					NecessaryDataRead[1] = true;
					this->ZMxZN = true;
				}
				else
					if (currentline == "false")
					{
						NecessaryDataRead[1] = true;
						this->ZMxZN = false;
					}
			}
			else
				// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				if (!reload && !NecessaryDataRead[10] && (currentline == "ZMxZNxZK"))
				{
					// get next line
					GetSaveLine(in, currentline);
					if (currentline == "true")
					{
						NecessaryDataRead[10] = true;
						this->ZMxZNxZK = true;
					}
					else
						if (currentline == "false")
						{
							NecessaryDataRead[10] = true;
							this->ZMxZNxZK = false;
						}
				}
				else
					// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					if (!reload && (currentline == "additional label"))
					{
						// get next line
						GetSaveLine(in, currentline);
						// and save the label
						this->additional_label = currentline.substr(0, currentline.find_first_not_of(AllCharacters));
					}
					else
						// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
						if (!reload && !NecessaryDataRead[2] && ((currentline == "begin twist") || (currentline == "begin twists")))
						{
							CTwistVector tmpTwist;

							j = 0;
							// get next lines until "end twist"
							while (GetSaveLine(in, currentline) && (currentline != "end twist") && (currentline != "end twists"))
							{
								convert_string_to_vector_of_rational(currentline, RationalVector);
								if (RationalVector.size() != 4)
								{
									cout << "\n  Warning in bool CSpaceGroup::LoadSpaceGroup(...): Check length of twist vector. Return false." << endl;
									return false;
								}
								for (i = 0; i < 4; ++i)
								{
									const rational<int> &tmp = RationalVector[i];
									tmpTwist[i] = ((double)tmp.numerator())/((double)tmp.denominator());
								}

								RationalTwistVectors.push_back(RationalVector);
								this->Twists[j] = tmpTwist;
								++j;
							}
							NecessaryDataRead[2] = true;
						}
						else
							// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
							if (!reload && (currentline == "lattice label"))
							{
								// get next line
								GetSaveLine(in, currentline);
								// and save the label
								this->lattice_label = currentline.substr(0, currentline.find_first_not_of(AllCharacters));
							}
							else
								// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
								if (!reload && !NecessaryDataRead[3] && (currentline == "begin lattice"))
								{
									double real = 0.0;
									double imag = 0.0;
									unsigned counter = 0;

									complexVector vec;

									// get next lines until "end lattice"
									while (GetSaveLine(in, currentline) && (currentline != "end lattice"))
									{
										if (currentline != "next lattice vector")
										{
											// read the real part of the lattice vector
											if (counter == 0)
											{
												if (currentline.find_first_not_of(AllRealNumbers) != string::npos)
												{
													cout << "\n  Warning in bool CSpaceGroup::LoadSpaceGroup(...): real part of lattice vector ill-defined. Return false." << endl;
													return false;
												}
												std::istringstream line(currentline);
												line >> real;
												counter = 1;
											}
											else
												// read the imaginary part of the lattice vector
												if (counter == 1)
												{
													if (currentline.find_first_not_of(AllRealNumbers) != string::npos)
													{
														cout << "\n  Warning in bool CSpaceGroup::LoadSpaceGroup(...): imaginary part of lattice vector ill-defined. Return false." << endl;
														return false;
													}
													std::istringstream line(currentline);
													line >> imag;
													counter = 0;
												}

											// save the component of the lattice vector
											if (counter == 0)
											{
												complex<double> tmp(real, imag);

												// begin: round the entries
												if (fabs(tmp.real()) < prec)
												{
													complex<double> tmp2(0, tmp.imag());
													tmp = tmp2;
												}
												if (fabs(tmp.imag()) < prec)
												{
													complex<double> tmp2(tmp.real(), 0);
													tmp = tmp2;
												}
												// end: round the entries

												vec.push_back(tmp);

												// save the lattice vector
												if (vec.size() == ComplexLatticeDim)
												{
													this->T6_Lattice.push_back(vec);
													vec.clear();
												}
											}
										}
									}
									if (this->T6_Lattice.size() == LatticeDim)
										NecessaryDataRead[3] = true;
								}
								else
									// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
									if (!reload && !NecessaryDataRead[4] && (currentline == "begin twist space group generators"))
									{
										// get next lines until "end space group generators"
										while (GetSaveLine(in, currentline) && (currentline != "end twist space group generators"))
										{
											if (!global_LoadSGElement(currentline, SGElement))
											{
												cout << "\n  Warning in bool CSpaceGroup::LoadSpaceGroup(...): Space group element ill specified. Return false." << endl;
												return false;
											}

											if (SGElement.NoTwist())
											{
												cout << "\n  Warning in bool CSpaceGroup::LoadSpaceGroup(...): Twist space group generator without twist is ill specified. Return false." << endl;
												return false;
											}
											this->SG_Generators_Twist.push_back(SGElement);
										}
										if (this->SG_Generators_Twist.size() != 0)
											NecessaryDataRead[4] = true;
									}
									else
										// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
										if (!reload && !NecessaryDataRead[5] && (currentline == "begin shift space group generators"))
										{
											// get next lines until "end space group generators"
											while (GetSaveLine(in, currentline) && (currentline != "end shift space group generators"))
											{
												if (!global_LoadSGElement(currentline, SGElement))
												{
													cout << "\n  Warning in bool CSpaceGroup::LoadSpaceGroup(...): Space group element ill specified. Return false." << endl;
													return false;
												}

												if (!SGElement.NoTwist())
												{
													cout << "\n  Warning in bool CSpaceGroup::LoadSpaceGroup(...): Shift space group generator with twist is ill specified. Return false." << endl;
													return false;
												}
												this->SG_Generators_Shift.push_back(SGElement);
											}
											if (this->SG_Generators_Shift.size() >= LatticeDim)
												NecessaryDataRead[5] = true;
										}
										else
											// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
											if (currentline == "find Wilson lines relations")
											{
												AutoFindRelations = true;
												NecessaryDataRead[6] = true;
												NecessaryDataRead[7] = true;
											}
											else
												// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
												if (!NecessaryDataRead[6] && (currentline == "allowed order of Wilson lines"))
												{
													if (AutoFindRelations)
														cout << "\n  Warning in bool CSpaceGroup::LoadSpaceGroup(...): using \"find Wilson lines relations\"." << endl;

													for (i = 0; i < LatticeDim; ++i)
													{
														// get next line
														GetSaveLine(in, currentline);
														if (currentline.find_first_not_of(AllPositiveDigits) != string::npos)
														{
															cout << "\n  Warning in bool CSpaceGroup::LoadSpaceGroup(...): order of Wilson line ill-defined. Return false." << endl;
															return false;
														}

														std::istringstream line1(currentline);
														if (!AutoFindRelations)
															line1 >> this->WL_AllowedOrders[i];
													}
													NecessaryDataRead[6] = true;
												}
												else
													// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
													if (!NecessaryDataRead[7] && (currentline == "begin identical WLs"))
													{
														if (AutoFindRelations)
															cout << "\n  Warning in bool CSpaceGroup::LoadSpaceGroup(...): using \"find Wilson lines relations\"." << endl;

														unsigned temp = 0;
														vector<unsigned> data;

														GetSaveLine(in, currentline);
														if (currentline == "none")
														{
															NecessaryDataRead[7] = true;
															AutoFindRelations = false;
															GetSaveLine(in, currentline); // read "end identical WLs"
														}
														else
														{
															// get next lines until "end identical WLs"
															while (currentline != "end identical WLs")
															{
																data.clear();

																std::istringstream line(currentline);
																while (line >> temp)
																{
																	if ((temp >= 1) && (temp <= LatticeDim))
																		data.push_back((unsigned)temp - 1);
																	else
																	{
																		cout << "\n  Warning in bool CSpaceGroup::LoadSpaceGroup(...): Cannot read identical Wilson lines: index of Wilson line must be 1,...," << LatticeDim << ".";
																		return false;
																	}
																}

																if (!AutoFindRelations)
																	this->WL_Relations.push_back(data);
																GetSaveLine(in, currentline);
															}
															NecessaryDataRead[7] = true;
														}
													}
													else
														// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
														if (!reload && (currentline == "begin discrete symmetries"))
														{
															bool case_gen  = false;
															bool case_R    = false;
															bool case_nonR = false;
															bool modular   = false;
															while (GetSaveLine(in, currentline) && (currentline != "end discrete symmetries"))
															{
																case_gen  = false;
																case_R    = false;
																case_nonR = false;
																if ((case_gen = (currentline == "new symmetry")) || (case_R = (currentline == "new R symmetry")) || (case_nonR = (currentline == "new non-R symmetry")) || (modular = (currentline == "new modular symmetry")))
																{
																	if (modular)
																	{
																		//struct SModularSymmetry
																		//{
																		//  rational<int>  Const;
																		//  unsigned       Index;
																		//  rationalVector ChargeOperator; // 3 components: Index-component of q_sh, N and \bar{N}
																		//};
																		SModularSymmetry NewSymmetry;

																		// get the label
																		GetSaveLine(in, currentline);
																		NewSymmetry.Label = currentline;

																		// get the constant
																		GetSaveLine(in, currentline);
																		convert_string_to_vector_of_rational(currentline, RationalVector);
																		if (RationalVector.size() != 1)
																		{
																			cout << "\n  Warning in bool CSpaceGroup::LoadSpaceGroup(...): \"Const\" of modular symmetry ill-defined. Return false." << endl;
																			return false;
																		}
																		NewSymmetry.Const = RationalVector[0];

																		// get the index
																		GetSaveLine(in, currentline);
																		if (currentline.find_first_not_of(AllPositiveDigits) != string::npos)
																		{
																			cout << "\n  Warning in bool CSpaceGroup::LoadSpaceGroup(...): \"Index\" of modular symmetry ill-defined. Return false." << endl;
																			return false;
																		}
																		std::istringstream line1(currentline);
																		line1 >> NewSymmetry.Index;
																		if ((NewSymmetry.Index != 1) && (NewSymmetry.Index != 2) && (NewSymmetry.Index != 3))
																		{
																			cout << "\n  Warning in bool CSpaceGroup::LoadSpaceGroup(...): \"Index\" of modular symmetry ill-defined. Return false." << endl;
																			return false;
																		}

																		// get the charge operator
																		GetSaveLine(in, currentline);
																		convert_string_to_vector_of_rational(currentline, RationalVector);
																		if (RationalVector.size() != 3)
																		{
																			cout << "\n  Warning in bool CSpaceGroup::LoadSpaceGroup(...): \"ChargeOperator\" of modular symmetry ill-defined. Return false." << endl;
																			return false;
																		}
																		NewSymmetry.ChargeOperator = RationalVector;
																		this->ModularSymmetries.push_back(NewSymmetry);
																	}
																	else
																	{
																		//struct SDiscreteSymmetry
																		//{
																		//  unsigned               Order;
																		//  rational<int>          SuperpotentialCharge;
																		//  rationalVector         ChargeOperator;
																		//  CSpaceGroupElement     SG_SymmetryGenerator;
																		//};
																		SDiscreteSymmetry NewSymmetry;

																		// order of the discrete symmetry
																		GetSaveLine(in, currentline);
																		if (currentline.find_first_not_of(AllPositiveDigits) != string::npos)
																		{
																			cout << "\n  Warning in bool CSpaceGroup::LoadSpaceGroup(...): order of discrete symmetry ill-defined. Return false." << endl;
																			return false;
																		}
																		std::istringstream line1(currentline);
																		line1 >> NewSymmetry.Order;

																		// (R) charge of the superpotential
																		GetSaveLine(in, currentline);
																		std::istringstream line(currentline);
																		line >> NewSymmetry.SuperpotentialCharge;

																		const bool WCharged = (NewSymmetry.SuperpotentialCharge.numerator() != 0);

																		if (!case_gen && ((case_R && !WCharged) || (case_nonR && WCharged)))
																		{
																			cout << "\n  Warning in bool CSpaceGroup::LoadSpaceGroup(...): R-charge of superpotential for discrete symmetry is ill-defined. Return false." << endl;
																			return false;
																		}

																		// the charge operator
																		GetSaveLine(in, currentline);
																		convert_string_to_vector_of_rational(currentline, RationalVector);
																		if (WCharged)
																		{
																			if (!case_gen && (RationalVector.size() != 4))
																			{
																				cout << "\n  Warning in bool CSpaceGroup::LoadSpaceGroup(...): \"ChargeOperator\" for R-type discrete symmetry ill-defined. Return false." << endl;
																				return false;
																			}
																		}
																		else
																		{
																			if (!case_gen && (RationalVector.size() != 8))
																			{
																				cout << "\n  Warning in bool CSpaceGroup::LoadSpaceGroup(...): \"ChargeOperator\" for non-R-type discrete symmetry ill-defined. Return false." << endl;
																				return false;
																			}
																		}
																		NewSymmetry.ChargeOperator = RationalVector;

																		// the space-group element generating the symmetry in the non-R case
																		/*if (case_nonR)			//hacking here!!! disabled for ZMxZNxZK orbifolds
																		{
																			GetSaveLine(in, currentline);

																			if (!global_LoadSGElement(currentline, SGElement))
																			{
																				cout << "\n  Warning in bool CSpaceGroup::LoadSpaceGroup(...): \"SG_SymmetryGenerator\" for non-R-type discrete symmetry ill specified. Return false." << endl;
																				return false;
																			}

																			NewSymmetry.SG_SymmetryGenerator = SGElement;
																		}*/

																		/*CPrint Print(Tstandard, &cout);
            Print.PrintDiscreteSymmetry(NewSymmetry);
            cout << endl;*/

																		if (WCharged)
																			this->DiscreteRSymmetries.push_back(NewSymmetry);
																		else
																			this->DiscreteNonRSymmetries.push_back(NewSymmetry);
																	}
																}
															}
														}
														else
															// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
															if (!NecessaryDataRead[8] && (currentline == "create constructing elements"))
															{
																AutoCreateConstructingElements = true;
																NecessaryDataRead[8] = true;
															}
															else
																// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
																if (!NecessaryDataRead[8] && (currentline == "begin constructing elements"))
																{
																	AutoCreateConstructingElements = false;

																	// next, read the constructing elements of the space group
																	while (GetSaveLine(in, currentline) && (currentline != "end constructing elements"))
																	{
																		if (currentline == "new sector")
																		{
																			if (New_Sector.size() != 0)
																			{
																				NecessaryDataRead[8] = true;
																				this->Sectors.push_back(New_Sector);
																				New_Sector.clear();
																			}
																		}
																		else
																		{
																			if (!global_LoadSGElement(currentline, SGElement))
																			{
																				cout << "\n  Warning in bool CSpaceGroup::LoadSpaceGroup(...): Space group element ill specified. Return false." << endl;
																				return false;
																			}

																			New_Sector.push_back(SGElement);
																		}
																	}
																	if (New_Sector.size() != 0)
																	{
																		NecessaryDataRead[8] = true;
																		this->Sectors.push_back(New_Sector);
																		New_Sector.clear();
																	}
																}
																else
																	// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
																	if (!NecessaryDataRead[9] && (currentline == "create centralizer elements"))
																	{
																		AutoCreateCentralizerElements = true;
																		NecessaryDataRead[9] = true;
																	}
																	else
																		// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
																		if (!NecessaryDataRead[9] && (currentline == "begin centralizer elements"))
																		{
																			AutoCreateCentralizerElements = false;

																			vector<CSpaceGroupElement> Current_Centralizer;
																			const CShiftVector NullShift;
																			const CTwistVector NullTwist;

																			s1 = this->Sectors.size();

																			unsigned counter = 0;
																			for (i = 0; i < s1; ++i)
																				counter += Sectors[i].size();

																			while (GetSaveLine(in, currentline) && (currentline != "end centralizer elements"))
																			{
																				if (currentline == "begin data")
																				{
																					Current_Centralizer.clear();

																					while (GetSaveLine(in, currentline) && (currentline != "end data"))
																					{
																						if (!global_LoadSGElement(currentline, SGElement))
																						{
																							cout << "\n  Warning in bool CSpaceGroup::LoadSpaceGroup(...): Space group element ill specified. Return false." << endl;
																							return false;
																						}
																						Current_Centralizer.push_back(SGElement);
																					}
																					if (Current_Centralizer.size() == 0)
																					{
																						cout << "\n  Warning in bool CSpaceGroup::LoadSpaceGroup(...): Centralizer is empty. Return false." << endl;
																						return false;
																					}

																					this->Centralizer.push_back(Current_Centralizer);
																				}
																			}
																			if (this->Centralizer.size() != counter)
																			{
																				cout << "\n  Warning in bool CSpaceGroup::LoadSpaceGroup(...): Not one centralizer per constructing element loaded. Return false." << endl;
																				return false;
																			}
																			NecessaryDataRead[9] = true;
																		}
		//else
	}
	in.close();

	// begin: error messages
	if ( (ZMxZN && (this->SG_Generators_Twist.size() != 2)) || (!ZMxZN && !ZMxZNxZK && (this->SG_Generators_Twist.size() != 1)) || (ZMxZNxZK && (this->SG_Generators_Twist.size() != 3)) )
	{
		cout << "\n  Warning in bool CSpaceGroup::LoadSpaceGroup(...): Check number of space group generators with twist. Return false." << endl;
		return false;
	}

	if (find(NecessaryDataRead.begin(), NecessaryDataRead.end(), false) != NecessaryDataRead.end())
	{
		vector<string> ErrorMessages(number_of_NecessaryData, "");
		ErrorMessages[0] = "Point group (M, N) is";
		ErrorMessages[1] = "ZMxZN is";
		ErrorMessages[2] = "The twist vector(s) is(are)";
		ErrorMessages[3] = "The six-dim. torus lattice is";
		ErrorMessages[4] = "The generators of the space group with non-trivial twist are";
		ErrorMessages[5] = "The generators of the space group without twist are";
		ErrorMessages[6] = "The allowed orders of the Wilsonlines are";
		ErrorMessages[7] = "The relations between the Wilsonlines are";
		ErrorMessages[8] = "The constructing elements are";
		ErrorMessages[9] = "The centralizer elements are";

		for (i = 0; i < number_of_NecessaryData; ++i)
		{
			if (!NecessaryDataRead[i])
				cout << "\n  Warning in bool CSpaceGroup::LoadSpaceGroup(...): " << ErrorMessages[i] << " not specified. Return false." << endl;
		}
		return false;
	}
	// end: error messages

	if (!reload)
	{
		this->CreateDualTorus();
		this->CreateTwistMatrices();
		this->CreateElementsForConjugacyClass();
	}

	if (AutoFindRelations)
		this->FindRelationsOfWilsonLines();

	if (AutoCreateConstructingElements || AutoCreateCentralizerElements)
	{
		if (this->SG_SetOfElements.size() == 0)
		{
			cout << "\n  create elements..." << flush;
			this->CreateSetOfElements();
			cout << "done." << endl;
		}
	}

	if (AutoCreateConstructingElements)
	{
		cout << "\n  Create constructing elements of the space group automatically." << endl;
		this->CreateConstructingElements();
	}

	if (AutoCreateCentralizerElements)
	{
		cout << "\n  Create centralizer elements automatically." << endl;
		this->CreateCentralizerElements();
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// begin: compute the prefactors for the modular invariance conditions
	// e.g. ShiftsWL_ScalarProductFactors[0][1] = 2  -> 2 (V_0 V_1 - v_0 v_1) = 0 mod 2
	// or   ShiftsWL_ScalarProductFactors[3][2] = 3  -> 3 W_1 V_2 = 0 mod 2
	vector<double> tmpSPF(8, 0.0);
	this->ShiftsWL_ScalarProductFactors.assign(9, tmpSPF);

	vector<unsigned> Orders;
	Orders.push_back(this->M);
	Orders.push_back(this->N);
	Orders.push_back(this->K);
	Orders.insert(Orders.end(), this->WL_AllowedOrders.begin(), this->WL_AllowedOrders.end());

	int ggt = 0;
	for (i = 0; i < 9; ++i)			//hacking here!!!
	{
		for (j = 0; j < 9; ++j)				//hacking here!!!
		{
			if (Orders[i] == Orders[j])
				ggt = Orders[i];
			else
				ggt = ggT(Orders[i], Orders[j]);

			this->ShiftsWL_ScalarProductFactors[i][j] = ggt;
			this->ShiftsWL_ScalarProductFactors[j][i] = ggt;
		}
	}
	// end: compute the prefactors for the modular invariance conditions
	////////////////////////////////////////////////////////////////////////////////////////////////////////

	if (AutoCreateConstructingElements || AutoCreateCentralizerElements || AutoFindRelations)
	{
		std::ofstream out(ifilename.data());

		this->PrintGeometryFile(out);

		if (!this->LoadSpaceGroup(ifilename, true))
		{
			cout << "\n  Warning in bool CSpaceGroup::LoadSpaceGroup(...): Space Group changed but cannot load new version. Return false." << endl;
			return false;
		}
	}

	return true;
}



/* ########################################################################################
######   CreateElementsForConjugacyClass()                                           ######
######                                                                               ######
######   Version: 20.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Creates "SG_AllNonStandardShifts" and "SG_AllRotations".                    ######
######################################################################################## */
bool CSpaceGroup::CreateElementsForConjugacyClass()
{
	unsigned i = 0;
	unsigned j = 0;
	unsigned k = 0;

	size_t s1 = 0;
	size_t s2 = 0;

	this->SG_AllNonStandardShifts.clear();
	this->SG_AllRotations.clear();

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// begin: create SG_AllNonStandardShifts
	vector<CSpaceGroupElement> BasisNonStandardTranslations;
	vector<unsigned> OrderOfNonStandardTranslations;

	bool standard_lattice_vector = true;

	unsigned order = 1;
	s1 = this->SG_Generators_Shift.size();
	for (i = 0; i < s1; ++i)
	{
		const CSpaceGroupElement &SGElement = this->SG_Generators_Shift[i];
		if (!SGElement.NoTwist())
		{
			cout << "\n  Warning in bool CSpaceGroup::CreateElementsForConjugacyClass: pure shift element failed. Return false." << endl;
			return false;
		}

		const CLatticeElement &n_alpha = SGElement.GetLatticeElement();

		order = 1;
		standard_lattice_vector = true;
		for (j = 0; j < LatticeDim; ++j)
		{
			if (n_alpha[j].denominator() != 1)
			{
				order = kgV(order, n_alpha[j].denominator());
				standard_lattice_vector = false;
			}
		}

		// add the non-standard lattice vector
		if (!standard_lattice_vector)
		{
			OrderOfNonStandardTranslations.push_back(order);
			BasisNonStandardTranslations.push_back(SGElement);
		}
	}

	const CSpaceGroupElement SG_Identity;
	CSpaceGroupElement currentElement;
	CSpaceGroupElement tmpElement;

	s2 = BasisNonStandardTranslations.size();
	vector<unsigned> currentCombination(s2, 0);

	if (s2 == 0)
		this->SG_AllNonStandardShifts.push_back(SG_Identity);
	else
	{
		currentCombination.assign(s2, 0);
		do
		{
			currentElement = SG_Identity;
			for (j = 0; j < s2; ++j)
			{
				if (currentCombination[j] != 0)
				{
					const CSpaceGroupElement &currentShift = BasisNonStandardTranslations[j];
					for (k = 0; k < currentCombination[j]; ++k)
					{
						this->SG_Multiply(currentElement, currentShift, tmpElement);
						currentElement = tmpElement;
					}
				}
			}
			this->SG_AllNonStandardShifts.push_back(currentElement);
		}
		while (NextNumber(currentCombination, OrderOfNonStandardTranslations, s2));
	}
	// end: create SG_AllNonStandardShifts
	////////////////////////////////////////////////////////////////////////////////////////////////////////


	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// begin: create SG_AllRotations
	s2 = this->SG_Generators_Twist.size()-1;			//hacking here!!! Witten twist is geometrically identity and does not to be considered here

	vector<unsigned> OrderOfRotations;
	if (this->ZMxZN or this->ZMxZNxZK) {
	OrderOfRotations.push_back(this->N);
	if (this->ZMxZNxZK)
		OrderOfRotations.push_back(this->K);
	}

	if (s2 == 0 or (!ZMxZN && !ZMxZNxZK))
		this->SG_AllRotations.push_back(SG_Identity);
	else
	{
		currentCombination.assign(s2, 0);
		do
		{
			currentElement = SG_Identity;
			for (j = 0; j < s2; ++j)
			{
				if (currentCombination[j] != 0)
				{
					const CSpaceGroupElement &currentRotation = SG_Generators_Twist[j+1];
					for (k = 0; k < currentCombination[j]; ++k)
					{
						this->SG_Multiply(currentRotation, currentElement, tmpElement);
						currentElement = tmpElement;
					}
				}
			}
			this->SG_AllRotations.push_back(currentElement);
		}
		while (NextNumber(currentCombination, OrderOfRotations, s2));
	}
	// end: create SG_AllRotations
	////////////////////////////////////////////////////////////////////////////////////////////////////////

	return true;
}



/* ########################################################################################
######   PrintPointGroup(out)                                           			 ######
######                                                                               ######
######   Version: 20.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) out : ostream object where the output is printed                         ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Prints the point group to "out".                                            ######
######################################################################################## */
void CSpaceGroup::PrintPointGroup(std::ostream &out) const
{
	if (this->ZMxZN)
		out << "Z" << this->M << " x Z" << this->N;
	else if (this->ZMxZNxZK)
		out << "Z" << this->M << " x Z" << this->N << " x Z" << this->K;
	else
		out << "Z" << this->M;

	if (this->additional_label != "")
		out << " - " << this->additional_label;

	out << " Orbifold";
}



/* ########################################################################################
######   CreateSubGroup(CPrint &Print, ...) const                                    ######
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
######   1) SubGroup            : the sub group of this orbifold group as specified  ######
######                            by "TwistFactors1" and "TwistFactors2"             ######
######   return value           : finished successfully?                             ######
###########################################################################################
######   description:                                                                ######
######   If this space group is Z_N: "TwistFactors1" contains one number k and the   ######
######   sub group has a twist k v_1. "TwistFactors2" contains one number l and the  ######
######   sub group has a second twist l v_2.                                         ######
######   If this space group is Z_N x Z_M: "TwistFactors1" contains two numbers      ######
######   (k_1, k_2) and the sub group has a twist k_1 v_1 + k_2 v_2. "TwistFactors2" ######
######   contains two numbers (l_1, l_2) and the sub group has a second twist        ######
######   l_1 v_1 + l_2 v_2.                                                          ######
######################################################################################## */
bool CSpaceGroup::CreateSubGroup(CPrint &Print, vector<unsigned> &TwistFactors1, vector<unsigned> &TwistFactors2, bool use_WL_in_FixedTori, CSpaceGroup &SubGroup) const
{
	// Set the precision
	const double prec = 0.001;

	unsigned i = 0;
	unsigned j = 0;
	unsigned k = 0;

	const bool newZMxZN = (TwistFactors2.size() != 0);

	CTwistVector Twist1, Twist2;
	bool ErrorMessage = false;

	if (this->ZMxZN)
	{
		if (TwistFactors1.size() != 2)
			ErrorMessage = true;
		else
		{
			TwistFactors1[0] = TwistFactors1[0] % this->M;
			TwistFactors1[1] = TwistFactors1[1] % this->N;
			Twist1 = (this->Twists[0] * TwistFactors1[0]) + (this->Twists[1] * TwistFactors1[1]);
		}

		if (newZMxZN)
		{
			if (TwistFactors2.size() != 2)
				ErrorMessage = true;
			else
			{
				TwistFactors2[0] = TwistFactors2[0] % this->M;
				TwistFactors2[1] = TwistFactors2[1] % this->N;
				Twist2 = (this->Twists[0] * TwistFactors2[0]) + (this->Twists[1] * TwistFactors2[1]);
			}
		}
	}
	else
	{
		if (TwistFactors1.size() != 1)
			ErrorMessage = true;
		else
		{
			TwistFactors1[0] = TwistFactors1[0] % this->M;
			Twist1 = this->Twists[0] * TwistFactors1[0];
		}

		if (newZMxZN)
		{
			if (TwistFactors2.size() != 1)
				ErrorMessage = true;
			else
			{
				TwistFactors2[0] = TwistFactors2[0] % this->M;
				Twist2 = this->Twists[0] * TwistFactors2[0];
			}
		}
	}

	// old orbifold was Z_M x Z_N, for example
	// 1) Z_4 x Z_4 with one factors (2,2) and ()
	//    sub-orbifold is Z_2
	// 2) Z_4 x Z_4 with one factors (0,2) and ()
	//    sub-orbifold is Z_2
	// 3) Z_6 x Z_6 with one factors (2,4) and ()
	//    sub-orbifold is Z_3
	// 4) Z_6 x Z_6 with TwistFactors1 = (2,4) and TwistFactors2 = (3,3)
	//    sub-orbifold is Z_3 x Z_2
	const unsigned new_M = Twist1.OrderOfTwist();
	const unsigned new_N = Twist2.OrderOfTwist();

	if ((newZMxZN && ((new_M == 1) || (new_N == 1))) || (!newZMxZN && (new_M == 1) && new_N == 1))
		ErrorMessage = true;

	if (ErrorMessage)
	{
		(*Print.out) << "\n  Warning in bool CSpaceGroup::CreateSubGroup(...) : check ZNxZM. Return false." << endl;
		return false;
	}

	// begin: check that Twist1 and Twist2 are linear independent
	if (newZMxZN)
	{
		bool LinearIndependent = true;

		for (i = 1; LinearIndependent && (i < new_M); ++i)
		{
			CVector i_Twist1 = Twist1 * i;
			for (j = 1; LinearIndependent && (j < new_N); ++j)
			{
				CVector j_Twist2 = Twist2 * j;

				LinearIndependent = false;
				for (k = 1; k < 4; ++k)
				{
					if (!is_integer(i_Twist1[k] - j_Twist2[k]))
						LinearIndependent = true;
				}
				if (!LinearIndependent)
				{
					(*Print.out) << "\n  Warning in bool CSpaceGroup::CreateSubGroup(...) : twist vectors of sub orbifold are linear dependent. Return false." << endl;
					return false;
				}
			}
		}
	}
	// end: check that Twist1 and Twist2 are linear independent

	SubGroup.M = new_M;
	SubGroup.N = new_N;
	SubGroup.ZMxZN = newZMxZN;

	string new_additional_label = "_Sub";

	std::ostringstream os1;
	os1 << TwistFactors1[0];
	new_additional_label += os1.str();

	if (this->ZMxZN)
	{
		std::ostringstream os2;
		os2 << TwistFactors1[1];
		new_additional_label += os2.str();
	}

	if (SubGroup.IsZMxZN())
	{
		new_additional_label += "_";

		std::ostringstream os1;
		os1 << TwistFactors2[0];
		new_additional_label += os1.str();

		if (this->ZMxZN)
		{
			std::ostringstream os2;
			os2 << TwistFactors2[1];
			new_additional_label += os2.str();
		}
	}

	string WLstring = "";
	if (use_WL_in_FixedTori)
		WLstring = "useWL";

	const string NewGeometryFilename = this->GeometryFilename.substr(0,this->GeometryFilename.find(".txt", 0)) + new_additional_label + WLstring + ".txt";

	std::ifstream TestFile;
	TestFile.open(NewGeometryFilename.data(), ifstream::in);
	if(!TestFile.is_open() || !TestFile.good())
	{
		TestFile.close();

		SubGroup.Twists.clear();
		SubGroup.Twists.push_back(Twist1);
		SubGroup.Twists.push_back(Twist2);

		SubGroup.additional_label              = new_additional_label;
		SubGroup.WL_AllowedOrders              = this->WL_AllowedOrders;
		SubGroup.WL_Relations                  = this->WL_Relations;
		SubGroup.lattice_label                 = this->lattice_label;
		SubGroup.ShiftsWL_ScalarProductFactors = this->ShiftsWL_ScalarProductFactors;
		SubGroup.T6_Lattice                    = this->T6_Lattice;
		SubGroup.DualT6_Lattice                = this->DualT6_Lattice;
		SubGroup.SG_Generators_Shift           = this->SG_Generators_Shift;

		SubGroup.SG_Generators_Twist.clear();

		const CSpaceGroupElement SGIdentity;
		CSpaceGroupElement tmpResult;
		CSpaceGroupElement NewSG_Twist(0,0, 0);		//hacking here!!!

		NewSG_Twist = SGIdentity;
		if (TwistFactors1[0] != 0)
		{
			const CSpaceGroupElement &SG_Twist1 = this->SG_Generators_Twist[0];
			for (i = 0; i < TwistFactors1[0]; ++i)
			{
				this->SG_Multiply(SG_Twist1, NewSG_Twist, tmpResult);
				NewSG_Twist = tmpResult;
			}
		}
		if (this->ZMxZN)
		{
			if (TwistFactors1[1] != 0)
			{
				const CSpaceGroupElement &SG_Twist2 = this->SG_Generators_Twist[1];
				for (i = 0; i < TwistFactors1[1]; ++i)
				{
					this->SG_Multiply(SG_Twist2, NewSG_Twist, tmpResult);
					NewSG_Twist = tmpResult;
				}
			}
		}
		NewSG_Twist.Set_m(1);  //hacking here!!! not optimized for ZMxZNxZK
		NewSG_Twist.Set_n(0);
		SubGroup.SG_Generators_Twist.push_back(NewSG_Twist);

		if (newZMxZN)
		{
			NewSG_Twist = SGIdentity;
			if (TwistFactors2[0] != 0)
			{
				const CSpaceGroupElement &SG_Twist1 = this->SG_Generators_Twist[0];
				for (i = 0; i < TwistFactors2[0]; ++i)
				{
					this->SG_Multiply(SG_Twist1, NewSG_Twist, tmpResult);
					NewSG_Twist = tmpResult;
				}
			}
			if (this->ZMxZN)
			{
				if (TwistFactors2[1] != 0)
				{
					const CSpaceGroupElement &SG_Twist2 = this->SG_Generators_Twist[1];
					for (i = 0; i < TwistFactors2[1]; ++i)
					{
						this->SG_Multiply(SG_Twist2, NewSG_Twist, tmpResult);
						NewSG_Twist = tmpResult;
					}
				}
			}
			NewSG_Twist.Set_m(0);
			NewSG_Twist.Set_n(1);
			SubGroup.SG_Generators_Twist.push_back(NewSG_Twist);
		}
		// end: create SG_Generators_Twist

		SubGroup.CreateTwistMatrices();
		SubGroup.CreateElementsForConjugacyClass();

		SubGroup.CreateSetOfElements();
		(*Print.out) << "  " << Print.cbegin << "Create new constructing elements." << Print.cend << endl;
		SubGroup.CreateConstructingElements();
		(*Print.out) << "  " << Print.cbegin << "Create new centralizer elements." << Print.cend << endl;
		SubGroup.CreateCentralizerElements();

		//begin: set order of WL in non-orbifolded directions to 1 and remove them from WL relations
		if (!use_WL_in_FixedTori)
		{
			vector<vector<unsigned> > &New_WL_Relations = SubGroup.WL_Relations;

			bool VectorRotated = false;
			complexVector RotatedTorusVector(ComplexLatticeDim,0);
			const double TwoPi = 2.0 * M_PI;
			size_t s1 = 0;
			vector<unsigned>::iterator it;

			// run through the six torus vectors
			for (i = 0; i < LatticeDim; ++i)
			{
				const complexVector &TorusVector = this->T6_Lattice[i];

				RotatedTorusVector.assign(ComplexLatticeDim,0);
				for (j = 0; j < ComplexLatticeDim; ++j)
					RotatedTorusVector[j] = polar(1.0, TwoPi * Twist1[j+1]) * TorusVector[j];

				VectorRotated = false;
				for (j = 0; j < ComplexLatticeDim; ++j)
				{
					if (norm(TorusVector[j] - RotatedTorusVector[j]) > prec)
						VectorRotated = true;
				}
				if (!VectorRotated && newZMxZN)
				{
					RotatedTorusVector.assign(ComplexLatticeDim,0);
					for (j = 0; j < ComplexLatticeDim; ++j)
						RotatedTorusVector[j] = polar(1.0, TwoPi * Twist2[j+1]) * TorusVector[j];

					for (j = 0; j < ComplexLatticeDim; ++j)
					{
						if (norm(TorusVector[j] - RotatedTorusVector[j]) > prec)
							VectorRotated = true;
					}
				}
				if (!VectorRotated)
				{
					SubGroup.WL_AllowedOrders[i] = 1;

					s1 = New_WL_Relations.size();
					for (j = 0; j < s1; ++j)
					{
						vector<unsigned> &New_WL_Relation = New_WL_Relations[j];

						it = find(New_WL_Relation.begin(), New_WL_Relation.end(), i);
						if (it != New_WL_Relation.end())
							New_WL_Relation.erase(it);
					}
				}
			}
			vector<vector<unsigned> >::iterator it2;
			for (it2 = New_WL_Relations.begin(); it2 < New_WL_Relations.end(); it2++)
			{
				if (it2->size() <= 1)
					New_WL_Relations.erase(it2);
			}
		}
		//end: set order of WL in non-orbifolded directions to 1 and remove them from WL relations

		std::ofstream File_NewSG;
		File_NewSG.open(NewGeometryFilename.data(), ofstream::out);
		SubGroup.PrintGeometryFile(File_NewSG);
		File_NewSG.close();
		(*Print.out) << "  " << Print.cbegin << "New space group saved to file \"" << NewGeometryFilename << "\"." << Print.cend << endl;
	}
	else
	{
		TestFile.close();

		if (!SubGroup.LoadSpaceGroup(NewGeometryFilename))
		{
			(*Print.out) << "\n  Warning in bool CSpaceGroup::CreateSubGroup(...) : geometry-file \"" << NewGeometryFilename << "\" is corrupt. Return false." << endl;
			return false;
		}

		//(*Print.out) << "  " << Print.cbegin << "Space group loaded from file \"" << NewGeometryFilename << "\"." << Print.cend << endl;

		if ((SubGroup.Twists[0] != Twist1) || (newZMxZN && (SubGroup.Twists[1] != Twist2)))
		{
			(*Print.out) << "\n  Warning in bool CSpaceGroup::CreateSubGroup(...) : twists loaded from the geometry-file differ from the ones computed. Return false." << endl;
			return false;
		}
		if ((SubGroup.GetM() != new_M) || (SubGroup.GetN() != new_N))
		{
			(*Print.out) << "\n  Warning in bool CSpaceGroup::CreateSubGroup(...) : ZMxZN loaded from the geometry-file differs from the one computed. Return false." << endl;
			return false;
		}
	}

	SubGroup.SpaceGroup_CheckStatus = NotChecked;
	return true;
}



/* ########################################################################################
######   SetSpaceGroup(const CSpaceGroup &SpaceGroup)                                ######
######                                                                               ######
######   Version: 09.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) SpaceGroup : a CSpaceGroup object                                        ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Sets the space group.                                                       ######
######################################################################################## */
void CSpaceGroup::SetSpaceGroup(const CSpaceGroup &SpaceGroup)
{
	*this = SpaceGroup;
	this->SpaceGroup_CheckStatus = NotChecked;
}



/* ########################################################################################
######   bool CSpaceGroup::PrintGeometryFile(...) const                              ######
######                                                                               ######
######   Version: 16.11.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) out       : ostream object containing the file where the space group     ######
######                  is saved                                                     ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Saves this space group to a file.                                           ######
######################################################################################## */
bool CSpaceGroup::PrintGeometryFile(std::ostream &out) const
{
	unsigned i = 0;
	unsigned j = 0;

	size_t s1 = 0;
	size_t s2 = 0;
	size_t s3 = 0;
	CPrint Print(Tstandard, &out);

	// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	out << "point group\n" << this->M << "\n" << this->N << "\n\n";

	// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	out << "ZMxZN\n";
	if (this->ZMxZN)
		out << "true\n\n";
	else
		out << "false\n\n";

	// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if (this->additional_label != "")
		out << "additional label\n" << this->additional_label << "\n\n";

	// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	out << "begin twist\n";

	Print.PrintRational(this->Twists[0]);
	out << "\n";
	if (this->ZMxZN)
	{
		Print.PrintRational(this->Twists[1]);
		out << "\n";
	}
	out << "end twist\n\n";

	// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if (this->lattice_label != "")
		out << "lattice label\n" << this->lattice_label << "\n\n";

	// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	out << "begin lattice\n";
	s1 = T6_Lattice.size();
	if (s1 != LatticeDim)
	{
		cout << "\n  Warning in bool CSpaceGroup::PrintGeometryFile(...) const: Number of T^6 basis vectors is not 6. Return false." << endl;
		return false;
	}

	for (i = 0; i < s1; ++i)
	{
		const complexVector &vec = this->T6_Lattice[i];
		if (vec.size() != ComplexLatticeDim)
		{
			cout << "\n  Warning in bool CSpaceGroup::PrintGeometryFile(...) const: T^6 basis vector is not " << ComplexLatticeDim << " (complex-) dimensional. Return false." << endl;
			return false;
		}
		else
		{
			out << "next lattice vector\n";
			for (j = 0; j < ComplexLatticeDim; ++j)
				out << setprecision(8) << real(vec[j]) << "\n" << imag(vec[j]) << "\n";
		}
	}
	out << "end lattice\n\n";

	// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	out << "begin twist space group generators\n";
	s1 = this->SG_Generators_Twist.size();
	for (i = 0; i < s1; ++i)
	{
		Print.PrintSGElement(this->SG_Generators_Twist[i], false);
		out << "\n";
	}
	out << "end twist space group generators\n\n";

	// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	out << "begin shift space group generators\n";
	s1 = this->SG_Generators_Shift.size();
	for (i = 0; i < s1; ++i)
	{
		Print.PrintSGElement(this->SG_Generators_Shift[i], false);
		out << "\n";
	}
	out << "end shift space group generators\n\n";

	// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if (this->WL_AllowedOrders.size() != LatticeDim)
	{
		cout << "\n  Warning in bool CSpaceGroup::PrintGeometryFile(...) const: Allowed Orders of Wilson lines not set. Return false." << endl;
		return false;
	}
	else
	{
		out << "allowed order of Wilson lines\n";
		for (i = 0; i < LatticeDim; ++i)
			out << this->WL_AllowedOrders[i] << "\n";
		out << "\n";
	}

	// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	s1 = WL_Relations.size();
	out << "begin identical WLs\n";
	if (s1 == 0)
		out << "none\n";
	else
	{
		for (i = 0; i < s1; ++i)
		{
			const vector<unsigned> &IdentifiedWLs = this->WL_Relations[i];
			s2 = IdentifiedWLs.size();
			for (j = 0; j < s2; ++j)
				out << IdentifiedWLs[j]+1 << " ";
			out << "\n";
		}
	}
	out << "end identical WLs\n\n" << endl;

	// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	s1 = this->DiscreteNonRSymmetries.size();
	s2 = this->DiscreteRSymmetries.size();
	s3 = this->ModularSymmetries.size();
	if ((s1 != 0) || (s2 != 0) || (s3 != 0))
	{
		out << "begin discrete symmetries\n";
		for (i = 0; i < s1; ++i)
		{
			const SDiscreteSymmetry &DiscreteNonRSymmetry = this->DiscreteNonRSymmetries[i];
			out << "new non-R symmetry\n";
			out << DiscreteNonRSymmetry.Order << "\n";
			out << DiscreteNonRSymmetry.SuperpotentialCharge << "\n";
			if (DiscreteNonRSymmetry.ChargeOperator.size() != 8)
			{
				cout << "\n  Warning in bool CSpaceGroup::PrintGeometryFile(...) const: Charge operator of discrete symmetry ill-defined. Return false." << endl;
				return false;
			}
			for (j = 0; j < 8; ++j)
			{
				out << " ";
				Print.PrintRational(DiscreteNonRSymmetry.ChargeOperator[j], false);
			}
			out << "\n";

			out << " " << DiscreteNonRSymmetry.SG_SymmetryGenerator.Get_k() << " " << DiscreteNonRSymmetry.SG_SymmetryGenerator.Get_n();			//hacking here!!! not optimized for ZMxZNxZK
			for (j = 0; j < 6; ++j)
			{
				out << " ";
				Print.PrintRational(DiscreteNonRSymmetry.SG_SymmetryGenerator.Get_n(j), false);
			}
			out << "\n";
		}

		for (i = 0; i < s2; ++i)
		{
			const SDiscreteSymmetry &DiscreteRSymmetry = this->DiscreteRSymmetries[i];
			out << "new R symmetry\n";
			out << DiscreteRSymmetry.Order << "\n";
			out << DiscreteRSymmetry.SuperpotentialCharge << "\n";
			if (DiscreteRSymmetry.ChargeOperator.size() != 4)
			{
				cout << "\n  Warning in bool CSpaceGroup::PrintGeometryFile(...) const: Charge operator of discrete R-symmetry ill-defined. Return false." << endl;
				return false;
			}
			for (j = 0; j < 4; ++j)
			{
				out << " ";
				Print.PrintRational(DiscreteRSymmetry.ChargeOperator[j], false);
			}
			out << "\n";
		}

		for (i = 0; i < s3; ++i)
		{
			const SModularSymmetry &ModularSymmetry = ModularSymmetries[i];
			out << "new modular symmetry\n";
			out << ModularSymmetry.Label << "\n";
			out << ModularSymmetry.Const << "\n";
			out << ModularSymmetry.Index << "\n";
			if (ModularSymmetry.ChargeOperator.size() != 3)
			{
				cout << "\n  Warning in bool CSpaceGroup::PrintGeometryFile(...) const: Charge operator of modular symmetry ill-defined. Return false." << endl;
				return false;
			}
			for (j = 0; j < 3; ++j)
				out << ModularSymmetry.ChargeOperator[j] << " ";
			out << "\n";
		}
		out << "end discrete symmetries\n\n";
	}

	// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	out << "begin constructing elements\n";
	s1 = this->Sectors.size();
	for (i = 0; i < s1; ++i)
	{
		const vector<CSpaceGroupElement> &Sector = Sectors[i];
		s2 = Sector.size();

		out << "new sector" << endl;
		for (j = 0; j < s2; ++j)
		{
			Print.PrintSGElement(Sector[j], false);
			out << "\n";
		}
	}
	out << "end constructing elements\n\n" << endl;

	// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	out << "begin centralizer elements";
	s1 = this->Centralizer.size();
	for (i = 0; i < s1; ++i)
	{
		const vector<CSpaceGroupElement> &Current_Centralizer = this->Centralizer[i];
		s2 = Current_Centralizer.size();

		out << "\nbegin data\n";
		for (j = 0; j < s2; ++j)
		{
			Print.PrintSGElement(Current_Centralizer[j], false);
			out << "\n";
		}
		out << "end data" << endl;
	}
	out << "end centralizer elements\n" << endl;

	return true;
}



/* ########################################################################################
######   ~CSpaceGroup()                                                              ######
######                                                                               ######
######   Version: 09.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard destructor of a CSpaceGroup object.                                ######
######################################################################################## */
CSpaceGroup::~CSpaceGroup()
{
}



/* ########################################################################################
######   SG_Inverse(const CSpaceGroupElement &SGElement) const                       ######
######                                                                               ######
######   Version: 13.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) SGElement : a CSpaceGroupElement object                                  ######
######   output:                                                                     ######
######   return value : the inverse element of "SGElement"                           ######
###########################################################################################
######   description:                                                                ######
######   Computes the inverse space group element of the CSpaceGroupElement object   ######
######   "SGElement".                                                                ######
######################################################################################## */
CSpaceGroupElement CSpaceGroup::SG_Inverse(const CSpaceGroupElement &SGElement) const
{
	CSpaceGroupElement result;
	result.Set_m(SGElement.Get_m());						//hacking here!!! For Witten twist g^(-1)=g since trivial twist by 2 \pi

	if ( (this->ZMxZN) or (this->ZMxZNxZK) ) {
		int tmp = -SGElement.Get_n();
		if (tmp < 0)
			tmp += this->N;
		result.Set_n(tmp);

		if (this->ZMxZNxZK)
		{
			tmp = -SGElement.Get_k();
			if (tmp < 0)
				tmp += this->K;
			result.Set_k(tmp);
		}
	}

	if (SGElement.GetLatticeElement().IsZero())
		return result;

	CLatticeElement RotatedVector;
	SGElement.GetLatticeElement().Rotate(this->TwistMatrices[result.Get_n()][result.Get_k()], RotatedVector);
	RotatedVector *= -1;
	result.SetLatticeElement(RotatedVector);

	return result;
}




/* ########################################################################################
######   bool CSpaceGroup::SG_Commute(...) const                                     ######
######                                                                               ######
######   Version: 05.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) SGElement1  : space group element                                        ######
######   2) SGElement2  : space group element                                        ######
######   output:                                                                     ######
######   return value  : do the elements commute?                                    ######
###########################################################################################
######   description:                                                                ######
######   Cheacks whether the two space group elements "SGElement1" and "SGElement2"  ######
######   commute, i.e. [g_1, g_2] = 0 where g_i denote the two space group elements. ######
######################################################################################## */
bool CSpaceGroup::SG_Commute(const CSpaceGroupElement &SGElement1, const CSpaceGroupElement &SGElement2) const
{
	const rationalMatrix &TwistMatrix1 = this->TwistMatrices[SGElement1.Get_n()][SGElement1.Get_n()]; //not tuned for ZMxZNxZK twists!!!
	const rationalMatrix &TwistMatrix2 = this->TwistMatrices[SGElement2.Get_k()][SGElement2.Get_k()];

	const CLatticeElement &LatticeVector1 = SGElement1.GetLatticeElement();
	const CLatticeElement &LatticeVector2 = SGElement2.GetLatticeElement();

	vector<rational<int> > n_alpha(LatticeDim, 0);

	unsigned j = 0;
	for (unsigned i = 0; i < LatticeDim; ++i)
	{
		rational<int> &n_i = n_alpha[i];
		const rationalVector &TwistVector1 = TwistMatrix1[i];
		const rationalVector &TwistVector2 = TwistMatrix2[i];

		n_i = LatticeVector1[i] - LatticeVector2[i];
		for (j = 0; j < LatticeDim; ++j)
			n_i += LatticeVector2[j] * TwistVector1[j] - LatticeVector1[j] * TwistVector2[j];

		if (n_i != 0)
			return false;
	}
	return true;
}



/* ########################################################################################
######   SG_FixedPoint(const CSpaceGroupElement &SGElement) const                    ######
######                                                                               ######
######   Version: 05.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) SGElement  : space group element                                         ######
######   output:                                                                     ######
######   return value  : is there a fixed point associated to "SGElement"?           ######
###########################################################################################
######   description:                                                                ######
######   Does the space group element "SGElement" have a dixed point or fixed brane? ######
######################################################################################## */
bool CSpaceGroup::SG_FixedPoint(const CSpaceGroupElement &SGElement) const
{
	// create the local twist "LocalTwist" of "SGElement"
	CVector LocalTwist = this->Twists[0] * SGElement.Get_n();
	if (this->ZMxZN)
		LocalTwist += this->Twists[1] * SGElement.Get_k();

	// create the reverse vector of "SGElement"
	complexVector ReverseVector;
	if (!this->SG_ReverseVector(SGElement, ReverseVector))
	{
		cout << "\n  Warning in bool CSpaceGroup::SG_FixedPoint(...) const : reverse vector failed. Return false." << endl;
		return false;
	}

	const double prec = 0.0001;

	// check whether "SGElement" has an associated fixed point
	for (unsigned i = 0; i < ComplexLatticeDim; ++i)
	{
		if (is_integer(LocalTwist[i+1]) && ((fabs(ReverseVector[i].real()) > prec) || (fabs(ReverseVector[i].imag()) > prec)))
			return false;
	}
	return true;
}



/* ########################################################################################
######   SG_FixedPoint(const CSpaceGroupElement &SGElement, ...) const               ######
######                                                                               ######
######   Version: 05.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) SGElement  : space group element                                         ######
######                                                                               ######
######   output:                                                                     ######
######   2) result     : the fixed point associated to the space group element       ######
######                   "SGElement" if it exists                                    ######
######   return value  : is there a fixed point associated to "SGElement"?           ######
###########################################################################################
######   description:                                                                ######
######   Computes the fixed point or fixed brane of "SGElement" and saves the result ######
######   in "result".                                                                ######
######################################################################################## */
bool CSpaceGroup::SG_FixedPoint(const CSpaceGroupElement &SGElement, FPCoordinates &result) const
{
	const double prec = 0.0001;

	result.Coordinates.assign(ComplexLatticeDim, complex<double>(0,0));
	result.FixedTorus.assign(ComplexLatticeDim,false);

	// create the local twist "LocalTwist" of "SGElement"
	CVector LocalTwist = this->Twists[0] * SGElement.Get_m();
	if (this->ZMxZN)
		LocalTwist += this->Twists[1] * SGElement.Get_n();
	if (this->ZMxZNxZK)
		LocalTwist += (this->Twists[1] * SGElement.Get_n() + this->Twists[2] * SGElement.Get_k());

	unsigned i = 0;

	// create the reverse vector of "SGElement"
	complexVector ReverseVector;
	if (!this->SG_ReverseVector(SGElement, ReverseVector))
	{
		cout << "\n  Warning in bool CSpaceGroup::SG_FixedPoint(...) const : reverse vector failed. Return false." << endl;
		return false;
	}

	// check whether "SGElement" has an associated fixed point
	for (i = 0; i < ComplexLatticeDim; ++i)
	{
		if (is_integer(LocalTwist[i+1]) && ((fabs(ReverseVector[i].real()) > prec) || (fabs(ReverseVector[i].imag()) > prec)))
			return false;
	}

	const complex<double> c_one(1,0);
	complex<double> tmp(0,0);

	// (theta^k, n_alpha e_alpha) f_g = f_g  <=>  f_g = (1 - theta^k)^(-1) n_alpha e_alpha
	for (i = 0; i < ComplexLatticeDim; ++i)
	{
		tmp = c_one - polar(1.0, 2.0 * M_PI * LocalTwist[i+1]);
		if ((fabs(real(tmp)) > prec) || (fabs(imag(tmp)) > prec))
			result.Coordinates[i] = (c_one/tmp) * ReverseVector[i];
		else
			result.FixedTorus[i] = true;
	}

	return true;
}




/* ########################################################################################
######   SG_DifferenceInLattice(...) const                                           ######
######                                                                               ######
######   Version: 28.02.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) fp1        : a fixed point coordinate                                    ######
######   2) fp2        : a fixed point coordinate                                    ######
######   output:                                                                     ######
######   return value  : is the difference of fp1 anf fp2 in the lattice spanned     ######
######                   by the six basis vectors?                                   ######
###########################################################################################
######   description:                                                                ######
######   Checks whether the difference of "fp1" and "fp2" in the lattice.            ######
######################################################################################## */
bool CSpaceGroup::SG_DifferenceInLattice(const FPCoordinates &fp1, const FPCoordinates &fp2) const
{
	CVector BasisVector(LatticeDim);
	complex<double> diff;

	unsigned i = 0;
	unsigned j = 0;
	for (i = 0; i < ComplexLatticeDim; ++i)
	{
		if (fp1.FixedTorus[i] || fp2.FixedTorus[i])
		{
			BasisVector[2 * i]     = 0.0;
			BasisVector[2 * i + 1] = 0.0;
		}
		else
		{
			diff = fp1.Coordinates[i] - fp2.Coordinates[i];
			BasisVector[2 * i]     = diff.real();
			BasisVector[2 * i + 1] = diff.imag();
		}
	}

	const double prec = 0.0001;
	bool FixedTorus = false;

	// is difference in the lattice?
	for (i = 0; i < LatticeDim; ++i)
	{
		const CVector &DualLatticeVector = this->DualT6_Lattice[i];

		if (!is_integer(BasisVector * DualLatticeVector))
		{
			FixedTorus = false;
			for (j = 0; j < ComplexLatticeDim; ++j)
			{
				if (fp1.FixedTorus[j] && fp2.FixedTorus[j])
				{
					FixedTorus = true;
					if ((fabs(DualLatticeVector[2 * j]) < prec) && (fabs(DualLatticeVector[2 * j + 1]) < prec))
						return false;
				}
			}
			if (!FixedTorus)
				return false;
		}
	}

	return true;
}



/* ########################################################################################
######   SG_FromSameConjugationClass(...) const                                      ######
######                                                                               ######
######   Version: 23.02.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Element1   : space group element g1                                      ######
######   2) Element2   : space group element g2                                      ######
######   output:                                                                     ######
######   return value  : does h in S exist such that g1 = h g2 h^-1 ?                ######
######                   or: can one rotate the fixed point of g2 such that its      ######
######                       difference to the one of g1 is in the six-lattice?      ######
###########################################################################################
######   description:                                                                ######
######   Checks whether "Element1" and "Element2" belong to the same conjugacy class.######
######################################################################################## */
bool CSpaceGroup::SG_FromSameConjugationClass(const CSpaceGroupElement &Element1, const CSpaceGroupElement &Element2) const
{
	if (this->DualT6_Lattice.size() != LatticeDim)
	{
		cout << "\n  Warning in bool CSpaceGroup::SG_FromSameConjugationClass(...) const : dual T6 lattice not defined. Return false." << endl;
		return false;
	}
	//hacking here!!!  both elements must be from the same Witten sector, i.e. have same w=0,1
	if (Element1.Get_m() != Element2.Get_m()) {
		return false;
	}
	const double prec = 0.0001;
	unsigned i = 0;

	if ( (this->ZMxZN) or (this->ZMxZNxZK) )
	{
		// create the local twist "LocalTwist" of "SGElement"
		CVector LocalTwist1 = this->Twists[1] * Element1.Get_n();
		CVector LocalTwist2 = this->Twists[1] * Element2.Get_n();
		if (this->ZMxZNxZK)
		{
			LocalTwist1 += this->Twists[2] * Element1.Get_k();
			LocalTwist2 += this->Twists[2] * Element2.Get_k();
		}
		// if the twists are different, the elements cannot belong to the same conjugacy class
		for (i = 1; i < 4; ++i)
		{
			if (fabs(LocalTwist1[i] - LocalTwist2[i]) > prec)
				return false;
		}
	}

	const complex<double> cnull(0,0);

	FPCoordinates fp1;
	FPCoordinates fp2;
	if (!this->SG_FixedPoint(Element1, fp1) || !this->SG_FixedPoint(Element2, fp2))
	{
		cout << "\n  Warning in bool CSpaceGroup::SG_FromSameConjugationClass(...) const: space group element has no associated fixed point. Return false." << endl;
		return false;
	}

	FPCoordinates rotated_fp2;
	rotated_fp2.Coordinates.assign(ComplexLatticeDim, 0);
	rotated_fp2.FixedTorus.assign(ComplexLatticeDim, false);

	unsigned j = 0;

	const size_t s1 = this->SG_AllRotations.size();
	const size_t s2 = this->SG_AllNonStandardShifts.size();

	const FPCoordinates original_fp2 = fp2;

	for (i = 0; i < s1; ++i)
	{
		const CSpaceGroupElement &SGElement = this->SG_AllRotations[i];

		for (j = 0; j < s2; ++j)
		{
			const CSpaceGroupElement &SGElementTranslation = this->SG_AllNonStandardShifts[j];

			fp2 = original_fp2;
			if (!this->SG_Multiply(SGElementTranslation, fp2, rotated_fp2))
			{
				cout << "\n  Warning in bool CSpaceGroup::SG_FromSameConjugationClass(...) const: multiplication failed. Return false." << endl;
				return false;
			}

			fp2 = rotated_fp2;
			if (!this->SG_Multiply(SGElement, fp2, rotated_fp2))
			{
				cout << "\n  Warning in bool CSpaceGroup::SG_FromSameConjugationClass(...) const: multiplication failed. Return false." << endl;
				return false;
			}

			if (this->SG_DifferenceInLattice(fp1, rotated_fp2))
				return true;
		}
	}
	return false;
}



/* ########################################################################################
######   SG_Multiply(...) const                                                      ######
######                                                                               ######
######   Version: 04.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) SGElement1 : space group element g1                                      ######
######   2) SGElement2 : space group element g2                                      ######
######   output:                                                                     ######
######   1) result     : the product g1 x g2                                         ######
######   return value  : no error?                                                   ######
###########################################################################################
######   description:                                                                ######
######   Multiplies "SGElement1" and "SGElement2" and saves the result in "result".  ######
######################################################################################## */
bool CSpaceGroup::SG_Multiply(const CSpaceGroupElement &SGElement1, const CSpaceGroupElement &SGElement2, CSpaceGroupElement &result) const
{
  if ((SGElement1.Get_n() >= this->TwistMatrices.size()) || (SGElement1.Get_k() >= this->TwistMatrices[SGElement1.Get_n()].size()))
  {
    cout << "\n  Warning in bool CSpaceGroup::SG_Multiply(...) const (1) : TwistMatrices not defined. Return false." << endl;
    return false;
  }

  result.Set_n((SGElement1.Get_n() + SGElement2.Get_n()) % (unsigned)this->N);
  result.Set_k((SGElement1.Get_k() + SGElement2.Get_k()) % (unsigned)this->K);

  CLatticeElement k_result;
  if (!SGElement2.GetLatticeElement().Rotate(this->TwistMatrices[SGElement1.Get_n()][SGElement1.Get_k()], k_result))
    return false;

  k_result += SGElement1.GetLatticeElement();
  result.SetLatticeElement(k_result);

  return true;
}



/* ########################################################################################
######   SG_Multiply(...) const                                                      ######
######                                                                               ######
######   Version: 13.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) SGElement     : space group element g                                    ######
######   2) LatticeVector : a complex 3-dim. lattice vector e                        ######
######   output:                                                                     ######
######   1) result     : the product g x e                                           ######
######   return value  : no error?                                                   ######
###########################################################################################
######   description:                                                                ######
######   Applies the space group element "SGElement" to the vector "LatticeVector"   ######
######   and saves the result in "result".                                           ######
######################################################################################## */
bool CSpaceGroup::SG_Multiply(const CSpaceGroupElement &SGElement, const complexVector &LatticeVector, complexVector &result) const
{
	if (LatticeVector.size() != ComplexLatticeDim)
	{
		cout << "\n  Warning in bool CSpaceGroup::SG_Multiply(...) const (2) : complexVector has incorrect length. Return false." << endl;
		return false;
	}
	result.assign(ComplexLatticeDim, complex<double>(0,0));
	//hacking here!!! not optimized for ZMxZNxZK
	// create the local twist "LocalTwist" of "SGElement"
	CVector LocalTwist = this->Twists[0] * SGElement.Get_n();
	if (this->ZMxZN)
		LocalTwist += this->Twists[1] * SGElement.Get_k();

	// create the reverse vector of "SGElement"
	complexVector ReverseVector;
	if (!this->SG_ReverseVector(SGElement, ReverseVector))
	{
		cout << "\n  Warning in bool CSpaceGroup::SG_Multiply(...) const (2) : reverse vector failed. Return false." << endl;
		return false;
	}

	// create the rotated lattice vector
	for (unsigned i = 0; i < ComplexLatticeDim; ++i)
		result[i] = (polar(1.0, 2.0 * M_PI * LocalTwist[i+1]) * LatticeVector[i]) + ReverseVector[i];

	return true;
}



/* ########################################################################################
######   SG_Multiply(...) const                                                      ######
######                                                                               ######
######   Version: 13.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) SGElement  : space group element g                                       ######
######   2) FixedPoint : coordinates of a fixed point f                              ######
######   output:                                                                     ######
######   1) result     : the product g x f                                           ######
######   return value  : no error?                                                   ######
###########################################################################################
######   description:                                                                ######
######   Applies the space group element "SGElement" to the point "FixedPoint"       ######
######   and saves the result in "result".                                           ######
######################################################################################## */
bool CSpaceGroup::SG_Multiply(const CSpaceGroupElement &SGElement, const FPCoordinates &FixedPoint, FPCoordinates &result) const
{
	if ((FixedPoint.Coordinates.size() != ComplexLatticeDim) || (FixedPoint.FixedTorus.size() != ComplexLatticeDim))
	{
		cout << "\n  Warning in bool CSpaceGroup::SG_Multiply(...) const (3) : FPCoordinates has incorrect length. Return false." << endl;
		return false;
	}
	result.Coordinates.assign(ComplexLatticeDim, complex<double>(0,0));
	result.FixedTorus = FixedPoint.FixedTorus;

	if ( (this->ZMxZN) or (this->ZMxZNxZK) )
	{
		// create the local twist "LocalTwist" of "SGElement"
		CVector LocalTwist = this->Twists[1] * SGElement.Get_n();
		if (this->ZMxZNxZK)
			LocalTwist += this->Twists[2] * SGElement.Get_k();

		// create the reverse vector of "SGElement"
		complexVector ReverseVector;
		if (!this->SG_ReverseVector(SGElement, ReverseVector))
		{
			cout << "\n  Warning in bool CSpaceGroup::SG_Multiply(...) const (3) : reverse vector failed. Return false." << endl;
			return false;
		}

		// create the rotated fixed point
		for (unsigned i = 0; i < ComplexLatticeDim; ++i)
		{
			if (FixedPoint.FixedTorus[i])
				result.Coordinates[i] = 0;
			else
				result.Coordinates[i] = (polar(1.0, 2.0 * M_PI * LocalTwist[i+1]) * FixedPoint.Coordinates[i]) + ReverseVector[i];
		}
	}
	else {			//if Z2W orbifold
		for (unsigned i = 0; i < ComplexLatticeDim; ++i) {
			result.Coordinates[i] = FixedPoint.Coordinates[i];
		}
	}

	return true;
}



/* ########################################################################################
######   bool CSpaceGroup::SG_ReverseVector(...) const                               ######
######                                                                               ######
######   Version: 13.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) SGElement     : an element of the space group                            ######
######                      (\theta^n\omega^k,n_\alpha e_\alpha)                     ######
######   2) ReverseVector : the vector n_\alpha e_\alpha                             ######
######   output:                                                                     ######
######   return value     : no error?                                                ######
###########################################################################################
######   description:                                                                ######
######   Computes the lattice vector n_\alpha e_\alpha.                              ######
######################################################################################## */
bool CSpaceGroup::SG_ReverseVector(const CSpaceGroupElement &SGElement, complexVector &ReverseVector) const
{
	if (this->T6_Lattice.size() != LatticeDim)
	{
		cout << "\n  Warning in bool CSpaceGroup::SG_ReverseVector(...) const: T" << LatticeDim << " lattice not defined. Return false." << endl;
		return false;
	}
	ReverseVector.assign(ComplexLatticeDim, 0.0);

	unsigned i = 0;
	unsigned j = 0;

	const CLatticeElement &n_alpha = SGElement.GetLatticeElement();

	for (j = 0; j < LatticeDim; ++j)
	{
		const rational<int> &r = n_alpha[j];
		if (r.numerator() != 0)
		{
			const complexVector &LatticeVector = this->T6_Lattice[j];
			for (i = 0; i < ComplexLatticeDim; ++i)
				ReverseVector[i] += ((complex<double>((double)r.numerator()/(double)r.denominator(), 0.0)) * LatticeVector[i]);
		}
	}
	return true;
}

