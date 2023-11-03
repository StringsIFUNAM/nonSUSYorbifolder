
#include "cstate.h"
#include "corbifold.h"
#include "coscillator.h"
#include "cprint.h"
#include "globalfunctions.h"
#include "clatticevector.h"

#define CHECKERROR true

using std::cout;
using std::setw;
using std::setprecision;
using std::endl;
using std::flush;
using std::vector;



/* ########################################################################################
######   CState()                                                                    ######
######                                                                               ######
######   Version: 26.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard constructor of a CField object. No content is specified.           ######
######################################################################################## */
CState::CState()
{
	this->internalIndex.assign(3,0);
	this->State_CheckStatus = CheckedAndFailed;
}



/* ########################################################################################
######   ~CState()                                                                   ######
######                                                                               ######
######   Version: 02.05.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard destructor of a CField object.                                     ######
######################################################################################## */
CState::~CState()
{
}



/* ########################################################################################
######   Create(const CHalfState &LeftMover, const CHalfState &RightMover)           ######
######                                                                               ######
######   Version: 17.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) LeftMover         : a CHalfState object describing the left-moving part  ######
######                          of the string                                        ######
######   2) RightMover        : a CHalfState object describing the right-moving part ######
######                          of the string                                        ######
######   3) Gamma_Centralizer : elements of the orbifold group used to define the    ######
######                          gamma phases                                         ######
######   output:                                                                     ######
######   return value         : finished successfully?                               ######
###########################################################################################
######   description:                                                                ######
######   Tensor left- and right-mover together to form a state and compute the gamma ######
######   phase of the state                                                          ######
######################################################################################## */
bool CState::Create(const CHalfState &LeftMover, const CHalfState &RightMover, const vector<COrbifoldGroupElement> &Gamma_Centralizer)
{
	this->qCharges.clear();
	this->Multiplets.clear();
	this->FieldIndices.clear();
	this->gamma_phases.clear();

	this->LeftMover  = LeftMover;
	this->RightMover = RightMover;

	// begin: gamma
	// The GAMMA phase
	//                      |f>        + exp(-2\pi i\gamma) |\theta f>   + exp(-4\pi i\gamma) |\theta^2 f>
	//               --h--> |\theta f> + exp(-2\pi i\gamma) |\theta^2 f> + exp(-4\pi i\gamma) |f>
	// = exp(2\pi i\gamma) (|f>        + exp(-2\pi i\gamma) |\theta f>   + exp(-4\pi i\gamma) |\theta^2 f>)
	//
	// |q_sh> x |p_sh> --h--> exp(2\pi i (-q_sh v + p_sh V)) |q_sh> x |p_sh>
	//
	// where h corresponds to v,V
	//
	// h - invariant:
	// |q_sh> x |p_sh> x (|f> + exp(-2\pi i\gamma) |\theta f> + exp(-4\pi i\gamma) |\theta^2 f>)
	//
	// with: -q_sh v + p_sh V + \gamma = 0 mod 1
	// thus: \gamma = q_sh v - p_sh V

	const vector<double> &LM_Eigenvalues = this->LeftMover.Eigenvalues;
	const vector<double> &RM_Eigenvalues = this->RightMover.Eigenvalues;

	double tmp = 0.0;

	const size_t s1 = Gamma_Centralizer.size();
	const size_t s2 = LM_Eigenvalues.size()-s1;

	this->gamma_phases.assign(s1, 0.0);

	for (unsigned i = 0; i < s1; ++i)
	{
		const CSpaceGroupElement &GammaElement = Gamma_Centralizer[i].SGElement;

		tmp = (GammaElement.Get_m() * RM_Eigenvalues[0]) + (GammaElement.Get_n() * RM_Eigenvalues[1]) + (GammaElement.Get_k() * RM_Eigenvalues[2]) - LM_Eigenvalues[i + s2];

		RoundDouble(tmp);
		tmp -= floor(tmp);
		this->gamma_phases[i] = tmp;
	}
	// end: gamma

	this->State_CheckStatus = CheckedAndGood;
	return true;
}



/* ########################################################################################
######   ContainsSUSYMultiplet(const vector<CField> &Fields, ...) const              ######
######                                                                               ######
######   Version: 11.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Fields    : vector of CField objects where to look for "Multiplet"       ######
######   2) Multiplet : the SUSY type (e.g. "LeftChiral") to look for                ######
######   output:                                                                     ######
######   return value : does "Fields" contain a field of SUSY type "Multiplet"?      ######
###########################################################################################
######   description:                                                                ######
######   Returns true if the set of fields specified by "Fields" contains at least   ######
######   one field of SUSY type "Multiplet", e.g. at least one is "LeftChiral".      ######
######################################################################################## */
bool CState::ContainsSUSYMultiplet(const vector<CField> &Fields, const SUSYMultiplet &Multiplet) const
{
	unsigned Index = 0;

	const size_t f1 = Fields.size();
	const size_t s1 = this->FieldIndices.size();
	for (unsigned i = 0; i < s1; ++i)
	{
		Index = this->FieldIndices[i];

		if (Index >= f1)
		{
			cout << "Warning in bool CState::ContainsSUSYMultiplet(...) const: \"Index\" out of range. Return false." << endl;
			return false;
		}
		if (Multiplet == Fields[Index].Multiplet)
			return true;
	}
	return false;
}



/* ########################################################################################
######   FindSUSYMultiplets(const CSector &Sector, ...)                              ######
######                                                                               ######
######   Version: 11.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Sectors                    : the sector this state belongs to, needed    ######
######                                   because the right-movers are stored there   ######
######   2) InvariantSupercharges      : a list of all orbifold invariant            ######
######                                   supercharges                                ######
######   3) AllSUSYChargesCombinations : for N=1 SUSY just + and - the supercharge   ######
######                                   for N>1 all +/- combinations of the         ######
######                                   supercharges; this is defined in            ######
######                                   "FindSUSYMultiplets(...)" in the class      ######
######                                   "CFixedBrane"                               ######
######   output:                                                                     ######
######   return value                  : finished successfully?                      ######
###########################################################################################
######   description:                                                                ######
######   Analyzes the right-movers to find out which SUSY multiplets are contained   ######
######   for this state.                                                             ######
######################################################################################## */
bool CState::FindSUSYMultiplets(const CSector &Sector, const vector<CVector> &InvariantSupercharges, const vector<vector<unsigned> > &AllSUSYChargesCombinations)
{
	// Set the precision
	const double prec = 0.0001;

	this->qCharges.clear();
	this->Multiplets.clear();

	const CHalfState       &RM                = this->GetRightMover();
	const vector<CVector>  &RM_Weights        = Sector.GetMasslessRightMover(RM.GetIndex()).Weights;
	const vector<unsigned> &RM_WeightsIndices = RM.Weights;

	const size_t NumberOfSUSY = InvariantSupercharges.size();

	unsigned i = 0;
	unsigned j = 0;
	unsigned k = 0;

	// begin: sort the q charges
	size_t s1 = RM_WeightsIndices.size(); //Run through the solutions of mass^2=0 in a particular chalfstate
	size_t s2 = AllSUSYChargesCombinations.size();
	size_t s3 = 0;

	CVector         tmp(4);
	vector<CVector> tmp_qCharges;
	vector<bool>    qChargesUnsorted(s1, true);
	bool SUSYtypeNotFixed = true;

	for (i = 0; i < s1; ++i)
	{
		if (NumberOfSUSY == 0) {										     //hacking here to include non-susy case  !!!!
			const CVector &qCharge = RM_Weights[RM_WeightsIndices[i]];  	//Save the  q_sh[i] = q + v_g  to a CVector

			if (fabs(qCharge[0]+1.) < prec)							// gauge bosons
			{
				Multiplets.push_back(Gauge);
			}
			else if (fabs(qCharge[0]-1.) < prec)
			{
				Multiplets.push_back(bGauge);
			}
			else if (fabs(qCharge[0]+0.5) < prec)					// fermions
			{
				Multiplets.push_back(LeftFermi);
			}
			else if (fabs(qCharge[0]-0.5) < prec)
			{
				Multiplets.push_back(RightFermi);
			}
			else if (fabs(qCharge[0]) < prec)						// scalars
			{
				int index;
				for (int l=1; l<4; l++)
				{
					if (fabs(qCharge[l]) > prec)
					{
						index = l;
						break;
					}
				}
				if (qCharge[index] > 0)
					Multiplets.push_back(bScalar);
				else
					Multiplets.push_back(Scalar);
			}
			else
			{
				cout << "\n  Warning in bool CState::FindSUSYMultiplets(...): Not all q charges have been identified. Return false." << endl;

				for (int r=0; r<4; r++)
					cout << qCharge[r]<<" ";
				cout << endl;

				return false;
			}
			SUSYtypeNotFixed = false;

			//Non-susy case alternative, slower!!
			/*   	  vector<double> Z2(4,1);									//Construct v_0=(0,1,1,1)
    	  Z2[0]=0;
    	  CVector CZ2;
    	  CZ2.operator=(Z2);
    	  CTwistVector Z2mod;

    	  if (Sector.Get_k() == 1) {									//Determine the sector w.r.t Z_2_mod to assign q to a clatticevector
    		  Z2mod.operator=(CZ2);
    	  }
    	  CLatticeVector q_weight(qCharge-Sector.GetTwist()+Z2mod);		//Change here if you want a chiral spectrum

    	  if (q_weight.From_SO8_Lattice()==SO8_V) {
    	  	  if (q_weight[0] != 0) {
    	  		  Multiplets.push_back(BOSE_mu);					//if q points in the uncompactified light-cone directions \mu
    	  	  }
    	  	  else {												//if q points in the complex compactified directions i
    	  		  int sum=0;
    	  		  for (int p=1; p<4; p++)	 {						//Find CPT conjugate
    	  			  cout <<q_weight[p]<<" ";
    	  			  sum=+q_weight[p];
    	  		  }
    	  		cout << sum<<endl;
    	  		  if (sum<0) {
    	  			  Multiplets.push_back(BOSE_bi);
    	  		  }
    	  		  else {
    	  			Multiplets.push_back(BOSE_i);
    	  		  }
    	  	  }
	  		  SUSYtypeNotFixed = false;
    	  }
    	  else if (q_weight.From_SO8_Lattice()==SO8_S) {
    		  if (fabs(qCharge[0] + 0.5) < prec) {
    			  Multiplets.push_back(LeftFermi_S);
    		  }
    		  else {
    			  Multiplets.push_back(RightFermi_S);
    		  }
			  SUSYtypeNotFixed = false;
    	  }
    	  else if (q_weight.From_SO8_Lattice()==SO8_C) {
    		  if (fabs(qCharge[0] + 0.5) < prec) {
    		     Multiplets.push_back(LeftFermi_C);
    		  }
    		  else {
    			  Multiplets.push_back(RightFermi_C);
    		  }
    		 SUSYtypeNotFixed = false;
       	  } */
			tmp_qCharges.clear();
			tmp_qCharges.push_back(qCharge);
			qChargesUnsorted[i] = false;
			this->qCharges.push_back(tmp_qCharges);
		}
		else {														//SUSY case
			if (qChargesUnsorted[i]) {
				const CVector &qCharge = RM_Weights[RM_WeightsIndices[i]];  	//Save the  q_sh[i] = q + v_g  to a CVector
				tmp_qCharges.clear();
				tmp_qCharges.push_back(qCharge);
				qChargesUnsorted[i] = false;

				// create all q charges that can be obtained from the current one by
				// adding and/or subtracting the supercharges
				for (j = 1; j < s2; ++j)  // (j = 0 is the zero-combination)
				{
					const vector<unsigned> &Combination = AllSUSYChargesCombinations[j];

					tmp = qCharge;
					for (k = 0; k < NumberOfSUSY; ++k)
					{
						if (Combination[k] == 1) // add the supercharge
							tmp += InvariantSupercharges[k];
						else
							if (Combination[k] == 2)  // subtract the supercharge
								tmp -= InvariantSupercharges[k];
					}
					// find the resulting q charge
					for (k = i+1; k < s1; ++k)
					{
						if (qChargesUnsorted[k])
						{
							const CVector &qCharge2 = RM_Weights[RM_WeightsIndices[k]];
							if (tmp == qCharge2)
							{
								tmp_qCharges.push_back(qCharge2);
								qChargesUnsorted[k] = false;
								break;
							}
						}
					}
				}
				this->qCharges.push_back(tmp_qCharges);
			}
			// end: sort the q charges

			if (find(qChargesUnsorted.begin(), qChargesUnsorted.end(), true) != qChargesUnsorted.end())
			{
				cout << "\n  Warning in bool CState::FindSUSYMultiplets(...): Not all q charges have been identified and sorted. Return false." << endl;
				return false;
			}

			s1 = this->qCharges.size();
			if (s1 == 0)
			{
				cout << "\n  Warning in bool CState::FindSUSYMultiplets(...): \"qCharges\" is empty. Return false." << endl;
				return false;
			}

			// begin: identify SUSY multiplet

			if (NumberOfSUSY == 1)
			{
				for (j = 0; j < s1; ++j)
				{
					const vector<CVector> &SetOfqCharges = this->qCharges[j];

					SUSYtypeNotFixed = true;

					s3 = SetOfqCharges.size();
					for (k = 0; SUSYtypeNotFixed && (k < s3); ++k)
					{
						if ((s3 == 2) && (fabs(SetOfqCharges[k][0] - 2.0) < prec))
						{
							Multiplets.push_back(GravityCC);
							SUSYtypeNotFixed = false;
						}
					}
					for (k = 0; SUSYtypeNotFixed && (k < s3); ++k)
					{
						if ((s3 == 2) && (fabs(SetOfqCharges[k][0] + 2.0) < prec))
						{
							Multiplets.push_back(Gravity);
							SUSYtypeNotFixed = false;
						}
					}
					for (k = 0; SUSYtypeNotFixed && (k < s3); ++k)
					{
						if ((s3 == 2) && (fabs(SetOfqCharges[k][0] - 1.0) < prec))
						{
							Multiplets.push_back(VectorCC);
							SUSYtypeNotFixed = false;
						}
					}
					for (k = 0; SUSYtypeNotFixed && (k < s3); ++k)
					{
						if ((s3 == 2) && (fabs(SetOfqCharges[k][0] + 1.0) < prec))
						{
							Multiplets.push_back(Vector);
							SUSYtypeNotFixed = false;
						}
					}
					for (k = 0; SUSYtypeNotFixed && (k < s3); ++k)
					{
						if ((s3 == 2) && (fabs(SetOfqCharges[k][0] - 0.5) < prec))
						{
							Multiplets.push_back(RightChiral);
							SUSYtypeNotFixed = false;
						}
					}
					for (k = 0; SUSYtypeNotFixed && (k < s3); ++k)
					{
						if ((s3 == 2) && (fabs(SetOfqCharges[k][0] + 0.5) < prec))
						{
							Multiplets.push_back(LeftChiral);
							SUSYtypeNotFixed = false;
						}
					}
				}
			}
			else if (NumberOfSUSY == 2)
			{
				for (j = 0; j < s1; ++j)
				{
					const vector<CVector> &SetOfqCharges = this->qCharges[j];

					SUSYtypeNotFixed = true;

					s3 = SetOfqCharges.size();
					for (k = 0; SUSYtypeNotFixed && (k < s3); ++k)
					{
						if ((s3 == 4) && (fabs(SetOfqCharges[k][0] - 2.0) < prec))
						{
							Multiplets.push_back(GravityCC);
							SUSYtypeNotFixed = false;
						}
					}
					for (k = 0; SUSYtypeNotFixed && (k < s3); ++k)
					{
						if ((s3 == 4) && (fabs(SetOfqCharges[k][0] + 2.0) < prec))
						{
							Multiplets.push_back(Gravity);
							SUSYtypeNotFixed = false;
						}
					}
					for (k = 0; SUSYtypeNotFixed && (k < s3); ++k)
					{
						if ((s3 == 4) && (fabs(SetOfqCharges[k][0] - 1.0) < prec))
						{
							Multiplets.push_back(VectorCC);
							SUSYtypeNotFixed = false;
						}
					}
					for (k = 0; SUSYtypeNotFixed && (k < s3); ++k)
					{
						if ((s3 == 4) && (fabs(SetOfqCharges[k][0] + 1.0) < prec))
						{
							Multiplets.push_back(Vector);
							SUSYtypeNotFixed = false;
						}
					}
					for (k = 0; SUSYtypeNotFixed && (k < s3); ++k)
					{
						if ((s3 == 4) && (fabs(SetOfqCharges[k][0] - 0.5) < prec))
						{
							Multiplets.push_back(Halfhyper);
							SUSYtypeNotFixed = false;
						}
					}
				}
			}
			else if (NumberOfSUSY>2)
			{
				cout << "\n  Warning in bool CState::FindSUSYMultiplets(...): Case N = " << NumberOfSUSY << " SUSY not defined. Return false." << endl;
				return false;
			}
		}
		// end: identify SUSY multiplet
	}		  			 //close the i-loop
	return true;
}




/* ########################################################################################
######   FindSUSYMultiplets(const CSector &Sector, ...)                              ######
######                                                                               ######
######   Version: 11.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Sectors                    : the sector this state belongs to, needed    ######
######                                   because the right-movers are stored there   ######
######   2) InvariantSupercharges      : a list of all orbifold invariant            ######
######                                   supercharges                                ######
######   3) AllSUSYChargesCombinations : for N=1 SUSY just + and - the supercharge   ######
######                                   for N>1 all +/- combinations of the         ######
######                                   supercharges; this is defined in            ######
######                                   "FindSUSYMultiplets(...)" in the class      ######
######                                   "CFixedBrane"                               ######
######   output:                                                                     ######
######   return value                  : finished successfully?                      ######
###########################################################################################
######   description:                                                                ######
######   Analyzes the right-movers to find out which SUSY multiplets are contained   ######
######   for this state.                                                             ######
######################################################################################## */
bool CState::TachyonicFindSUSYMultiplets(const CSector &Sector, const vector<CVector> &InvariantSupercharges, const vector<vector<unsigned> > &AllSUSYChargesCombinations)
{
	// Set the precision
	const double prec = 0.0001;

	this->qCharges.clear();
	this->Multiplets.clear();

	const CHalfState       &RM                = this->GetRightMover();
	const vector<CVector>  &RM_Weights        = Sector.GetRTachyons()[0].Weights;			//CTachyonHalfState object containing the tachyonic q, the excited tachyonic state has the same q
	const vector<unsigned> &RM_WeightsIndices = RM.Weights;

	//const size_t NumberOfSUSY = InvariantSupercharges.size();

	unsigned i = 0;
	unsigned j = 0;
	unsigned k = 0;

	// begin: sort the q charges
	size_t s1 = RM_WeightsIndices.size(); 				// Run through the solutions of mass^2<0 in a particular chalfstate
	//size_t s2 = AllSUSYChargesCombinations.size();
	//size_t s3 = 0;

	CVector         tmp(4);
	vector<CVector> tmp_qCharges;
	vector<bool>    qChargesUnsorted(s1, true);
	bool SUSYtypeNotFixed = true;

	for (i = 0; i < s1; ++i)
	{
		const CVector &qCharge = RM_Weights[RM_WeightsIndices[i]];  	//Save the  q_sh[i] = q + v_g  to a CVector
		if (fabs(qCharge[0]) < prec)						// scalars
		{
			int index;
			for (int l=1; l<4; l++)
			{
				if (fabs(qCharge[l]) > prec)
				{
					index = l;
					break;
				}
			}
			if (RM.Excited == 0)									//This switch is not tested yet
			{
				if (qCharge[index] > 0)
					Multiplets.push_back(bTachyon);
				else
					Multiplets.push_back(Tachyon);
			}
			else
			{
				if (qCharge[index] > 0)
					Multiplets.push_back(Excited_bTachyon);
				else
					Multiplets.push_back(Excited_Tachyon);
			}
		}
		else
		{
			cout << "\n  Warning in bool CState::TachyonicFindSUSYMultiplets(...): Cannot identify tachyon type. Return false." << endl;

			for (int r=0; r<4; r++)
				cout << qCharge[r]<<" ";
			cout << endl;

			return false;
		}

		tmp_qCharges.clear();
		tmp_qCharges.push_back(qCharge);
		qChargesUnsorted[i] = false;
		this->qCharges.push_back(tmp_qCharges);
	}

	return true;
}



/* ########################################################################################
######   SetInternalIndex(const unsigned &i, const unsigned &j, ...)                 ######
######                                                                               ######
######   Version: 26.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) i : index of this states' sector in the member variable "Sectors" in the ######
######          class "COrbifold"                                                    ######
######   2) j : index of this states' fixed brane in the member variable             ######
######          "FixedBranes" in the class "Sector"                                  ######
######   3) k : index of this state in the member variable "InvariantStates" in the  ######
######          class "FixedBrane"                                                   ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Set the internal index of this state (specifies the sector, the fixed brane ######
######   and the index of the state).                                                ######
######################################################################################## */
void CState::SetInternalIndex(const unsigned &i, const unsigned &j, const unsigned &k)
{
	this->internalIndex[0] = i;
	this->internalIndex[1] = j;
	this->internalIndex[2] = k;
}



/* ########################################################################################
######   TachyonicCreateRepresentations(...)                                                  ######
######                                                                               ######
######   Version: 29.02.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) FixedBrane    : the CFixedBrane object corresponding to the localization ######
######                      of this state                                            ######
######   2) Vacuum        : the SConfig object "Vacuum" contains the gauge group and ######
######                      the identified fields will be saved here                 ######
######   3) FieldCounters : counters of the field-label index                        ######
######                      (e.g. singlet n with field-counter 17 yields label n_17) ######
######   output:                                                                     ######
######   return value     : finished successfully?                                   ######
###########################################################################################
######   description:                                                                ######
######   Computes the U(1) charges and the non-Abelian representations from the left-######
######   moving shifted momenta p_sh of this state (p_sh are stored in "FixedBrane"),######
######   creates the corresponding CField objects and gives them standard labels     ######
######   using "FieldCounters" to enumerate them and save the result in the vacuum   ######
######   "Vacuum".                                                                   ######
######################################################################################## */
bool CState::TachyonicCreateRepresentations(const CFixedBrane &FixedBrane, SConfig &Vacuum, vector<unsigned> &FieldCounters)
{
	const CGaugeGroup &GaugeGroup = Vacuum.SymmetryGroup.GaugeGroup;
	vector<CField>    &Fields     = Vacuum.Fields;

	// Set the precision
	const double prec = 0.0001;

	const size_t number_of_U1s     = GaugeGroup.u1directions.size();
	const size_t number_of_factors = GaugeGroup.factor.size();

	long double Charge = 0.0;

	size_t s1 = 0;
	size_t s2 = 0;
	unsigned i = 0;
	unsigned j = 0;
	unsigned k = 0;
	unsigned l = 0;

	// get the set of weights for this state
	const vector<unsigned> &unsigned_Weights = this->LeftMover.Weights;
	s1 = unsigned_Weights.size();

	vector<CVector>                  List_U1Charges;
	vector<vector<vector<double> > > List_Weights;
	vector<vector<unsigned> >        List_Indices;
	vector<vector<bool> >            List_NonAbelianRep;

	vector<unsigned>        tmp_Indices;
	vector<vector<double> > tmp_Weights;
	vector<CVector>::iterator pos;

	CVector      U1Charges(number_of_U1s);
	vector<bool> NonAbelianRep(number_of_factors, false);

	const vector<CVector> &all_Weights = FixedBrane.TachyonicGetMassLeftMover(this->LeftMover.GetIndex()).Weights;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// begin: sort all weights according to their "U1Charges" and non-Abelian charges "NonAbelianRep"
	for (i = 0; i < s1; ++i) // run through all weights
	{
		const CVector &weight = all_Weights.at(unsigned_Weights[i]);

		// begin: compute the U(1) charges
		for (j = 0; j < number_of_U1s; ++j)
		{
			const doubleVector &U1Direction = GaugeGroup.u1directions[j];

			// compute the U(1) charge
			Charge = 0.0;
			for (k = 0; k < 16; ++k)
				Charge += U1Direction[k] * weight[k];

			RoundCharge(Charge);
			U1Charges[j] = (double)Charge;
		}
		// end: compute the U(1) charges

		// begin: is the weight charged under the j-th gauge group factor
		NonAbelianRep.assign(number_of_factors, false);
		for (j = 0; j < number_of_factors; ++j)
		{
			const vector<vector<double> > &SimpleRoots = GaugeGroup.factor[j].simpleroots;
			s2 = SimpleRoots.size();
			for (k = 0; !NonAbelianRep[j] && (k < s2); ++k)
			{
				const vector<double> &SimpleRoot = SimpleRoots[k];

				Charge = 0.0;
				for (l = 0; l < 16; ++l)
					Charge += SimpleRoot[l] * weight[l];

				if (fabs(Charge) > prec)
					NonAbelianRep[j] = true;
			}
		}
		// end: is the weight charged under the j-th gauge group factor

		// begin: sort
		bool unknown = true;
		s2 = List_U1Charges.size();
		for (j = 0; unknown && (j < s2); ++j)
		{
			if ((List_U1Charges[j] == U1Charges) && (List_NonAbelianRep[j] == NonAbelianRep))
			{
				List_Weights[j].push_back(weight);
				List_Indices[j].push_back(unsigned_Weights[i]);
				unknown = false;
			}
		}
		if (unknown)
		{
			List_U1Charges.push_back(U1Charges);
			List_NonAbelianRep.push_back(NonAbelianRep);

			tmp_Weights.clear();
			tmp_Weights.push_back(weight);
			List_Weights.push_back(tmp_Weights);

			tmp_Indices.clear();
			tmp_Indices.push_back(unsigned_Weights[i]);
			List_Indices.push_back(tmp_Indices);
		}
		// end: sort
	}
	// end: sort all weights according to their "U1Charges" and non-Abelian charges "NonAbelianRep"
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////


	cirrep<double> Irrep;
	vector<int> HW;
	SDimension DimOfRep;

	if (FieldCounters.size() != 15)			//hacking here!!!
	{
		cout << "Warning in bool CState::CreateRepresentations(...): \"FieldCounters\" ill-defined. Return false." << endl;
		return false;
	}
	size_t number_of_SUSYMultiplets = this->Multiplets.size();

	CVector Vector_m;
	CVector Vector_p;
	Vector_m.Assign(4,-0.5);
	Vector_p.Assign(4, 0.5);

	// begin: create new fields
	s1 = List_Weights.size();
	for (i = 0; i < s1; ++i)
	{
		const vector<vector<double> > &Weights = List_Weights[i];

		Irrep = findHighestWeight<double>(GaugeGroup, Weights);

		CField NewField;
		NewField.VEVs.Assign(Weights.size(), 0.0);

		// "NewField.internalIndex" is set in bool COrbifold::CreateRepresentations(...)

		NewField.AccU1Charges.clear();
		NewField.BmLCharge        = 0.0;
		NewField.WeightIndices    = List_Indices[i];
		NewField.U1Charges        = List_U1Charges[i];
		NewField.SGElement        = FixedBrane.GetSGElement();
		NewField.gamma_phases     = this->gamma_phases;
		NewField.OsciContribution = this->OsciContribution;

		for (j = 0; j < number_of_factors; ++j)
		{
			const gaugeGroupFactor<double> &ggf = GaugeGroup.factor[j];
			HW = findDynkinLabels(ggf, Irrep.highestweight);
			if (!this->DetermineDimension(ggf, HW, DimOfRep))
				return false;

			NewField.Dimensions.push_back(DimOfRep);
			NewField.HighestWeights_DL.push_back(HW);
		}


		for (j = 0; j < number_of_SUSYMultiplets; ++j)
		{
			this->FieldIndices.push_back(Fields.size());

			NewField.Multiplet = this->Multiplets[j];
			//hacking here!!!
			// begin: untwisted modulus
			/*if (NewField.SGElement.IsZero() && IsSinglet(NewField.Dimensions) && NewField.U1Charges.IsZero())
			{
				if (NewField.Multiplet == LeftChiral) {
					NewField.Multiplet = LCModulus;
				}
				else if (NewField.Multiplet == RightChiral) {
					NewField.Multiplet = RCModulus;
				}
				else if (NewField.Multiplet == Gauge) {						//hacking here for non-susy gravity and...
					NewField.Multiplet = Gravity;
					NewField.Labels[0]  = "G";
					NewField.Numbers[0] = FieldCounters[1];
					++FieldCounters[1];
				}
				else if (NewField.Multiplet == bGauge) {						//hacking here for non-susy gravity and...
					NewField.Multiplet = GravityCC;
					NewField.Labels[0]  = "bG";
					NewField.Numbers[0] = FieldCounters[2];
					++FieldCounters[2];
				}
				else if (NewField.Multiplet == Scalar) {						//...untwisted moduli
					NewField.Multiplet = moduli;
					NewField.Labels[0]  = "m";
					NewField.Numbers[0] = FieldCounters[3];
					++FieldCounters[3];
				}
				else if (NewField.Multiplet == bScalar) {						//untwisted moduli bar
					NewField.Multiplet = bmoduli;
					NewField.Labels[0]  = "bm";
					NewField.Numbers[0] = FieldCounters[4];
					++FieldCounters[4];
				}
			}*/
			// end: untwisted modulus

			//     if (NumberOfSUSY==0) {						//enable check for susy - non-susy

			switch (NewField.Multiplet)
			{
			/*case Gauge:									//N=0 states   hacking here!!!
			{
				NewField.Labels[0]  = "g";							//vector bosons
				NewField.Numbers[0] = FieldCounters[5];
				++FieldCounters[5];
				break;
			}
			case bGauge:
			{
				NewField.Labels[0]  = "bg";				//bar vector bosons
				NewField.Numbers[0] = FieldCounters[6];
				++FieldCounters[6];
				break;
			}*/
			case Scalar:
			{
				NewField.Labels[0]  = "s";				//scalar since index i points in the internal space osc.No.=0
				NewField.Numbers[0] = FieldCounters[7];
				++FieldCounters[7];
				break;
			}
			case bScalar:
			{
				NewField.Labels[0]  = "bs";				//bar scalar
				NewField.Numbers[0] = FieldCounters[8];
				++FieldCounters[8];
				break;
			}
			/*case LeftFermi:
			{
				NewField.Labels[0]  = "f";
				NewField.Numbers[0] = FieldCounters[11];
				++FieldCounters[11];
				break;
			}
			case RightFermi:
			{
				NewField.Labels[0]  = "bf";
				NewField.Numbers[0] = FieldCounters[12];
				++FieldCounters[12];
				break;
			}*/
			}
			if (NewField.Multiplet == NOT_DEF_SUSY) {
				cout << "Warning in bool CState::TachyonicCreateRepresentations(...): SUSY type not defined. Return false.!" << endl;
				return false;
			}
			Fields.push_back(NewField);
		}
	}
	// end: create new fields

	return true;
}



/* ########################################################################################
######   ExcitedTachyonicCreateRepresentations(...)                                                  ######
######                                                                               ######
######   Version: 29.02.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) FixedBrane    : the CFixedBrane object corresponding to the localization ######
######                      of this state                                            ######
######   2) Vacuum        : the SConfig object "Vacuum" contains the gauge group and ######
######                      the identified fields will be saved here                 ######
######   3) FieldCounters : counters of the field-label index                        ######
######                      (e.g. singlet n with field-counter 17 yields label n_17) ######
######   output:                                                                     ######
######   return value     : finished successfully?                                   ######
###########################################################################################
######   description:                                                                ######
######   Computes the U(1) charges and the non-Abelian representations from the left-######
######   moving shifted momenta p_sh of this state (p_sh are stored in "FixedBrane"),######
######   creates the corresponding CField objects and gives them standard labels     ######
######   using "FieldCounters" to enumerate them and save the result in the vacuum   ######
######   "Vacuum".                                                                   ######
######################################################################################## */
bool CState::ExcitedTachyonicCreateRepresentations(const CFixedBrane &FixedBrane, SConfig &Vacuum, vector<unsigned> &FieldCounters)
{
	const CGaugeGroup &GaugeGroup = Vacuum.SymmetryGroup.GaugeGroup;
	vector<CField>    &Fields     = Vacuum.Fields;

	// Set the precision
	const double prec = 0.0001;

	const size_t number_of_U1s     = GaugeGroup.u1directions.size();
	const size_t number_of_factors = GaugeGroup.factor.size();

	long double Charge = 0.0;

	size_t s1 = 0;
	size_t s2 = 0;
	unsigned i = 0;
	unsigned j = 0;
	unsigned k = 0;
	unsigned l = 0;

	// get the set of weights for this state
	const vector<unsigned> &unsigned_Weights = this->LeftMover.Weights;
	s1 = unsigned_Weights.size();

	vector<CVector>                  List_U1Charges;
	vector<vector<vector<double> > > List_Weights;
	vector<vector<unsigned> >        List_Indices;
	vector<vector<bool> >            List_NonAbelianRep;

	vector<unsigned>        tmp_Indices;
	vector<vector<double> > tmp_Weights;
	vector<CVector>::iterator pos;

	CVector      U1Charges(number_of_U1s);
	vector<bool> NonAbelianRep(number_of_factors, false);

	const vector<CVector> &all_Weights = FixedBrane.ExcitedTachyonicGetMassLeftMover(this->LeftMover.GetIndex()).Weights;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// begin: sort all weights according to their "U1Charges" and non-Abelian charges "NonAbelianRep"
	for (i = 0; i < s1; ++i) // run through all weights
	{
		const CVector &weight = all_Weights.at(unsigned_Weights[i]);

		// begin: compute the U(1) charges
		for (j = 0; j < number_of_U1s; ++j)
		{
			const doubleVector &U1Direction = GaugeGroup.u1directions[j];

			// compute the U(1) charge
			Charge = 0.0;
			for (k = 0; k < 16; ++k)
				Charge += U1Direction[k] * weight[k];

			RoundCharge(Charge);
			U1Charges[j] = (double)Charge;
		}
		// end: compute the U(1) charges

		// begin: is the weight charged under the j-th gauge group factor
		NonAbelianRep.assign(number_of_factors, false);
		for (j = 0; j < number_of_factors; ++j)
		{
			const vector<vector<double> > &SimpleRoots = GaugeGroup.factor[j].simpleroots;
			s2 = SimpleRoots.size();
			for (k = 0; !NonAbelianRep[j] && (k < s2); ++k)
			{
				const vector<double> &SimpleRoot = SimpleRoots[k];

				Charge = 0.0;
				for (l = 0; l < 16; ++l)
					Charge += SimpleRoot[l] * weight[l];

				if (fabs(Charge) > prec)
					NonAbelianRep[j] = true;
			}
		}
		// end: is the weight charged under the j-th gauge group factor

		// begin: sort
		bool unknown = true;
		s2 = List_U1Charges.size();
		for (j = 0; unknown && (j < s2); ++j)
		{
			if ((List_U1Charges[j] == U1Charges) && (List_NonAbelianRep[j] == NonAbelianRep))
			{
				List_Weights[j].push_back(weight);
				List_Indices[j].push_back(unsigned_Weights[i]);
				unknown = false;
			}
		}
		if (unknown)
		{
			List_U1Charges.push_back(U1Charges);
			List_NonAbelianRep.push_back(NonAbelianRep);

			tmp_Weights.clear();
			tmp_Weights.push_back(weight);
			List_Weights.push_back(tmp_Weights);

			tmp_Indices.clear();
			tmp_Indices.push_back(unsigned_Weights[i]);
			List_Indices.push_back(tmp_Indices);
		}
		// end: sort
	}
	// end: sort all weights according to their "U1Charges" and non-Abelian charges "NonAbelianRep"
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////


	cirrep<double> Irrep;
	vector<int> HW;
	SDimension DimOfRep;

	if (FieldCounters.size() != 15)			//hacking here!!!
	{
		cout << "Warning in bool CState::CreateRepresentations(...): \"FieldCounters\" ill-defined. Return false." << endl;
		return false;
	}
	size_t number_of_SUSYMultiplets = this->Multiplets.size();

	CVector Vector_m;
	CVector Vector_p;
	Vector_m.Assign(4,-0.5);
	Vector_p.Assign(4, 0.5);

	// begin: create new fields
	s1 = List_Weights.size();
	for (i = 0; i < s1; ++i)
	{
		const vector<vector<double> > &Weights = List_Weights[i];

		Irrep = findHighestWeight<double>(GaugeGroup, Weights);

		CField NewField;
		NewField.VEVs.Assign(Weights.size(), 0.0);

		// "NewField.internalIndex" is set in bool COrbifold::CreateRepresentations(...)

		NewField.AccU1Charges.clear();
		NewField.BmLCharge        = 0.0;
		NewField.WeightIndices    = List_Indices[i];
		NewField.U1Charges        = List_U1Charges[i];
		NewField.SGElement        = FixedBrane.GetSGElement();
		NewField.gamma_phases     = this->gamma_phases;
		NewField.OsciContribution = this->OsciContribution;

		for (j = 0; j < number_of_factors; ++j)
		{
			const gaugeGroupFactor<double> &ggf = GaugeGroup.factor[j];
			HW = findDynkinLabels(ggf, Irrep.highestweight);
			if (!this->DetermineDimension(ggf, HW, DimOfRep))
				return false;

			NewField.Dimensions.push_back(DimOfRep);
			NewField.HighestWeights_DL.push_back(HW);
		}


		for (j = 0; j < number_of_SUSYMultiplets; ++j)
		{
			this->FieldIndices.push_back(Fields.size());

			NewField.Multiplet = this->Multiplets[j];
			//hacking here!!!
			// begin: untwisted modulus
			/*if (NewField.SGElement.IsZero() && IsSinglet(NewField.Dimensions) && NewField.U1Charges.IsZero())
			{
				if (NewField.Multiplet == LeftChiral) {
					NewField.Multiplet = LCModulus;
				}
				else if (NewField.Multiplet == RightChiral) {
					NewField.Multiplet = RCModulus;
				}
				else if (NewField.Multiplet == Gauge) {						//hacking here for non-susy gravity and...
					NewField.Multiplet = Gravity;
					NewField.Labels[0]  = "G";
					NewField.Numbers[0] = FieldCounters[1];
					++FieldCounters[1];
				}
				else if (NewField.Multiplet == bGauge) {						//hacking here for non-susy gravity and...
					NewField.Multiplet = GravityCC;
					NewField.Labels[0]  = "bG";
					NewField.Numbers[0] = FieldCounters[2];
					++FieldCounters[2];
				}
				else if (NewField.Multiplet == Scalar) {						//...untwisted moduli
					NewField.Multiplet = moduli;
					NewField.Labels[0]  = "m";
					NewField.Numbers[0] = FieldCounters[3];
					++FieldCounters[3];
				}
				else if (NewField.Multiplet == bScalar) {						//untwisted moduli bar
					NewField.Multiplet = bmoduli;
					NewField.Labels[0]  = "bm";
					NewField.Numbers[0] = FieldCounters[4];
					++FieldCounters[4];
				}
			}*/
			// end: untwisted modulus

			//     if (NumberOfSUSY==0) {						//enable check for susy - non-susy

			switch (NewField.Multiplet)
			{
			/*case Gauge:									//N=0 states   hacking here!!!
			{
				NewField.Labels[0]  = "g";							//vector bosons
				NewField.Numbers[0] = FieldCounters[5];
				++FieldCounters[5];
				break;
			}
			case bGauge:
			{
				NewField.Labels[0]  = "bg";				//bar vector bosons
				NewField.Numbers[0] = FieldCounters[6];
				++FieldCounters[6];
				break;
			}*/
			case Scalar:
			{
				NewField.Labels[0]  = "s";				//scalar since index i points in the internal space osc.No.=0
				NewField.Numbers[0] = FieldCounters[7];
				++FieldCounters[7];
				break;
			}
			case bScalar:
			{
				NewField.Labels[0]  = "bs";				//bar scalar
				NewField.Numbers[0] = FieldCounters[8];
				++FieldCounters[8];
				break;
			}
			/*case rScalar:
			{
				NewField.Labels[0]  = "r";				//r scalar
				NewField.Numbers[0] = FieldCounters[9];
				++FieldCounters[9];
				break;
			}
			case brScalar:
			{
				NewField.Labels[0]  = "br";				//bar r scalar
				NewField.Numbers[0] = FieldCounters[10];
				++FieldCounters[10];
				break;
			}
			case LeftFermi_S:
			{
				NewField.Labels[0]  = "f_S";
				NewField.Numbers[0] = FieldCounters[11];
				++FieldCounters[11];
				break;
			}
			case RightFermi_S:
			{
				NewField.Labels[0]  = "bf_S";
				NewField.Numbers[0] = FieldCounters[12];
				++FieldCounters[12];
				break;
			}
			case LeftFermi_C:
			{
				NewField.Labels[0]  = "f_C";
				NewField.Numbers[0] = FieldCounters[13];
				++FieldCounters[13];
				break;
			}
			case RightFermi_C:
			{
				NewField.Labels[0]  = "bf_C";
				NewField.Numbers[0] = FieldCounters[14];
				++FieldCounters[14];
				break;
			}*/
			}
			if (NewField.Multiplet == NOT_DEF_SUSY) {
				cout << "Warning in bool CState::CreateRepresentations(...): SUSY type not defined. Return false.!" << endl;
				return false;
			}
			Fields.push_back(NewField);
		}
	}
	// end: create new fields

	return true;
}



/* ########################################################################################
######   CreateRepresentations(...)                                                  ######
######                                                                               ######
######   Version: 29.02.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) FixedBrane    : the CFixedBrane object corresponding to the localization ######
######                      of this state                                            ######
######   2) Vacuum        : the SConfig object "Vacuum" contains the gauge group and ######
######                      the identified fields will be saved here                 ######
######   3) FieldCounters : counters of the field-label index                        ######
######                      (e.g. singlet n with field-counter 17 yields label n_17) ######
######   output:                                                                     ######
######   return value     : finished successfully?                                   ######
###########################################################################################
######   description:                                                                ######
######   Computes the U(1) charges and the non-Abelian representations from the left-######
######   moving shifted momenta p_sh of this state (p_sh are stored in "FixedBrane"),######
######   creates the corresponding CField objects and gives them standard labels     ######
######   using "FieldCounters" to enumerate them and save the result in the vacuum   ######
######   "Vacuum".                                                                   ######
######################################################################################## */
bool CState::CreateRepresentations(const CFixedBrane &FixedBrane, SConfig &Vacuum, vector<unsigned> &FieldCounters)
{
	const CGaugeGroup &GaugeGroup = Vacuum.SymmetryGroup.GaugeGroup;
	vector<CField>    &Fields     = Vacuum.Fields;

	// Set the precision
	const double prec = 0.0001;

	const size_t number_of_U1s     = GaugeGroup.u1directions.size();
	const size_t number_of_factors = GaugeGroup.factor.size();

	long double Charge = 0.0;

	size_t s1 = 0;
	size_t s2 = 0;
	unsigned i = 0;
	unsigned j = 0;
	unsigned k = 0;
	unsigned l = 0;

	// get the set of weights for this state
	const vector<unsigned> &unsigned_Weights = this->LeftMover.Weights;
	s1 = unsigned_Weights.size();

	vector<CVector>                  List_U1Charges;
	vector<vector<vector<double> > > List_Weights;
	vector<vector<unsigned> >        List_Indices;
	vector<vector<bool> >            List_NonAbelianRep;

	vector<unsigned>        tmp_Indices;
	vector<vector<double> > tmp_Weights;
	vector<CVector>::iterator pos;

	CVector      U1Charges(number_of_U1s);
	vector<bool> NonAbelianRep(number_of_factors, false);

	const vector<CVector> &all_Weights = FixedBrane.GetMasslessLeftMover(this->LeftMover.GetIndex()).Weights;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// begin: sort all weights according to their "U1Charges" and non-Abelian charges "NonAbelianRep"
	for (i = 0; i < s1; ++i) // run through all weights
	{
		const CVector &weight = all_Weights.at(unsigned_Weights[i]);

		// begin: compute the U(1) charges
		for (j = 0; j < number_of_U1s; ++j)
		{
			const doubleVector &U1Direction = GaugeGroup.u1directions[j];

			// compute the U(1) charge
			Charge = 0.0;
			for (k = 0; k < 16; ++k)
				Charge += U1Direction[k] * weight[k];

			RoundCharge(Charge);
			U1Charges[j] = (double)Charge;
		}
		// end: compute the U(1) charges

		// begin: is the weight charged under the j-th gauge group factor
		NonAbelianRep.assign(number_of_factors, false);
		for (j = 0; j < number_of_factors; ++j)
		{
			const vector<vector<double> > &SimpleRoots = GaugeGroup.factor[j].simpleroots;
			s2 = SimpleRoots.size();
			for (k = 0; !NonAbelianRep[j] && (k < s2); ++k)
			{
				const vector<double> &SimpleRoot = SimpleRoots[k];

				Charge = 0.0;
				for (l = 0; l < 16; ++l)
					Charge += SimpleRoot[l] * weight[l];

				if (fabs(Charge) > prec)
					NonAbelianRep[j] = true;
			}
		}
		// end: is the weight charged under the j-th gauge group factor

		// begin: sort
		bool unknown = true;
		s2 = List_U1Charges.size();
		for (j = 0; unknown && (j < s2); ++j)
		{
			if ((List_U1Charges[j] == U1Charges) && (List_NonAbelianRep[j] == NonAbelianRep))
			{
				List_Weights[j].push_back(weight);
				List_Indices[j].push_back(unsigned_Weights[i]);
				unknown = false;
			}
		}
		if (unknown)
		{
			List_U1Charges.push_back(U1Charges);
			List_NonAbelianRep.push_back(NonAbelianRep);

			tmp_Weights.clear();
			tmp_Weights.push_back(weight);
			List_Weights.push_back(tmp_Weights);

			tmp_Indices.clear();
			tmp_Indices.push_back(unsigned_Weights[i]);
			List_Indices.push_back(tmp_Indices);
		}
		// end: sort
	}
	// end: sort all weights according to their "U1Charges" and non-Abelian charges "NonAbelianRep"
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////


	cirrep<double> Irrep;
	vector<int> HW;
	SDimension DimOfRep;

	if (FieldCounters.size() != 15)			//hacking here!!!
	{
		cout << "Warning in bool CState::CreateRepresentations(...): \"FieldCounters\" ill-defined. Return false." << endl;
		return false;
	}
	size_t number_of_SUSYMultiplets = this->Multiplets.size();

	CVector Vector_m;
	CVector Vector_p;
	Vector_m.Assign(4,-0.5);
	Vector_p.Assign(4, 0.5);

	// begin: create new fields
	s1 = List_Weights.size();
	for (i = 0; i < s1; ++i)
	{
		const vector<vector<double> > &Weights = List_Weights[i];

		Irrep = findHighestWeight<double>(GaugeGroup, Weights);

		CField NewField;
		NewField.VEVs.Assign(Weights.size(), 0.0);

		// "NewField.internalIndex" is set in bool COrbifold::CreateRepresentations(...)

		NewField.AccU1Charges.clear();
		NewField.BmLCharge        = 0.0;
		NewField.WeightIndices    = List_Indices[i];
		NewField.U1Charges        = List_U1Charges[i];
		NewField.SGElement        = FixedBrane.GetSGElement();
		NewField.gamma_phases     = this->gamma_phases;
		NewField.OsciContribution = this->OsciContribution;

		for (j = 0; j < number_of_factors; ++j)
		{
			const gaugeGroupFactor<double> &ggf = GaugeGroup.factor[j];
			HW = findDynkinLabels(ggf, Irrep.highestweight);
			if (!this->DetermineDimension(ggf, HW, DimOfRep))
				return false;

			NewField.Dimensions.push_back(DimOfRep);
			NewField.HighestWeights_DL.push_back(HW);
		}


		for (j = 0; j < number_of_SUSYMultiplets; ++j)
		{
			this->FieldIndices.push_back(Fields.size());

			NewField.Multiplet = this->Multiplets[j];
			//hacking here!!!
			// begin: untwisted modulus
			if (NewField.SGElement.IsZero() && IsSinglet(NewField.Dimensions) && NewField.U1Charges.IsZero())
			{
				if (NewField.Multiplet == LeftChiral) {
					NewField.Multiplet = LCModulus;
				}
				else if (NewField.Multiplet == RightChiral) {
					NewField.Multiplet = RCModulus;
				}
				else if (NewField.Multiplet == Gauge) {						//hacking here for non-susy gravity and...
					NewField.Multiplet = Gravity;
					NewField.Labels[0]  = "G";
					NewField.Numbers[0] = FieldCounters[1];
					++FieldCounters[1];
				}
				else if (NewField.Multiplet == bGauge) {						//hacking here for non-susy gravity and...
					NewField.Multiplet = GravityCC;
					NewField.Labels[0]  = "bG";
					NewField.Numbers[0] = FieldCounters[2];
					++FieldCounters[2];
				}
				else if (NewField.Multiplet == Scalar) {						//...untwisted moduli
					NewField.Multiplet = moduli;
					NewField.Labels[0]  = "m";
					NewField.Numbers[0] = FieldCounters[3];
					++FieldCounters[3];
				}
				else if (NewField.Multiplet == bScalar) {						//untwisted moduli bar
					NewField.Multiplet = bmoduli;
					NewField.Labels[0]  = "bm";
					NewField.Numbers[0] = FieldCounters[4];
					++FieldCounters[4];
				}
			}
			// end: untwisted modulus

			//     if (NumberOfSUSY==0) {						//enable check for susy - non-susy

			switch (NewField.Multiplet)
			{
			case Gauge:									//N=0 states   hacking here!!!
			{
				NewField.Labels[0]  = "g";							//vector bosons
				NewField.Numbers[0] = FieldCounters[5];
				++FieldCounters[5];
				break;
			}
			case bGauge:
			{
				NewField.Labels[0]  = "bg";				//bar vector bosons
				NewField.Numbers[0] = FieldCounters[6];
				++FieldCounters[6];
				break;
			}
			case Scalar:
			{
				NewField.Labels[0]  = "s";				//scalar since index i points in the internal space osc.No.=0
				NewField.Numbers[0] = FieldCounters[7];
				++FieldCounters[7];
				break;
			}
			case bScalar:
			{
				NewField.Labels[0]  = "bs";				//bar scalar
				NewField.Numbers[0] = FieldCounters[8];
				++FieldCounters[8];
				break;
			}
			case LeftFermi:
			{
				NewField.Labels[0]  = "f";
				NewField.Numbers[0] = FieldCounters[11];
				++FieldCounters[11];
				break;
			}
			case RightFermi:
			{
				NewField.Labels[0]  = "bf";
				NewField.Numbers[0] = FieldCounters[12];
				++FieldCounters[12];
				break;
			}
			}
			//      }
			/*      else { 											//susy states, disabled, hacking here!!!

      const vector<CVector> &SetOfqCharges = this->qCharges[j];
      s2 =  SetOfqCharges.size();

      if (Vacuum.InvariantSupercharges.size() == 1)						//susy relevant
      {
        for (k = 0; k < s2; ++k)
        {
          // vector
          if (SetOfqCharges[k] == Vector_m)
            NewField.q_sh = SetOfqCharges[k] - Vacuum.InvariantSupercharges[0];
          else
          if (SetOfqCharges[k] == Vector_p)
            NewField.q_sh = SetOfqCharges[k] + Vacuum.InvariantSupercharges[0];
          else
          // left-chiral
          if (fabs(SetOfqCharges[k][0] + 0.5) < prec)
            NewField.q_sh = SetOfqCharges[k] + Vacuum.InvariantSupercharges[0];
          else
          // right-chiral
          if (fabs(SetOfqCharges[k][0] - 0.5) < prec)
            NewField.q_sh = SetOfqCharges[k] - Vacuum.InvariantSupercharges[0];
        }
      }
      else
      {
        for (k = 0; k < s2; ++k)
        {
          if (fabs(SetOfqCharges[k][0]) < prec)
            NewField.q_sh = SetOfqCharges[k];
        }
      }

      const bool LC = (NewField.Multiplet == LCModulus);

      // begin: untwisted modulus
      if (LC || (NewField.Multiplet == RCModulus))
      {
        double checkT = -1.0;
        if (LC)
          checkT = 1.0;

        unsigned Index_q = 0;
        unsigned Index_osci = 0;

        bool TModulus = true;
        for (k = 1; k < 4; ++k)
        {
          if (fabs(NewField.q_sh[k]) > prec)
            Index_q = k;
          if (fabs(NewField.OsciContribution[k]) > prec)
          {
            if (fabs(NewField.OsciContribution[k] + checkT) < prec)
              TModulus = false;
            Index_osci = k;
          }
        }
        if (Index_q == Index_osci)
          NewField.Numbers[0] = Index_q;
        else
          NewField.Numbers[0] = 10 * Index_q + Index_osci;

        if (LC)
        {
          if (TModulus)
          {
            NewField.Labels[0]  = "T";
            ++FieldCounters[8];
          }
          else
          {
            NewField.Labels[0]  = "U";
            ++FieldCounters[9];
          }
        }
        else
        {
          if (TModulus)
          {
            NewField.Labels[0]  = "bT";
            ++FieldCounters[10];
          }
          else
          {
            NewField.Labels[0]  = "bU";
            ++FieldCounters[11];
          }
        }
      }
      // end: untwisted modulus
      else {
      // begin: other multiplets
    	switch (NewField.Multiplet) {
          case LeftChiral:
          {
            NewField.Labels[0]  = "F";
            NewField.Numbers[0] = FieldCounters[0];
            ++FieldCounters[0];
            break;
          }
          case RightChiral:
          {
            NewField.Labels[0]  = "bF";
            NewField.Numbers[0] = FieldCounters[1];
            ++FieldCounters[1];
            break;
          }
          case Vector:
          {
            NewField.Labels[0]  = "V";
            NewField.Numbers[0] = FieldCounters[2];
            ++FieldCounters[2];
            break;
          }
          case VectorCC:
          {
            NewField.Labels[0]  = "bV";
            NewField.Numbers[0] = FieldCounters[3];
            ++FieldCounters[3];
            break;
          }
          case Hyper:
          {
            NewField.Labels[0]  = "H";
            NewField.Numbers[0] = FieldCounters[4];
            ++FieldCounters[4];
            break;
          }
          case Halfhyper:
          {
            NewField.Labels[0]  = "HH";
            NewField.Numbers[0] = FieldCounters[5];
            ++FieldCounters[5];
            break;
          }
        case Gravity:
          {
            NewField.Labels[0]  = "G";
            NewField.Numbers[0] = FieldCounters[6];
            ++FieldCounters[6];
            break;
          }
          case GravityCC:
          {
            NewField.Labels[0]  = "bG";
            NewField.Numbers[0] = FieldCounters[7];
            ++FieldCounters[7];
            break;
          }

        default :					//hacking here to put check for undefined multiplet type outside loop
          {
            cout << "Warning in bool CState::CreateRepresentations(...): SUSY type not defined. Return false." << endl;
           // return true;
          }
        }
      }
      // end: other multiplets
      }*/

			if (NewField.Multiplet == NOT_DEF_SUSY) {
				cout << "Warning in bool CState::CreateRepresentations(...): SUSY type not defined. Return false.!" << endl;
				return false;
			}
			Fields.push_back(NewField);
		}
	}
	// end: create new fields

	return true;
}



/* ########################################################################################
######   DetermineDimension(const gaugeGroupFactor<double> &ggf, ...) const          ######
######                                                                               ######
######   Version: 26.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) ggf              : current gauge group factor under consideration        ######
######   2) HighestWeight_DL : the highets weight of a representation of "ggf"       ######
######                         in Dynkin labels                                      ######
######   3) DimOfRep         : reference to an SDimension object where the result    ######
######                         will be saved                                         ######
######   output:                                                                     ######
######   return value        : representation identified and result saved in         ######
######                         "DimOfRep"?                                           ######
###########################################################################################
######   description:                                                                ######
######   Identifies the dimension of the representation specified by the highest     ######
######   weight "HighestWeight_DL" of the gauge group "ggf" and saves the result to  ######
######   "DimOfRep".                                                                 ######
######################################################################################## */
bool CState::DetermineDimension(const gaugeGroupFactor<double> &ggf, const intVector &HighestWeight_DL, SDimension &DimOfRep) const
{
	const string adj = "adj";

	const unsigned rank = ggf.rank;

	if (HighestWeight_DL.size() != rank)
	{
		cout << "Warning in bool CState::DetermineDimension(...) const: Length of the highest vector is not equal to the rank. Return false." << endl;
		return false;
	}

	intVector Singlet(rank, 0);
	if (HighestWeight_DL == Singlet)
	{
		DimOfRep.Dimension = 1;
		DimOfRep.AdditionalLabel = "";
		return true;
	}

	switch (ggf.algebra[0])
	{
	////////////////////////////////////////////////////////////////////
	// A
	////////////////////////////////////////////////////////////////////
	case 'A':
	{
		switch (rank)
		{
		////////////////////////////////////////////////////////////////
		// A1
		////////////////////////////////////////////////////////////////
		case 1:
		{
			if (HighestWeight_DL[0] == 1)
			{
				DimOfRep.Dimension = 2;
				DimOfRep.AdditionalLabel = "";
				return true;
			}
			if (HighestWeight_DL[0] == 2)
			{
				DimOfRep.Dimension = 3;
				DimOfRep.AdditionalLabel = adj;
				return true;
			}
		}
		////////////////////////////////////////////////////////////////
		// A2
		////////////////////////////////////////////////////////////////
		case 2:
		{
			intVector A2_3(2,0);     A2_3[0]     = 1;
			intVector A2_3bar(2,0);  A2_3bar[1]  = 1;
			intVector A2_8(2,0);     A2_8[0]     = 1; A2_8.at(1) = 1;

			if (HighestWeight_DL == A2_3)
			{
				DimOfRep.Dimension = 3;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A2_3bar)
			{
				DimOfRep.Dimension = -3;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A2_8)
			{
				DimOfRep.Dimension = 8;
				DimOfRep.AdditionalLabel = adj;
				return true;
			}
		}
		////////////////////////////////////////////////////////////////
		// A3
		////////////////////////////////////////////////////////////////
		case 3:
		{
			intVector A3_4(3,0);     A3_4[0]     = 1;
			intVector A3_4bar(3,0);  A3_4bar[2]  = 1;
			intVector A3_6(3,0);     A3_6[1]     = 1;
			intVector A3_15(3,0);    A3_15[0]    = 1; A3_15[2] = 1;

			if (HighestWeight_DL == A3_4)
			{
				DimOfRep.Dimension = 4;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A3_4bar)
			{
				DimOfRep.Dimension = -4;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A3_6)
			{
				DimOfRep.Dimension = 6;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A3_15)
			{
				DimOfRep.Dimension = 15;
				DimOfRep.AdditionalLabel = adj;
				return true;
			}
		}
		////////////////////////////////////////////////////////////////
		// A4
		////////////////////////////////////////////////////////////////
		case 4:
		{
			intVector A4_5(4,0);     A4_5[0]     = 1;
			intVector A4_5bar(4,0);  A4_5bar[3]  = 1;
			intVector A4_10(4,0);    A4_10[1]    = 1;
			intVector A4_10bar(4,0); A4_10bar[2] = 1;
			intVector A4_24(4,0);    A4_24[0]    = 1; A4_24[3] = 1;

			if (HighestWeight_DL == A4_5)
			{
				DimOfRep.Dimension = 5;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A4_5bar)
			{
				DimOfRep.Dimension = -5;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A4_10)
			{
				DimOfRep.Dimension = 10;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A4_10bar)
			{
				DimOfRep.Dimension = -10;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A4_24)
			{
				DimOfRep.Dimension = 24;
				DimOfRep.AdditionalLabel = adj;
				return true;
			}
		}
		////////////////////////////////////////////////////////////////
		// A5
		////////////////////////////////////////////////////////////////
		case 5:
		{
			intVector A5_6(5,0);     A5_6[0]     = 1;
			intVector A5_6bar(5,0);  A5_6bar[4]  = 1;
			intVector A5_20(5,0);    A5_20[2]    = 1;
			intVector A5_15(5,0);    A5_15[1]    = 1;
			intVector A5_15bar(5,0); A5_15bar[3] = 1;
			intVector A5_35(5,0);    A5_35[0]    = 1; A5_35[4] = 1;

			if (HighestWeight_DL == A5_6)
			{
				DimOfRep.Dimension = 6;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A5_6bar)
			{
				DimOfRep.Dimension = -6;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A5_15)
			{
				DimOfRep.Dimension = 15;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A5_15bar)
			{
				DimOfRep.Dimension = -15;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A5_20)
			{
				DimOfRep.Dimension = 20;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A5_35)
			{
				DimOfRep.Dimension = 35;
				DimOfRep.AdditionalLabel = adj;
				return true;
			}
		}
		////////////////////////////////////////////////////////////////
		// A6
		////////////////////////////////////////////////////////////////
		case 6:
		{
			intVector A6_7(6,0);     A6_7[0]     = 1;
			intVector A6_7bar(6,0);  A6_7bar[5]  = 1;
			intVector A6_21(6,0);    A6_21[1]    = 1;
			intVector A6_21bar(6,0); A6_21bar[4] = 1;
			intVector A6_35(6,0);    A6_35[2]    = 1;
			intVector A6_35bar(6,0); A6_35bar[3] = 1;
			intVector A6_48(6,0);    A6_48[0]    = 1; A6_48[5] = 1;

			if (HighestWeight_DL == A6_7)
			{
				DimOfRep.Dimension = 7;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A6_7bar)
			{
				DimOfRep.Dimension = -7;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A6_21)
			{
				DimOfRep.Dimension = 21;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A6_21bar)
			{
				DimOfRep.Dimension = -21;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A6_35)
			{
				DimOfRep.Dimension = 35;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A6_35bar)
			{
				DimOfRep.Dimension = -35;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A6_48)
			{
				DimOfRep.Dimension = 48;
				DimOfRep.AdditionalLabel = adj;
				return true;
			}
		}
		////////////////////////////////////////////////////////////////
		// A7
		////////////////////////////////////////////////////////////////
		case 7:
		{
			intVector A7_8(7,0);     A7_8[0]     = 1;
			intVector A7_8bar(7,0);  A7_8bar[6]  = 1;
			intVector A7_28(7,0);    A7_28[1]    = 1;
			intVector A7_28bar(7,0); A7_28bar[5] = 1;
			intVector A7_56(7,0);    A7_56[2]    = 1;
			intVector A7_56bar(7,0); A7_56bar[4] = 1;
			intVector A7_63(7,0);    A7_63[0]    = 1; A7_63[6] = 1;
			intVector A7_70(7,0);    A7_70[3]    = 1;

			if (HighestWeight_DL == A7_8)
			{
				DimOfRep.Dimension = 8;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A7_8bar)
			{
				DimOfRep.Dimension = -8;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A7_28)
			{
				DimOfRep.Dimension = 28;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A7_28bar)
			{
				DimOfRep.Dimension = -28;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A7_56)
			{
				DimOfRep.Dimension = 56;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A7_56bar)
			{
				DimOfRep.Dimension = -56;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A7_63) // gauge group
			{
				DimOfRep.Dimension = 63;
				DimOfRep.AdditionalLabel = adj;
				return true;
			}

			if (HighestWeight_DL == A7_70)
			{
				DimOfRep.Dimension = 70;
				DimOfRep.AdditionalLabel = "";
				return true;
			}
		}
		////////////////////////////////////////////////////////////////
		// A8
		////////////////////////////////////////////////////////////////
		case 8:
		{
			intVector A8_9(8,0);     A8_9[0]     = 1;
			intVector A8_9bar(8,0);  A8_9bar[7]  = 1;
			intVector A8_36(8,0);    A8_36[1]    = 1;
			intVector A8_36bar(8,0); A8_36bar[6] = 1;
			intVector A8_80(8,0);    A8_80[0]    = 1; A8_80[7]    = 1;
			intVector A8_84(8,0);    A8_84[2]    = 1;
			intVector A8_84bar(8,0); A8_84bar[5] = 1;

			if (HighestWeight_DL == A8_9)
			{
				DimOfRep.Dimension = 9;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A8_9bar)
			{
				DimOfRep.Dimension = -9;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A8_36)
			{
				DimOfRep.Dimension = 36;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A8_36bar)
			{
				DimOfRep.Dimension = -36;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A8_80)
			{
				DimOfRep.Dimension = 80;
				DimOfRep.AdditionalLabel = adj;
				return true;
			}

			if (HighestWeight_DL == A8_84)
			{
				DimOfRep.Dimension = 84;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A8_84bar)
			{
				DimOfRep.Dimension = -84;
				DimOfRep.AdditionalLabel = "";
				return true;
			}
		}
		////////////////////////////////////////////////////////////////
		// A9
		////////////////////////////////////////////////////////////////
		case 9:
		{
			intVector A9_10(9,0);    A9_10[0]     = 1;
			intVector A9_10bar(9,0); A9_10bar[8]  = 1;
			intVector A9_45(9,0);    A9_45[1]    = 1;
			intVector A9_45bar(9,0); A9_45bar[7] = 1;
			intVector A9_99(9,0);    A9_99[0]    = 1; A9_99[8]    = 1;

			if (HighestWeight_DL == A9_10)
			{
				DimOfRep.Dimension = 10;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A9_10bar)
			{
				DimOfRep.Dimension = -10;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A9_45)
			{
				DimOfRep.Dimension = 45;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A9_45bar)
			{
				DimOfRep.Dimension = -45;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A9_99)
			{
				DimOfRep.Dimension = 99;
				DimOfRep.AdditionalLabel = adj;
				return true;
			}
		}
		////////////////////////////////////////////////////////////////
		// A10
		////////////////////////////////////////////////////////////////
		case 10:
		{
			intVector A10_11(10,0);    A10_11[0]    = 1;
			intVector A10_11bar(10,0); A10_11bar[9] = 1;
			intVector A10_55(10,0);    A10_55[1]    = 1;
			intVector A10_55bar(10,0); A10_55bar[8] = 1;
			intVector A10_120(10,0);   A10_120[0]   = 1; A10_120[9] = 1;

			if (HighestWeight_DL == A10_11)
			{
				DimOfRep.Dimension = 11;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A10_11bar)
			{
				DimOfRep.Dimension = -11;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A10_55)
			{
				DimOfRep.Dimension = 55;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A10_55bar)
			{
				DimOfRep.Dimension = -55;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A10_120)
			{
				DimOfRep.Dimension = 120;
				DimOfRep.AdditionalLabel = adj;
				return true;
			}
		}
		////////////////////////////////////////////////////////////////
		// A11
		////////////////////////////////////////////////////////////////
		case 11:
		{
			intVector A11_12(11,0);    A11_12[0]    = 1;
			intVector A11_12bar(11,0); A11_12bar[10] = 1;
			intVector A11_66(11,0);    A11_66[1]    = 1;
			intVector A11_66bar(11,0); A11_66bar[9] = 1;
			intVector A11_143(11,0);   A11_143[0]   = 1; A11_143[10]= 1;

			if (HighestWeight_DL == A11_12)
			{
				DimOfRep.Dimension = 12;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A11_12bar)
			{
				DimOfRep.Dimension = -12;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A11_66)
			{
				DimOfRep.Dimension = 66;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A11_66bar)
			{
				DimOfRep.Dimension = -66;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A11_143)
			{
				DimOfRep.Dimension = 143;
				DimOfRep.AdditionalLabel = adj;
				return true;
			}
		}
		////////////////////////////////////////////////////////////////
		// A12
		////////////////////////////////////////////////////////////////
		case 12:
		{
			intVector A12_13(12,0);    A12_13[0]    = 1;
			intVector A12_13bar(12,0); A12_13bar[11] = 1;
			intVector A12_78(12,0);    A12_78[1]    = 1;
			intVector A12_78bar(12,0); A12_78bar[10] = 1;
			intVector A12_168(12,0);   A12_168[0]   = 1; A12_168[11]= 1;

			if (HighestWeight_DL == A12_13)
			{
				DimOfRep.Dimension = 13;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A12_13bar)
			{
				DimOfRep.Dimension = -13;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A12_78)
			{
				DimOfRep.Dimension = 78;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A12_78bar)
			{
				DimOfRep.Dimension = -78;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A12_168)
			{
				DimOfRep.Dimension = 168;
				DimOfRep.AdditionalLabel = adj;
				return true;
			}
		}
		////////////////////////////////////////////////////////////////
		// A13
		////////////////////////////////////////////////////////////////
		case 13:
		{
			intVector A13_14(13,0);    A13_14[0]    = 1;
			intVector A13_14bar(13,0); A13_14bar[12] = 1;
			intVector A13_91(13,0);    A13_91[1]    = 1;
			intVector A13_91bar(13,0); A13_91bar[11] = 1;
			intVector A13_195(13,0);   A13_195[0]   = 1; A13_195[12]= 1;

			if (HighestWeight_DL == A13_14)
			{
				DimOfRep.Dimension = 14;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A13_14bar)
			{
				DimOfRep.Dimension = -14;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A13_91)
			{
				DimOfRep.Dimension = 91;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A13_91bar)
			{
				DimOfRep.Dimension = -91;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A13_195)
			{
				DimOfRep.Dimension = 195;
				DimOfRep.AdditionalLabel = adj;
				return true;
			}
		}
		////////////////////////////////////////////////////////////////
		// A14
		////////////////////////////////////////////////////////////////
		case 14:
		{
			intVector A14_15(14,0);    A14_15[0]    = 1;
			intVector A14_15bar(14,0); A14_15bar[13] = 1;
			intVector A14_105(14,0);    A14_105[1]    = 1;
			intVector A14_105bar(14,0); A14_105bar[12] = 1;
			intVector A14_224(14,0);   A14_224[0]   = 1; A14_224[13]= 1;

			if (HighestWeight_DL == A14_15)
			{
				DimOfRep.Dimension = 15;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A14_15bar)
			{
				DimOfRep.Dimension = -15;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A14_105)
			{
				DimOfRep.Dimension = 105;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A14_105bar)
			{
				DimOfRep.Dimension = -105;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A14_224)
			{
				DimOfRep.Dimension = 224;
				DimOfRep.AdditionalLabel = adj;
				return true;
			}
		}
		////////////////////////////////////////////////////////////////
		// A15
		////////////////////////////////////////////////////////////////
		case 15:
		{
			intVector A15_16(15,0);    A15_16[0]    = 1;
			intVector A15_16bar(15,0); A15_16bar[14] = 1;
			intVector A15_120(15,0);    A15_120[1]    = 1;
			intVector A15_120bar(15,0); A15_120bar[13] = 1;
			intVector A15_255(15,0);   A15_255[0]   = 1; A15_255[14]= 1;

			if (HighestWeight_DL == A15_16)
			{
				DimOfRep.Dimension = 16;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A15_16bar)
			{
				DimOfRep.Dimension = -16;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A15_120)
			{
				DimOfRep.Dimension = 120;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A15_120bar)
			{
				DimOfRep.Dimension = -120;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == A15_255)
			{
				DimOfRep.Dimension = 255;
				DimOfRep.AdditionalLabel = adj;
				return true;
			}
		}

		}
		break;
	}
	////////////////////////////////////////////////////////////////////
	// D
	////////////////////////////////////////////////////////////////////
	case 'D':
	{
		intVector DN_vector(rank,0);
		DN_vector[0] = 1;
		if (HighestWeight_DL == DN_vector)
		{
			DimOfRep.Dimension = 2*rank;

			if (rank == 4)
				DimOfRep.AdditionalLabel = "v";
			else
				DimOfRep.AdditionalLabel = "";

			return true;
		}

		intVector DN_adjoint(rank,0);
		DN_adjoint[1] = 1;
		if (HighestWeight_DL == DN_adjoint)
		{
			DimOfRep.Dimension = (2*rank*rank) - rank;
			DimOfRep.AdditionalLabel = adj;
			return true;
		}

		if (rank == 4)
		{
			intVector D4_8s(4,0);  D4_8s[3] = 1;
			intVector D4_8c(4,0);  D4_8c[2] = 1;

			if (HighestWeight_DL == D4_8s)
			{
				DimOfRep.Dimension = 8;
				DimOfRep.AdditionalLabel = "s";
				return true;
			}

			if (HighestWeight_DL == D4_8c)
			{
				DimOfRep.Dimension = 8;
				DimOfRep.AdditionalLabel = "c";
				return true;
			}
		}

		if (rank == 5)
		{
			intVector D5_16(5,0);     D5_16[4]     = 1;
			intVector D5_16bar(5,0);  D5_16bar[3]  = 1;

			if (HighestWeight_DL == D5_16)
			{
				DimOfRep.Dimension = 16;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == D5_16bar)
			{
				DimOfRep.Dimension =-16;
				DimOfRep.AdditionalLabel = "";
				return true;
			}
		}
		if (rank == 6)
		{
			intVector D6_32(6,0);     D6_32[5]     = 1;
			intVector D6_32bar(6,0);  D6_32bar[4]  = 1;

			if (HighestWeight_DL == D6_32)
			{
				DimOfRep.Dimension = 32;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == D6_32bar)
			{
				DimOfRep.Dimension = -32;
				DimOfRep.AdditionalLabel = "";
				return true;
			}
		}
		if (rank == 7)
		{
			intVector D7_64(7,0);     D7_64[6]     = 1;
			intVector D7_64bar(7,0);  D7_64bar[5]  = 1;

			if (HighestWeight_DL == D7_64)
			{
				DimOfRep.Dimension = 64;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == D7_64bar)
			{
				DimOfRep.Dimension = -64;
				DimOfRep.AdditionalLabel = "";
				return true;
			}
		}
		if (rank == 8)
		{
			intVector D8_128(8,0);    D8_128[7]    = 1;
			intVector D8_128bar(8,0); D8_128bar[6] = 1;

			if (HighestWeight_DL == D8_128)
			{
				DimOfRep.Dimension = 128;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == D8_128bar)
			{
				DimOfRep.Dimension = -128;
				DimOfRep.AdditionalLabel = "";
				return true;
			}
		}
		break;
	}
	////////////////////////////////////////////////////////////////////
	// E
	////////////////////////////////////////////////////////////////////
	case 'E':
	{
		//////////////////////////////////////////////////////////////////
		// E6
		//////////////////////////////////////////////////////////////////
		if (rank == 6)
		{
			intVector E6_27(6,0);    E6_27[0]    = 1;
			intVector E6_27bar(6,0); E6_27bar[4] = 1;
			intVector E6_78(6,0);    E6_78[5]    = 1;

			if (HighestWeight_DL == E6_27)
			{
				DimOfRep.Dimension = 27;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == E6_27bar)
			{
				DimOfRep.Dimension = -27;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == E6_78)
			{
				DimOfRep.Dimension = 78;
				DimOfRep.AdditionalLabel = adj;
				return true;
			}
		}
		//////////////////////////////////////////////////////////////////
		// E7
		//////////////////////////////////////////////////////////////////
		if (rank == 7)
		{
			intVector E7_56(7,0);    E7_56[5]    = 1;
			intVector E7_133(7,0);   E7_133[0]   = 1;

			if (HighestWeight_DL == E7_56)
			{
				DimOfRep.Dimension = 56;
				DimOfRep.AdditionalLabel = "";
				return true;
			}

			if (HighestWeight_DL == E7_133)
			{
				DimOfRep.Dimension = 133;
				DimOfRep.AdditionalLabel = adj;
				return true;
			}
		}
		//////////////////////////////////////////////////////////////////
		// E8
		//////////////////////////////////////////////////////////////////
		if (rank == 8)
		{
			intVector E8_248(8,0);   E8_248[6]   = 1;

			if (HighestWeight_DL == E8_248)
			{
				DimOfRep.Dimension = 248;
				DimOfRep.AdditionalLabel = adj;
				return true;
			}
		}
		break;
	}
	}
	cout << "Warning in bool CState::DetermineDimension(...) const: Highest weight " << HighestWeight_DL << " of " << ggf.algebra << " unknown. Return false." << endl;
	return false;
}



/* ########################################################################################
######   GetFieldIndices(const vector<CField> &Fields, ...) const                    ######
######                                                                               ######
######   Version: 26.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Fields       : vector of CField objects where to look for "Multiplet"    ######
######   2) Multiplet    : the SUSY type (e.g. left-chiral) to look for              ######
######   3) FieldIndices : list of internal field-indices from "Fields" of SUSY type ######
######                     "Multiplet"                                               ######
######   output:                                                                     ######
######   return value     : finished successfully?                                   ######
###########################################################################################
######   description:                                                                ######
######   Finds all fields from "Fields" with SUSY type "Multiplet" and stores the    ######
######   corresponding field indices in "FieldIndices". For example, find all        ######
######   "LeftChiral" fields of this state.                                          ######
######################################################################################## */
bool CState::GetFieldIndices(const vector<CField> &Fields, const SUSYMultiplet &Multiplet, vector<unsigned> &FieldIndices) const
{
	const bool UseAnySUSYKind = (Multiplet == AnyKind);

	const size_t f1 = Fields.size();
	for (unsigned i = 0; i < f1; ++i)
	{
		const CField &Field = Fields[i];

		if ((Field.GetInternalIndex() == this->internalIndex) && (UseAnySUSYKind || (Multiplet == Field.Multiplet)))
			FieldIndices.push_back(i);
	}
	return true;
}



/* ########################################################################################
######   RecalculateU1Charges(...) const                                             ######
######                                                                               ######
######   Version: 11.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) AllWeightsFromFixedBrane : list of all weights resulting from the local  ######
######                                 shift V_loc of the current fixed point        ######
######   2) Vacuum                   : the SConfig object contains the new U(1)      ######
######                                 generators and the list of fields whose U(1)  ######
######                                 charges have to be recalculated               ######
######   output:                                                                     ######
######   return value                : finished successfully?                        ######
###########################################################################################
######   description:                                                                ######
######   Recalculate the U(1) charges of the fields from "Vacuum" using the U(1)     ######
######   generators, also specified in "Vacuum". Used by "Config_ComputeNewU1Charges"######
######   in the class "COrbifold".                                                   ######
######################################################################################## */
bool CState::RecalculateU1Charges(const vector<CVector> &AllWeightsFromFixedBrane, SConfig &Vacuum) const
{
	const CGaugeGroup &GaugeGroup = Vacuum.SymmetryGroup.GaugeGroup;
	vector<CField>    &Fields     = Vacuum.Fields;

	const size_t number_of_U1s = GaugeGroup.u1directions.size();

	long double Charge = 0.0;

	size_t s1 = 0;
	unsigned i = 0;
	unsigned j = 0;
	unsigned k = 0;

	// run through all representations (and corresponding fields)
	s1 = this->FieldIndices.size();
	for (i = 0; i < s1; ++i)
	{
		CField                 &Field   = Fields[this->FieldIndices[i]];
		const vector<unsigned> &Weights = Field.WeightIndices;

		CVector &U1Charges = Field.U1Charges;
		U1Charges.Assign(number_of_U1s, 0.0);

		if (Weights.size() == 0)
		{
			cout << "Warning in bool CState::RecalculateU1Charges(...) const: Weights of the representation are unknown. Return false." << endl;
			return false;
		}
		// take the first weight and compute the U(1) charges
		const CVector &weight = AllWeightsFromFixedBrane.at(Weights[0]);
		for (j = 0; j < number_of_U1s; ++j)
		{
			const doubleVector &U1Direction = GaugeGroup.u1directions[j];

			// compute the U(1) charge
			Charge = 0.0;
			for (k = 0; k < 16; ++k)
				Charge += U1Direction[k] * weight[k];

			RoundCharge(Charge);
			U1Charges[j] = (double)Charge;
		}
	}
	return true;
}



/* ########################################################################################
######   &CState::GetLeftMover() const                                               ######
######                                                                               ######
######   Version: 17.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : reference to the left-moving CHalfState object.              ######
###########################################################################################
######   description:                                                                ######
######   Allows access to the content of the private member variable "LeftMover".    ######
######################################################################################## */
const CHalfState &CState::GetLeftMover() const
{
#ifdef CHECKERROR
	if (this->State_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const CHalfState &CState::GetLeftMover() const : State ill-defined." << endl;
#endif

	return this->LeftMover;
}



/* ########################################################################################
######   &CState::GetRightMover() const                                              ######
######                                                                               ######
######   Version: 17.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : reference to the right-moving CHalfState object.             ######
###########################################################################################
######   description:                                                                ######
######   Allows access to the content of the private member variable "RightMover".   ######
######################################################################################## */
const CHalfState &CState::GetRightMover() const
{
#ifdef CHECKERROR
	if (this->State_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const CHalfState &CState::GetRightMover() const : State ill-defined." << endl;
#endif

	return this->RightMover;
}



/* ########################################################################################
######   &GetGammaPhase(const unsigned &i) const                                     ######
######                                                                               ######
######   Version: 17.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) i         : index in the private member variable "gamma_phases"          ######
######   output:                                                                     ######
######   return value : reference to the gamma phase in the i-th position of         ######
######                  "gamma_phases"                                               ######
###########################################################################################
######   description:                                                                ######
######   Allows access to the content of the private member variable "gamma_phases". ######
######################################################################################## */
const double &CState::GetGammaPhase(const unsigned &i) const
{
#ifdef CHECKERROR
	if (this->State_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const double &CState::GetGammaPhase(...) const : State ill-defined." << endl;

	if (i >= this->gamma_phases.size())
	{
		cout << "\n  Warning in const double &CState::GetGammaPhase(...) const : Index i out of range. Set i = 0." << endl;
		return this->gamma_phases[0];
	}
#endif

	return this->gamma_phases[i];
}



/* ########################################################################################
######   &GetFieldIndex(const unsigned &i) const                                     ######
######                                                                               ######
######   Version: 11.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) i         : index in the private member variable "FieldIndices"          ######
######   output:                                                                     ######
######   return value : reference to the field index in the i-th position of         ######
######                  "FieldIndices"                                               ######
###########################################################################################
######   description:                                                                ######
######   Allows access to the content of the private member variable "FieldIndices". ######
######################################################################################## */
const unsigned &CState::GetFieldIndex(const unsigned &i) const
{
#ifdef CHECKERROR
	if (this->State_CheckStatus != CheckedAndGood)
		cout << "\n  Warning in const unsigned &CState::GetFieldIndex(...) const : State ill-defined." << endl;

	if (i >= this->FieldIndices.size())
	{
		cout << "\n  Warning in const unsigned &CState::GetFieldIndex(...) const : Index i out of range. Set i = 0." << endl;
		return this->FieldIndices[0];
	}
#endif

	return this->FieldIndices[i];
}
