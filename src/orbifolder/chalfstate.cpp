#include <iostream>
#include <math.h>

#include "chalfstate.h"
#include "corbifoldgroupelement.h"
#include "clatticevector.h"
#include "globalfunctions.h"

using std::cout;
using std::endl;
using std::flush;



/* ########################################################################################
######   CHalfState(MoversType Type, unsigned Index)                                 ######
######                                                                               ######
######   Version: 26.10.2010                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Type  : "LeftMover" or "RightMover"                                      ######
######   2) Index : refers to the corresponding CMasslessHalfState object            ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Constructor of a CHalfState object.                                         ######
######################################################################################## */
CHalfState::CHalfState(MoversType Type, unsigned Index)
: Index(Index), Type(Type)
{
}



/* ########################################################################################
######   CHalfState()                                                                ######
######                                                                               ######
######   Version: 26.10.2010                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard constructor of a CHalfState object.                                ######
######################################################################################## */
CHalfState::CHalfState()
: Index(0), Type(NOT_DEF_MOVTYPE), Excited(0)
{
}



/* ########################################################################################
######   ~CHalfState()                                                               ######
######                                                                               ######
######   Version: 18.10.2010                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard destructor of a CHalfState object.                                 ######
######################################################################################## */
CHalfState::~CHalfState()
{
}



/* ########################################################################################
######   CHalfState()                                                                ######
######                                                                               ######
######   Version: 26.10.2010                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard constructor of a CTachyonHalfState object.                         ######
######################################################################################## */
CTachyonHalfState::CTachyonHalfState(MoversType Type, const S_OscillatorExcitation &Excitation)
{
	this->Type       = Type;
	this->Excitation = Excitation;
}



/* ########################################################################################
######   ~CHalfState()                                                               ######
######                                                                               ######
######   Version: 18.10.2010                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard destructor of a CTachyonHalfState object.                          ######
######################################################################################## */
CTachyonHalfState::~CTachyonHalfState()
{
}



/* ########################################################################################
######   CMasslessHalfState(MoversType Type, ...)                                    ######
######                                                                               ######
######   Version: 19.10.2010                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Type       : "LeftMover" or "RightMover"                                 ######
######   2) Excitation : specifies possible oscillator excitations of the ground     ######
######                   state                                                       ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Constructor of a CMasslessHalfState object. Call "SolveMassEquation" to     ######
######   create the massless left- or right-movers (half-states).                    ######
######################################################################################## */
CMasslessHalfState::CMasslessHalfState(MoversType Type, const S_OscillatorExcitation &Excitation)
{
	this->Type       = Type;
	this->Excitation = Excitation;
}



/* ########################################################################################
######   ~CMasslessHalfState()                                                       ######
######                                                                               ######
######   Version: 18.10.2010                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard destructor of a CMasslessHalfState object.                         ######
######################################################################################## */
CMasslessHalfState::~CMasslessHalfState()
{
}



/* ########################################################################################
######   SolveMassEquation(const CVector &constructing_Element, ...)                 ######
######                                                                               ######
######   Version: 20.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) constructing_Element : the COrbifoldGroupElement object containing the   ######
######                             local shift V_loc                                 ######
######   2) Lattice              : E8xE8, Spin32 or SO8                              ######
######   output:                                                                     ######
######   return value            : finished successfully?                            ######
###########################################################################################
######   description:                                                                ######
######   Solves the equation for massless left- or right-movers (specified by the    ######
######   member variable "Type") using the oscillator excitation (specified by the   ######
######   member variable "Excitation") and using the "constructing_Element" V_g or   ######
######   v_g. The solutions p_sh = p + V_g or q_sh = q + v_g are stored in the       ######
######   member variable "Weights", where p (or q) from "Lattice".                   ######
######################################################################################## */
bool CMasslessHalfState::SolveMassEquation(const CVector &constructing_Element, const SelfDualLattice &Lattice)
{
	if (this->Weights.size() != 0)
	{
		cout << "\n  Warning in bool CMasslessHalfState::SolveMassEquation(...) : set of weights is not empty - now cleared." << endl;
		this->Weights.clear();
	}

	const size_t dim = constructing_Element.GetSize();

	if (this->Type == LeftMover)
	{
		if (dim != 16)
		{
			cout << "\n  Warning in bool CMasslessHalfState::SolveMassEquation(...) : local shift does not have 16 components. Return false." << endl;
			return false;
		}
		if ((Lattice != E8xE8) && (Lattice != Spin32))
		{
			cout << "\n  Warning in bool CMasslessHalfState::SolveMassEquation(...) : \"Lattice\" must be \"E8xE8\" or \"Spin32\". Return false." << endl;
			return false;
		}
	}
	else
		if (this->Type == RightMover)
		{
			if (dim != 4)
			{
				cout << "\n  Warning in bool CMasslessHalfState::SolveMassEquation(...) : local shift does not have 4 components. Return false." << endl;
				return false;
			}
			if (Lattice != SO8)
			{
				cout << "\n  Warning in bool CMasslessHalfState::SolveMassEquation(...) : \"Lattice\" must be \"SO8\". Return false." << endl;
				return false;
			}
		}

	double c = this->Excitation.ZeroPointEnergy;
	if ((c > -0.00001) && (c < 0.00001))
		c = 0.0;

	if (c < 0)
	{
		cout << "\n  Warning in bool CMasslessHalfState::SolveMassEquation(...) : local shift squared is negative: " << c << ". Return false." << endl;
		return false;
	}

	const double sqrt_c = sqrt(c);

	vector<double> upper_bound(dim), lower_bound(dim);
	double bound = 0;
	double tmp   = 0;
	unsigned i = 0;

	// compute upper and lower bound for the root vectors
	// for all i: (p_i + V_i)^2 <= c <= 1 for right-movers <= for left movers

	for (i = 0; i < dim; ++i)
	{
		//upper_bound: p_i <= sqrt(c) - V_i
		bound = sqrt_c - constructing_Element[i];
		tmp   = round_double_to_int(bound);
		if (fabs(bound - tmp) > 0.00001)
		{
			tmp += 0.5;
			if (tmp-0.00001 > bound)
				tmp -= 0.5;
			if (tmp-0.00001 > bound)
				tmp -= 0.5;
		}
		upper_bound[i] = tmp;

		bound = -sqrt_c - constructing_Element[i];
		tmp   = round_double_to_int(bound);
		if (fabs(bound - tmp) > 0.00001)
		{
			tmp -= 0.5;
			if (tmp+0.00001 < bound)
				tmp += 0.5;
			if (tmp+0.00001 < bound)
				tmp += 0.5;
		}
		lower_bound[i] = tmp;

		if (lower_bound[i] > upper_bound[i])
			return false;
	}

	vector<double> lower_bound_int(dim);
	vector<double> lower_bound_halfint(dim);

	for (i = 0; i < dim; ++i)
	{
		if (fmod(lower_bound[i],1.0) != 0)
		{
			lower_bound_halfint[i] = lower_bound[i];
			lower_bound_int[i]     = lower_bound[i] + 0.5;
		}
		else
		{
			lower_bound_halfint[i] = lower_bound[i] + 0.5;
			lower_bound_int[i]     = lower_bound[i];
		}
	}

	const vector<double> *use_lower_bound = NULL;
	const vector<double> *use_lower_bound2 = NULL;

	// create new root vectors tmp which lie between the boundaries
	vector<double> p(dim, 0.0);
	CLatticeVector TestVector(dim);

	// seperate (p + V)^2 = p^2 + 2pV + V^2
	//                    = p_1^2 + p_2^2 + ...
	//                    = 2p_1 V_1 + 2 p_2 V_2 + ...
	//                    = V_1^2 + V_2^2 + ...
	vector<double> all_p_sqr(dim-1, 0.0);
	vector<double> all_pV(dim-1, 0.0);

	vector<double> all_VV(dim-1, 0.0);
	for (i = 0; i < dim-1; ++i)
		all_VV[i] = constructing_Element.GetSqrTo(i);

	vector<double> all_c_minus_VV(dim-1, 0.0);
	for (i = 0; i < dim-1; ++i)
		all_c_minus_VV[i] = c + 0.0001 - all_VV[i];

	const double VV_minus_c = constructing_Element * constructing_Element - c;

	bool go_on = true;

	for (p[0] = lower_bound[0]; p[0] <= upper_bound[0]; p[0] += 0.5)
	{
		if (fmod(p[0],1.0) != 0)
			use_lower_bound = &lower_bound_halfint;
		else
			use_lower_bound = &lower_bound_int;

		all_p_sqr[0] = p[0] * p[0];
		all_pV[0]    = 2.0 * p[0] * constructing_Element[0];

		for (p[1] = (*use_lower_bound)[1]; p[1] <= upper_bound[1]; ++p[1])
		{
			all_p_sqr[1] = all_p_sqr[0] + (p[1] * p[1]);
			all_pV[1]    = all_pV[0]    + (2.0 * p[1] * constructing_Element[1]);

			if (all_p_sqr[1] + all_pV[1] <= all_c_minus_VV[1])
			{
				for (p[2] = (*use_lower_bound)[2]; p[2] <= upper_bound[2]; ++p[2])
				{
					all_p_sqr[2] = all_p_sqr[1] + (p[2] * p[2]);
					all_pV[2]    = all_pV[1]    + (2.0 * p[2] * constructing_Element[2]);

					if (all_p_sqr[2] + all_pV[2] <= all_c_minus_VV[2])
					{
						for (p[3] = (*use_lower_bound)[3]; p[3] <= upper_bound[3]; ++p[3])
						{
							if (dim == 4)								//If a right-mover charge is computed
							{
								if (fabs(all_p_sqr[2] + (p[3] * p[3]) + all_pV[2] + (2.0 * p[3] * constructing_Element[3]) + VV_minus_c) < 0.00001)
								{
									TestVector = p;
									if (TestVector.From_SO8S_Lattice() != 0) {			//hacking here!!!! //to input only the cospinor lattice to get the spinorial of SO8 after moding by Z_2 change to SO8C
										this->Weights.push_back(TestVector + constructing_Element);
									}
								}
							}
							else
							{
								all_p_sqr[3] = all_p_sqr[2] + (p[3] * p[3]);
								all_pV[3]    = all_pV[2]    + (2.0 * p[3] * constructing_Element[3]);

								if (all_p_sqr[3] + all_pV[3] <= all_c_minus_VV[3])
								{
									for (p[4] = (*use_lower_bound)[4]; p[4] <= upper_bound[4]; ++p[4])
									{
										all_p_sqr[4] = all_p_sqr[3] + (p[4] * p[4]);
										all_pV[4]    = all_pV[3]    + (2.0 * p[4] * constructing_Element[4]);

										if (all_p_sqr[4] + all_pV[4] <= all_c_minus_VV[4])
										{
											for (p[5] = (*use_lower_bound)[5]; p[5] <= upper_bound[5]; ++p[5])
											{
												all_p_sqr[5] = all_p_sqr[4] + (p[5] * p[5]);
												all_pV[5]    = all_pV[4]    + (2.0 * p[5] * constructing_Element[5]);

												if (all_p_sqr[5] + all_pV[5] <= all_c_minus_VV[5])
												{
													for (p[6] = (*use_lower_bound)[6]; p[6] <= upper_bound[6]; ++p[6])
													{
														all_p_sqr[6] = all_p_sqr[5] + (p[6] * p[6]);
														all_pV[6]    = all_pV[5]    + (2.0 * p[6] * constructing_Element[6]);

														if (all_p_sqr[6] + all_pV[6] <= all_c_minus_VV[6])
														{
															for (p[7] = (*use_lower_bound)[7]; p[7] <= upper_bound[7]; ++p[7])
															{
																all_p_sqr[7] = all_p_sqr[6] + (p[7] * p[7]);
																all_pV[7]    = all_pV[6]    + (2.0 * p[7] * constructing_Element[7]);

																// is first half in the E8 lattice?
																go_on = true;
																if (Lattice == E8xE8)
																{
																	TestVector = p;
																	if (!TestVector.From_E8_Lattice(1))												//hacking here!!!!
																		go_on = false;
																}
																if( (all_p_sqr[7] + all_pV[7] <= all_c_minus_VV[7]) && go_on )
																{
																	// second E8 is independent of the first one
																	for(p[8] = lower_bound[8]; p[8] <= upper_bound[8]; p[8] += 0.5)
																	{
																		if (Lattice == E8xE8)
																		{
																			if (fmod(p[8],1.0) != 0)
																				use_lower_bound2 = &lower_bound_halfint;
																			else
																				use_lower_bound2 = &lower_bound_int;
																		}
																		else
																			use_lower_bound2 = use_lower_bound;

																		all_p_sqr[8] = all_p_sqr[7] + (p[8] * p[8]);
																		all_pV[8]    = all_pV[7]    + (2.0 * p[8] * constructing_Element[8]);

																		if (all_p_sqr[8] + all_pV[8] <= all_c_minus_VV[8])
																		{
																			for (p[9] = (*use_lower_bound2)[9]; p[9] <= upper_bound[9]; ++p[9])
																			{
																				all_p_sqr[9] = all_p_sqr[8] + (p[9] * p[9]);
																				all_pV[9]    = all_pV[8]    + (2.0 * p[9] * constructing_Element[9]);

																				if (all_p_sqr[9] + all_pV[9] <= all_c_minus_VV[9])
																				{
																					for (p[10] = (*use_lower_bound2)[10]; p[10] <= upper_bound[10]; ++p[10])
																					{
																						all_p_sqr[10] = all_p_sqr[9] + (p[10] * p[10]);
																						all_pV[10]    = all_pV[9]    + (2.0 * p[10] * constructing_Element[10]);

																						if (all_p_sqr[10] + all_pV[10] <= all_c_minus_VV[10])
																						{
																							for (p[11] = (*use_lower_bound2)[11]; p[11] <= upper_bound[11]; ++p[11])
																							{
																								all_p_sqr[11] = all_p_sqr[10] + (p[11] * p[11]);
																								all_pV[11]    = all_pV[10]    + (2.0 * p[11] * constructing_Element[11]);

																								if (all_p_sqr[11] + all_pV[11] <= all_c_minus_VV[11])
																								{
																									for (p[12] = (*use_lower_bound2)[12]; p[12] <= upper_bound[12]; ++p[12])
																									{
																										all_p_sqr[12] = all_p_sqr[11] + (p[12] * p[12]);
																										all_pV[12]    = all_pV[11]    + (2.0 * p[12] * constructing_Element[12]);

																										if (all_p_sqr[12] + all_pV[12] <= all_c_minus_VV[12])
																										{
																											for (p[13] = (*use_lower_bound2)[13]; p[13] <= upper_bound[13]; ++p[13])
																											{
																												all_p_sqr[13] = all_p_sqr[12] + (p[13] * p[13]);
																												all_pV[13]    = all_pV[12]    + (2.0 * p[13] * constructing_Element[13]);

																												if (all_p_sqr[13] + all_pV[13] <= all_c_minus_VV[13])
																												{
																													for (p[14] = (*use_lower_bound2)[14]; p[14] <= upper_bound[14]; ++p[14])
																													{
																														all_p_sqr[14] = all_p_sqr[13] + (p[14] * p[14]);
																														all_pV[14]    = all_pV[13]    + (2.0 * p[14] * constructing_Element[14]);

																														if (all_p_sqr[14] + all_pV[14] <= all_c_minus_VV[14])
																														{
																															for (p[15] = (*use_lower_bound2)[15]; p[15] <= upper_bound[15]; ++p[15])
																															{
																																if (fabs(all_p_sqr[14] + (p[15] * p[15]) + all_pV[14] + (2.0 * p[15] * constructing_Element[15]) + VV_minus_c) < 0.00001)
																																{

																																	TestVector = p;							//hacking here!!!!
																																	if (((Lattice == E8xE8) && TestVector.From_E8_Lattice(2)) || ((Lattice == Spin32) && TestVector.From_Spin32_Lattice()))
																																		this->Weights.push_back(TestVector + constructing_Element);
																																}
																															}
																														}
																													}
																												}
																											}
																										}
																									}
																								}
																							}
																						}
																					}
																				}
																			}
																		}
																	}
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	if (this->Weights.size() != 0)
		return true;

	return false;
}


/* ########################################################################################			//hacking here!!!
######  TachyonSolver(const CVector &constructing_Element, ...)                      ######
######                                                                               ######
######   Version: 20.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) constructing_Element : the COrbifoldGroupElement object containing the   ######
######                             local shift V_loc                                 ######
######   2) Lattice              : E8xE8, Spin32 or SO8                              ######
######   output:                                                                     ######
######   return value            : finished successfully?                            ######
###########################################################################################
######   description:                                                                ######
######   Solves the equation for massless left- or right-movers (specified by the    ######
######   member variable "Type") using the oscillator excitation (specified by the   ######
######   member variable "Excitation") and using the "constructing_Element" V_g or   ######
######   v_g. The solutions p_sh = p + V_g or q_sh = q + v_g are stored in the       ######
######   member variable "Weights", where p (or q) from "Lattice".                   ######
######################################################################################## */
bool CTachyonHalfState::TachyonSolver(const unsigned &nindex, const unsigned &kindex, const CVector &constructing_Element, const SelfDualLattice &Lattice)
{
	bool verbose = false; 
	
	if (this->Weights.size() != 0)
	{
		cout << "\n  Warning in bool CMasslessHalfState::SolveMassEquation(...) : set of weights is not empty - now cleared." << endl;
		this->Weights.clear();
	}
	if (Lattice != SO8) {
		{
			cout << "\n  Warning in bool CMasslessHalfState::SolveMassEquation(...) : \"Lattice\" must be \"SO8\". Return false." << endl;
			return false;
		}
	}

	double prec = 0.00001;

	double c = (0.5*this->Excitation.ZeroPointEnergy);		//ZeroPointEnergy=2*c
	if ((c > -prec) && (c < prec))
		c = 0.0;

	if (c < 0)
	{
			cout << "\n  Warning in bool CMasslessHalfState::SolveMassEquation(...) : local shift squared is negative: " << c << ". Return false." << endl;
			return false;
	}


	//compute q charges for m_R^2<0

	double sum, m_R;
	vector<double> q(4, 0.0);
	CLatticeVector TestVector(4);

	for (q[0]=-7.; q[0]<1.; q[0]++) {                     //for q_V    //-1<q_0<1 to for m_R<0
		for (q[1]=-7.; q[1]<6.; q[1]++) {
			for (q[2]=-7.; q[2]<6.; q[2]++) {
				for (q[3]=-7.; q[3]<6.; q[3]++) {

					m_R=0.;
					TestVector = q;

					if (TestVector.From_SO8S_Lattice() != 0) {
						for (int i=0; i<4; i++)               												//Compute (q+v_g)^2/2
							m_R += 0.5*(q[i]+constructing_Element[i])*(q[i]+constructing_Element[i]);
						m_R -= c;

						if ( (-m_R > prec) && (m_R < prec) ) {

							this->Weights.push_back(TestVector + constructing_Element);
							this->tachyonmass.push_back(m_R);

						   if(verbose)
						   {
							cout << " potential tachyons in sector (1, " << nindex << ", " << kindex << ") with weight q_i = (";
							for (int i=0; i<4; i++) {
								cout << q[i]+constructing_Element[i] <<" ";
							}
							cout << ")\n possible tacyonic Excitation N = " << this->Excitation.NumberOperator << endl;
							cout << " and mass M_R^2 = " << m_R << endl << endl;
						   }
						}
					}
				}
			}
		}
	}
	//we do not need to check for the spinors, no tachyons from there!
	/*for (q[0]=-0.5; q[0]<0.5; q[0]++) {                     //for q_S		//-1<q_0<1 to for m_R<0
		for (q[1]=-4.5; q[1]<6.; q[1]++) {
			for (q[2]=-4.5; q[2]<6.; q[2]++) {
				for (q[3]=-4.5; q[3]<6.; q[3]++) {
					m_R=0.;

					TestVector = q;
					if (TestVector.From_SO8S_Lattice() != 0) {
						for (int i=0; i<4; i++) {               //Compute m_R^2/8 = (q+v_g)^2/2 -1/2 + \delta c
							m_R += 0.5*(q[i]+constructing_Element[i])*(q[i]+constructing_Element[i]);					//formula assumes standard c==-1/2+ \delta c
						}
						if (m_R < -c) {
							this->Weights.push_back(TestVector + constructing_Element);
							tachyonmass.push_back(m_R+c);
							cout << "potential tachyons in sector (1, " << nindex  << ", " << kindex << ") with weight q_i\n";
							for (int i=0; i<4; i++) {
								cout << q[i]+constructing_Element[i] <<" ";
							}
							cout <<endl;
						}
					}
				}
			}
		}
	}*/

	if (this->Weights.size() != 0)
		return true;

	return false;

}
