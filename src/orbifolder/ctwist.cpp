#include <iostream>
#include <cstdlib>

#include "ctwist.h"
#include "globalfunctions.h"

#define CHECKERROR true

using std::cout;
using std::endl;
using std::exit;



/* ########################################################################################
######   CTwistVector()                                                              ######
######                                                                               ######
######   Version: 19.03.2010                                                         ######
######   Check-Level: 2                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard constructor of a CTwistVector object. Creates a 4-dim. null vector.######
######################################################################################## */
CTwistVector::CTwistVector()
{
	this->SetToZero();
}



/* ########################################################################################
######   CTwistVector(double v_0, double v_1, double v_2, double v_3)                ######
######                                                                               ######
######   Version: 11.08.2010                                                         ######
######   Check-Level: 2                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) - 4) v_i              : v_0 should be zero and v_i, i=1,2,3 are the      ######
######                              rotation angles in the three complex internal    ######
######                              directions                                       ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Constructor of a CTwistVector object. Creates a 4-dim. vector with entries  ###### 
######   "v_i" for i = 0,..,3 and calls "UpdateData()" to check and update the       ######
######   content.                                                                    ######
######################################################################################## */
CTwistVector::CTwistVector(double v_0, double v_1, double v_2, double v_3)
{
	this->resize(4,0);
	this->Size = 4;

	this->at(0)  = v_0;
	this->at(1)  = v_1;
	this->at(2)  = v_2;
	this->at(3)  = v_3;

	this->UpdateData();
}



/* ########################################################################################
######   ~CTwistVector()                                                             ######
######                                                                               ######
######   Version: 11.08.2010                                                         ######
######   Check-Level: 2                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard destructor of a CTwistVector object.                               ######
######################################################################################## */
CTwistVector::~CTwistVector()
{
}



/* ########################################################################################
######   operator=(const CVector &Vector)                                            ######
######                                                                               ######
######   Version: 11.08.2010                                                         ######
######   Check-Level: 2                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Vector : a CVector object with 4 components                              ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Specify the content of this CTwistVector object by a CVector object. Then   ######
######   call "UpdateData()" to check and update the content.                        ######
######################################################################################## */
void CTwistVector::operator=(const CVector &Vector)
{
#ifdef CHECKERROR
	if (Vector.Size != 4)
	{
		cout << "\n  Warning in void CTwistVector::operator=(...): Vector has length " << Vector.size() << " != 4. Hence, it is set to zero!" << endl;
		this->SetToZero();
		return;
	}
#endif

	this->at(0) = Vector[0];
	this->at(1) = Vector[1];
	this->at(2) = Vector[2];
	this->at(3) = Vector[3];
	this->UpdateData();
}



/* ########################################################################################
######   operator+=(const CVector &Vector2)                                          ######
######                                                                               ######
######   Version: 01.04.2010                                                         ######
######   Check-Level: 2                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Vector2 : a CVector object with 4 components                             ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Uses the standard addition of two vectors. Stores the result in this        ######
######   CTwistVector object and calls "UpdateData()" to check and update the        ######
######   content.                                                                    ######
######################################################################################## */
void CTwistVector::operator+=(const CVector &Vector2)
				{
#ifdef CHECKERROR
	if (Vector2.Size != 4)
	{
		cout << "\n  Warning in void CTwistVector::operator+=(...): Vector has length " << Vector2.size() << " != 4. Hence, nothing added!" << endl;
		return;
	}
#endif

	CVector::operator +=(Vector2);
	this->UpdateData();
				}



/* ########################################################################################
######   bool CTwistVector::SetToZero()                                              ######
######                                                                               ######
######   Version: 19.03.2010                                                         ######
######   Check-Level: 2                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Set this CTwistVector object to the null vector.                            ######
######################################################################################## */
void CTwistVector::SetToZero()
{
	this->resize(4,0);
	this->Size  = 4;
	this->Order = 1;
	this->a_L   = -1.0;
	this->a_R   = -0.5;
}



/* ########################################################################################
######   UpdateData()                                                                ######
######                                                                               ######
######   Version: 16.05.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : twist vector updated successfully?                           ######
###########################################################################################
######   description:                                                                ######
######   Checks this CTwistVector object and determines the order and the left- and  ######
######   right-moving zero-point energies associated to the twisted sector of this   ######
######   twist vector. Saves the result in the member variables "Order", "a_L" and   ######
######   "a_R", respectively.                                                        ######
######################################################################################## */
bool CTwistVector::UpdateData()
{
	// Set the precision
	const double prec = 0.0001;

	unsigned i = 0;

	// compute the zero point energy
	double tmp = 0;
	double entry = 0.0;
	for (i = 0; i < 4; ++i)
	{
		entry = fmod(this->at(i), 1.0);
		if ((fabs(entry - 1.0) < prec) || (fabs(entry) < prec))
			entry = 0;
		else
			while(entry < -prec) ++entry;

		// now, "entry" should be in the interval [0,1)
		if ((entry < 0) || (entry >= 1))
		{
			cout << "\n  Warning in bool CTwistVector::UpdateData(): Cannot compute the zero point energy. Return false." << endl;
			return false;
		}

		tmp += entry * (1.0 - entry);
	}

	this->a_L = -1.0 + (0.5 * tmp);
	this->a_R = -0.5 + (0.5 * tmp);

	// compute the order of the twist
	const unsigned min_Order =  1;
	const unsigned max_Order = 20;

	bool Order_found = false;

	double sum=0.;
	for (int i=0; i<4; i++) {
		entry = (*this)[i];
		sum += entry;
	}

	double test_sum;

	if ( sum > prec )							//for SuSy-breaking twists, assuming sum>0
	{
		for (unsigned n = min_Order; n < max_Order; ++n)
		{
			test_sum = n*sum;
			if ( is_even(test_sum) ) {

				Order_found  = true;

				for (i = 0; i < 4; ++i)
				{
					entry = (*this)[i] * n;

					if (!is_integer(entry))
					{
						Order_found = false;
						break;
					}
				}
				if (Order_found) {
					this->Order = n;
					return true;
				}
			}
		}
		for (i = 0; i < 4; ++i)
			cout << (*this)[i] << " ";
		cout << "\nSum v_i = " << sum;

		cout << "\nError in determining the order of SuSy-breaking twist\n";
		return false;
	}

	for (unsigned n = min_Order; n < max_Order; ++n)
	{
		Order_found  = true;

		for (i = 0; i < 4; ++i)						//first entry is always zero
		{
			entry = (*this)[i] * n;

			if (!is_integer(entry))
			{
				Order_found = false;
				break;
			}
		}
		if (Order_found)
		{
			this->Order = n;
			return true;
		}
	}

	cout << "\n  Warning in bool CTwistVector::UpdateData(): Cannot find the order of the twist. Return false." << endl;
	return false;
}
