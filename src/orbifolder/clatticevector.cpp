#include <iostream>
#include <cstdlib>

#include "clatticevector.h"
#include "globalfunctions.h"

#define CHECKERROR true

using std::cout;
using std::endl;



/* ########################################################################################
######   CLatticeVector()                                                            ######
######                                                                               ######
######   Version: 25.06.2009                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard constructor of a CLatticeVector object. No content is specified.   ######
######################################################################################## */
CLatticeVector::CLatticeVector()
{
}



/* ########################################################################################
######   CLatticeVector(const unsigned length)                                       ######
######                                                                               ######
######   Version: 02.05.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) dim : the dimension of the lattice vector                                ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Constructor of a CLatticeVector object. An empty "dim"-dimensional vector   ######
######   is created.                                                                 ######
######################################################################################## */
CLatticeVector::CLatticeVector(const unsigned dim)
{
	this->assign(dim, 0.0);
	this->Size = dim;
}



/* ########################################################################################
######   CLatticeVector(const CVector &Vector2)                                      ######
######                                                                               ######
######   Version: 25.06.2009                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Vector2 : a CVector object                                               ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Constructor of a CLatticeVector object. The content of "Vector2" is         ######
######   assigned to this CLatticeVector.                                            ######
######################################################################################## */
CLatticeVector::CLatticeVector(const CVector &Vector2)
{
	this->clear();
	this->Size = Vector2.size();

	this->insert(this->end(), Vector2.begin(), Vector2.end());
}



/* ########################################################################################
######   ~CLatticeVector()                                                           ######
######                                                                               ######
######   Version: 25.06.2009                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard destructor of a CLatticeVector object.                             ######
######################################################################################## */
CLatticeVector::~CLatticeVector()
{
}



/* ########################################################################################
######   operator=(const CVector &Vector)                                            ######
######                                                                               ######
######   Version: 25.06.2009                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Vector2 : a CVector object                                               ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Specify the content of this CLatticeVector object from a CVector object.    ######
######################################################################################## */
void CLatticeVector::operator=(const CVector &Vector)
{
	this->clear();
	this->Size = Vector.size();

	this->insert(this->end(), Vector.begin(), Vector.end());
}



/* ########################################################################################
######   operator=(const vector<double> &DVector)                                    ######
######                                                                               ######
######   Version: 25.06.2009                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) DVector : a vector of double                                             ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Specify the content of this CLatticeVector object from a vector of doubles. ######
######################################################################################## */
void CLatticeVector::operator=(const vector<double> &DVector)
{
	this->clear();
	this->Size = DVector.size();

	this->insert(this->end(), DVector.begin(), DVector.end());
}



/* ########################################################################################
######   From_E8_Lattice(unsigned part) const                                        ######
######                                                                               ######
######   Version: 16.05.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) part      : "part" specifies first (1) or second (2) E8 factor           ######
######   output:                                                                     ######
######   return value : is the first/second part of the vector from E8?              ######
###########################################################################################
######   description:                                                                ######
######   Checks whether the first or second part of this CLatticeVector is from the  ######
######   E_8 root lattice.                                                           ######
######################################################################################## */
bool CLatticeVector::From_E8_Lattice(unsigned part) const
{
	unsigned from = 0;
	unsigned to   = 8;

	if (part == 2)
	{
		from = 8;
		to   = 16;
	}

#ifdef CHECKERROR
	if (this->Size < to)
	{
		cout << "Warning in bool CLatticeVector::From_E8_Lattice(...) const : vector has length != " << to << ".  Return false." << endl;
		return false;
	}
#endif

	const double prec = 0.0001;

	// begin: all entries integer or all entries half-integer
	double tmp = (*this)[from];
	const double first_entry = fabs(round_double_to_int(tmp) - tmp);

	if ((first_entry > prec) && (fabs(first_entry - 0.5) > prec))
		return false;

	unsigned i = from + 1;
	for (; i < to; ++i)
	{
		tmp = (*this)[i];

		if (fabs(first_entry - fabs(round_double_to_int(tmp) - tmp)) > prec)
			return false;
	}
	// end: all entries integer or all entries half-integer

	// sum over all entries
	double sum = 0.0;
	for (i = from; i < to; ++i)
		sum += (*this)[i];

	tmp = round_double_to_int(sum);

	// sum over all entries must be even
	if ((fabs(tmp - sum) < prec) && (fmod(tmp, 2.0) == 0))
		return true;

	return false;
}


/* ########################################################################################		//hacking here!!!!!
######   From_G8_Lattice(unsigned part) const                                        ######
######                                                                               ######
######   Version: 16.05.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) part      : "part" specifies first (1) or second (2) G8 factor           ######
######   output:                                                                     ######
######   return value : is the first/second part of the vector from G8?              ######
###########################################################################################
######   description:                                                                ######
######   Checks whether the first or second part of this CLatticeVector is from the  ######
######   Gamma_8 root lattice or cospinor of Gamma_8.                                ######
######################################################################################## */
bool CLatticeVector::From_G8_Lattice(unsigned part) const
{
	unsigned from = 0;
	unsigned to   = 8;

	if (part == 2)
	{
		from = 8;
		to   = 16;
	}

#ifdef CHECKERROR
	if (this->Size < to)
	{
		cout << "Warning in bool CLatticeVector::From_G8_Lattice(...) const : vector has length != " << to << ".  Return false." << endl;
		return false;
	}
#endif

	const double prec = 0.0001;

	// begin: all entries integer or all entries half-integer
	double tmp = (*this)[from];
	unsigned i = from + 1;

	// assume INT type
	bool INTType = true;

	if (fabs(fmod(tmp,1)) < prec) {					//if integer
		for (; i < to; ++i) {
			tmp = (*this)[i];
			if (fabs(fmod(tmp,1)) > prec) {
				return false;
			}
		}
	}
	else if (fabs(fmod((tmp-0.5),1)) < prec) {		//if half-integer
		for (; i < to; ++i) {
			tmp = (*this)[i];
			if (fabs(fmod((tmp-0.5),1)) > prec) {
				return false;
			}
		}
		INTType = false;
	}
	// end: all entries integer or all entries half-integer

	// sum over all entries must be integer
	double sum = 0.0;
	for (i = from; i < to; ++i) {
		sum += (*this)[i];
	}
	tmp = round_double_to_int(sum);

	if (INTType==true) {													// sum over all entries must be even
		if ( (fabs(tmp - sum) < prec ) && ( fabs(fmod(tmp, 2.0)) < prec ) ) {
			return true;
		}
	}
	else {																	// sum over all entries must be odd
		if ((fabs(tmp - sum) < prec) && ( fabs(fmod(tmp-1., 2.0)) < prec ) ) {
			return true;
		}
	}

	return false;
}

/* ########################################################################################
######   From_SO8C_Lattice() const                                                    ######
######                                                                               ######
######   Version: 16.05.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : the SO(8) lattice type: SO8_V, SO8_C, or NO_LATTICE         ######
###########################################################################################
######   description:                                                                ######
######   Checks whether this CLatticeVector is from the SO(8) spinorial or vectorial ######
######   weight lattice.                                                             ######
######################################################################################## */
LatticeType CLatticeVector::From_SO8C_Lattice() const
{
	if (this->Size != 4)
	{
		cout << "Warning in LatticeType CLatticeVector::From_SO8_Lattice() const : vector has length != 4. Return NO_LATTICE." << endl;
		return NO_LATTICE;
	}

	// begin: all entries integer or all entries half-integer
	double tmp = (*this)[0];
	const double first_entry = fabs(round_double_to_int(tmp) - tmp);

	const double prec = 0.0001;

	// assume INT type
	bool INTType = true;
	if (first_entry > prec)
	{
		if (fabs(first_entry - 0.5) < prec)
			INTType = false;
		else
			return NO_LATTICE;
	}

	unsigned i = 1;
	for (; i < 4; ++i)
	{
		tmp = (*this)[i];

		if (fabs(first_entry - fabs(round_double_to_int(tmp) - tmp)) > prec)
			return NO_LATTICE;
	}
	// end: all entries integer or all entries half-integer


	// sum over all entries
	double sum = 0.0;
	for (i = 0; i < 4; ++i) {
		sum +=(*this)[i];
	}

	tmp = round_double_to_int(sum);

	// sum over all entries must be an integer
	if (fabs(tmp - sum) > prec)
		return NO_LATTICE;

	tmp = fabs(fmod(tmp, 2.0));

	// odd integer
	if (fabs(tmp - 1.0) < prec) { //integral lattice  or cospinor lattice
		if (INTType) {
			return SO8_V;
		}
		else
		{
			return SO8_C;
		}
	}
	return NO_LATTICE;
}




/* ########################################################################################
######   From_SO8_Lattice() const                                                    ######
######                                                                               ######
######   Version: 16.05.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : the SO(8) lattice type: SO8_V, SO8_S, or NO_LATTICE         ######
###########################################################################################
######   description:                                                                ######
######   Checks whether this CLatticeVector is from the SO(8) spinorial or vectorial ######
######   weight lattice.                                                             ######
######################################################################################## */
LatticeType CLatticeVector::From_SO8S_Lattice() const
{
  if (this->Size != 4)
  {
    cout << "Warning in LatticeType CLatticeVector::From_SO8_Lattice() const : vector has length != 4. Return NO_LATTICE." << endl;
    return NO_LATTICE;
  }

  // begin: all entries integer or all entries half-integer
  double tmp = (*this)[0];
  const double first_entry = fabs(round_double_to_int(tmp) - tmp);

  const double prec = 0.0001;

  // assume INT type
  bool INTType = true;
  if (first_entry > prec)
  {
    if (fabs(first_entry - 0.5) < prec)
       INTType = false;
    else
      return NO_LATTICE;
  }

  unsigned i = 1;
  for (; i < 4; ++i)
  {
    tmp = (*this)[i];

    if (fabs(first_entry - fabs(round_double_to_int(tmp) - tmp)) > prec)
      return NO_LATTICE;
  }
  // end: all entries integer or all entries half-integer


  // sum over all entries
  double sum = 0.0;
  for (i = 0; i < 4; ++i)
    sum += (*this)[i];

  tmp = round_double_to_int(sum);

  // sum over all entries must be an integer
  if (fabs(tmp - sum) > prec)
    return NO_LATTICE;

  tmp = fabs(fmod(tmp, 2.0));

  // odd integer
  if (INTType)
  {
    if (fabs(tmp - 1.0) < prec) //integral lattice
      return SO8_V;
  }
  // even integer
  else
  {
    if (tmp < prec)             //half-integral lattice
      return SO8_S;
  }
  return NO_LATTICE;
}




/* ########################################################################################
######   From_Spin32_Lattice() const                                                 ######
######                                                                               ######
######   Version: 16.05.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : is the vector from the spin32/Z2 lattice?                    ######
###########################################################################################
######   description:                                                                ######
######   Checks whether this CLatticeVector is from the spin32/Z_2 weight lattice.   ######
######################################################################################## */
bool CLatticeVector::From_Spin32_Lattice() const
{
	if (this->Size != 16)
	{
		cout << "Warning in bool CLatticeVector::From_Spin32_Lattice() const : vector has length != 16. Return false." << endl;
		return false;
	}

	// begin: all entries integer or all entries half-integer
	double tmp = (*this)[0];
	const double first_entry = fabs(round_double_to_int(tmp) - tmp);

	const double prec = 0.0001;

	if ((first_entry > prec) && (fabs(first_entry - 0.5) > prec))
		return false;

	unsigned i = 1;
	for (; i < 16; ++i)
	{
		tmp = (*this)[i];

		if (fabs(first_entry - fabs(round_double_to_int(tmp) - tmp)) > prec)
			return false;
	}
	// end: all entries integer or all entries half-integer

	// sum over all entries
	double sum = 0.0;
	for (i = 0; i < 16; ++i)
		sum += (*this)[i];

	tmp = round_double_to_int(sum);

	// sum over all entries must be even
	if ((fabs(tmp - sum) < prec) && (fmod(tmp, 2.0) == 0))
		return true;

	return false;
}

