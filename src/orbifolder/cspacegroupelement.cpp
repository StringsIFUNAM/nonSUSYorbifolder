#include <iostream>
#include <stdexcept>

#include "cspacegroupelement.h"

using std::cout;
using std::endl;
using std::out_of_range;


/* ########################################################################################
######   CLatticeElement()                                                           ######
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
######   Standard constructor of a CLatticeElement object. Creates a 6-dim. null     ######
######   vector.                                                                     ######
######################################################################################## */
CLatticeElement::CLatticeElement()
{
  this->assign(LatticeDim, 0);
}



/* ########################################################################################
######   ~CLatticeElement()                                                          ######
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
######   Standard destructor of a CLatticeElement object.                            ######
######################################################################################## */
CLatticeElement::~CLatticeElement()
{
}



/* ########################################################################################
######   CLatticeElement(const rationalVector &n_alpha)                              ######
######                                                                               ######
######   Version: 11.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Constructor of a CLatticeElement object. Creates a 6-dim. vector with       ###### 
######   entries "n_alpha".                                                          ######
######################################################################################## */
CLatticeElement::CLatticeElement(const rationalVector &n_alpha)
{
  if (n_alpha.size() != LatticeDim)
  {
    cout << "Warning in CLatticeElement::CLatticeElement(...) : n_alpha does not contain " << LatticeDim << " components. Hence, set to zero." << endl;
    this->assign(LatticeDim, 0);
    return;
  }
  for (unsigned i = 0; i < LatticeDim; ++i)
    this->at(i) = n_alpha[i];
}



/* ########################################################################################
######   operator=(const rationalVector &n_alpha)                                    ######
######                                                                               ######
######   Version: 11.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) n_alpha : a six-dim vector of rational numbers                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Sets the content of this CLatticeElement object to "n_alpha".               ######
######################################################################################## */
void CLatticeElement::operator=(const rationalVector &n_alpha)
{
  if (n_alpha.size() != LatticeDim)
  {
    cout << "Warning in void CLatticeElement::operator=(...) : n_alpha does not contain " << LatticeDim << " components. Hence, set to zero." << endl;
    this->assign(LatticeDim, 0);
    return;
  }
  for (unsigned i = 0; i < LatticeDim; ++i)
    this->at(i) = n_alpha[i];
}



/* ########################################################################################
######   operator==(const CLatticeElement &Element2) const                           ######
######                                                                               ######
######   Version: 05.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Element2  : a CLatticeElement object                                     ######
######   output:                                                                     ######
######   return value : are the two CLatticeElement objects equal?                   ######
###########################################################################################
######   description:                                                                ######
######   Compares two CLatticeElement objects. Returns "true" if they are equal and  ######
######   "false" if they are not.                                                    ######
######################################################################################## */
bool CLatticeElement::operator==(const CLatticeElement &Element2) const
{
  for (unsigned i = 0; i < LatticeDim; ++i)
  {
    if (this->at(i) != Element2[i])
      return false;
  }

  return true;
}



/* ########################################################################################
######   operator!=(const CLatticeElement &Element2) const                           ######
######                                                                               ######
######   Version: 05.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Element2  : a CLatticeElement object                                     ######
######   output:                                                                     ######
######   return value : are the two CLatticeElement objects not equal?               ######
###########################################################################################
######   description:                                                                ######
######   Compares two CLatticeElement objects. Returns "false" if they are equal     ######
######   and "true" if they are not.                                                 ######
######################################################################################## */
bool CLatticeElement::operator!=(const CLatticeElement &Element2) const
{
  for (unsigned i = 0; i < LatticeDim; ++i)
  {
    if (this->at(i) != Element2[i])
      return true;
  }

  return false;
}



/* ########################################################################################
######   operator+=(const CLatticeElement &b)                                        ######
######                                                                               ######
######   Version: 05.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) b : the CLatticeElement object to be added                               ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Defines the standard addition of two CLatticeElement objects. Stores the    ######
######   result in this CLatticeElement object.                                      ######
######################################################################################## */
void CLatticeElement::operator+=(const CLatticeElement &b)
{
  for (unsigned i = 0; i < LatticeDim; ++i)
    this->at(i) += b[i];
}



/* ########################################################################################
######   operator-=(const CLatticeElement &b)                                        ######
######                                                                               ######
######   Version: 05.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) b : the CLatticeElement object to be added                               ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Defines the standard subtraction of two CLatticeElement objects. Stores the ######
######   result in this CLatticeElement object.                                      ######
######################################################################################## */
void CLatticeElement::operator-=(const CLatticeElement &b)
{
  for (unsigned i = 0; i < LatticeDim; ++i)
    this->at(i) -= b[i];
}



/* ########################################################################################
######   operator*=(const int &b)                                                    ######
######                                                                               ######
######   Version: 05.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) b : an integer scale factor                                              ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Defines the standard multiplication of a CLatticeElement object with an     ######
######   integer. Stores the result in this CLatticeElement object.                  ######
######################################################################################## */
void CLatticeElement::operator*=(const int &b)
{
  for (unsigned i = 0; i < LatticeDim; ++i)
    this->at(i) *= b;
}



/* ########################################################################################
######   operator*=(const rational<int> &b)                                          ######
######                                                                               ######
######   Version: 05.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) b : a rational scale factor                                              ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Defines the standard multiplication of a CLatticeElement object with a      ######
######   rational number. Stores the result in this CLatticeElement object.          ######
######################################################################################## */
void CLatticeElement::operator*=(const rational<int> &b)
{
  for (unsigned i = 0; i < LatticeDim; ++i)
    this->at(i) *= b;
}



/* ########################################################################################
######   operator*(const rationalVector &ChargeOperator) const                       ######
######                                                                               ######
######   Version: 05.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) ChargeOperator : a vector of six rational numbers                        ######
######   output:                                                                     ######
######   return value      : the scalar product of this CLatticeElement object and   ######
######                       the rational vector "ChargeOperator"                    ######
###########################################################################################
######   description:                                                                ######
######   Defines the standard scalar product. Used to compute the charges with       ######
######   respect to a discrete symmetry.                                             ######
######################################################################################## */
rational<int> CLatticeElement::operator*(const rationalVector &ChargeOperator) const
{
  if (ChargeOperator.size() != (LatticeDim + 2))
  {
    cout << "Warning in rational<int> CLatticeElement::operator*(...) const : charge operator of the discrete symmetry is ill-defined. Return charge 0." << endl;
    return 0;
  }

  rational<int> DiscreteCharge = 0;

  for (unsigned i = 0; i < LatticeDim; ++i)
    DiscreteCharge += (this->at(i) * ChargeOperator[i + 2]);

  return DiscreteCharge;
}



/* ########################################################################################
######   operator*(const int &b) const                                               ######
######                                                                               ######
######   Version: 04.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) b         : an integer scale factor                                      ######
######   output:                                                                     ######
######   return value : the rescaled vector                                          ######
###########################################################################################
######   description:                                                                ######
######   Defines the standard multiplication of a CLatticeElement object with an     ######
######   integer. Returns the result.                                                ######
######################################################################################## */
CLatticeElement CLatticeElement::operator*(const int &b) const
{
  CLatticeElement result;
  for (unsigned i = 0; i < LatticeDim; ++i)
    result.at(i) = this->at(i) * b;
  return result;
}



/* ########################################################################################
######   &operator[](const unsigned &alpha) const                                    ######
######                                                                               ######
######   Version: 05.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) alpha     : index from 0 to 5                                            ######
######   output:                                                                     ######
######   return value : the content of n_alpha                                       ######
###########################################################################################
######   description:                                                                ######
######   Defines the operator [] to gain access to the alpha-th component.           ######
######################################################################################## */
const rational<int> &CLatticeElement::operator[](const unsigned &alpha) const
{
  try {
    return this->at(alpha);
  }
  catch (out_of_range& oor) {
    cout << "Out of Range error: " << oor.what() << endl;
  }
  return this->at(alpha);
}



/* ########################################################################################
######   IsZero() const                                                              ######
######                                                                               ######
######   Version: 13.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : is this vector empty?                                        ######
###########################################################################################
######   description:                                                                ######
######   Checks whether this vector is empty. Returns "true" if it is empty and      ######
######   "false" otherwise.                                                          ######
######################################################################################## */
bool CLatticeElement::IsZero() const
{
  for (unsigned i = 0; i < LatticeDim; ++i)
  {
    if (this->at(i) != 0)
      return false;
  }
  return true;
}



/* ########################################################################################
######   Rotate(const rationalMatrix &TwistMatrix, CLatticeElement &result) const    ######
######                                                                               ######
######   Version: 05.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) TwistMatrix : the rotation matrix                                        ######
######   output:                                                                     ######
######   1) result      : the rotated vector                                         ######
######   return value   : finished successfully?                                     ######
###########################################################################################
######   description:                                                                ######
######   Rotates this CLatticeElement object using the rotation matrix "TwistMatrix" ######
######   and save the result in "result".                                            ######
######################################################################################## */
bool CLatticeElement::Rotate(const rationalMatrix &TwistMatrix, CLatticeElement &result) const
{
  bool NoError = true;
  if (TwistMatrix.size() != LatticeDim)
    NoError = false;

  unsigned j = 0;
  for (j = 0; NoError && (j < LatticeDim); ++j)
  {
    if (TwistMatrix[j].size() != LatticeDim)
      NoError = false;
  }

  if (!NoError)
  {
    cout << "Warning in bool CLatticeElement::Rotate(...) const: TwistMatrix ill-defined. Return false." << endl;
    return false;
  }

  for (unsigned i = 0; i < LatticeDim; ++i)
  {
    rational<int> &n_i = result.at(i);
    n_i = 0;
    const rationalVector &TwistVector = TwistMatrix[i];
    
    for (j = 0; j < LatticeDim; ++j)
      n_i += this->at(j) * TwistVector[j];
  }
  return true;
}



/* ########################################################################################
######   Set_n_alpha(const unsigned &alpha, const rational<int> &n_alpha)            ######
######                                                                               ######
######   Version: 05.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) alpha     : index from 0 to 5                                            ######
######   2) n_alpha   : set the alpha-th component of this CLatticeElement vector to ######
######                  "n_alpha"                                                    ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Sets the alpha-th component of this CLatticeElement vector to "n_alpha".    ######
######################################################################################## */
bool CLatticeElement::Set_n_alpha(const unsigned &alpha, const rational<int> &n_alpha)
{
  if (alpha >= LatticeDim)
  {
    cout << "Warning in bool CLatticeElement::Set_n_alpha(...) : alpha out of range. Return false." << endl;
    return false;
  }
  this->at(alpha) = n_alpha;
  return true;
}







/* ########################################################################################
######   ~CSpaceGroupElement()                                                       ######
######                                                                               ######
######   Version: 11.12.2006                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard destructor of a CSpaceGroupElement object.                         ######
######################################################################################## */
CSpaceGroupElement::~CSpaceGroupElement()
{
}



/* ########################################################################################
######   CSpaceGroupElement()                                                        ######
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
######   Standard constructor of a CSpaceGroupElement object. Creates the identity.  ######
######################################################################################## */
CSpaceGroupElement::CSpaceGroupElement()
{
	this->m = 0;
	this->n = 0;
	this->k = 0;
}



/* ########################################################################################
######   CSpaceGroupElement(unsigned m, unsigned n, ...)                             ######
######                                                                               ######
######   Version: 04.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) k       : labels the twisted sector                                      ######
######   2) l       : labels the twisted sector                                      ######
######   3) n_alpha : a vector of six rational numbers                               ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Constructor of a CSpaceGroupElement object. Creates the element             ######
######     (\theta^k \omega^l, n_\alpha e_\alpha).                                   ######
######################################################################################## */
CSpaceGroupElement::CSpaceGroupElement(unsigned m, unsigned n, unsigned k, const rationalVector &n_alpha)
{
    this->m = m;
	this->n = n;
	this->k = k;
  this->LatticeElement = n_alpha;
}



/* ########################################################################################
######   CSpaceGroupElement(unsigned k, unsigned l)                                  ######
######                                                                               ######
######   Version: 19.01.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) k       : labels the twisted sector                                      ######
######   2) l       : labels the twisted sector                                      ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Constructor of a CSpaceGroupElement object. Creates the element             ######
######     (\theta^k \omega^l, 0).                                                   ######
######################################################################################## */
CSpaceGroupElement::CSpaceGroupElement(unsigned m, unsigned n, unsigned k)
{
    this->m = m;
	this->n = n;
	this->k = k;
}





/* ########################################################################################
######   operator==(const CSpaceGroupElement &Element2) const                        ######
######                                                                               ######
######   Version: 04.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Element2  : a CSpaceGroupElement object                                  ######
######   output:                                                                     ######
######   return value : are the two CSpaceGroupElement objects equal?                ######
###########################################################################################
######   description:                                                                ######
######   Compares two CSpaceGroupElement objects. Returns "true" if they are equal   ######
######   and "false" if they are not.                                                ######
######################################################################################## */
bool CSpaceGroupElement::operator==(const CSpaceGroupElement &Element2) const
{
  if (this->m != Element2.m)
    return false;

  if (this->n != Element2.n)
    return false;

  if (this->k != Element2.k)
    return false;

  if (this->LatticeElement != Element2.LatticeElement)
    return false;

  return true;
}



/* ########################################################################################
######   operator!=(const CSpaceGroupElement &Element2) const                        ######
######                                                                               ######
######   Version: 04.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Element2  : a CSpaceGroupElement object                                  ######
######   output:                                                                     ######
######   return value : are the two CSpaceGroupElement objects not equal?            ######
###########################################################################################
######   description:                                                                ######
######   Compares two CSpaceGroupElement objects. Returns "false" if they are equal  ######
######   and "true" if they are not.                                                 ######
######################################################################################## */
bool CSpaceGroupElement::operator!=(const CSpaceGroupElement &Element2) const
{
	  if (this->m != Element2.m)
	    return true;

	  if (this->n != Element2.n)
	    return true;

	  if (this->k != Element2.k)
	    return true;

  if (this->LatticeElement != Element2.LatticeElement)
    return true;

  return false;
}



/* ########################################################################################
######   operator*(const rationalVector &ChargeOperator) const                       ######
######                                                                               ######
######   Version: 04.08.2011                                                         ######
######   Check-Level: 2                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) ChargeOperator : a vector of eight rational numbers                      ######
######   output:                                                                     ######
######   return value      : the scalar product of this CSpaceGroupElement object    ######
######                       and the rational vector "ChargeOperator"                ######
###########################################################################################
######   description:                                                                ######
######   Defines the standard scalar product. Used to compute the charges with       ######
######   respect to a discrete symmetry, i.e. kc_1 + lc_2 + n_alpha c_{2+alpha}.     ######
######################################################################################## */
rational<int> CSpaceGroupElement::operator*(const rationalVector &ChargeOperator) const
{
  if (ChargeOperator.size() != 8)
  {
    cout << "Warning in rational<int> CSpaceGroupElement::operator*(...) const : charge operator of the discrete symmetry is ill-defined. Return charge 0." << endl;
    return 0;
  }

  return ((int(this->n) * ChargeOperator[0]) + (int(this->k) * ChargeOperator[1]) + (this->LatticeElement * ChargeOperator)); //return ((this->n * ChargeOperator[0]) + (this->k * ChargeOperator[1]) + (this->LatticeElement * ChargeOperator));
}



/* ########################################################################################
######   IsZero() const                                                              ######
######                                                                               ######
######   Version: 13.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : is this space group element the identity?                    ######
###########################################################################################
######   description:                                                                ######
######   Checks whether this space group element is the identity. Returns "true" if  ######
######   it is empty and "false" otherwise.                                          ######
######################################################################################## */
bool CSpaceGroupElement::IsZero() const
{
  return ((this->m == 0) && (this->n == 0) && (this->k == 0) && this->LatticeElement.IsZero());
}



/* ########################################################################################
######   NoTwist() const                                                             ######
######                                                                               ######
######   Version: 24.01.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : is the twist zero?                                           ######
###########################################################################################
######   description:                                                                ######
######   Checks weather the rotation of this CSpaceGroupElement object is trivial.    ######
######################################################################################## */
bool CSpaceGroupElement::NoTwist() const
{
  if ((this->m == 0) && (this->n == 0) && (this->k == 0))
    return true;
  
  return false;
}



/* ########################################################################################
######   Set_m(const unsigned &m)                                                    ######
######                                                                               ######
######   Version: 04.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) k         : set "k" of this CSpaceGroupElement object to "k"             ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Sets the twisted sector "k" of this CSpaceGroupElement object to "k"        ######
######################################################################################## */
bool CSpaceGroupElement::Set_m(const unsigned &m)
{
  this->m = m;
  return true;
}



/* ########################################################################################
######   Set_n(const unsigned &n)                                                    ######
######                                                                               ######
######   Version: 04.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) l         : set "l" of this CSpaceGroupElement object to "l"             ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Sets the twisted sector "l" of this CSpaceGroupElement object to "l"        ######
######################################################################################## */
bool CSpaceGroupElement::Set_n(const unsigned &n)
{
  this->n = n;
  return true;
}



/* ########################################################################################
######   Set_k(const unsigned &k)                                                    ######
######                                                                               ######
######   Version: 04.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) l         : set "l" of this CSpaceGroupElement object to "l"             ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Sets the twisted sector "l" of this CSpaceGroupElement object to "l"        ######
######################################################################################## */
bool CSpaceGroupElement::Set_k(const unsigned &k)
{
  this->k = k;
  return true;
}



/* ########################################################################################
######   SetLatticeElement(const CLatticeElement &LatticeElement)                    ######
######                                                                               ######
######   Version: 04.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) LatticeElement : a six-dim. lattice vector                               ######
######   output:                                                                     ######
######   return value      : finished successfully?                                  ######
###########################################################################################
######   description:                                                                ######
######   Sets the lattice vector n_alpha (stored in the private member variable      ######
######   "LatticeElement") to "LatticeElement".                                      ######
######################################################################################## */
bool CSpaceGroupElement::SetLatticeElement(const CLatticeElement &LatticeElement)
{
  this->LatticeElement = LatticeElement;
  return true;
}



/* ########################################################################################
######   Set_n_alpha(const rationalVector &n_alpha)                                  ######
######                                                                               ######
######   Version: 04.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   3) n_alpha   : a vector of six rational numbers                             ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Sets the lattice vector n_alpha (stored in the private member variable      ######
######   "LatticeElement") to "n_alpha".                                             ######
######################################################################################## */
bool CSpaceGroupElement::Set_n_alpha(const rationalVector &n_alpha)
{
  this->LatticeElement = n_alpha;
  return true;
}



/* ########################################################################################
######   Set_n_alpha(const unsigned &alpha, const rational<int> &n_alpha)            ######
######                                                                               ######
######   Version: 04.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) alpha     : index of the component to be set                             ######
######   1) n_alpha   : new value of component alpha                                 ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Sets the alpha-th component of this CLatticeElement vector to "n_alpha".    ######
######################################################################################## */
bool CSpaceGroupElement::Set_n_alpha(const unsigned &alpha, const rational<int> &n_alpha)
{
  return this->LatticeElement.Set_n_alpha(alpha, n_alpha);
}
