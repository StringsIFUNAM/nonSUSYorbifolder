
#include <iostream>
#include <iomanip>
#include <cstdlib>

#include "cvector.h"
#include "globalfunctions.h"

#define CHECKERROR true

using std::cout;
using std::endl;
using std::exit;
using std::setw;
using std::setprecision;


namespace H
{
  double tmp(0.0);
  unsigned i(0);
  const double prec(0.001);
}



/* ########################################################################################
######   CVector()                                                                   ######
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
######   Standard constructor of a CVector object. Creates a 0-dim. vector.          ######
######################################################################################## */
CVector::CVector()
 : Size(0)
{
}



/* ########################################################################################
######   CVector(const unsigned dim)                                                 ######
######                                                                               ######
######   Version: 07.12.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) dim : the dimension of the vector                                        ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Constructor of a CVector object. Creates a vector of dimension "dim".       ######
######################################################################################## */
CVector::CVector(const unsigned dim)
{
  this->assign(dim, 0.0);
  this->Size = dim;
}



/* ########################################################################################
######   ~CVector()                                                                  ######
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
######   Standard destructor of a CVector object.                                    ######
######################################################################################## */
CVector::~CVector()
{
}



/* ########################################################################################
######   operator=(const doubleVector &Vector)                                       ######
######                                                                               ######
######   Version: 21.01.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Vector : a vector of double                                              ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Specify the content of this CVector object by a vector of doubles.          ######
######################################################################################## */
void CVector::operator=(const doubleVector &Vector)
{
  this->clear();
  this->Size = Vector.size();

  this->insert(this->end(), Vector.begin(), Vector.end());
}



/* ########################################################################################
######   operator=(const intVector &Vector)                                          ######
######                                                                               ######
######   Version: 21.01.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Vector : a vector of integer                                             ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Specify the content of this CVector object by a vector of integers.         ######
######################################################################################## */
void CVector::operator=(const intVector &Vector)
{
  this->clear();
  this->Size = Vector.size();

  this->insert(this->end(), Vector.begin(), Vector.end());
}



/* ########################################################################################
######   operator=(const vector<rational<int> > &RationalVector)                     ######
######                                                                               ######
######   Version: 13.09.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) RationalVector : a vector of rational numbers                            ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Specify the content of this CVector object by a vector of rational numbers. ######
######################################################################################## */
void CVector::operator=(const vector<rational<int> > &RationalVector)
{
  this->Size = RationalVector.size();

  this->assign(this->Size, 0);

  for (H::i = 0; H::i < this->Size; ++H::i)
  {
    const rational<int> &tmp = RationalVector[H::i];
    this->at(H::i) = ((double)tmp.numerator())/((double)tmp.denominator());
  }
}



/* ########################################################################################
######   operator+(const CVector &Vector2) const                                     ######
######                                                                               ######
######   Version: 21.01.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Vector2   : a CVector object                                             ######
######   output:                                                                     ######
######   return value : the sum of the two vectors                                   ######
###########################################################################################
######   description:                                                                ######
######   Defines the standard addition of two vectors and returns the result.        ######
######################################################################################## */
CVector CVector::operator+(const CVector &Vector2) const
{
  CVector result(this->Size);

  #ifdef CHECKERROR
  if (this->Size != Vector2.Size)
  {
    cout << "\n  Warning in CVector CVector::operator+(...) const : the vectors have different lengths. Return empty vector." << endl;
    return result;
  }
  #endif

  for (H::i = 0; H::i < this->Size; ++H::i)
    result[H::i] = this->at(H::i) + Vector2[H::i];

  return result;
}



/* ########################################################################################
######   operator-(const CVector &Vector2) const                                     ######
######                                                                               ######
######   Version: 21.01.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Vector2   : a CVector object                                             ######
######   output:                                                                     ######
######   return value : the difference of the two vectors                            ######
###########################################################################################
######   description:                                                                ######
######   Defines the standard subtraction of two vectors and returns the result.     ######
######################################################################################## */
CVector CVector::operator-(const CVector &Vector2) const
{
  CVector result(this->Size);
  #ifdef CHECKERROR
  if (this->Size != Vector2.Size)
  {
    cout << "\n  Warning in CVector CVector::operator-(...) const : the vectors have different lengths. Return empty vector." << endl;
    return result;
  }
  #endif

  for (H::i = 0; H::i < this->Size; ++H::i)
    result[H::i] = this->at(H::i) - Vector2[H::i];

  return result;
}



/* ########################################################################################
######   operator*(double a) const                                                   ######
######                                                                               ######
######   Version: 25.06.2009                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) a         : a double number                                              ######
######   output:                                                                     ######
######   return value : the vector itself rescaled by the factor "a"                 ######
###########################################################################################
######   description:                                                                ######
######   Defines the standard multiplication of a vector and a double.               ######
######   Returns the result.                                                         ######
######################################################################################## */
CVector CVector::operator*(double a) const
{
  CVector result(this->Size);

  if (a != 0)
  {
    for (H::i = 0; H::i < this->Size; ++H::i)
      result[H::i] = this->at(H::i) * a;
  }

  return result;
}



/* ########################################################################################
######   operator*(const vector<CVector>::const_iterator &it_Vector2) const          ######
######                                                                               ######
######   Version: 21.01.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) it_Vector2 : a const_iterator to a CVector object                        ######
######   output:                                                                     ######
######   return value  : the scalar product of the two vectors                       ######
###########################################################################################
######   description:                                                                ######
######   Defines the standard scalar product of a vector and a iterator to a vector. ######
######   Returns the result.                                                         ######
######################################################################################## */
double CVector::operator*(const vector<CVector>::const_iterator &it_Vector2) const
{
  #ifdef CHECKERROR
  if (this->Size != it_Vector2->Size)
  {
    cout << "\n  Warning in double CVector::operator*(...) const : the vectors have different lengths. Return -1.0." << endl;
    return -1.0;
  }
  #endif

  H::tmp = 0.0;
  for (H::i = 0; H::i < this->Size; ++H::i)
    H::tmp += this->at(H::i) * it_Vector2->at(H::i);

  return H::tmp;
}



/* ########################################################################################
######   operator*(const CVector &Vector2) const                                     ######
######                                                                               ######
######   Version: 21.01.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Vector2   : a CVector object                                             ######
######   output:                                                                     ######
######   return value : the scalar product of the two vectors                        ######
###########################################################################################
######   description:                                                                ######
######   Defines the standard scalar product of two CVector objects.                 ######
######   Returns the result.                                                         ######
######################################################################################## */
double CVector::operator*(const CVector &Vector2) const
{
  #ifdef CHECKERROR
  if (this->Size != Vector2.Size)
  {
    cout << "\n  Warning in double CVector::operator*(...) const : the vectors have different lengths. Return -1.0." << endl;
    return -1.0;
  }
  #endif

  H::tmp = 0.0;
  for (H::i = 0; H::i < this->Size; ++H::i)
    H::tmp += this->at(H::i) * Vector2[H::i];

  return H::tmp;
}



/* ########################################################################################
######   operator+=(const CVector &Vector2)                                          ######
######                                                                               ######
######   Version: 21.01.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Vector2 : a CVector object                                               ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Defines the standard addition of two CVector objects. Stores the result in  ######
######   this CVector object.                                                        ######
######################################################################################## */
void CVector::operator+=(const CVector &Vector2)
{
  #ifdef CHECKERROR
  if (this->Size != Vector2.Size)
  {
    cout << "\n  Warning in void CVector::operator+=(...) : the vectors have different lengths. Return." << endl;
    return;
  }
  #endif

  for (H::i = 0; H::i < this->Size; ++H::i)
    this->at(H::i) += Vector2[H::i];
}



/* ########################################################################################
######   operator-=(const CVector &Vector2)                                          ######
######                                                                               ######
######   Version: 02.02.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Vector2 : a CVector object                                               ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Defines the standard subtraction of two CVector objects. Stores the result  ######
######   in this CVector object.                                                     ######
######################################################################################## */
void CVector::operator-=(const CVector &Vector2)
{
  #ifdef CHECKERROR
  if (this->Size != Vector2.Size)
  {
    cout << "\n  Warning in void CVector::operator-=(...) : the vectors have different lengths. Return." << endl;
    return;
  }
  #endif

  for (H::i = 0; H::i < this->Size; ++H::i)
    this->at(H::i) -= Vector2[H::i];
}



/* ########################################################################################
######   operator==(const CVector &Vector2) const                                    ######
######                                                                               ######
######   Version: 21.01.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Vector2   : a CVector object                                             ######
######   output:                                                                     ######
######   return value : are the vectors equal?                                       ######
###########################################################################################
######   description:                                                                ######
######   Compares two CVector objects. Returns "true" if they are equal and "false"  ######
######   if they are not.                                                            ######
######################################################################################## */
bool CVector::operator==(const CVector &Vector2) const
{
  #ifdef CHECKERROR
  if (this->Size != Vector2.Size)
  {
    cout << "\n  Warning in bool CVector::operator==(...) const : the vectors have different lengths. Return false." << endl;
    return false;
  }
  #endif

  for (H::i = 0; H::i < this->Size; ++H::i)
  {
    if (fabs(this->at(H::i) - Vector2[H::i]) > H::prec)
      return false;
  }

  return true;
}



/* ########################################################################################
######   operator!=(const CVector &Vector2) const                                    ######
######                                                                               ######
######   Version: 21.01.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Vector2   : a CVector object                                             ######
######   output:                                                                     ######
######   return value : are the vectors unequal?                                     ######
###########################################################################################
######   description:                                                                ######
######   Compares two CVector objects. Returns "true" if they are not equal and      ######
######   "false" if they are.                                                        ######
######################################################################################## */
bool CVector::operator!=(const CVector &Vector2) const
{
  #ifdef CHECKERROR
  if (this->Size != Vector2.Size)
  {
    cout << "\n  Warning in bool CVector::operator!=(...) const : the vectors have different lengths. Return false." << endl;
    return false;
  }
  #endif
    
  for (H::i = 0; H::i < this->Size; ++H::i)
  {
    if (fabs(this->at(H::i) - Vector2[H::i]) > H::prec)
      return true;
  }

  return false;
}



/* ########################################################################################
######   Assign(const unsigned dim)                                                  ######
######                                                                               ######
######   Version: 02.02.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) dim : dimension of the vector                                            ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Changes this CVector object to a null vector of dimension "dim".            ######
######################################################################################## */
void CVector::Assign(const unsigned dim)
{
  this->assign(dim, 0);
  this->Size = dim;
}



/* ########################################################################################
######   Assign(const unsigned dim, const double value)                              ######
######                                                                               ######
######   Version: 02.02.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) dim   : dimension of the vector                                          ######
######   2) value : assigns "dim" copies of "value" to the vector                    ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Changes this CVector object to a vector of dimension "dim" and sets the     ######
######   values of all components to "value".                                        ######
######################################################################################## */
void CVector::Assign(const unsigned dim, const double value)
{
  this->assign(dim, value);
  this->Size = dim;
}



/* ########################################################################################
######   CreateRandomVector(...)                                                     ######
######                                                                               ######
######   Version: 21.01.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Lattice        : E8xE8 or Spin32                                         ######
######   2) SymmetryBlocks : the symmetry structure w.r.t. permutations of the       ######
######                       vector entries                                          ######
######   3) Part1_Blocks   : for E8xE8, pieces of various lengths that might be      ######
######                       assigned to the first part of the vector (first E8)     ######
######   4) Part2_Blocks   : as above, for second part (second E8)                   ######
######   output:                                                                     ######
######   return value      : random vector created successfully?                     ######
###########################################################################################
######   description:                                                                ######
######   Creates a random vector of lattice type "Lattice" respecting the symmetry   ###### 
######   structure (as specified in "SymmetryBlocks") and using the contents of      ######
######   "Part1_Blocks" and "Part2_Blocks". Returns true if the vector has been      ######
######   created successfully. This routine is used by "CreateRandom(...)" in the    ###### 
######   class CShiftVector and by "CreateRandom(...) in the class "CWilsonLine".    ######
######################################################################################## */
bool CVector::CreateRandomVector(SelfDualLattice Lattice, const vector<unsigned> &SymmetryBlocks, const vector<vector<vector<double> > > *Part1_Blocks, const vector<vector<vector<double> > > *Part2_Blocks)
{
  double max = RAND_MAX+0.000001;

  unsigned length = 0;
  unsigned total_length = 0;

  this->Size = 16;
  this->clear();

  const vector<vector<vector<double> > > *Blocks = NULL;
  // begin: fill Vector with pieces from Blocks respecting the SymmetryBlocks
  const size_t s1 = SymmetryBlocks.size();
  for (unsigned i = 0; i < s1; ++i)
  {
    length = SymmetryBlocks[i];
    total_length += length;

    if (Lattice == E8xE8)
    {
      if (total_length < 9)
        Blocks = Part1_Blocks;
      else
        Blocks = Part2_Blocks;
    }
    else
      Blocks = Part1_Blocks;

    #ifdef CHECKERROR
    if (Blocks->size() < length)
    {
      cout << "\n  Warning in bool CVector::CreateRandomVector(...) : \"Blocks\" does not contain the symmetry length " << length << ". Return false." << endl;
      this->assign(16, 0.0);
      return false;
    }
    #endif
      
    const vector<vector<double> > &CurrentBlock = Blocks->at(length-1);
    const vector<double>          &PartOfVector = CurrentBlock[(unsigned)(rand() * (double)CurrentBlock.size()/max)];
    this->insert(this->end(), PartOfVector.begin(), PartOfVector.end());
  }
  // end: fill Vector with pieces from Blocks respecting the SymmetryBlocks

  return true;
}



/* ########################################################################################
######   GetLength() const                                                           ######
######                                                                               ######
######   Version: 02.05.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : the length of the vector                                     ######
###########################################################################################
######   description:                                                                ######
######   Computes the length of a vector using the standard definition and returns   ######
######   the result.                                                                 ######
######################################################################################## */
double CVector::GetLength() const
{
  if (this->Size == 1)
    return fabs(this->at(0));

  H::tmp = 0.0;
  for (H::i = 0; H::i < this->Size; ++H::i)
    H::tmp += this->at(H::i) * this->at(H::i);

  return sqrt(H::tmp);
}



/* ########################################################################################
######   GetSqrTo(unsigned cut) const                                                ######
######                                                                               ######
######   Version: 21.01.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) cut       : vector index where to crop the vector, from 0,..,this->Size  ######
######   output:                                                                     ######
######   return value : the length-square of the cropped vector at index "cut"       ######
###########################################################################################
######   description:                                                                ######
######   Computes the squared length of this CVector object from the first component ######
######   to the position "cut" and returns the result, i.e.                          ######
######   result = \sum_{i=0}^{cut} v_i^2                                             ######
######################################################################################## */
double CVector::GetSqrTo(unsigned cut) const
{
  #ifdef CHECKERROR
  if (cut > this->Size)
  {
    cout << "\n  Warning in double CVector::GetSqrTo(...) const : variable \"cut\" out of range. Return -1.0." << endl;
    return -1.0;
  }
  #endif

  H::tmp = 0.0;
  for (H::i = 0; H::i < cut; ++H::i)
    H::tmp += this->at(H::i) * this->at(H::i);

  return H::tmp;
}



/* ########################################################################################
######   IsScalarProductInteger(const CVector &Vector2) const                        ######
######                                                                               ######
######   Version: 16.05.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Vector2   : a CVector object                                             ######
######   output:                                                                     ######
######   return value : is the scalar product of the two vectors integer?            ######
###########################################################################################
######   description:                                                                ######
######   Checks whether the scalar product of the two vectors is integer.            ######
######################################################################################## */
bool CVector::IsScalarProductInteger(const CVector &Vector2) const
{
  #ifdef CHECKERROR
  if (this->Size != Vector2.Size)
  {
    cout << "\n  Warning in bool CVector::IsScalarProductInteger(...) const : the vectors have different lengths. Return false." << endl;
    return false;
  }
  #endif

  H::tmp = 0.0;
  for (H::i = 0; H::i < this->Size; ++H::i)
    H::tmp += this->at(H::i) * Vector2[H::i];

  return is_integer( H::tmp);
}



/* ########################################################################################
######   Push_back(double value)                                                     ######
######                                                                               ######
######   Version: 25.06.2009                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) value : adds a new component with value "value" to the vector            ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Adds a component with value "value" to this CVector object.                 ######
######################################################################################## */
void CVector::Push_back(double value)
{
  this->push_back(value);
  this->Size = this->size();
}



/* ########################################################################################
######   IsZero() const                                                              ######
######                                                                               ######
######   Version: 18.04.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : are all components of the vector zero?                       ######
###########################################################################################
######   description:                                                                ######
######   Checks whether this CVector object is a null vector or not.                 ######
######################################################################################## */
bool CVector::IsZero() const
{
  for (H::i = 0; H::i < this->Size; ++H::i)
  {
    if (fabs(this->at(H::i)) > H::prec)
      return false;
  }
  return true;
}



/* ########################################################################################
######   SetToZero()                                                                 ######
######                                                                               ######
######   Version: 25.02.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Set this CVector object to the null vector.                                 ######
######################################################################################## */
void CVector::SetToZero()
{
  this->assign(this->Size, 0.0);
}


