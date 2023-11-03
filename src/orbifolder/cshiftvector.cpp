#include <iostream>
#include <cstdlib>

#include "cshiftvector.h"
#include "clatticevector.h"
#include "globalfunctions.h"

#define CHECKERROR true

using std::cout;
using std::endl;



/* ########################################################################################
######   CShiftVector()                                                              ######
######                                                                               ######
######   Version: 19.03.2010                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard constructor of a CShiftVector object. Creates a 16-dim. null       ######
######   vector.                                                                     ######
######################################################################################## */
CShiftVector::CShiftVector()
{
  this->SetToZero();
  this->Lattice = UNSPECIFIED_LATTICE;
}



/* ########################################################################################
######   CShiftVector(SelfDualLattice Lattice)                                       ######
######                                                                               ######
######   Version: 19.03.2010                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Lattice : E8xE8 or Spin32                                                ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Constructor of a CShiftVector object of lattice type "Lattice". Creates a   ###### 
######   16-dim. null vector for either the E8xE8 or the Spin32 gauge lattice.       ######
######################################################################################## */
CShiftVector::CShiftVector(SelfDualLattice Lattice)
{
  if (Lattice == UNSPECIFIED_LATTICE)
  {
    cout << "\n  Warning in CShiftVector::CShiftVector(...): The even and self-dual lattice was not defined. Hence, it is set to E8xE8." << endl;
    this->Lattice = E8xE8;
  }
  else
    this->Lattice = Lattice;

  this->SetToZero();
}



/* ########################################################################################
######   CShiftVector(SelfDualLattice Lattice, ...)                                  ######
######                                                                               ######
######   Version: 19.03.2010                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Lattice   : E8xE8 or Spin32                                              ######
######   2) - 17) V_i : the 16 components of the shift                               ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Constructor of a CShiftVector object of lattice type "Lattice". Creates a   ###### 
######   16-dim. vector for either the E8xE8 or the Spin32 gauge lattice with        ######
######   entries "V_i" for i = 0,..,15.                                              ######
######################################################################################## */
CShiftVector::CShiftVector(SelfDualLattice Lattice, double V_0, double V_1, double  V_2, double  V_3, double  V_4, double  V_5, double  V_6, double  V_7,
                           double V_8, double V_9, double V_10, double V_11, double V_12, double V_13, double V_14, double V_15)
{
  if (Lattice == UNSPECIFIED_LATTICE)
  {
    cout << "\n  Warning in CShiftVector::CShiftVector(...): The even and self-dual lattice was not defined. Hence, it is set to E8xE8." << endl;
    this->Lattice = E8xE8;
  }
  else
    this->Lattice = Lattice;

  this->resize(16,0);
  this->Size = 16;

  (*this)[0]  = V_0;
  (*this)[1]  = V_1;
  (*this)[2]  = V_2;
  (*this)[3]  = V_3;
  (*this)[4]  = V_4;
  (*this)[5]  = V_5;
  (*this)[6]  = V_6;
  (*this)[7]  = V_7;
  (*this)[8]  = V_8;
  (*this)[9]  = V_9;
  (*this)[10] = V_10;
  (*this)[11] = V_11;
  (*this)[12] = V_12;
  (*this)[13] = V_13;
  (*this)[14] = V_14;
  (*this)[15] = V_15;


  if (!this->UpdateData())
  {
    cout << "\n  Warning in CShiftVector::CShiftVector(...): Shift vector set to zero." << endl;
    this->SetToZero();
    return;
  }
}



/* ########################################################################################
######   ~CShiftVector()                                                             ######
######                                                                               ######
######   Version: 19.03.2010                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard destructor of a CShiftVector object.                               ######
######################################################################################## */
CShiftVector::~CShiftVector()
{
}



/* ########################################################################################
######   operator=(const CVector &Vector)                                            ######
######                                                                               ######
######   Version: 19.03.2010                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Vector : a CVector obejct of dimension 16 and Order = 1, ..., 24         ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Specify the content of this CShiftVector object by a CVector object.        ######
######################################################################################## */
void CShiftVector::operator=(const CVector &Vector)
{
  #ifdef CHECKERROR
  if (Vector.Size != 16)
  {
    cout << "\n  Warning in void CShiftVector::operator=(...): Shift vector is not a 16-dim. vector. Hence shift vector is set to zero." << endl;
    this->SetToZero();
    return;
  }
  #endif

  this->clear();
  this->insert(this->end(), Vector.begin(), Vector.end());

  if (!this->UpdateData())
  {
    cout << "\n  Warning in void CShiftVector::operator=(...): Shift vector set to zero." << endl;
    this->SetToZero();
    return;
  }
}



/* ########################################################################################
######   operator=(const vector<rational<int> > &RationalVector)                     ######
######                                                                               ######
######   Version: 19.03.2010                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) RationalVector : a vector of rantional<int> with 16 entries              ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Specify the content of this CShiftVector object by a vector of rational     ######
######   numbers.                                                                    ######
######################################################################################## */
void CShiftVector::operator=(const vector<rational<int> > &RationalVector)
{
  #ifdef CHECKERROR
  if (RationalVector.size() != 16)
  {
    cout << "\n  Warning in void CShiftVector::operator=(...): Shift vector is not a 16-dim. vector. Hence shift vector is set to zero." << endl;
    this->SetToZero();
    return;
  }
  #endif

  for (unsigned i = 0; i < 16; ++i)
  {
    const rational<int> &tmp = RationalVector[i];
    (*this)[i] = ((double)tmp.numerator())/((double)tmp.denominator());
  }

  if (!this->UpdateData())
  {
    cout << "\n  Warning in void CShiftVector::operator=(...): Shift vector set to zero." << endl;
    this->SetToZero();
    return;
  }
}



/* ########################################################################################
######   CreateRandom(unsigned TwistOrder, double TwistSqr, ...)                     ######
######                                                                               ######
######   Version: 16.05.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) TwistOrder      : order N of the twist vector v                          ######
######   2) TwistSqr        : the twist vector squared: v^2                          ######
######   3) SymmetryBlocks  : the symmetry structure w.r.t. permutations of the      ######
######                        vector entries                                         ######
######   4) VectorialBlocks : pieces of vectorial lattice vectors of various lengths ######
######                        that might be assigned to the shift                    ######
######   5) SpinorialBlocks : pieces of spinorial lattice vectors of various lengths ######
######                        that might be assigned to the shift                    ######
######   output:                                                                     ######
######   return value       : finished successfully?                                 ######
###########################################################################################
######   description:                                                                ######
######   Creates a random shift vector of order "TwistOrder" respecting the symmetry ###### 
######   structure (as specified in "SymmetryBlocks") and using the contents of      ######
######   "VectorialBlocks" and "SpinorialBlocks". Ensures that the modular           ######
######   invariance condition                                                        ######
######     N (V^2 - v^2) = 0 mod 2                                                   ######
######   is fulfilled. Returns true if the shift vector has been created             ###### 
######   successfully.                                                               ######
######################################################################################## */
bool CShiftVector::CreateRandom(unsigned TwistOrder, double TwistSqr, const vector<unsigned> &SymmetryBlocks, const vector<vector<vector<double> > > &VectorialBlocks, const vector<vector<vector<double> > > &SpinorialBlocks)
{
  unsigned i = 0;

  double sqr = 0.0;
  double sum = 0.0;

  const vector<vector<vector<double> > > *Part1_Blocks = NULL;
  const vector<vector<vector<double> > > *Part2_Blocks = NULL;
  if (this->Lattice == E8xE8)
  {
    if (rand()/(double)RAND_MAX < 0.5)
      Part1_Blocks = &VectorialBlocks;
    else
      Part1_Blocks = &SpinorialBlocks;

    if (rand()/(double)RAND_MAX < 0.5)
      Part2_Blocks = &VectorialBlocks;
    else
      Part2_Blocks = &SpinorialBlocks;
  }
  else
  {
    if (rand()/(double)RAND_MAX < 0.5)
    {
      Part1_Blocks = &VectorialBlocks;
      Part2_Blocks = &VectorialBlocks;
    }
    else
    {
      Part1_Blocks = &SpinorialBlocks;
      Part2_Blocks = &SpinorialBlocks;
    }
  }

  const unsigned MAX_Emergency_Exit = 250000;
  unsigned Emergency_Exit = 0;

  while (Emergency_Exit < MAX_Emergency_Exit)
  {
    ++Emergency_Exit;

    this->CreateRandomVector(this->Lattice, SymmetryBlocks, Part1_Blocks, Part2_Blocks);

    // begin: check that Vector is really from the lattice
    sqr = 0.0;
    sum = 0.0;
    for (i = 0; i < 8; ++i)
    {
      sqr += this->at(i) * this->at(i);
      sum += this->at(i);
    }
    if (this->Lattice == E8xE8)
      sum *= TwistOrder;

    // is first part of N V in the E8xE8 lattice?
    if (((this->Lattice == E8xE8) && is_even(sum)) || (this->Lattice == Spin32))
    {
      if (this->Lattice == E8xE8)
        sum = 0.0;
      for (i = 8; i < 16; ++i)
      {
        sqr += this->at(i) * this->at(i);
        sum += this->at(i);
      }
      sum *= TwistOrder;

      // is second part of N V in the E8xE8 lattice?
      if (is_even(sum))
      {
        sqr = TwistOrder * (sqr - TwistSqr);

        // Check modular invariance condition: N (V^2 - v^2) = 0 mod 2
        if (is_even(sqr))
        {
          if (this->UpdateData()) {
            return true;
          }
          else
          {
            this->SetToZero();
            return false;
          }
        }
      }
    }
    // end: check that Vector is really from the lattice
  }
  this->SetToZero();
  return false;
}



/* ########################################################################################
######   CreateRandomBrother(const vector<CVector> &LatticeVectors, ...)             ######
######                                                                               ######
######   Version: 16.05.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) LatticeVectors : list of (E8xE8 or Spin32/Z2) lattice vectors that can   ######
######                       be added to the current shift in order to constrcut a   ######
######                       brother model                                           ######
######   2) TwistSqr       : the twist vector squared: v^2                           ######
######   3) TwistOrder     : order N of the twist vector v                           ######
######   output:                                                                     ######
######   return value      : finished successfully?                                  ######
###########################################################################################
######   description:                                                                ######
######   Creates a random brother shift vector (i.e. a shift vector that differs     ######
######   from the original one by adding lattice vectors). Adds lattice vectors from ######
######   the set "LatticeVectors" at most two times each vector. Ensures that the    ######
######   modular invariance condition                                                ######
######     N (V^2 - v^2) = 0 mod 2                                                   ######
######   is fulfilled by the randomly chosen brother shift vector.                   ######  
######################################################################################## */
bool CShiftVector::CreateRandomBrother(const vector<CVector> &LatticeVectors, const double &TwistSqr, const unsigned TwistOrder)
{
  const double prob = 0.003;

  const size_t s1 = LatticeVectors.size();

  const double dRAND_MAX = (double)RAND_MAX;

  CVector tmp(16);

  unsigned i = 0;
  double   sqr = 0.0;

  while (true)
  {
    for (i = 0; i < 16; ++i)
      tmp[i] = this->at(i);

    for (i = 0; i < s1; ++i)
    {
      if (rand()/dRAND_MAX < prob)
      {
        if (rand()/dRAND_MAX < 0.5)
        {
          if (rand()/dRAND_MAX < 0.8)
            tmp = tmp + LatticeVectors[i];
          else
            tmp = tmp + (LatticeVectors[i] * 2.0);
        }
        else
        {
          if (rand()/dRAND_MAX < 0.8)
            tmp = tmp - LatticeVectors[i];
          else
            tmp = tmp - (LatticeVectors[i] * 2.0);
        }
      }
    }
    // Check modular invariance condition: N (V^2 - v^2) = 0 mod 2
    sqr = 0.0;
    for (i = 0; i < 16; ++i)
      sqr += tmp[i] * tmp[i];
    sqr -= TwistSqr;
    sqr *= (double)TwistOrder;

    if (is_even(sqr))
    {
      this->clear();
      this->insert(this->end(), tmp.begin(), tmp.end());
      break;
    }
  }
  return true;
}



/* ########################################################################################
######   LoadShiftVector(SelfDualLattice Lattice, ifstream &in)                      ######
######                                                                               ######
######   Version: 24.08.2010                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Lattice   : E8xE8 or Spin32                                              ######
######   1) in        : ifstream object from which the shift will be loaded          ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Loads the shift vector from an ifstream object (e.g. from a file).          ######
######   The lattice type is specified by "Lattice".                                 ######
######################################################################################## */
bool CShiftVector::LoadShiftVector(SelfDualLattice Lattice, ifstream &in)
{
  if (Lattice == UNSPECIFIED_LATTICE)
  {
    cout << "\n  Warning in void CShiftVector::LoadShiftVector(...): The even and self-dual lattice was not defined. Hence, it is set to E8xE8." << endl;
    this->Lattice = E8xE8;
  }
  else
    this->Lattice = Lattice;

  this->resize(16,0);
  rationalVector RationalVector(16, rational<int>(0));

  // begin: load from file
  string currentline = "";

  if (!GetSaveLine(in, currentline))
  {
    this->SetToZero();
    return false;
  }
  convert_string_to_vector_of_rational(currentline, RationalVector);
  // end: load from file

  if (RationalVector.size() != 16)
  {
    cout << "\n  Warning in void CShiftVector::LoadShiftVector(...): Shift has length " << RationalVector.size() << " != 16. Hence, it is set to zero!" << endl;
    this->SetToZero();
    return false;
}

  this->resize(16,0);
  this->Size = 16;

  for (unsigned i = 0; i < 16; ++i)
  {
    const rational<int> &tmp = RationalVector[i];
    (*this)[i] = ((double)tmp.numerator())/((double)tmp.denominator());
  }

  if (!this->UpdateData())
  {
    cout << "\n  Warning in void CShiftVector::LoadShiftVector(...): Shift vector set to zero." << endl;
    this->SetToZero();
    return false;
  }
  return true;
}



/* ########################################################################################
######   SetToZero()                                                                 ######
######                                                                               ######
######   Version: 19.03.2010                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Set this CShiftVector object to the null vector.                            ######
######################################################################################## */
void CShiftVector::SetToZero()
{
  this->assign(16,0);
  this->Size = 16;

  this->Order = 1;
}



/* ########################################################################################
######   UpdateData()                                                                ######
######                                                                               ######
######   Version: 19.03.2010                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Checks this CWilsonLine object and determines the order of the shift        ######
######   vector. Saves the result to the member variable "Order".                    ######
######################################################################################## */
bool CShiftVector::UpdateData()
{
  #ifdef CHECKERROR
  if (this->Size != 16)
  {
    cout << "\n  Warning in bool CShiftVector::UpdateData(): ShiftVector is not a 16-dim. vector." << endl;
    return false;
  }
  if (this->Lattice == UNSPECIFIED_LATTICE)
  {
    cout << "\n  Warning in bool CShiftVector::UpdateData(): Even and self-dual lattice is not specified." << endl;
    return false;
  }
  #endif

  // compute the order of the shift
  const unsigned min_Order =  1;
  const unsigned max_Order = 24;

  CLatticeVector TestVector(16);
  // E8 x E8' lattice
  if (this->Lattice == E8xE8)
  {
    for (unsigned n = min_Order; n < max_Order; ++n)
    {
      TestVector = this->operator*(n);

      if (TestVector.From_E8_Lattice(1) && TestVector.From_E8_Lattice(2))
      {
        this->Order = n;
        return true;
      }
    }
  }
  // Spin32/Z2 lattice
  else
  {
    for (unsigned n = min_Order; n < max_Order; ++n)
    {
      TestVector = this->operator*(n);

      if (TestVector.From_Spin32_Lattice())
      {
        this->Order = n;
        return true;
      }
    }
  }

  cout << "\n  Warning in bool CShiftVector::UpdateData(): Order of ShiftVector not found." << endl;
  return false;
}
