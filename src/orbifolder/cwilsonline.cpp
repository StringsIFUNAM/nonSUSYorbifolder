#include <iostream>

#include "clatticevector.h"
#include "cwilsonline.h"

#include "globalfunctions.h"

using std::cout;
using std::endl;
using std::exit;



/* ########################################################################################
######   CWilsonLine()                                                               ######
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
######   Standard constructor of a CWilsonLine object. Creates a 16-dim. null vector.######
######################################################################################## */
CWilsonLine::CWilsonLine()
{
  this->SetToZero();
  this->Lattice = UNSPECIFIED_LATTICE;
}



/* ########################################################################################
######   CWilsonLine(const SelfDualLattice &Lattice)                                 ######
######                                                                               ######
######   Version: 04.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Lattice : E8xE8 or Spin32                                                ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Constructor of a CWilsonLine object of lattice type "Lattice". Creates a    ###### 
######   16-dim. null vector for either the E8xE8 or the Spin32 gauge lattice.       ######
######################################################################################## */
CWilsonLine::CWilsonLine(const SelfDualLattice &Lattice)
{
  if (Lattice == UNSPECIFIED_LATTICE)
    cout << "Warning in CWilsonLine::CWilsonLine(...) : Even and self-dual lattice is not specified." << endl;

  this->assign(16,0);
  this->Size = 16;

  this->Is_Zero = true;
  this->Lattice = Lattice;
  this->Order   = 1;
}



/* ########################################################################################
######   CWilsonLine(SelfDualLattice Lattice, double A_0, double A_1,...)            ######
######                                                                               ######
######   Version: 19.03.2010                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Lattice      : E8xE8 or Spin32                                           ######
######   2)-17)  A_i     : the 16 components of the wilson line                      ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Constructor of a CWilsonLine object of lattice type "Lattice". Creates a    ###### 
######   16-dim. vector for either the E8xE8 or the Spin32 gauge lattice with        ######
######   entries "A_i" for i = 0,..,15.                                              ######
######################################################################################## */
CWilsonLine::CWilsonLine(SelfDualLattice Lattice, double A_0, double A_1, double  A_2, double  A_3, double  A_4, double  A_5, double  A_6, double  A_7,
                         double A_8, double A_9, double A_10, double A_11, double A_12, double A_13, double A_14, double A_15)
{
  if (Lattice == UNSPECIFIED_LATTICE)
    cout << "Warning in CWilsonLine::CWilsonLine(...) : Even and self-dual lattice is not specified." << endl;

  this->resize(16,0);
  this->Size = 16;

  this->at(0)  = A_0;
  this->at(1)  = A_1;
  this->at(2)  = A_2;
  this->at(3)  = A_3;
  this->at(4)  = A_4;
  this->at(5)  = A_5;
  this->at(6)  = A_6;
  this->at(7)  = A_7;
  this->at(8)  = A_8;
  this->at(9)  = A_9;
  this->at(10) = A_10;
  this->at(11) = A_11;
  this->at(12) = A_12;
  this->at(13) = A_13;
  this->at(14) = A_14;
  this->at(15) = A_15;

  this->Lattice = Lattice;
  if (!this->UpdateData())
  {
    this->SetToZero();
    cout << "Warning in void CWilsonLine::CWilsonLine(...) : Wilson line set to zero." << endl;
    return;
  }
}



/* ########################################################################################
######   ~CWilsonLine()                                                              ######
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
######   Standard destructor of a CWilsonLine object.                                ######
######################################################################################## */
CWilsonLine::~CWilsonLine()
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
######   1) Vector : a CVector with 16 entries                                       ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Specify the content of this CWilsonLine object by a CVector object.         ######
######################################################################################## */
void CWilsonLine::operator=(const CVector &Vector)
{
  if (Vector.Size != 16)
  {
    cout << "Warning in void CWilsonLine::operator=(...) : Vector has length " << Vector.size() << " != 16. Hence Wilson line set to zero." << endl;
    this->SetToZero();
    return;
  }

  this->assign(16,0);
  for (unsigned i = 0; i < 16; ++i)
    this->at(i) = Vector[i];

  if (!this->UpdateData())
  {
    cout << "Warning in void CWilsonLine::operator=(...) : Wilson line set to zero." << endl;
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
######   Specify the content of this CWilsonLine object by a vector of rational      ######
######   numbers.                                                                    ######
######################################################################################## */
void CWilsonLine::operator=(const vector<rational<int> > &RationalVector)
{
  if (RationalVector.size() != 16)
  {
    cout << "Warning in void CWilsonLine::operator=(...) : Vector has length " << RationalVector.size() << " != 16. Hence Wilson line set to zero." << endl;
    this->SetToZero();
    return;
  }

  this->assign(16,0);
  for (unsigned i = 0; i < 16; ++i)
  {
    const rational<int> &tmp = RationalVector[i];
    this->at(i) = ((double)tmp.numerator())/((double)tmp.denominator());
  }

  if (!this->UpdateData())
  {
    cout << "Warning in void CWilsonLine::operator=(...) : Wilson line set to zero." << endl;
    this->SetToZero();
    return;
  }
}



/* ########################################################################################
######   CreateRandom(...)                                                           ######
######                                                                               ######
######   Version: 16.05.2010                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) WLOrder         : the order of the Wilson line                           ######
######   2) SymmetryBlocks  : the symmetry structure w.r.t. permutations of the      ######
######                        vector entries                                         ######
######   3) VectorialBlocks : possible parts of the random vector from the vectorial ######
######                        part of the lattice                                    ######
######   4) SpinorialBlocks : possible parts of the random vector from the spinorial ######
######                        part of the lattice                                    ######
######   output:                                                                     ######
######   return value       : random vector created successfully?                    ######
###########################################################################################
######   description:                                                                ######
######   Creates a random Wilson line of order "WLOrder" respecting the symmetry     ###### 
######   structure (as specified in "SymmetryBlocks") and using the contents of      ######
######   "VectorialBlocks" and "SpinorialBlocks". Returns true if the Wilson line    ######
######   has been created successfully.                                              ###### 
######################################################################################## */
bool CWilsonLine::CreateRandom(unsigned WLOrder, const vector<unsigned> &SymmetryBlocks, const vector<vector<vector<double> > > &VectorialBlocks, const vector<vector<vector<double> > > &SpinorialBlocks)
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

  while (true)
  {
    this->CreateRandomVector(this->Lattice, SymmetryBlocks, Part1_Blocks, Part2_Blocks);

    if (this->Lattice == E8xE8)
    {
      // set first 8 entries to zero
      if (rand()/(double)RAND_MAX < 0.02)
      {
        for (i = 0; i < 8; ++i)
          this->at(i) = 0.0;
      }
      // set last 8 entries to zero
      if (rand()/(double)RAND_MAX < 0.02)
      {
        for (i = 8; i < 16; ++i)
          this->at(i) = 0.0;
      }
    }
    else
    {
      if (rand()/(double)RAND_MAX < 0.02)
      {
        this->SetToZero();
        return true;
      }
    }

    // begin: check that Vector is really from the lattice
    sqr = 0.0;
    sum = 0.0;
    for (i = 0; i < 8; ++i)
    {
      sqr += this->at(i) * this->at(i);
      sum += this->at(i);
    }
    if (this->Lattice == E8xE8)
      sum *= WLOrder;

    // is first part of N_A A in the E8xE8 lattice or is the lattice Spin32?
    if (((this->Lattice == E8xE8) && is_even(sum)) || (this->Lattice == Spin32))
    {
      if (this->Lattice == E8xE8)
        sum = 0.0;
      for (i = 8; i < 16; ++i)
      {
        sqr += this->at(i) * this->at(i);
        sum += this->at(i);
      }
      sum *= WLOrder;

      // is second part of N_A A in the E8xE8 or Spin32 lattice?
      if (is_even(sum))
      {
        sqr *= WLOrder;

        // Check modular invariance condition: N_A A^2 = 0 mod 2
        if (is_even(sqr))
        {
          if (this->UpdateData())
            return true;
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
######   CreateRandomBrother(const vector<CVector> &LatticeVectors)                  ######
######                                                                               ######
######   Version: 16.05.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) LatticeVectors : the lattice vector that can be added                    ######
######   output:                                                                     ######
######   return value      : random brother created successfully?                    ######
###########################################################################################
######   description:                                                                ######
######   Creates a random brother Wilson line (i.e. a Wilson line that differs from  ######
######   the original one by adding lattice vectors). Adds lattice vectors from the  ######
######   set "LatticeVectors" at most three times each vector.                       ######  
######################################################################################## */
bool CWilsonLine::CreateRandomBrother(const vector<CVector> &LatticeVectors)
{
  const double prob = 0.003;

  const size_t s1 = LatticeVectors.size();

  const double dOrder    = (double)this->Order;
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
          {
            if (rand()/dRAND_MAX < 0.5)
              tmp = tmp + (LatticeVectors[i] * 2.0);
            else
              tmp = tmp + (LatticeVectors[i] * 3.0);
          }
        }
        else
        {
          if (rand()/dRAND_MAX < 0.8)
            tmp = tmp - LatticeVectors[i];
          else
          {
            if (rand()/dRAND_MAX < 0.5)
              tmp = tmp - (LatticeVectors[i] * 2.0);
            else
              tmp = tmp - (LatticeVectors[i] * 3.0);
          }
        }
      }
    }
    // Check modular invariance condition: N_A A^2 = 0 mod 2
    sqr = 0.0;
    for (i = 0; i < 16; ++i)
      sqr += tmp[i] * tmp[i];

    sqr *= dOrder;

    if (is_even(sqr))
    {
      this->clear();
      this->insert(this->end(), tmp.begin(), tmp.end());
      break;
    }
  }
  this->Is_Zero = false;
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
######   Set this CWilsonLine object to the null vector.                             ######
######################################################################################## */
void CWilsonLine::SetToZero()
{
  this->assign(16,0);
  this->Size = 16;

  this->Is_Zero = true;
  this->Order   = 1;
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
######   Checks this CWilsonLine object and determines the order of the Wilson line. ######
######   Saves the result to the member variable "Order". If this CWilsonLine object ######
######   is a null vector set the member variable "Is_Zero" to true.                 ######
######################################################################################## */
bool CWilsonLine::UpdateData()
{
  if (this->Size != 16)
  {
    cout << "Warning in bool CWilsonLine::UpdateData() : WilsonLine is not a 16-dim. vector." << endl;
    return false;
  }
  if (this->Lattice == UNSPECIFIED_LATTICE)
  {
    cout << "Warning in bool CWilsonLine::UpdateData() : Even and self-dual lattice is not specified." << endl;
    return false;
  }

  // Set the precision
  const double prec = 0.00001;

  unsigned i = 0;

  // check if the Wilson line is zero
  this->Is_Zero = true;
  for (i = 0; this->Is_Zero && (i < 16); ++i)
  {
    if (fabs(this->at(i)) > prec)
      this->Is_Zero = false;
  }

  if (this->Is_Zero)
  {
    this->Order = 1;
    return true;
  }

  // compute the order of the wilson line
  const unsigned min_Order =  1;
  const unsigned max_Order = 20;

  // E8 x E8' lattice
  if (this->Lattice == E8xE8)
  {
    for (i = min_Order; i < max_Order; ++i)
    {
      CLatticeVector TestVector = (*this) * i;

      if (TestVector.From_E8_Lattice(1) && TestVector.From_E8_Lattice(2))
      {
        this->Order = i;
        return true;
      }
    }
  }
  // Spin32/Z2 lattice
  else
  {
    for (i = min_Order; i < max_Order; ++i)
    {
      CLatticeVector TestVector = (*this) * i;

      if (TestVector.From_Spin32_Lattice())
      {
        this->Order = i;
        return true;
      }
    }
  }

  cout << "Warning in bool CWilsonLine::UpdateData() : Order of WilsonLine not found." << endl;
  return false;
}
