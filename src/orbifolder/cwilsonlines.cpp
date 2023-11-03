#include <iostream>
#include <fstream>
#include <sstream>

#include "cwilsonlines.h"
#include "globalfunctions.h"

using std::cout;
using std::endl;
using std::vector;



/* ########################################################################################
######   operator*                                                                   ######
######                                                                               ######
######   Version: 05.08.2011                                                         ######
######   Check-Level: 2                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   LatticeElement : a six. dimensional vector of rational numbers              ######
######   output:                                                                     ######
######   return         : a CVector object                                           ######
###########################################################################################
######   description:                                                                ######
######   Computes the vector                                                         ######
######     \sum_{\alpha=1}^6 n_alpha W_alpha                                         ######
######   and returns the result as a CVector object.                                 ######
######################################################################################## */
CVector CWilsonLines::operator*(const CLatticeElement &LatticeElement) const
{
  CVector result(16);
  if (this->WL_CheckStatus == NotChecked)
  {
    cout << "\n  Warning in CVector CWilsonLines::operator*(...) const : data not checked. Return zero." << endl;
    return result;
  }

  double n = 0.0;

  for (unsigned i = 0; i < LatticeDim; ++i)
  {
    if ((LatticeElement[i].numerator() != 0) && !this->Set[i].GetIs_Zero())
    {
      n = ((double)LatticeElement[i].numerator())/((double)LatticeElement[i].denominator());
      result += (this->Set[i] * n);
    }
  }
  return result;
}



/* ########################################################################################
######   CWilsonLines(const SelfDualLattice &Lattice)                                ######
######                                                                               ######
######   Version: 05.08.2011                                                         ######
######   Check-Level: 2                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Lattice : E8xE8 or Spin32                                                ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Constructor of a CWilsonLines object. Creates six Wilson lines of lattice   ######
######   type "Lattice", each being a null vector.                                   ######
######################################################################################## */
CWilsonLines::CWilsonLines(const SelfDualLattice &Lattice)
  : WL_CheckStatus(NotChecked)
{
  if (Lattice == UNSPECIFIED_LATTICE)
    cout << "Warning in CWilsonLines::CWilsonLines(...) : Even and self-dual lattice is not specified." << endl;

  this->SetToZero(Lattice);
}



/* ########################################################################################
######   CWilsonLines()                                                              ######
######                                                                               ######
######   Version: 05.08.2011                                                         ######
######   Check-Level: 2                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard constructor of a CWilsonLines object. Creates six Wilson lines,    ######
######   each being a null vector.                                                   ######
######################################################################################## */
CWilsonLines::CWilsonLines()
  : WL_CheckStatus(NotChecked)
{
  this->SetToZero(UNSPECIFIED_LATTICE);
}



/* ########################################################################################
######   LoadWilsonLines(const SelfDualLattice &Lattice, ifstream &in)               ######
######                                                                               ######
######   Version: 24.08.2011                                                         ######
######   Check-Level: 2                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Lattice   : E8xE8 or Spin32                                              ######
######   2) in        : ifstream object that contains the file from which the Wilson ######
######                  lines are read                                               ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Loads the six Wilson lines from an ifstream object (e.g. from a file).      ######
######   The lattice type is specified by "Lattice".                                 ######
######################################################################################## */
bool CWilsonLines::LoadWilsonLines(const SelfDualLattice &Lattice, ifstream &in)
{
  if (Lattice == UNSPECIFIED_LATTICE)
  {
    cout << "Warning in bool CWilsonLines::LoadWilsonLines(...) : Even and self-dual lattice is not specified. Return false." << endl;
    return false;
  }

  string currentline = "";

  rationalVector RationalVector;
  CVector Vector(16);
  unsigned j = 0;

  // convert the data to 16-dim. vectors of rationals
  for (unsigned i = 0; i < LatticeDim; ++i)
  {
    if (!GetSaveLine(in, currentline))
    {
      this->WL_CheckStatus = NotChecked;
      return false;
    }
    convert_string_to_vector_of_rational(currentline, RationalVector);

    if (RationalVector.size() != 16)
    {
      cout << "Warning in bool CWilsonLines::LoadWilsonLines(...) : WilsonLine is not a 16-dim. vector. Hence set to zero." << endl;
      Vector = CVector(16);
    }
    else
    {
      for (j = 0; j < 16; ++j)
      {
        const rational<int> &tmp = RationalVector[j];
        Vector[j] = ((double)tmp.numerator())/((double)tmp.denominator());
      }
    }
    this->Set[i] = Vector;
  }
  this->WL_CheckStatus = NotChecked;
  return true;
}



/* ########################################################################################
######   SetLattice(const SelfDualLattice &Lattice)                                  ######
######                                                                               ######
######   Version: 05.08.2011                                                         ######
######   Check-Level: 2                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Lattice   : E8xE8 or Spin32                                              ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Change the lattice type of the six Wilson lines to "Lattice".               ######
######################################################################################## */
bool CWilsonLines::SetLattice(const SelfDualLattice &Lattice)
{
  if (Lattice == UNSPECIFIED_LATTICE)
  {
    cout << "Warning in bool CWilsonLines::SetLattice(...) : lattice not defined. Return false." << endl;
    return false;
  }

  for (unsigned i = 0; i < LatticeDim; ++i)
    this->Set[i].SetLattice(Lattice);

  this->WL_CheckStatus = NotChecked;
  return true;
}



/* ########################################################################################
######   SetToZero(const SelfDualLattice &Lattice)                                   ######
######                                                                               ######
######   Version: 10.08.2011                                                         ######
######   Check-Level: 2                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Sets the six CWilsonLine objects to null vectors of lattice type "Lattice". ######
######################################################################################## */
void CWilsonLines::SetToZero(const SelfDualLattice &Lattice)
{
  if (Lattice == UNSPECIFIED_LATTICE)
  {
    const CWilsonLine null;
    this->Set.assign(LatticeDim, null);

    this->WL_CheckStatus = NotChecked;
  }
  else
  {
    const CWilsonLine null(Lattice);
    this->Set.assign(LatticeDim, null);

    this->WL_CheckStatus = CheckedAndGood;
  }
}



/* ########################################################################################
######   SetWilsonLines(const CWilsonLines &W)                                       ######
######                                                                               ######
######   Version: 08.08.2011                                                         ######
######   Check-Level: 2                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Copies the content of "W" to this CWilsonLines object.                      ######
######################################################################################## */
void CWilsonLines::SetWilsonLines(const CWilsonLines &W)
{
  this->WL_CheckStatus = NotChecked;
  this->Set = W.Set;
}



/* ########################################################################################
######   SetWilsonLines(...)                                                         ######
######                                                                               ######
######   Version: 05.08.2011                                                         ######
######   Check-Level: 2                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) W_1       : Wilson line corresponding to e_1                             ######
######   1) W_2       : Wilson line corresponding to e_2                             ######
######   1) W_3       : Wilson line corresponding to e_3                             ######
######   1) W_4       : Wilson line corresponding to e_4                             ######
######   1) W_5       : Wilson line corresponding to e_5                             ######
######   1) W_6       : Wilson line corresponding to e_6                             ######
######                                                                               ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Sets the six CWilsonLine objects to "W_i" for i = 1,...,6.                  ######
######################################################################################## */
bool CWilsonLines::SetWilsonLines(const CWilsonLine &W_1, const CWilsonLine &W_2, const CWilsonLine &W_3,
                                  const CWilsonLine &W_4, const CWilsonLine &W_5, const CWilsonLine &W_6)
{
  if (LatticeDim != 6)
  {
    cout << "Warning in bool CWilsonLines::SetWilsonLines(...) : \"LatticeDim\" not 6. Return false." << endl;
    return false;
  }

  const SelfDualLattice Lattice_W2 = W_2.GetLattice();
  const SelfDualLattice Lattice_W3 = W_3.GetLattice();
  const SelfDualLattice Lattice_W4 = W_4.GetLattice();
  const SelfDualLattice Lattice_W5 = W_5.GetLattice();

  if ((W_1.GetLattice() != Lattice_W2) || (Lattice_W2 != Lattice_W3) || (Lattice_W3 != Lattice_W4) || (Lattice_W4 != Lattice_W5) || (Lattice_W5 != W_6.GetLattice()))
  {
    cout << "Warning in bool CWilsonLines::SetWilsonLines(...) : the Wilson lines belong to different 16-dim. lattices. Return false." << endl;
    return false;
  }

  this->Set[0] = W_1;
  this->Set[1] = W_2;
  this->Set[2] = W_3;
  this->Set[3] = W_4;
  this->Set[4] = W_5;
  this->Set[5] = W_6;

  this->WL_CheckStatus = NotChecked;
  return true;
}



/* ########################################################################################
######   SetWilsonLine(...)                                                          ######
######                                                                               ######
######   Version: 05.08.2011                                                         ######
######   Check-Level: 2                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) i         : i-th dircetion e_i                                           ######
######   2) W         : Wilson line corresponding to e_i                             ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Sets the i-th CWilsonLine object to "W".                                    ######
######################################################################################## */
bool CWilsonLines::SetWilsonLine(unsigned i, const CWilsonLine &W)
{
  if (i >= LatticeDim)
  {
    cout << "Warning in bool CWilsonLines::SetWilsonLine(...) : index \"i\" out of range. Return false." << endl;
    return false;
  }

  this->Set[i] = W;
  this->WL_CheckStatus = NotChecked;
  return true;
}


  
/* ########################################################################################
######   ~CWilsonLines()                                                             ######
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
######   Standard destructor of a CWilsonLines object.                               ######
######################################################################################## */
CWilsonLines::~CWilsonLines()
{
}



/* ########################################################################################
######   GetWilsonLine(unsigned i) const                                             ######
######                                                                               ######
######   Version: 04.08.2011                                                         ######
######   Check-Level: 2                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) i         : index i from 0,..,5                                          ######
######                                                                               ######
######   output:                                                                     ######
######   return value : if i = 0,..,5: i-th Wilson line; else first Wilson line      ######
###########################################################################################
######   description:                                                                ######
######   Get (constant) access to the i-th CWilsonLine object for i = 0,..,5.        ######
######################################################################################## */
const CWilsonLine &CWilsonLines::GetWilsonLine(unsigned i) const
{
  if (this->WL_CheckStatus == NotChecked)
    cout << "\n  Warning in const CWilsonLine &CWilsonLines::GetWilsonLine(...) const : data not checked." << endl;

  if (i < LatticeDim)
    return this->Set[i];

  cout << "Warning in GetWilsonLine(CWilsonLine) : index i = " << i << " out of range. Return first Wilson line." << endl;
  return this->Set[0];
}



/* ########################################################################################
######   GetWilsonLines() const                                                      ######
######                                                                               ######
######   Version: 09.08.2011                                                         ######
######   Check-Level: 2                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : the set of Wilson lines                                      ######
###########################################################################################
######   description:                                                                ######
######   Get (constant) access to the set of six CWilsonLine objects.                ######
######################################################################################## */
const vector<CWilsonLine> &CWilsonLines::GetWilsonLines() const
{
  if (this->WL_CheckStatus == NotChecked)
    cout << "\n  Warning in const vector<CWilsonLine> &CWilsonLines::GetWilsonLines() const : data not checked." << endl;

  return this->Set;
}



/* ########################################################################################
######   Check(...)                                                                  ######
######                                                                               ######
######   Version: 10.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Check that the six CWilsonLine objects are of the correct order and that    ######
######   they fulfill the relations among each other as imposed by the geometry of   ######
######   the orbifold. For example for Z_3 the six Wilson lines W_i must fulfill:    ######
######     1) 3W_i must be in the gauge lattice and                                  ######
######     2) W_1 = W_2, W_3 = W_4 and W_5 = W_6.                                    ######
######################################################################################## */
bool CWilsonLines::Check(const vector<vector<unsigned> > &WL_Relations, const vector<unsigned> &WL_AllowedOrders)
{
  if (this->WL_CheckStatus == NotChecked)
  {
    if (WL_AllowedOrders.size() != LatticeDim)
    {
      cout << "\n  Warning in bool CWilsonLines::Check(...) : Allowed orders of Wilson lines not known. Check cancelled and return false." << endl;
      this->WL_CheckStatus = CheckedAndFailed;
      return false;
    }

    size_t s2 = 0;
    unsigned i = 0;
    unsigned j = 0;

    for (i = 0; i < LatticeDim; ++i)
    {
      j = WL_AllowedOrders[i];
      // check only if the allowed order of the wilson line is known
      if ((j != 0) && ((j % this->Set[i].GetOrder()) != 0))
      {
        this->WL_CheckStatus = CheckedAndFailed;
        return false;
      }
    }

    const size_t s1 = WL_Relations.size();
    for (i = 0; i < s1; ++i)
    {
      const vector<unsigned> &Relation = WL_Relations[i];

      s2 = Relation.size();
      if (s2 == 0)
      {
        cout << "\n  Warning in bool CWilsonLines::Check(...) : Relation between Wilson lines not well-defined. Check cancelled and return false." << endl;
        this->WL_CheckStatus = CheckedAndFailed;
        return false;
      }

      const CWilsonLine &WL = this->Set[Relation[0]];

      for (j = 1; j < s2; ++j)
      {
        if (WL != this->Set[Relation[j]])
        {
          this->WL_CheckStatus = CheckedAndFailed;
          return false;
        }
      }
    }

    // if all checks passed
    this->WL_CheckStatus = CheckedAndGood;
  }
  return true;
}
