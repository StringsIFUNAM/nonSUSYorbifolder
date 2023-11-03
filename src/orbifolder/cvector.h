
#ifndef CVECTOR_H
#define CVECTOR_H

#include <iostream>
#include <iomanip>
#include <cstdlib>

#include <vector>
#include <boost/config.hpp>
#include <boost/rational.hpp>

//! CVector.
/*!
A CVector object is a vector of dimension Size.
 */

using boost::rational;
using std::vector;

typedef vector<int>            intVector;
typedef vector<double>         doubleVector;
typedef vector<rational<int> > rationalVector;

#ifndef ENUM_SELFDUALLATTICE
#define ENUM_SELFDUALLATTICE
enum SelfDualLattice {UNSPECIFIED_LATTICE = 0, E8xE8, Spin32, SO8};
#endif

class CVector : public doubleVector{
  friend class CWilsonLine;
  friend class CTwistVector;
  friend class CShiftVector;
  friend class CLatticeVector;

public: 
// member functions
  CVector();
  CVector(const unsigned dim);
  ~CVector();

  void    operator=(const doubleVector &Vector);
  void    operator=(const intVector &Vector);
  void    operator=(const vector<rational<int> > &RationalVector);
  CVector operator+(const CVector &Vector2) const;
  CVector operator-(const CVector &Vector2) const;
  CVector operator*(double a) const;
  double  operator*(const vector<CVector>::const_iterator &it_Vector2) const;
  double  operator*(const CVector &Vector2) const;
  void    operator+=(const CVector &Vector2);
  void    operator-=(const CVector &Vector2);
  bool    operator==(const CVector &Vector2) const;
  bool    operator!=(const CVector &Vector2) const;

  void            Assign(const unsigned dim);
  void            Assign(const unsigned dim, const double value);
  bool            CreateRandomVector(SelfDualLattice Lattice, const vector<unsigned> &SymmetryBlocks, const vector<vector<vector<double> > > *Part1_Blocks, const vector<vector<vector<double> > > *Part2_Blocks);
  double          GetLength() const;
  const unsigned &GetSize() const {return this->Size;};
  double          GetSqrTo(unsigned cut) const;
  bool            IsScalarProductInteger(const CVector &Vector2) const;
  bool            IsZero() const;
  void            Push_back(double value);
  void            SetToZero();

private:
// member variables
  unsigned Size;
};

#endif
