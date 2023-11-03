#ifndef CMONOMIAL_H
#define CMONOMIAL_H

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <vector>

#include "cvector.h"
#include "corbifold.h"

using std::string;
using std::vector;

//! CMonomial.
/*!
  A CMonomial object M describes a gauge invariant monomial in the fields \phi_i, i.e.
    M = \phi_i^(n_i) ... \phi_j^(n_j)
  with n_i being positive integers. They correspond to D-flat directions of the theory.
 */

class COrbifold;

class CMonomial{
public:
// member functions
  //! A constructor.
  CMonomial(unsigned NumberOfU1s);
  //! A destructor.
  ~CMonomial();

  void operator*=(const CMonomial &Monomial2);

  bool CheckGaugeInvariance(const COrbifold &Orbifold, const SConfig &Vacuum, bool demand_D0 = true, bool check_non_Abelian_invariance = true) const;
  bool ContainsField(unsigned FieldIndex) const;
  bool ReadMonomial(string input, const vector<CField> &Fields, unsigned use_Labels);
  void SetVEVs(double VEV_Scale, vector<CField> &Fields);

// member variables
  CVector                   U1Charges;
  vector<vector<unsigned> > GaugeEquivalentFields;
  vector<unsigned>          Exponents;
  bool                      VEV_TurnedOn;
};

#endif
