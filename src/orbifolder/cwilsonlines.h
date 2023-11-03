
#ifndef CWILSONLINES_H
#define CWILSONLINES_H
#include <vector>

#include "cwilsonline.h"
#include "cspacegroupelement.h"

//! CWilsonLines.
/*!
  A CWilsonLines object is a set of six CWilsonLine obejcts. Each one is associated to one of the six torus translations.
 */

using std::ifstream;

#ifndef ENUM_CHECKSTATUS
#define ENUM_CHECKSTATUS
enum CheckStatus {UNSPECIFIED_CHECKSTATUS = 0, NotChecked, CheckedAndGood, CheckedAndFailed};
#endif


class CWilsonLines {
public:
  CWilsonLines();
  CWilsonLines(const SelfDualLattice &Lattice);
  ~CWilsonLines();

  CVector                    operator*(const CLatticeElement &LatticeElement) const;

  bool                       Check(const vector<vector<unsigned> > &WL_Relations, const vector<unsigned> &WL_AllowedOrders);
  void                       WL_NotChecked() {this->WL_CheckStatus = NotChecked;};

  bool                       LoadWilsonLines(const SelfDualLattice &Lattice, ifstream &in);
  bool                       SetLattice(const SelfDualLattice &Lattice);
  void                       SetToZero(const SelfDualLattice &Lattice);
  void                       SetWilsonLines(const CWilsonLines &W);
  bool                       SetWilsonLines(const CWilsonLine &W_1, const CWilsonLine &W_2, const CWilsonLine &W_3,
                                           const CWilsonLine &W_4, const CWilsonLine &W_5, const CWilsonLine &W_6);
  bool                       SetWilsonLine(unsigned i, const CWilsonLine &W);

  const CheckStatus         &GetCheckStatus() const {return this->WL_CheckStatus;};
  const CWilsonLine         &GetWilsonLine(unsigned i) const;
  const vector<CWilsonLine> &GetWilsonLines() const;

private:
  vector<CWilsonLine>        Set;
  CheckStatus                WL_CheckStatus;
};

#endif
