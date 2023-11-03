#ifndef CSPACEGROUPELEMENT_H
#define CSPACEGROUPELEMENT_H

#include <vector>
#include <boost/rational.hpp>

using boost::rational;
using std::vector;

typedef vector<rational<int> >          rationalVector;
typedef vector<vector<rational<int> > > rationalMatrix;

const unsigned LatticeDim = 6;
const unsigned ComplexLatticeDim = 3;

class CLatticeElement : private rationalVector{
public:
  ~CLatticeElement();
  CLatticeElement();
  CLatticeElement(const rationalVector &n_alpha);

  void                          operator=(const rationalVector &n_alpha);
  bool                          operator==(const CLatticeElement &Element2) const;
  bool                          operator!=(const CLatticeElement &Element2) const;
  void                          operator+=(const CLatticeElement &b);
  void                          operator-=(const CLatticeElement &b);
  void                          operator*=(const int &b);
  void                          operator*=(const rational<int> &b);
  rational<int>                 operator*(const rationalVector &ChargeOperator) const;
  CLatticeElement               operator*(const int &b) const;
  const rational<int>          &operator[](const unsigned &alpha) const;

  bool                          IsZero() const;
  bool                          Rotate(const rationalMatrix &TwistMatrix, CLatticeElement &result) const;
  bool                          Set_n_alpha(const unsigned &alpha, const rational<int> &n_alpha);
  const rationalVector         &Vector() const {return *this;};
};


class CSpaceGroupElement{
public:
  ~CSpaceGroupElement();
  CSpaceGroupElement();
  CSpaceGroupElement(unsigned m, unsigned n, unsigned k);
  CSpaceGroupElement(unsigned m, unsigned n, unsigned k, const rationalVector &n_alpha);

  bool                          operator==(const CSpaceGroupElement &Element2) const;
  bool                          operator!=(const CSpaceGroupElement &Element2) const;
  rational<int>                 operator*(const rationalVector &ChargeOperator) const;
  
  bool                          NoTwist() const;
  bool                          IsZero() const;

  bool                          Set_m(const unsigned &m);
  bool                          Set_n(const unsigned &n);
  bool                          Set_k(const unsigned &k);
  bool                          SetLatticeElement(const CLatticeElement &LatticeElement);
  bool                          Set_n_alpha(const rationalVector &n_alpha);
  bool                          Set_n_alpha(const unsigned &alpha, const rational<int> &n_alpha);
  
  const unsigned               &Get_m() const {return this->m;};
  const unsigned               &Get_n() const {return this->n;};
  const unsigned               &Get_k() const {return this->k;};
  const CLatticeElement        &GetLatticeElement() const {return this->LatticeElement;};
  const rational<int>          &Get_n(const unsigned &alpha) const {return this->LatticeElement[alpha];};
  
private:
  unsigned                      m;
  unsigned                      n;
  unsigned                      k;
  CLatticeElement               LatticeElement;
};

#endif
