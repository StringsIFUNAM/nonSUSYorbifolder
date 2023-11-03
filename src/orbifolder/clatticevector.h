
#ifndef CLATTICEVECTOR_H
#define CLATTICEVECTOR_H

#include "cvector.h"

//enum LatticeType {NO_LATTICE=0, SO8_V, SO8_S, SO8_C}; //hacking here to include cospinor lattice
enum LatticeType {NO_LATTICE=0, SO8_V, SO8_S, SO8_C};

class CLatticeVector : public CVector{
public: 
// member functions
  CLatticeVector();
  CLatticeVector(const unsigned dim);
  CLatticeVector(const CVector &Vector2);
  ~CLatticeVector();

  void        operator=(const CVector &Vector);
  void        operator=(const vector<double> &DVector);

  bool        From_E8_Lattice(unsigned part) const;
  bool        From_G8_Lattice(unsigned part) const;
  LatticeType From_SO8S_Lattice() const;					//hacking here!!!
  LatticeType From_SO8C_Lattice() const;					//hacking here!!!
  bool        From_Spin32_Lattice() const;
};

#endif
