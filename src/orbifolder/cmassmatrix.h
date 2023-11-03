#ifndef CMASSMATRIX_H
#define CMASSMATRIX_H

#include "corbifold.h"

struct YukawaCoupling;

class CMassMatrix {
public: 
  CMassMatrix();
  CMassMatrix(const SConfig &Vacuum, const string &Label_Row, const string &Label_Column, bool AutoTranspose);

  bool IsEmpty() const;
  bool IsRankMaximal(unsigned &row_rank, unsigned &column_rank) const;
  void SetRandom_vev();
  void Update(const SConfig &Vacuum);
  const vector<vector<vector<YukawaCoupling> > > &GetMatrix() const {return this->Matrix;};
  
  ~CMassMatrix();

  vector<unsigned> FieldIndices_Column;
  vector<unsigned> FieldIndices_Row;

  string Label_Column;
  string Label_Row;

private:
  vector<vector<vector<YukawaCoupling> > > Matrix;
};

#endif
