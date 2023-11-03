#ifndef CGAUGEINVARIANCE_H
#define CGAUGEINVARIANCE_H

#include "groupTheory.hpp"

#ifndef TYPEDEF_CGAUGEGROUP
#define TYPEDEF_CGAUGEGROUP
typedef gaugeGroup<double> CGaugeGroup;
#endif

class COrbifold;
struct SConfig;

class CGaugeInvariance{
public:
  CGaugeInvariance();

  ~CGaugeInvariance();

  bool CheckInvariance(const COrbifold &Orbifold, const SConfig &VEVConfig, const vector<unsigned> &FieldCoupling);
  bool CheckUnknownInvariance(const COrbifold &Orbifold, const SConfig &VEVConfig, const vector<unsigned> &FieldCoupling, unsigned ggf) const;
  bool CheckUnknownInvarianceEasyCases(const gaugeGroupFactor<double> &GaugeGroupFactor, const vector<vector<int> > &Coupling_HWs, bool &CurrentCaseSupported) const;
  bool LoadTensorProducts(const CGaugeGroup &GaugeGroup);

  bool TensorProduct(const vector<vector<int> > &DL_weights1, const vector<vector<int> > &DL_weights2, const gaugeGroupFactor<double> &ggf, vector<vector<vector<int> > > &Decomposed_DL_Weights) const;
  bool TensorProductWithItself(const vector<vector<int> > &DL_weights, const gaugeGroupFactor<double> &ggf, vector<vector<vector<int> > > &Decomposed_DL_Weights_sym, vector<vector<vector<int> > > &Decomposed_DL_Weights_asym) const;
  bool DecomposeIntoIrreps(vector<vector<int> > &DL_weights, const gaugeGroupFactor<double> &ggf, vector<vector<vector<int> > > &result_weights, vector<int> &result_dimensions) const;

  vector<vector<vector<vector<int> > > > AllInvariantTensorProducts;
  vector<vector<vector<vector<int> > > > AllNonInvariantTensorProducts;

  vector<string> Algebras;
  vector<string> Filenames;
};


class CGaugeIndices{
public:
  CGaugeIndices();
  ~CGaugeIndices();

  const vector<vector<int> >          &GetCubicIndices() const {return this->CubicIndices;};
  const vector<vector<vector<int> > > &GetQuadraticIndices() const {return this->QuadraticIndices;};
  const vector<vector<int> >          &GetQuadraticIndicesAdj() const {return this->QuadraticIndicesAdj;};

  bool GetQuadraticIndex(const gaugeGroupFactor<double> &ggf, const unsigned &Dimension, unsigned &Index) const;
  bool GetQuadraticIndexAdj(const gaugeGroupFactor<double> &ggf, unsigned &Index) const;
  
private:
  vector<vector<int> >          CubicIndices;
  vector<vector<vector<int> > > QuadraticIndices;
  vector<vector<int> >          QuadraticIndicesAdj;
};

#endif
