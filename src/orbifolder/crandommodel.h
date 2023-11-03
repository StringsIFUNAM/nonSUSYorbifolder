#ifndef CRANDOMMODEL_H
#define CRANDOMMODEL_H

#include "corbifoldgroup.h"

//! CRandomModel.
/*!
The class CRandomModel is used by the member function CreateRandom(..) of the class COrbifoldGroup in order to create random shifts and Wilson lines.
 */


class CRandomModel{
public:
  CRandomModel(const SelfDualLattice &Lattice);
  ~CRandomModel();

  bool                                      Initiate(const COrbifoldGroup &OrbifoldGroup, const vector<bool> &UseOrigShiftsAndWilsonLines, const vector<CVector> &UnbrokenRoots);
  vector<unsigned>                          Obtain_Symmetry_Factors(const vector<unsigned> &oldfactors, const CVector &Shift_or_WL) const;
  bool                                      CheckUnbrokenRoots(const CVector &Vector) const;

  const vector<unsigned>                   &GetOrigFactors() const {return this->OrigFactors;};
  const vector<CVector>                    &GetLatticeVectors() const {return this->LatticeVectors;};

  const vector<vector<vector<double> > >   &GetVectorialBlocksOfOrder(const unsigned &Order) const;
  const vector<vector<vector<double> > >   &GetSpinorialBlocksOfOrder(const unsigned &Order) const;

  const int                                &GetCreateOrIdentifyWL(const unsigned &i) const {return this->CreateOrIdentifyWL[i];};
  const vector<bool>                       &GetUseOriginalVectors() const {return this->UseOriginalVectors;};

  const unsigned                           &GetMAX_Emergency_Exit() const {return this->MAX_Emergency_Exit;};
  const double                             &Get_v1v2() const {return this->v1v2;};
  const double                             &Get_v0v1() const {return this->v0v1;};
  const double                             &Get_v0v2() const {return this->v0v2;};

private:
  vector<unsigned>                          OrigFactors;
  vector<CVector>                           UnbrokenRoots;

  vector<CVector>                           LatticeVectors;

  vector<vector<vector<vector<double> > > > VectorialBlocks;
  vector<vector<vector<vector<double> > > > SpinorialBlocks;

  vector<int>                               CreateOrIdentifyWL;
  vector<bool>                              UseOriginalVectors;

  unsigned                                  MAX_Emergency_Exit;
  double                                    v1v2;
  double                                    v0v1;
  double                                    v0v2;
};

#endif
