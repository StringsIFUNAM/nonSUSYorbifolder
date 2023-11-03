#include "cgaugeinvariance.h"

#include "cfixedpoint.h"
#include "cstate.h"
#include "corbifold.h"
#include "cprint.h"
#include "globalfunctions.h"




/* ########################################################################################
######   bool CGaugeIndices::CGaugeIndices()                                         ######
######                                                                               ######
######   Version: 09.02.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######                                                                               ######
######   output:                                                                     ######
######                                                                               ######
######################################################################################## */
CGaugeIndices::CGaugeIndices()
{
  vector<int>          max_dim_of_rep(221, -1);
  vector<vector<int> > max_rank(16, max_dim_of_rep);
  
  unsigned dim   = 0;
  int index = 1;
  unsigned i = 1;
  unsigned j = 1;
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////
  // begin: define cubic indices of reps
  this->CubicIndices = max_rank;

  // A_N, highest weight L = (1,0,...,0)_DL
  for (i = 2; i < 16; ++i)
  {
    dim   = i + 1;
    index = 1;
    this->CubicIndices.at(i).at(dim) = index;
  }
  // for example       A     A1    2    is 1
  // for example       A     A4    5    is 1

  // A_N, highest weight L = (0,1,0,...,0)_DL
  for (i = 3; i < 16; ++i)
  {
    dim   = (unsigned)round_double_to_int(1.0/2.0 * (i + 1.0) * i);
    index = i - 3;
    this->CubicIndices.at(i).at(dim) = index;
  }
  //for example        A     A3     6   is 0
  //for example        A     A4    10   is 1
  //for example        A     A5    15   is 2

  // A_N, highest weight L = (0,0,1,0,...,0)_DL
  for (i = 5; i < 12; ++i)
  {
    dim   = (unsigned)round_double_to_int(1.0/6.0 * (i + 1.0) * i * (i - 1.0));
    index = round_double_to_int(0.5 * (i - 2.0) * (i - 5.0));
    this->CubicIndices.at(i).at(dim) = index;
  }
  //for example        A     A5    20   is 0
  //for example        A     A7    56   is 5
  //for example        A     A8    84   is 9

  // A_N, highest weight L = (0,0,0,1,0,...,0)_DL
  for (i = 7; i < 10; ++i)
  {
    dim   = (unsigned)round_double_to_int(1.0/24.0 * (i + 1.0) * i * (i - 1.0) * (i - 2.0));
    index = round_double_to_int(1.0/6.0 * (i - 2.0) * (i - 3.0) * (i - 7.0));
    this->CubicIndices.at(i).at(dim) = index;
  }
  //for example        A     A7    70   is  0
  //for example        A     A8   126   is  5
  //for example        A     A9   210   is 14

  //singlets
  for (i = 2; i < 16; ++i)
    this->CubicIndices[i][1] = 0;
  // end: define cubic indices of reps
  /////////////////////////////////////////////////////////////////////////////////////////////////////
    

  /////////////////////////////////////////////////////////////////////////////////////////////////////
  // begin: define indices of reps
  vector<int> Indices_of_AN(17, -1);
  for (i = 1; i < 17; ++i)
    Indices_of_AN[i] = 2*(i+1);

  vector<int> Indices_of_DN(17, -1);
  for (i = 1; i < 17; ++i)
    Indices_of_DN[i] = 4*(i-1);

  vector<int> Indices_of_EN(9, -1);
  Indices_of_EN[6] = 24;
  Indices_of_EN[7] = 36;
  Indices_of_EN[8] = 60;
  
  this->QuadraticIndicesAdj.push_back(Indices_of_AN);
  this->QuadraticIndicesAdj.push_back(Indices_of_DN);
  this->QuadraticIndicesAdj.push_back(Indices_of_EN);


  this->QuadraticIndices.assign(3,max_rank);
  // A_N, highest weight L = (1,0,...,0)_DL
  for (i = 1; i < 16; ++i)
  {
    dim   = i + 1;
    index = 1;
    this->QuadraticIndices.at(0).at(i).at(dim) = index;
  }
  // for example       A     A1    2    is 1
  // for example       A     A4    5    is 1

  // A_N, highest weight L = (0,1,0,...,0)_DL
  for (i = 3; i < 16; ++i)
  {
    dim   = (unsigned)round_double_to_int(1.0/2.0 * (i + 1.0) * i);
    index = i - 1;
    this->QuadraticIndices.at(0).at(i).at(dim) = index;
  }
  //for example        A     A3     6   is 2
  //for example        A     A4    10   is 3
  //for example        A     A5    15   is 4

  // A_N, highest weight L = (0,0,1,0,...,0)_DL
  for (i = 5; i < 12; ++i)
  {
    dim   = (unsigned)round_double_to_int(1.0/6.0 * (i + 1.0) * i * (i - 1.0));
    index = round_double_to_int(0.5 * (i - 1.0) * (i - 2.0));
    this->QuadraticIndices.at(0).at(i).at(dim) = index;
  }
  //for example        A     A5    20   is  6
  //for example        A     A7    56   is 15
  //for example        A     A8    84   is 21

  // A_N, highest weight L = (0,0,0,1,0,...,0)_DL
  for (i = 7; i < 10; ++i)
  {
    dim   = (unsigned)round_double_to_int(1.0/24.0 * (i + 1.0) * i * (i - 1.0) * (i - 2.0));
    index = round_double_to_int(1.0/6.0 * (i - 1.0) * (i - 2.0) * (i - 3.0));
    this->QuadraticIndices.at(0).at(i).at(dim) = index;
  }
  //for example        A     A7    70   is 20
  //for example        A     A8   126   is 35
  //for example        A     A9   210   is 56

  // D_N, highest weight L = (1,0,...,0)_DL
  for (i = 4; i < 16; ++i)
  {
    dim   = 2 * i;
    index = 2;
    this->QuadraticIndices.at(1).at(i).at(dim) = index;
  }
  //for example          D     D4    8  is  2
  //for example          D     D5   10  is  2
  //for example          D     D6   12  is  2

  // D_N, highest weight L = (0,...,0,1)_DL
  for (i = 5; i < 9; ++i)
  {
    dim   = (unsigned)round_double_to_int(pow(2.0,double(i-1)));
    index = round_double_to_int(pow(2.0,double(i-3)));
    this->QuadraticIndices.at(1).at(i).at(dim) = index;
  }
  //for example          D     D5   16  is   4
  //for example          D     D6   32  is   8
  //for example          D     D7   64  is  16
  //for example          D     D8  128  is  32

  //singlets
  for (i = 0; i < 3; ++i)
  {
    for (j = 1; j < 16; ++j)
      this->QuadraticIndices[i][j][1] = 0;
  }

  // E6
  this->QuadraticIndices[2][6][27] = 6;

  // E7
  this->QuadraticIndices[2][7][56] = 12;
  // end: define indices of reps
  /////////////////////////////////////////////////////////////////////////////////////////////////////
}


/* ########################################################################################
######   bool CGaugeIndices::~CGaugeIndices()                                        ######
######                                                                               ######
######   Version: 09.02.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######                                                                               ######
######   output:                                                                     ######
######                                                                               ######
######################################################################################## */
CGaugeIndices::~CGaugeIndices()
{
}



/* ########################################################################################
######   bool CGaugeIndices::GetQuadraticIndex(...) const                            ######
######                                                                               ######
######   Version: 09.02.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######                                                                               ######
######   output:                                                                     ######
######                                                                               ######
######################################################################################## */
bool CGaugeIndices::GetQuadraticIndex(const gaugeGroupFactor<double> &ggf, const unsigned &Dimension, unsigned &Index) const
{
  const string &algebra = ggf.algebra;

  unsigned gg_type = 0;
  if (algebra[0] == 'A')
    gg_type = 0;
  else
  if (algebra[0] == 'D')
    gg_type = 1;
  else
  if (algebra[0] == 'E')
    gg_type = 2;

  Index = this->QuadraticIndices[gg_type][ggf.rank][Dimension];
  if (Index == -1)
  {
    cout << "Warning in bool CGaugeIndices::GetQuadraticIndex(...) const: Index of rep " << Dimension << " of " << ggf.algebra << " is not known. Return false." << endl;
    return false;
  }
  return true;
}



/* ########################################################################################
######   bool CGaugeIndices::GetQuadraticIndexAdj(...) const                         ######
######                                                                               ######
######   Version: 09.02.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######                                                                               ######
######   output:                                                                     ######
######                                                                               ######
######################################################################################## */
bool CGaugeIndices::GetQuadraticIndexAdj(const gaugeGroupFactor<double> &ggf, unsigned &Index) const
{
  const string &algebra = ggf.algebra;

  unsigned gg_type = 0;
  if (algebra[0] == 'A')
    gg_type = 0;
  else
  if (algebra[0] == 'D')
    gg_type = 1;
  else
  if (algebra[0] == 'E')
    gg_type = 2;

  Index = this->QuadraticIndicesAdj[gg_type][ggf.rank];
  if (Index == -1)
  {
    cout << "Warning in bool CGaugeIndices::GetQuadraticIndexAdj(...) const: Index of adjoint of " << ggf.algebra << " is not known. Return false." << endl;
    return false;
  }
  return true;
}



/* ##########################################################################
######   CGaugeInvariance::CGaugeInvariance()                          ######
######                                                                 ######
######   Version: 17.01.2007                                           ######
########################################################################## */
CGaugeInvariance::CGaugeInvariance()
{
}


/* ##########################################################################
######   CGaugeInvariance::~CGaugeInvariance()                         ######
######                                                                 ######
######   Version: 17.01.2007                                           ######
########################################################################## */
CGaugeInvariance::~CGaugeInvariance()
{
}



/* ##########################################################################
######   HW_criterion(const vector<int> &HW1, const vector<int> &HW2)  ######
######                                                                 ######
######   Version: 16.01.2007                                           ######
########################################################################## */
bool HW_criterion(const vector<int> &HW1, const vector<int> &HW2)
{
  const size_t s1 = HW1.size();
  if (s1 != HW2.size())
  {
    cout << "Error in bool HWcriterion(const vector<int> &HW1, const vector<int> &HW2): Return false." << endl;
    return false;
  }
  for (unsigned i = 0; i < s1; ++i)
  {
    if (HW1[i] > HW2[i])
      return true;
    if (HW1[i] < HW2[i])
      return false;
  }
  
  return false;
}



/* ########################################################################################
######   bool CGaugeInvariance::CheckInvariance(...)                                 ######
######                                                                               ######
######   Version: 05.04.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######                                                                               ######
######   output:                                                                     ######
######                                                                               ######
######################################################################################## */
bool CGaugeInvariance::CheckInvariance(const COrbifold &Orbifold, const SConfig &VEVConfig, const vector<unsigned> &FieldCoupling)
{
  const size_t OriginalOrder = FieldCoupling.size();
  const vector<CField> &Fields = VEVConfig.Fields;

  const vector<gaugeGroupFactor<double> > &factors = VEVConfig.SymmetryGroup.GaugeGroup.factor;
  const size_t number_of_factors = factors.size();
  
  size_t Order = 0;
  unsigned j = 0;
  unsigned k = 0;
  unsigned algebra_no = 0;
  vector<vector<int> > Coupling_HWs;
  vector<unsigned>     NonTrivialFieldCoupling;

  // check invariance for each group i
  for (unsigned i = 0; i < number_of_factors; ++i)
  {
    const string &Algebra = factors[i].algebra;

    // if allowed couplings for this gauge group factor have not been loaded before
    vector<string>::iterator pos = find(this->Algebras.begin(), this->Algebras.end(), Algebra);
    if (pos == this->Algebras.end())
    {
      cout << "Warning in bool CGaugeInvariance::CheckInvariance(...): Algebra not known. return false." << endl;
      return false;
    }

    algebra_no = distance(this->Algebras.begin(), pos);

    // begin: get the highest weights (ignoring singlets) corresponding to the i-th gauge group factor
    Coupling_HWs.clear();
    NonTrivialFieldCoupling.clear();
    for (j = 0; j < OriginalOrder; ++j)
    {
      const vector<int> &HW = Fields[FieldCoupling[j]].HighestWeights_DL[i];
      if (count(HW.begin(), HW.end(), 0) != HW.size())
      {
        Coupling_HWs.push_back(HW);
        NonTrivialFieldCoupling.push_back(FieldCoupling[j]);
      }
    }
    Order = Coupling_HWs.size();
    // end: get the highest weights (ignoring singlets) corresponding to the i-th gauge group factor
    
    // if tensor product is not the trivial product
    if (Order != 0)
    {
      // unique sorting of the highest weights in "Coupling_HWs"
      stable_sort(Coupling_HWs.begin(), Coupling_HWs.end(), HW_criterion);

      const vector<vector<vector<int> > > &current_Invariants = this->AllInvariantTensorProducts[algebra_no];

      // is current coupling NOT in the list of gauge invariant couplings?
      if (find(current_Invariants.begin(), current_Invariants.end(), Coupling_HWs) == current_Invariants.end())
      {
        const vector<vector<vector<int> > > &current_NonInvariants = this->AllNonInvariantTensorProducts[algebra_no];

        // is current coupling in the list of non gauge invariant couplings?
        if (find(current_NonInvariants.begin(), current_NonInvariants.end(), Coupling_HWs) != current_NonInvariants.end())
          return false;
        else
        {
          bool CurrentCaseSupported = true;
          bool invariant = this->CheckUnknownInvarianceEasyCases(factors[i], Coupling_HWs, CurrentCaseSupported);

          if (!CurrentCaseSupported)
            invariant = this->CheckUnknownInvariance(Orbifold, VEVConfig, NonTrivialFieldCoupling, i);

          if (!invariant)
          {
            // coupling is NOT gauge invariant and was not known before
            // save new coupling
            this->AllNonInvariantTensorProducts[algebra_no].push_back(Coupling_HWs);

            string nongi_filename = this->Filenames.at(algebra_no);
            nongi_filename.insert(16, "non_");
            
            std::ofstream app_nongi_file;
            app_nongi_file.open(nongi_filename.data(), ofstream::app);

            // begin: print result to file and to screen
            cout << "\nWarning: Forbidden coupling: ( ";

            for (j = 0; j < Order-1; ++j)
            {
              const vector<int> &Coupling_HW = Coupling_HWs[j];
              for (k = 0; k < Coupling_HW.size(); ++k)
              {
                cout << Coupling_HW[k] << " ";
                app_nongi_file << Coupling_HW[k] << " ";
              }
              cout << "x";
              app_nongi_file << "x";
            }
            const vector<int> &Coupling_HW = Coupling_HWs[Order-1];
            for (k = 0; k < Coupling_HW.size(); ++k)
            {
              cout << Coupling_HW[k] << " ";
              app_nongi_file << Coupling_HW[k] << " ";
            }
            app_nongi_file << endl;
            cout << ") of gauge group " << Algebra << " saved to file " << nongi_filename << "." << endl;
            // end: print result to file and to screen

            app_nongi_file.close();

            return false;
          }
          else
          {
            // coupling is gauge invariant and was not known before
            // save new coupling
            this->AllInvariantTensorProducts[algebra_no].push_back(Coupling_HWs);

            string gi_filename = this->Filenames.at(algebra_no);

            std::ofstream app_gi_file;
            app_gi_file.open(gi_filename.data(), ofstream::app);

            // begin: print result to file and to screen
            cout << "\nWarning: Allowed coupling: ( ";

            for (j = 0; j < Order-1; ++j)
            {
              const vector<int> &Coupling_HW = Coupling_HWs[j];
              for (k = 0; k < Coupling_HW.size(); ++k)
              {
                cout << Coupling_HW[k] << " ";
                app_gi_file << Coupling_HW[k] << " ";
              }
              cout << "x";
              app_gi_file << "x";
            }
            const vector<int> &Coupling_HW = Coupling_HWs[Order-1];
            for (k = 0; k < Coupling_HW.size(); ++k)
            {
              cout << Coupling_HW[k] << " ";
              app_gi_file << Coupling_HW[k] << " ";
            }
            app_gi_file << endl;
            cout << ") of gauge group " << Algebra << " saved to file " << gi_filename << "." << endl;
            // end: print result to file and to screen

            app_gi_file.close();
          }
        }
      }
    }
  }
  return true;
}



/* ########################################################################################
######   bool CGaugeInvariance::CheckUnknownInvarianceEasyCases(...) const           ######
######                                                                               ######
######   Version: 01.04.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######                                                                               ######
######   output:                                                                     ######
######                                                                               ######
######################################################################################## */
bool CGaugeInvariance::CheckUnknownInvarianceEasyCases(const gaugeGroupFactor<double> &GaugeGroupFactor, const vector<vector<int> > &Coupling_HWs, bool &CurrentCaseSupported) const
{
  const size_t s1 = Coupling_HWs.size();
  unsigned i = 0;

  if (GaugeGroupFactor.algebra[0] == 'A')
  {
    if (GaugeGroupFactor.rank == 1)
    {
      CurrentCaseSupported = true;
      vector<int> Doublet(1,1);

      unsigned number_of_doublets = 0;
      for (i = 0; i < s1; ++i)
      {
        if (Coupling_HWs[i] == Doublet)
          ++number_of_doublets;
      }
      if (number_of_doublets % 2 == 0)
      {
        //cout << "number of doublets is even: " << number_of_doublets << endl;
        return true;
      }
      else
      {
        //cout << "number of doublets is odd: " << number_of_doublets << endl;
        return false;
      }
    }
    else
    {
      unsigned number_of_fundamental = 0;
      unsigned number_of_antifundamental = 0;

      const vector<int> Singlet(GaugeGroupFactor.rank, 0);
      vector<int> Fundamental(GaugeGroupFactor.rank, 0);
      vector<int> Antifundamental(GaugeGroupFactor.rank, 0);
      Fundamental[0] = 1;
      Antifundamental[GaugeGroupFactor.rank-1] = 1;

      for (i = 0; i < s1; ++i)
      {
        if (Coupling_HWs[i] == Fundamental)
          ++number_of_fundamental;
        else
        if (Coupling_HWs[i] == Antifundamental)
          ++number_of_antifundamental;
        else
        {
          if (Coupling_HWs[i] != Singlet)
          {
            //cout << "this case is not supported" << endl;
            CurrentCaseSupported = false;
            return false;
          }
        }
      }
      CurrentCaseSupported = true;

      int chiral = number_of_fundamental - number_of_antifundamental;
      double rest = fmod(fabs((double)chiral), (GaugeGroupFactor.rank+1));
      if (fabs(rest) < 0.0001)
      {
        //cout << "rank = " << GaugeGroupFactor.rank << ". Invariant: number of chiral fundamentals: " << chiral << endl;
        return true;
      }
      else
      {
        //cout << "rank = " << GaugeGroupFactor.rank << ". Not invariant: number of chiral fundamentals: " << chiral << endl;
        return false;
      }
    }
  }
  CurrentCaseSupported = false;
  return false;
}



/* ########################################################################################
######   bool CGaugeInvariance::CheckUnknownInvariance(...) const                    ######
######                                                                               ######
######   Version: 17.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######                                                                               ######
######   output:                                                                     ######
######                                                                               ######
######################################################################################## */
bool CGaugeInvariance::CheckUnknownInvariance(const COrbifold &Orbifold, const SConfig &VEVConfig, const vector<unsigned> &FieldCoupling, unsigned ggf) const
{
  bool info = false;
  
  const size_t Order = FieldCoupling.size();
  
  const vector<CField> &Fields = VEVConfig.Fields;
  
  const gaugeGroupFactor<double> &GaugeGroupFactor = Orbifold.StandardConfig.SymmetryGroup.GaugeGroup.factor[ggf];
        
  // begin: define some constants and variables
  unsigned i = 0;
  unsigned j = 0;
  size_t s1 = 0;
  
  intMatrix               Weights_DL;
  intMatrix::iterator     pos;
  intMatrix               different_Weights_DL;
  vector<intMatrix>       all_Weights_DL;
  vector<vector<double> > Weights_CW_Copy;
  vector<unsigned>        Number_Of_Weights(Order, 0);
  // end: define some constants and variables

  if (info)
    cout << "Check Invariance of:\n";
  // run through the coupling

  for (i = 0; i < Order; ++i)
  {
    const CField &Field = Fields[FieldCoupling[i]];
    
    if (info)
    {
      CPrint Print(Tstandard, &cout);
      Print.PrintRep(Field.Dimensions);
      cout << " ";
    }

    s1 = Field.GetNumberOfLMWeights();
    Weights_CW_Copy.clear();
    for (j = 0; j < s1; ++j)
      Weights_CW_Copy.push_back(Field.GetLMWeight(j, Orbifold.GetSectors()));

    // save only the different weights in Dynkin labels
    different_Weights_DL.clear();
    Weights_DL = findDynkinLabels(GaugeGroupFactor, Weights_CW_Copy);
    
    for (j = 0; j < s1; ++j)
    {
      const intVector &Weight_DL = Weights_DL[j];

      pos = find(different_Weights_DL.begin(), different_Weights_DL.end(), Weight_DL);
      if (pos == different_Weights_DL.end())
        different_Weights_DL.push_back(Weight_DL);
    }
    
    Number_Of_Weights[i] = different_Weights_DL.size();
    all_Weights_DL.push_back(different_Weights_DL);
  }
        
  if (info)
  {
    cout << endl;
    string tmp;
    cin >> tmp;
    if (tmp == "y")
    {
      cout << "gauge invariant" << endl;
      return true;
    }
    if (tmp == "n")
    {
      cout << "not gauge invariant" << endl;
      return false;
    }
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // begin: open files for saving the combinations
  const string gi_combi_filename = "GaugeInvariance/tmp_gi_combinations";

  vector<std::ofstream*> gi_Combinations;
  std::ofstream gi_Combination01, gi_Combination02, gi_Combination03;
  std::ofstream gi_Combination04, gi_Combination05, gi_Combination06;
  std::ofstream gi_Combination07, gi_Combination08, gi_Combination09;
  std::ofstream gi_Combination10, gi_Combination11, gi_Combination12;
  std::ofstream gi_Combination13, gi_Combination14, gi_Combination15;
  std::ofstream gi_Combination16, gi_Combination17, gi_Combination18;
  std::ofstream gi_Combination19, gi_Combination20, gi_Combination21;
  std::ofstream gi_Combination22, gi_Combination23, gi_Combination24;
  std::ofstream gi_Combination25, gi_Combination26, gi_Combination27;
  std::ofstream gi_Combination28, gi_Combination29, gi_Combination30;
  std::ofstream gi_Combination31, gi_Combination32, gi_Combination33;
  std::ofstream gi_Combination34, gi_Combination35, gi_Combination36;
  std::ofstream gi_Combination37, gi_Combination38, gi_Combination39;
  std::ofstream gi_Combination40, gi_Combination41, gi_Combination42;

  gi_Combinations.push_back(&gi_Combination01);   gi_Combinations.push_back(&gi_Combination02);
  gi_Combinations.push_back(&gi_Combination03);   gi_Combinations.push_back(&gi_Combination04);
  gi_Combinations.push_back(&gi_Combination05);   gi_Combinations.push_back(&gi_Combination06);
  gi_Combinations.push_back(&gi_Combination07);   gi_Combinations.push_back(&gi_Combination08);
  gi_Combinations.push_back(&gi_Combination09);   gi_Combinations.push_back(&gi_Combination10);
  gi_Combinations.push_back(&gi_Combination11);   gi_Combinations.push_back(&gi_Combination12);
  gi_Combinations.push_back(&gi_Combination13);   gi_Combinations.push_back(&gi_Combination14);
  gi_Combinations.push_back(&gi_Combination15);   gi_Combinations.push_back(&gi_Combination16);
  gi_Combinations.push_back(&gi_Combination17);   gi_Combinations.push_back(&gi_Combination18);
  gi_Combinations.push_back(&gi_Combination19);   gi_Combinations.push_back(&gi_Combination20);
  gi_Combinations.push_back(&gi_Combination21);   gi_Combinations.push_back(&gi_Combination22);
  gi_Combinations.push_back(&gi_Combination23);   gi_Combinations.push_back(&gi_Combination24);
  gi_Combinations.push_back(&gi_Combination25);   gi_Combinations.push_back(&gi_Combination26);
  gi_Combinations.push_back(&gi_Combination27);   gi_Combinations.push_back(&gi_Combination28);
  gi_Combinations.push_back(&gi_Combination29);   gi_Combinations.push_back(&gi_Combination30);
  gi_Combinations.push_back(&gi_Combination31);   gi_Combinations.push_back(&gi_Combination32);
  gi_Combinations.push_back(&gi_Combination33);   gi_Combinations.push_back(&gi_Combination34);
  gi_Combinations.push_back(&gi_Combination35);   gi_Combinations.push_back(&gi_Combination36);
  gi_Combinations.push_back(&gi_Combination37);   gi_Combinations.push_back(&gi_Combination38);
  gi_Combinations.push_back(&gi_Combination39);   gi_Combinations.push_back(&gi_Combination40);
  gi_Combinations.push_back(&gi_Combination41);   gi_Combinations.push_back(&gi_Combination42);

  const size_t number_of_gi = gi_Combinations.size();
  
  string tmp_filename = "";
  for (i = 0; i < number_of_gi; ++i)
  {
    tmp_filename = gi_combi_filename;

    std::ostringstream os1;
    os1 << i;

    tmp_filename += os1.str();
    tmp_filename += ".txt";
    gi_Combinations[i]->open(tmp_filename.data(), ofstream::out);
  }
  // end: open files for saving the combinations
  ////////////////////////////////////////////////////////////////////////////////////////////////////////

  // begin: create the combinations
  vector<vector<unsigned> > vec_adjoining_Positions;
  vector<unsigned>          tmp_Combination(Order, 0);
  unsigned long             counter = 0;
  RecursiveCounting_to_files(tmp_Combination, 0, Number_Of_Weights[0], Number_Of_Weights, gi_Combinations, counter);
  // end: create the combinations

  // close the files now in order to open them later read-only
  for (i = 0; i < number_of_gi; ++i)
    gi_Combinations[i]->close();

    
  // InvariantCombinations:
  // i-th line corresponds to the i-th representation
  // j-th column corresponds to the j-th weight of this representation
  // InvariantCombinations[i][j] = false  ->  combination of weights does not give a singlet
  // InvariantCombinations[i][j] = true   ->  combination of weights gives a singlet
  //
  // in order to be invariant, InvariantCombinations[i][j] must be true for all i,j
  vector<vector<bool> > InvariantCombinations;
  {
    vector<bool> tmp_InvariantCombinations;
      
    for (i = 0; i < Order; ++i)
    {
      tmp_InvariantCombinations.assign(Number_Of_Weights[i], false);
      InvariantCombinations.push_back(tmp_InvariantCombinations);
    }
  }

  const CVector NullVector_x(all_Weights_DL[0][0].size());
  
  CVector Sum_Weights;
  CVector weight;

  string tmp_string = "";
  unsigned temp = 0 ;
  vector<unsigned> Combination;
  
  // begin: run through all files containing the combinations
  for (i = 0; i < number_of_gi; ++i)
  {
    // begin: open the file
    std::ostringstream os1;
    os1 << i;

    tmp_filename = gi_combi_filename;
    tmp_filename += os1.str();
    tmp_filename += ".txt";

    std::ifstream in_gi_Combinations;
    in_gi_Combinations.open(tmp_filename.data(), ifstream::in);
    // end: open the file

    tmp_string = "";
   
    // read the lines of this file
    while (getline(in_gi_Combinations, tmp_string) )
    {
      // begin: convert the line to a combination
      Combination.clear();
      temp = 0;

      std::istringstream line(tmp_string);
      while (line >> temp)
        Combination.push_back(temp);
      // end: convert the line to a combination

      // begin: add the weights
      Sum_Weights = NullVector_x;
      for (j = 0; j < Order; ++j)
      {
        weight = all_Weights_DL[j][Combination[j]];

        Sum_Weights = Sum_Weights + weight;
      }
      // end: add the weights
      
      // begin: check invariance
      if (Sum_Weights == NullVector_x)
      {
        for (j = 0; j < Order; ++j)
          InvariantCombinations[j][Combination[j]] = true;
      }
      // end: check invariance
    }
    in_gi_Combinations.close();

    // delete the content of the file
    std::ofstream gi_Combinations;
    gi_Combinations.open(tmp_filename.data(), ofstream::out);
    gi_Combinations.close();
  }
  // end: run through all files containing the combinations

  size_t s2 = 0;
  for (i = 0; i < Order; ++i)
  {
    const vector<bool> &tmp_InvariantCombinations = InvariantCombinations[i];
    s2 = tmp_InvariantCombinations.size();
    for (j = 0; j < s2; ++j)
    {
      if (!tmp_InvariantCombinations[j])
        return false;
    }
  }
  return true;
}



/* ########################################################################################
######   bool CGaugeInvariance::LoadTensorProducts(...)                              ######
######                                                                               ######
######   Version: 01.04.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######                                                                               ######
######   output:                                                                     ######
######                                                                               ######
######################################################################################## */
bool CGaugeInvariance::LoadTensorProducts(const CGaugeGroup &GaugeGroup)
{
  this->Algebras.clear();
  this->AllInvariantTensorProducts.clear();
  this->AllNonInvariantTensorProducts.clear();
  this->Filenames.clear();

  const vector<gaugeGroupFactor<double> > &factors = GaugeGroup.factor;
  const size_t number_of_factors        = factors.size();

  bool new_file = false;
  
  string gi_filename = "";
  string nongi_filename = "";
  string currentline = "";
  vector<vector<vector<int> > > InvariantTensorProducts;
  vector<vector<vector<int> > > NonInvariantTensorProducts;
  vector<vector<int> >          tmp_InvariantTensorProduct;
  vector<vector<int> >          tmp_NonInvariantTensorProduct;
  
  // run through all gauge group factors
  for (unsigned i = 0; i < number_of_factors; ++i)
  {
    const string &Algebra = factors[i].algebra;

    // if allowed couplings for this gauge group factor have not been loaded before
    if (find(this->Algebras.begin(), this->Algebras.end(), Algebra) == this->Algebras.end())
    {
      // check, whether a new file needs to be created
      new_file = false;

      // begin: create the filenames
      gi_filename = "gauge_invariance_";
      gi_filename += Algebra;
      gi_filename += ".txt";

      nongi_filename = "GaugeInvariance/non_";
      nongi_filename += gi_filename;
      gi_filename = "GaugeInvariance/" + gi_filename;
      // end: create the filenames

      // begin: gauge invariance
      std::ifstream input_file1;
      input_file1.open(gi_filename.data(), ifstream::in);

      // begin: create new file
      if((!input_file1.is_open()) || (!input_file1.good()))
      {
        #ifdef DEBUG
          cout << "Warning! Gauge invariance written to file: " << gi_filename << endl;
        #endif
        input_file1.close();
      
        std::ofstream app_gi_file;
        app_gi_file.open(gi_filename.data(), ofstream::out);
        app_gi_file << Algebra << " gauge invariance" << endl;
        app_gi_file.close();
      }
      // end: create new file
      
      // begin: load existing file
      else
      {
        #ifdef DEBUG
          cout << "Gauge invariance loaded from file: " << gi_filename << endl;
        #endif

        // for the line "...gauge invariance"
        getline(input_file1, currentline);

        int         Dim = 0;
        vector<int> HW;

        InvariantTensorProducts.clear();

        while (getline(input_file1, currentline))
        {
          tmp_InvariantTensorProduct.clear();

          string::size_type curr_pos = 0;
          string::size_type to_pos = 0;
          string tmp_string = "";

          while (curr_pos != string::npos)
          {
            to_pos = currentline.find("x", curr_pos);

            if (to_pos == string::npos)
              tmp_string = currentline.substr(curr_pos, string::npos);
            else
              tmp_string = currentline.substr(curr_pos, to_pos - curr_pos);

            HW.clear();
            std::istringstream Line_Sectors(tmp_string);
            while (Line_Sectors >> Dim)
              HW.push_back(Dim);

            tmp_InvariantTensorProduct.push_back(HW);

            if (to_pos == string::npos)
              curr_pos = to_pos;
            else
              curr_pos = to_pos+1;
          }

          InvariantTensorProducts.push_back(tmp_InvariantTensorProduct);
        }
        input_file1.close();
      }
      // end: load existing file
      // end: gauge invariance

      // begin: non gauge invariance
      std::ifstream input_file2;
      input_file2.open(nongi_filename.data(), ifstream::in);

      // begin: create new file
      if((!input_file2.is_open()) || (!input_file2.good()))
      {
        #ifdef DEBUG
          cout << "Warning! Non gauge invariance written to file: " << nongi_filename << endl;
        #endif
        input_file2.close();
      
        std::ofstream app_gi_file;
        app_gi_file.open(nongi_filename.data(), ofstream::out);
        app_gi_file << Algebra << " non gauge invariance" << endl;
        app_gi_file.close();
      }
      // end: create new file
      
      // begin: load existing file
      else
      {
        #ifdef DEBUG
          cout << "Non gauge invariance loaded from file: " << nongi_filename << endl;
        #endif

        // for the line "...gauge invariance"
        getline(input_file2, currentline);

        int         Dim = 0;
        vector<int> HW;

        NonInvariantTensorProducts.clear();

        while (getline(input_file2, currentline))
        {
          tmp_NonInvariantTensorProduct.clear();

          string::size_type curr_pos = 0;
          string::size_type to_pos = 0;
          string tmp_string = "";

          while (curr_pos != string::npos)
          {
            to_pos = currentline.find("x", curr_pos);

            if (to_pos == string::npos)
              tmp_string = currentline.substr(curr_pos, string::npos);
            else
              tmp_string = currentline.substr(curr_pos, to_pos - curr_pos);

            HW.clear();
            std::istringstream Line_Sectors(tmp_string);
            while (Line_Sectors >> Dim)
              HW.push_back(Dim);

            tmp_NonInvariantTensorProduct.push_back(HW);

            if (to_pos == string::npos)
              curr_pos = to_pos;
            else
              curr_pos = to_pos+1;
          }

          NonInvariantTensorProducts.push_back(tmp_NonInvariantTensorProduct);
        }
        input_file2.close();
      }
      // end: load existing file
      // end: gauge invariance
      
      this->Algebras.push_back(Algebra);
      this->Filenames.push_back(gi_filename);
      this->AllInvariantTensorProducts.push_back(InvariantTensorProducts);
      this->AllNonInvariantTensorProducts.push_back(NonInvariantTensorProducts);
    }
  }
  return true;
}



/* ##########################################################################
######   bool DecomposeIntoIrreps(...)                                 ######
######                                                                 ######
######   Version: 01.03.2007                                           ######
######                                                                 ######
######   Does not work!!!                                              ######
########################################################################## */
bool CGaugeInvariance::DecomposeIntoIrreps(vector<vector<int> > &DL_weights, const gaugeGroupFactor<double> &ggf, vector<vector<vector<int> > > &result_weights, vector<int> &result_dimensions) const
{
  /*
  cout << "weights " << DL_weights.size() << endl;
  
  for (unsigned i = 0; i < DL_weights.size(); ++i)
  {
    for (unsigned j = 0; j < DL_weights[i].size(); ++j)
      cout << DL_weights[i][j] << " ";
    cout << endl;
  }
  cout << endl;
  */
  const size_t s1 = DL_weights.size();
  if (s1 == 0)
    return true;
  
  const size_t s2 = DL_weights[0].size();
  /*
  CGaugeGroup TmpGroup;
  TmpGroup.factor.push_back(ggf);
  
  cirrep<double> tmp = findHighestWeight<double>(TmpGroup, vector<vector<T> > myWeights);
  */
  vector<vector<int> > HighestWeights_DL = findHighestWeights(DL_weights);
  
  if (HighestWeights_DL.size() == 0)
  {
    cout << "\nError in bool CGaugeInvariance::DecomposeIntoIrreps(...) const: Highest weight not found. Return false." << endl;
    return false;
  }
  //cout << "Highest weights" << endl;
  //for (unsigned j = 0; j < HighestWeights_DL[0].size(); ++j)
  //  cout << HighestWeights_DL[0][j] << " ";
  //cout << endl;
  
  //cout << "ggf.algebra " << ggf.algebra.substr(0,1) << endl;
  
   
  freudenthal Freud(ggf.simpleroots, highestRoot(ggf.algebra.substr(0,1), ggf.rank), HighestWeights_DL[0]);
  
  map<intVector, int> WeightsFromFreudenthal = Freud.getWeightMap();
  
  CState tmp;
  SDimension dim;
  dim.Dimension = 1;
 // tmp.DetermineDimension(ggf, HighestWeights_DL[0], dim);
  int dimension = dim.Dimension;
  //cout << HighestWeights_DL[0] << " dim = " << dimension << endl;
  //Freud.print();
  
  vector<intVector> WeightsFromFreudenthalCopy;
  
  unsigned i = 0;
  bool equal = true;
  unsigned mult = 1;
  unsigned counter = 0;
  vector<int> FreudWeight;
  
  for( map<intVector, int>::iterator FreudWeightMap = WeightsFromFreudenthal.begin();
       FreudWeightMap != WeightsFromFreudenthal.end(); ++FreudWeightMap)
  {
    FreudWeight = FreudWeightMap->first;
    mult        = FreudWeightMap->second;
    
    for (i = 0; i < mult; ++i)
      WeightsFromFreudenthalCopy.push_back(FreudWeight);
  
    counter = 0;
    
    for( vector<vector<int> >::iterator DL_weight = DL_weights.begin();
         ((DL_weight != DL_weights.end()) && (counter != mult)); ++DL_weight)
    {
      equal = true;
      for (i = 0; equal && (i < s2); ++i)
      {
        if (DL_weight->at(i) != FreudWeight.at(i))
          equal = false;
      }

      if (equal)
      {
        ++counter;
        DL_weights.erase(DL_weight--);
      }
    }
    if (counter != mult)
    {
      cout << "Error in bool CGaugeInvariance::DecomposeIntoIrreps(...) const: Weight from Freudenthal algorithm not found." << endl;
      return false;
    }
  }  
  
  result_weights.push_back(WeightsFromFreudenthalCopy);
  result_dimensions.push_back(dimension);
  this->DecomposeIntoIrreps(DL_weights, ggf, result_weights, result_dimensions);
  
  return true;
}

  
/* ##########################################################################
######   bool TensorProduct(...)                                       ######
######                                                                 ######
######   Version: 05.02.2007                                           ######
########################################################################## */
bool CGaugeInvariance::TensorProduct(const vector<vector<int> > &DL_weights1, const vector<vector<int> > &DL_weights2, const gaugeGroupFactor<double> &ggf, vector<vector<vector<int> > > &Decomposed_DL_Weights) const
{
  const size_t s1 = DL_weights1.size();
  const size_t s2 = DL_weights2.size();
  
  if ((s1 == 0) || (s2 == 0))
  {
    cout << "Warning in bool CGaugeInvariance::TensorProduct(...) const: DL_weights1 or DL_weights2 is empty!" << endl;
    return false;
  }
  const size_t s3 = DL_weights1[0].size();
  if (s3 != DL_weights2[0].size())
  {
    cout << "Error in bool CGaugeInvariance::TensorProduct(...) const: Weights have different lengths!" << endl;
    return false;
  }
  
  if (Decomposed_DL_Weights.size() != 0)
  {
    cout << "Warning in bool CGaugeInvariance::TensorProduct(...) const: Vectors Decomposed_DL_Weights is not empty! Now cleared!" << endl;
    Decomposed_DL_Weights.clear();
  }
  
  unsigned i = 0;
  unsigned j = 0;
  unsigned k = 0;
  vector<int> sum(s3,0);
  
  vector<vector<int> > weightspace;
  
  for (i = 0; i < s1; ++i)
  {
    const vector<int> &DL_weight1 = DL_weights1[i];
    for (j = 0; j < s2; ++j)
    {
      const vector<int> &DL_weight2 = DL_weights2[j];
        
      sum.assign(s3,0);
      for (k = 0; k < s3; ++k)
        sum[k] = DL_weight1[k] + DL_weight2[k];
      
      weightspace.push_back(sum);
    }
  }
  vector<int> result_dimensions;
  return this->DecomposeIntoIrreps(weightspace, ggf, Decomposed_DL_Weights, result_dimensions);
}



/* ##########################################################################
######   bool TensorProductWithItself(...)                             ######
######                                                                 ######
######   Version: 02.02.2007                                           ######
########################################################################## */
bool CGaugeInvariance::TensorProductWithItself(const vector<vector<int> > &DL_weights, const gaugeGroupFactor<double> &ggf, vector<vector<vector<int> > > &Decomposed_DL_Weights_sym, vector<vector<vector<int> > > &Decomposed_DL_Weights_asym) const
{
  const size_t s1 = DL_weights.size();
  if (s1 == 0)
  {
    cout << "Warning in bool CGaugeInvariance::TensorProductWithItself(...) const: DL_weights is empty!" << endl;
    return false;
  }
  const size_t s2 = DL_weights[0].size();
  
  if ((Decomposed_DL_Weights_sym.size() != 0) || (Decomposed_DL_Weights_asym.size() != 0))
  {
    cout << "Warning in bool CGaugeInvariance::TensorProductWithItself(...) const: Vectors Decomposed_DL_Weights_sym and Decomposed_DL_Weights_asym are not empty! Now cleared!" << endl;
    Decomposed_DL_Weights_sym.clear();
    Decomposed_DL_Weights_asym.clear();
  }
  
  unsigned i = 0;
  unsigned j = 0;
  unsigned k = 0;
  vector<int> sum(s2,0);
  
  vector<vector<int> > sym_weightspace;
  vector<vector<int> > antisym_weightspace;
  
  for (i = 0; i < s1; ++i)
  {
    const vector<int> &DL_weight1 = DL_weights[i];
    for (j = i; j < s1; ++j)
    {
      const vector<int> &DL_weight2 = DL_weights[j];
        
      sum.assign(s2,0);
      for (k = 0; k < s2; ++k)
        sum[k] = DL_weight1[k] + DL_weight2[k];
      
      // for j == i and j != i, include the sum to the symmetric weight space
      sym_weightspace.push_back(sum);
      
      // only if j != i, include the sum to the anti-symmetric weight space
      if (j != i)
        antisym_weightspace.push_back(sum);
    }
  }
  vector<int> result_dimensions_sym;
  vector<int> result_dimensions_asym;
  return (
    this->DecomposeIntoIrreps(sym_weightspace, ggf, Decomposed_DL_Weights_sym, result_dimensions_sym) 
 && this->DecomposeIntoIrreps(antisym_weightspace, ggf, Decomposed_DL_Weights_asym, result_dimensions_asym));
}
