
#include <stdio.h>

#include "cprompt.h"



// Mis inclusiones
#include <iostream>
#include <cstring>
#include <readline/readline.h>





using namespace std;

unsigned SELFDUALLATTICE;

vector<int>              current_folder_global;








int main(int argc, char *argv[])
{
/*  SELFDUALLATTICE = 1;
  cout << "Starting..." << endl;
  std::ifstream in("modelBl.txt");
  if((!in.is_open()) || (!in.good()))
  {
    cout << "\n  File \"modelBl.txt\" not found.\n" << endl;
    exit(1);
  }

  string ProgramFilename = "";
  unsigned i = 0;
  unsigned j = 0;
  unsigned FP = 0;
  bool use_freelyWL = false;
  cout << "argc = " << argc << endl;
  if (argc == 3)
  {
    string tmp = "";

    cout << "\n  use local shift: ";
    tmp = argv[1];
    if (tmp == "0")
    {
      FP = 0;
      cout << "-";
    }
    else
    if (tmp == "1")
    {
      FP = 1;
      cout << "V_1 ";
    }
    else
    if (tmp == "2")
    {
      FP = 2;
      cout << "V_1 + W_5";
    }
    else
    if (tmp == "3")
    {
      FP = 3;
      cout << "V_1 + W_6";
    }
    else
    if (tmp == "4")
    {
      FP = 4;
      cout << "V_1 + W_5 + W_6";
    }

    tmp = argv[2];
    if (tmp == "1")
    {
      cout << " and freely acting WL";
      use_freelyWL = true;
    }
    cout << endl;
    COrbifoldGroup OrbifoldGroup;
    if (OrbifoldGroup.LoadOrbifoldGroup(in, ProgramFilename))
    {
      // begin: create the E8 x E8' or SO(32) gauge group
      const CVector Null_Shift(16);
      S_OscillatorExcitation Excitation;
      Excitation.NumberOperator = 0.0;
      Excitation.ZeroPointEnergy = 2.0;
      CMasslessHalfState Roots10D(LeftMover, Excitation);
      Roots10D.SolveMassEquation(Null_Shift, E8xE8);
      // end: create the E8 x E8' or SO(32) gauge group

      CPrint Print(Tstandard, &cout);

      bool invariant = true;
      vector<vector<double> > Roots;
      for (i = 0; i < 480; ++i)
      {
        const CVector &Root = Roots10D.Weights[i];
        invariant = true;
      
        if (invariant)
        {
          CVector tmp = OrbifoldGroup.GetShift(0) + OrbifoldGroup.GetShift(1);
          if (!is_integer(tmp * Root))
            invariant = false;
        }

        for (j = 0; invariant && (j < 4); ++j)
        {
          if (!is_integer(OrbifoldGroup.GetWilsonLines().GetWilsonLine(j) * Root))
            invariant = false;
        }

        CVector tmp = OrbifoldGroup.GetWilsonLines().GetWilsonLine(1) * 1.5;
        if (invariant && use_freelyWL && !is_integer(tmp * Root))
          invariant = false;

        if (invariant)
        {
          CVector tmp(16);
          if (FP == 1)
            tmp = OrbifoldGroup.GetShift(0);
          if (FP == 2)
            tmp = OrbifoldGroup.GetShift(0) + OrbifoldGroup.GetWilsonLines().GetWilsonLine(4);
          if (FP == 3)
            tmp = OrbifoldGroup.GetShift(0) + OrbifoldGroup.GetWilsonLines().GetWilsonLine(5);
          if (FP == 4)
            tmp = OrbifoldGroup.GetShift(0) + OrbifoldGroup.GetWilsonLines().GetWilsonLine(4) + OrbifoldGroup.GetWilsonLines().GetWilsonLine(5);
          
          if (!is_integer(tmp * Root))
            invariant = false;
        }

        if (invariant)
          Roots.push_back(Root);
      }

      CGaugeGroup GaugeGroup = determineAlgebra(Roots);
    
      const vector<gaugeGroupFactor<double> > copy_factors = GaugeGroup.factor;
      vector<vector<double> > SimpleRoots;
      const double prec = 0.0001;

      // collect all simple roots
      size_t s1 = copy_factors.size();
      for (i = 0; i < s1; ++i)
      {
        const vector<vector<double> > &ggf_SimpleRoots = copy_factors[i].simpleroots;
        SimpleRoots.insert(SimpleRoots.end(), ggf_SimpleRoots.begin(), ggf_SimpleRoots.end());
      }

      GaugeGroup.factor.clear();

      unsigned PosAnd = 0;

      bool only_from_first_E8  = true;
      bool only_from_second_E8 = true;

      // from first E_8
      for (i = 0; i < s1; ++i)
      {
        const gaugeGroupFactor<double> &ggf            = copy_factors[i];
        const vector<double>           &ggf_SimpleRoot = ggf.simpleroots[0];

        only_from_first_E8 = true;
        for (j = 8; only_from_first_E8 && (j < 16); ++j)
        {
          if (fabs(ggf_SimpleRoot[j]) > prec)
            only_from_first_E8 = false;
        }
        if (only_from_first_E8)
        {
          GaugeGroup.factor.push_back(ggf);
          ++PosAnd;
        }
      }
      // from second E_8
      for (i = 0; i < s1; ++i)
      {
        const gaugeGroupFactor<double> &ggf            = copy_factors[i];
        const vector<double>           &ggf_SimpleRoot = ggf.simpleroots[0];

        only_from_second_E8 = true;
        for (j = 0; only_from_second_E8 && (j < 8); ++j)
        {
          if (fabs(ggf_SimpleRoot[j]) > prec)
            only_from_second_E8 = false;
        }
        if (only_from_second_E8)
          GaugeGroup.factor.push_back(ggf);
      }
      if (s1 != GaugeGroup.factor.size())
      {
        cout << "\n  Warning in bool COrbifold::FindGaugeGroup(...) const: Cannot sort the gauge group factors. Number of factors before/after sorting: " << s1 << "/" << GaugeGroup.factor.size() << ". Return false." << endl;
        return false;
      }
      // reorder the string "GaugeGroup.algebra" and insert the word "and"
      if ((PosAnd != 0)  && (PosAnd != s1))
      {
        GaugeGroup.algebra = "";
        for (i = 0; i < PosAnd; ++i)
        {
          GaugeGroup.algebra += GaugeGroup.factor[i].algebra;
          if (i + 1 < PosAnd)
            GaugeGroup.algebra += " + ";
        }
        GaugeGroup.algebra += " and ";
        for (i = PosAnd; i < s1; ++i)
        {
          GaugeGroup.algebra += GaugeGroup.factor[i].algebra;
          if (i + 1 < s1)
            GaugeGroup.algebra += " + ";
        }
      }

      cout << "\n  unbroken gauge group " << GaugeGroup.algebra << "\n" << endl;
      Print.PrintSimpleRoots(GaugeGroup, -1);
      cout << endl;
    }
}*/
  cout << "\n  ###########################################################################\n";
  cout << "  #  The C++ Orbifolder                                                     #\n";
  cout << "  #  Version: 1.2                                                           #\n";
  cout << "  #  by H.P. Nilles, S. Ramos-Sanchez, P.K.S. Vaudrevange and A. Wingerter  #\n";
  cout << "  ###########################################################################\n" << endl;




  if (argc == 1)
  {
    CPrompt Prompt;
    Prompt.StartPrompt();
  }
  else
  if (argc == 2)
  {
    const string ModelFilename = argv[1];
    if (ModelFilename == "web")
    {
      cout << "  Initiate web interface...\n" << endl;
      CPrompt Prompt;
      Prompt.StartPrompt("", false, true);
    }
    else
    {
      CPrompt Prompt(ModelFilename.data());
      Prompt.StartPrompt();
    }
  }
  return EXIT_SUCCESS;
}
