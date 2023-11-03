#include <stdio.h>
#include "cprompt.h"
using namespace std;

int main(int argc, char *argv[])
{
  ifstream in("modelKRZ_A1.txt");
  if((!in.is_open()) || (!in.good()))
    exit(1);

  CPrint Print(Tstandard, &cout);
  string ProgramFilename = "";

  COrbifoldGroup OrbifoldGroup;
  if (OrbifoldGroup.LoadOrbifoldGroup(in, ProgramFilename))// load from file
  {
    cout << "\n-> Model file \"modelKRZ_A1.txt\" loaded." << endl;
    COrbifold KRZ_A1(OrbifoldGroup);                       // create the orbifold
    cout << "-> Orbifold \"KRZ_A1\" created.\n" << endl;

    cout << "-> Print shift and Wilson lines:" << endl;
    Print.PrintShift(OrbifoldGroup.GetShift(0));
    Print.PrintWilsonLines(OrbifoldGroup.GetWilsonLines(), true);

    cout << "\n-> Print spectrum, first with and then without U(1) charges:" << endl;
    Print.PrintSummaryOfVEVConfig(KRZ_A1.StandardConfig);  // print with U(1)s
    SConfig VEVConfig = KRZ_A1.StandardConfig;             // create new vev-config.
    VEVConfig.ConfigLabel = "TestConfig";                  // rename new vev-config.
    VEVConfig.SymmetryGroup.observable_sector_U1s.clear(); // change obs. sector
    Print.PrintSummaryOfVEVConfig(VEVConfig);              // print without U(1)s

    cout << "-> Analyze model:" << endl;
    vector<SConfig> AllVEVConfigs;
    bool SM = true;                                        // look for SM
    bool PS = true;                                        // look for PS models
    bool SU5 = true;                                       // look for SU(5) models
    // analyze the configuration "KRZ_A1.StandardConfig" of "KRZ_A1"
    // and save the result to "AllVEVConfigs"
    CAnalyseModel Analyze;
    Analyze.AnalyseModel(KRZ_A1, KRZ_A1.StandardConfig, SM, PS, SU5, 
    AllVEVConfigs, Print);
    if (SM || PS || SU5) // if one of the three possibilities is true 
    {
      cout << "-> Model has 3 generations plus vector-like exotics:" << endl;
      const size_t s1 = AllVEVConfigs.size();              // print all new configs.
      for (unsigned i = 0; i < s1; ++i)
        Print.PrintSummaryOfVEVConfig(AllVEVConfigs[i], LeftChiral, true);
    }
  }
  return EXIT_SUCCESS;
}
