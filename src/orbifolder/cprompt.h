#ifndef CPROMPT_H
#define CPROMPT_H

#include <iostream>
#include <iomanip>
#include <cstdlib>

#include <vector>
#include <boost/rational.hpp>

#include "corbifold.h"
#include "cprint.h"

using std::string;
using boost::rational;

// ConditionTypeId    variable          ComparisonTypeId    equality          ValueTypeId         value
//----------------------------          ----------------------------          ----------------------------
//     0              length                 0              ==                     0              even
//     1              vev                    1              !=                     1              odd
//     2              B-L                    2              >                      2              even/odd
//     3              Q_i                    3              >=                     3              rational number
//     4              acc. Q_i               4              <                      4              string
//     5              rep_i                  5              <=
//     6              p_sh_i                 6              involves
//     7              R_i                    7              !involves
//     8              #osci.
//     9              label
// where i = index
#ifndef STRUCT_SCONDITION
#define STRUCT_SCONDITION
struct SCondition
{
  unsigned ConditionTypeId;
  unsigned ComparisonTypeId;
  unsigned ValueTypeId;

  unsigned      ith_entry;
  rational<int> value;
  string        field_label;
};
#endif

#ifndef STRUCT_PROMPVARIABLES
#define STRUCT_PROMPVARIABLES
struct PromptVariables
{
  vector<vector<string> > Known_GG_string;
  vector<vector<vector<vector<RepVector> > > > Known_Zero_SortedCoupling;
  vector<vector<vector<vector<RepVector> > > > Known_NonZero_SortedCoupling;

  vector<string> AvailableLatticesFilenames;
  vector<string> AvailableLatticesLabels;
  vector<string> AvailableAdditionalLabels;
};
#endif


class CPrompt{
public:
// member functions
  CPrompt();
  CPrompt(const COrbifold &Orbifold);
  CPrompt(const string &Filename);
 ~CPrompt();

  bool              StartPrompt(string ifilename = "", bool stop_when_file_done = false, bool online_mode = false);
  bool              ExecuteCommand(string command);
  bool              ExecuteOrbifoldCommand(string command);
  bool              LoadProgram(const string &Filename, vector<string> &Commands);
  bool              LoadOrbifolds(const string &Filename, bool inequivalent = false, unsigned compare_couplings_up_to_order = 0);

  bool              FindCommandType0(const string &inputstring, const string &command) const;
  bool              FindCommandType1(const string &inputstring, const string &command, string &parameters) const;
  bool              FindCommandType2(const string &inputstring, const string &command, string &command_parameter, string &other_parameters) const;
  bool              FindParameterType1(string &inputstring, const string &parameter) const;
  bool              FindParameterType2(string &inputstring, const string &parameter, string &parameter_value) const;
  bool              FindParameterType3(string &inputstring, const string &parameter, unsigned &index) const;
  bool              FindParameterType4(string &inputstring, const string &parameter, unsigned &number) const;

  void              MessageHelpCreateNewOrbifold(unsigned StartWithLine = 0) const;
  bool              MessageLabelError(const string &Label) const;
  bool              MessageOrbifoldAlreadyExists(const string &OrbifoldLabel, bool PrintOutput = true) const;
  bool              MessageOrbifoldNotKnown(const string &OrbifoldLabel, unsigned &index) const;
  bool              MessageParameterNotKnown(const string &parameter_string) const;
  bool              MessageVEVConfigAlreadyExists(const vector<SConfig> &VEVConfigs, const string &VEVConfigLabel, const unsigned &VEVConfigNumber) const;
  bool              MessageVEVConfigNotKnown(const vector<SConfig> &VEVConfigs, const string &VEVConfigLabel, const unsigned &VEVConfigNumber, unsigned &index) const;
  bool              MessageXAlreadyExists(const vector<string> &NamesOfSetsOfX, const string &X, const string &X_Label) const;
  bool              MessageXNotKnown(const vector<string> &NamesOfSetsOfX, const string &X, const string &X_Label, unsigned &index) const;
    
  bool              PrintCurrentDirectory(string &output) const;
  void              PrintCommandsConditions() const;
  void              PrintCommandsMonomials() const;
  void              PrintCommandsProcesses() const;
  void              PrintCommandsSets() const;
  void              PrintFor(unsigned number_of_Type, const string &Type, const string &Var) const;
    
  bool              FindSpaceGroupsInDirectory(const unsigned &M, const unsigned &N, const string &directory);//orig N1
  //bool              FindSpaceGroupsInDirectory(const unsigned &N, const unsigned &K, const string &directory); //added J22

  bool              SplitVEVConfigLabel(string &VEVConfigLabel, unsigned &VEVConfigNumber) const;
//  void              ExtractLabels(const SUSYMultiplet &Multiplet, string input, vector<string> &FieldLabels);  //original
  void              ExtractLabels(const vector<SUSYMultiplet>  &Multiplet, string input, vector<string> &FieldLabels);  //added June20,2020

  

  vector<unsigned>  GetIndicesOnlyFieldWithNumber(const vector<string> &FieldLabels) const;
  vector<unsigned>  GetIndices(const vector<string> &FieldLabels) const;
  bool              GetLocalization(const string &Localization, CSpaceGroupElement &result) const;
  CFixedBrane      &AccessFixedBrane(const string &input, COrbifold &Orbifold, bool &FixedBraneFound) const;
  bool              FindConditions(string &input_string, vector<SCondition> &Conditions) const;
  bool              FindSUSYType(string &input_string, int NumberOfSupersymmetry, SUSYMultiplet &Multiplet) const;
  bool              ApplyConditions(const vector<SCondition> &Conditions, vector<unsigned> &FieldIndices) const;
  bool              ApplyConditions(vector<string> &Commands, unsigned &exec_command);
  bool              FindConditionsAndFilterFieldIndices(string &input_string, vector<unsigned> &FieldIndices) const;

// member variables
  CPrint                   Print;
  PromptVariables          PV;

  CGaugeIndices            GaugeIndices;
  vector<COrbifold>        Orbifolds;
  int                      OrbifoldIndex;
  vector<vector<SConfig> > AllVEVConfigs;
  vector<unsigned>         AllVEVConfigsIndex;

private:
  bool                     wait_for_processes;
  bool                     exit;           // Exit the prompt if true
  bool                     pre_exit;       // The question "Do you really want to quit?" was asked, ask again until the answer is "yes" or "no"

  bool                     print_output;   // can be set to true using the paramter "no output"
  bool                     online_mode;    // for the web interface

  bool                     keep_output_to_file;
  string                   output_filename;

  vector<int>              current_folder;  // current_folder[0] = -1 is the main directory, else the corresponding orbifold directory
                                            // current_folder[1] gives the subdirectory in an orbifold directory: 
                                            //   0 for no subdirectory
                                            //   1 for "model"
                                            //   2 for "gauge group"
                                            //   3 for "spectrum"
                                            //   4 for "couplings"
                                            //   5 for "vev-config"
                                            // current_folder[2] gives the further subdirectory of "vev-config" 
                                            //   0 for no subdirectory
                                            //   1 for "labels"
};

#endif
