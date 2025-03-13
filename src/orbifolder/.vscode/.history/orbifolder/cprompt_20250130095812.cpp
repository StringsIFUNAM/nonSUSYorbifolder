#include "cprompt.h"
#include "canalysemodel.h"
#include "globalfunctions.h"
#include "crandommodel.h"
#include "cinequivalentspectra.h"

#include <unistd.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <vector>

#include <cstdlib>
#include <iostream>
#include <cstring>
#include <readline/readline.h>
#include <readline/history.h>


extern unsigned SELFDUALLATTICE;
extern vector<int> current_folder_global;

using std::vector;
using std::cout;
using std::endl;

const char* command_names_main[] = {
    "create",
    "orbifold(",
    "with point group(",
    "cd",
    "random orbifold from(",
    "if(",
    "generations",
    "inequivalent",
    "save",
    "orbifolds",
    "use(",
    "#models(",
    "delete orbifold",
    "dir",
    "load",
    "orbifolds(",
    "rename orbifold(",
    "to(",
    "help",
    "print info",
    "when done",
    "program(",
    "from(",
    //"do not check anomalies",
    nullptr
};

const char* command_names_orbifold_model[] = {
    "dir",
    "help",
    "cd",
    "print",
    "orbifold label",
    "heterotic string type",
    "available space groups",
    "point group",
    "space group",
    "twist",
    "#SUSY",
    "Wilson lines",
    "use space group(",
    "set",
    "shift",
    "standard embedding",
    "WL W(",
    nullptr
};

const char* command_names_orbifold[] = {
    "m",
    "gg",
    "s",
    "v",
    "l",
    "cd",
    "model",
    "gauge group",
    "spectrum",
    "vev-config",
    "vev-config/labels",
    "dir",
    "help",
    "print",
    nullptr
};

const char* command_names_orbifold_gauge[] = {
    "dir",
    "help",
    "cd",
    "print",
    "set",
    "U1(",
    //"B-L",
    "print",
    "gauge group",
    "beta coefficients",
    "simple root(",
    "simple roots",
    //"FI term",
    "anomaly info",
    //"B-L generator",
    "U1",
    "generator(",
    "generators",
    nullptr
};


const char* command_names_orbifold_spectrum[] = {
    "dir",
    "help",
    "cd",
    "print",
    "all states",
    "summary",
    "print(",
    "list of charges(",
    "tex table(",
    "of",
    "sectors",
    "fixed"
    ,"points",
    "sector T(",
    "point(",
    "no U1s",
    "with labels",
    nullptr
};


const char* command_names_orbifold_couplings[] = {
    "dir",
    nullptr
};


const char* command_names_orbifold_vev_main[] = {
    "dir",
    "help",
    "cd",
    "use config(",
    "create config(",
    "rename config(",
    "to(",
    "delete config(",
    "print",
    "configs",
    "gauge group",
    "analyze config",
    "labels",
    "select observable sector:",
    "gauge group(",
    "full gauge group",
    "no gauge groups",
    "U1s(",
    "all U1s",
    "no U1s",
    nullptr
};


const char* command_names_orbifold_vev_labels[] = {
    "dir",
    "help",
    "cd",
    "change label(",
    "to(",
    "create labels",
    "assign label(",
    "to fixed point(",
    "print labels",
    "use label(",
    "load labels(",
    "save labels(",
    nullptr
};


const char** getCommandNames(vector<int>&current_folder_global) {
    if (current_folder_global[0]==-1)
    {
      return command_names_main;
    }
    else
    {
    switch (current_folder_global[1]) {
        case 0:
        {
            return command_names_orbifold;
        }
        case 1:
        {
            return command_names_orbifold_model;
        }
        case 2:
        {
            return command_names_orbifold_gauge;
        }
        case 3:
        {
            return command_names_orbifold_spectrum;
        }
        case 4:
        {
            return command_names_orbifold_couplings;
        }
        case 5:
        {
            switch (current_folder_global[2])
            {
              case 0:
              {
                return command_names_orbifold_vev_main;
              }
              case 1:
              {
                return command_names_orbifold_vev_labels;
              }
            }


        }
        default:
            return nullptr;
    }
    }
}






char** command_name_completion(const char*, int, int);
char* command_name_generator(const char*, int);

char** command_name_completion(const char* text, int start, int end) {
    rl_attempted_completion_over = 1;
    return rl_completion_matches(text, command_name_generator);
}

char* command_name_generator(const char* text, int state) {
    static int list_index, len;
    const char* name;

    if (!state) {
        list_index = 0;
        len = strlen(text);
    }

    const char** command_names = getCommandNames(current_folder_global);

    while ((name = command_names[list_index++])) {
        rl_completion_append_character = '\0';
        if (strncmp(name, text, len) == 0) {
            return strdup(name);
        }
    }

    return nullptr;
}





string CPrompt::Trim(const string &command) const {
    size_t start = command.find_first_not_of(" \t\n\r\f\v");
    size_t end = command.find_last_not_of(" \t\n\r\f\v");
    if (start == string::npos) {
        return "";
    }
    return command.substr(start, end - start + 1);
}




bool containsErrorKeywords(const string &error_message, const vector<string> &keywords) {
    for (const auto &keyword : keywords) {
        if (error_message.find(keyword) != string::npos) {
            return true;
        }
    }
    return false;
}




vector<string> error_keywords = {
    "No existe el archivo o el directorio",  // Español
    "Ninguna entrada del manual para",       // Español
    "No manual entry for",                   // Inglés
    "No such file or directory",             // Inglés
    "Aucune entrée de manuel pour",          // Francés
    "Aucun fichier ou dossier de ce type",   // Francés
    "Nessuna voce di manuale per",           // Italiano
    "Nessun file o directory"                // Italiano
    // Puedes agregar más palabras clave según sea necesario
};




/* ########################################################################################
######   CPrompt(COrbifold &Orbifold)                                                ######
######                                                                               ######
######   Version: 06.01.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard constructor of a CPrompt object. No content is specified.          ######
######################################################################################## */
CPrompt::CPrompt()
  : Print(Tstandard, &cout)
{
  // start the prompt in its main directory
  this->current_folder.assign(10,0);
  this->current_folder[0]  = -1;
  this->OrbifoldIndex = -1;

  this->output_filename     = "";
  this->keep_output_to_file = false;
  this->pre_exit            = false;
  this->exit                = false;
  this->online_mode         = false;
  this->print_output        = true;
  this->wait_for_processes  = false;


}



/* ########################################################################################
######   CPrompt(const COrbifold &Orbifold)                                          ######
######                                                                               ######
######   Version: 06.01.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Orbifold : COrbifold object to be loaded into the prompt                 ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Constructor of a CPrompt object. Loads "Orbifold" as a directory into the   ######
######   prompt.                                                                     ######
######################################################################################## */
CPrompt::CPrompt(const COrbifold &Orbifold)
  : Print(Tstandard, &cout)
{
  // begin: add orbifold
  this->Orbifolds.push_back(Orbifold);
  this->OrbifoldIndex = 0;

  SConfig TestConfig = Orbifold.StandardConfig;
  TestConfig.ConfigLabel = "TestConfig";
  TestConfig.ConfigNumber = 1;

  vector<SConfig> Configs;
  Configs.push_back(Orbifold.StandardConfig);
  Configs.push_back(TestConfig);
  this->AllVEVConfigs.push_back(Configs);
  this->AllVEVConfigsIndex.push_back(1);
  // end: add orbifold

  // start the prompt in the directory of "Orbifold"
  this->current_folder.assign(10,0);
  this->current_folder[0]   = 0;

  this->output_filename     = "";
  this->keep_output_to_file = false;
  this->pre_exit            = false;
  this->exit                = false;
  this->online_mode         = false;
  this->print_output        = true;
  this->wait_for_processes  = false;

}



/* ########################################################################################
######   CPrompt(const string &Filename                                              ######
######                                                                               ######
######   Version: 06.01.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Filename : file name of a model file                                     ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Constructor of a CPrompt object. Load orbifolds from model file "Filename". ######
######################################################################################## */
CPrompt::CPrompt(const string &Filename)
  : Print(Tstandard, &cout)
{
  // start the prompt in its main directory
  this->current_folder.assign(10,0);
  this->current_folder[0]  = -1;
  this->OrbifoldIndex = -1;

  this->output_filename     = "";
  this->keep_output_to_file = false;
  this->pre_exit            = false;
  this->exit                = false;
  this->online_mode         = false;
  this->print_output        = true;
  this->wait_for_processes  = false;




  (*this->Print.out) << "\n";
  this->LoadOrbifolds(Filename, false, 0);
}



/* ########################################################################################
######   ~CPrompt()                                                                  ######
######                                                                               ######
######   Version: 06.01.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Standard destructor of a CPrompt object.                                    ######
######################################################################################## */
CPrompt::~CPrompt()
{
}


/* ########################################################################################
######   StartPrompt(string ifilename, bool stop_when_file_done, bool online_mode)   ######
######                                                                               ######
######   Version: 26.09.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) ifilename           : load commands from this file if specified          ######
######   2) stop_when_file_done : only execute commands from file then stop prompt   ######
######   3) online_mode         : for the web interface                              ######
######   output:                                                                     ######
######   return value           : prompt stoppt without problems?                    ######
###########################################################################################
######   description:                                                                ######
######   Starts the linux-style command line called the prompt.                      ######
######################################################################################## */
bool CPrompt::StartPrompt(string ifilename, bool stop_when_file_done, bool online_mode)
{
  this->online_mode        = online_mode;
  this->wait_for_processes = false;
  this->exit               = false;
  this->pre_exit           = false;

  this->Print.out  = &cout;
  this->Print.SetOutputType(Tstandard);

  // if "online_mode" than write all output to file "result.txt"
  // otherwise "keep_output_to_file" might be turned on using "@begin print to file"
  this->keep_output_to_file = this->online_mode;

  unsigned exec_command = 0;
  vector<string> Commands;

  string tmp_string1 = "";
  string command     = "";





  if (this->online_mode)
  {
    //ifilename = "program.txt";
    size_t pos_point = ifilename.find('.');
    string name_without_extension = "";
    if (pos_point != std::string::npos)
    {
        name_without_extension = ifilename.substr(0, pos_point);
    }


    this->output_filename  = "result_"+name_without_extension+".txt";
  }
  else
  {
    if (ifilename != "")
    {
      this->LoadProgram(ifilename, Commands);
      if ((Commands.size() == 0) && stop_when_file_done)
        return true;

      command   = "";
      ifilename = "";
    }
  }

  string parameter_string1 = "";
  string parameter_string2 = "";

  std::ofstream output_file;
  std::ofstream output_file_aux;

  size_t old_number_of_commands = 0;
  bool output_file_open = false;
  bool use_cin          = false;









  while (true)
  {
    // STEP 0 //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // begin: online mode
    if (this->online_mode)
    {
      struct stat stFileInfo;

      old_number_of_commands = Commands.size();
      while ((exec_command == Commands.size()) )//|| (stat(ifilename.c_str(),&stFileInfo) == 0))
      {
        // begin: wait for the file "ifilename" to exist
        if (stat(ifilename.c_str(), &stFileInfo) != 0)
        {
          cout << "  Waiting for file \"" << ifilename << "\"." << endl;

          while (stat(ifilename.c_str(), &stFileInfo) != 0)
            usleep(1000);
        }
        // end: wait for the file "ifilename" to exist

        // begin: load the commands
        usleep(5000);
        std::ifstream input(ifilename.data());
        while (GetSaveLine(input, command))
          Commands.push_back(command);
        input.close();
        // end: load the commands

        // begin: delete the file "ifilename"
        //usleep(500);
        //tmp_string1 = "rm ";
        //tmp_string1 += ifilename;
        //if (system(tmp_string1.data()) != 0)
        //  cout << "  System(" << tmp_string1 << ") failed." << endl;

        // end: delete the file "ifilename"
      }
      if (old_number_of_commands != Commands.size())
      {
        if (output_file_open)
        {
          (*this->Print.out) << endl << "eof" << endl;
          usleep(1000);
          output_file.close();
          usleep(7500);
        }
        output_file.open(this->output_filename.data());
        output_file_open = true;

        this->Print.out = &output_file;
      }
    }
    // end: online mode
    // STEP 0 //////////////////////////////////////////////////////////////////////////////////////////////////////////

    /*cout << endl << endl << "current commmand: " << Commands[exec_command] << endl;
    cout << "begin: all commands before ApplyConditions" << endl;
    for (unsigned q= 0; q < Commands.size(); ++q)
    cout << Commands[q] << endl;
    cout << "end: all commands before ApplyConditions" << endl;
    cout << "exec_command = " << exec_command << endl;*/

    this->ApplyConditions(Commands, exec_command);

    /*cout << "begin: all commands after ApplyConditions" << endl;
    for (unsigned q= 0; q < Commands.size(); ++q)
    cout << Commands[q] << endl;
    cout << "end: all commands after ApplyConditions" << endl;
    cout << "exec_command = " << exec_command << endl;*/

    // STEP 1 //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // begin: read commands
    use_cin = false;
    if (!this->online_mode && (exec_command == Commands.size()))
    {
      if (stop_when_file_done)
      {
        if (output_file_open)
        {
          (*this->Print.out) << flush;

          usleep(1000);
          output_file.close();
          output_file_open = false;
        }
      }
      // It's comment, for elimine double prompt
      /*if (!this->keep_output_to_file)
      {
        this->PrintCurrentDirectory(tmp_string1);
        (*this->Print.out) << tmp_string1 << flush;
      }*/

      use_cin = true;

      // just replace this

      char* user_input;

      string my_folder;

      string my_prompt = ""; //"/> "

      current_folder_global = current_folder;

      rl_attempted_completion_function = command_name_completion;
      const char* rl_completer_word_break_characters = " ";

      if (PrintCurrentDirectory(my_folder)){
        //CPrompt::command_name_completion

        const char* rl_completer_word_break_characters = " ";
        user_input = readline((my_folder+my_prompt).c_str());
        command = user_input;
        add_history(user_input);
      }



      free(user_input);
      /*
      while (() != nullptr) {

      }
      */
      //GetSaveLine(cin, command);


      // just replace this

      cin.clear();
      Commands.push_back(Trim(command));
    }

    command = Commands[exec_command];
    ++exec_command;
    // end: read commands
    // STEP 1 //////////////////////////////////////////////////////////////////////////////////////////////////////////

    // STEP 2 //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // begin: replace variables
    if (command.find("$OrbifoldLabel$", 0) != string::npos)
    {
      if ((this->Orbifolds.size() == 0) || (this->OrbifoldIndex == -1))
        (*this->Print.out) << "\n  " << this->Print.cbegin << "Cannot replace variable \"$OrbifoldLabel$\"." << this->Print.cend << "\n" << endl;
      else
        global_ReplaceString(command, "$OrbifoldLabel$", this->Orbifolds[this->OrbifoldIndex].OrbifoldGroup.Label);
    }
    if (command.find("$VEVConfigLabel$", 0) != string::npos)
    {
      if ((this->Orbifolds.size() == 0) || (this->OrbifoldIndex == -1))
        (*this->Print.out) << "\n  " << this->Print.cbegin << "Cannot replace variable \"$VEVConfigLabel$\"." << this->Print.cend << "\n" << endl;
      else
      {
        const SConfig &VEVConfig = this->AllVEVConfigs[this->OrbifoldIndex][this->AllVEVConfigsIndex[this->OrbifoldIndex]];
        global_ReplaceString(command, "$VEVConfigLabel$", VEVConfig.ConfigLabel);
      }
    }
    if (command.find("$Directory$", 0) != string::npos)
    {
      string tmp = "";
      this->PrintCurrentDirectory(tmp);
      global_ReplaceString(command, "$Directory$", tmp);
    }
    // end: replace variables
    // STEP 2 //////////////////////////////////////////////////////////////////////////////////////////////////////////

    this->print_output = !this->FindParameterType1(command, "no output");

    // STEP 3 //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // begin: write output to file?
    //if (this->online_mode)
    //{
    //  if (this->FindParameterType1(command, "@begin print"))
    //  {
    //    if (this->print_output)
    //      (*this->Print.out) << "\n  " << this->Print.cbegin << "Command disabled in the web interface." << this->Print.cend << "\n" << endl;
    //    command = "";
    //  }
    //}
    //else
    {
      if (this->FindParameterType2(command, "@begin print to file(", this->output_filename))
      {
        this->keep_output_to_file = true;

        output_file.open(this->output_filename.data(), ofstream::app | ios::ate);
        output_file_open = true;

        this->Print.out = &output_file;
      }
      else
        if (this->FindParameterType2(command, " to file(", tmp_string1))
      {
        if (this->keep_output_to_file)
        {
          if (this->print_output)
          {
            if (this->online_mode)
            {
                this->output_filename_aux = tmp_string1;
                (*this->Print.out) << "\n  " << this->Print.cbegin << "Result written to file \"" << this->output_filename_aux << "\"." << this->Print.cend << "\n" << endl;

                output_file_aux.open(this->output_filename_aux.data(), ofstream::out | ios::ate);


                if (output_file_aux.is_open()&& output_file_aux.good())
                {
                    this->Print.out = &output_file_aux;
                    (*this->Print.out) << flush;
                }

            }
            else
                (*this->Print.out) << "\n  " << this->Print.cbegin << "Output permanently written to file \"" << this->output_filename << "\". Ignoring parameter \"to file(" << tmp_string1 << ")\"." << this->Print.cend << "\n" << endl;
          }
        }
        else
        {
          this->output_filename = tmp_string1;
          if (this->print_output)
            (*this->Print.out) << "\n  " << this->Print.cbegin << "Result written to file \"" << this->output_filename << "\"." << this->Print.cend << "\n" << endl;

          output_file.open(this->output_filename.data(), ofstream::app | ios::ate);
          output_file_open = true;

          this->Print.out = &output_file;
        }
      }
    }


    if (output_file_open && (!output_file.is_open() || !output_file.good()))
    {
      if (this->print_output)
        (*this->Print.out) << "\n  " << this->Print.cbegin << "Cannot write to file \"" << this->output_filename << "\"." << this->Print.cend << "\n" << endl;
      this->Print.out = &cout;

      output_file.close();
      output_file_open = false;
    }
    // end: write output to file?
    // STEP 3 //////////////////////////////////////////////////////////////////////////////////////////////////////////


    // STEP 4 //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // begin: execute commands
    if (!this->online_mode && this->FindParameterType1(command, "@end print to file"))
    {
      (*this->Print.out) << flush;
      this->keep_output_to_file = false;

      usleep(1000);
      output_file.close();
      output_file_open = false;

      this->Print.out = &cout;
      usleep(500);
    }

    if (this->online_mode || (!use_cin && !output_file_open))
    {
      this->PrintCurrentDirectory(tmp_string1);
      (*this->Print.out) << tmp_string1 << command << endl;
    }

    if (this->FindParameterType2(command, "load program(", ifilename))
      this->LoadProgram(ifilename, Commands);

    if (this->pre_exit)
    {
      if (FindCommandType0(command, "yes") || FindCommandType0(command, "Yes"))
        this->exit = true;
      else
      {
        if (FindCommandType0(command, "no") || FindCommandType0(command, "No"))
        {
          if (this->print_output)
            (*this->Print.out) << "\n  " << this->Print.cbegin << "Program not stopped." << this->Print.cend << "\n" << endl;

          this->pre_exit = false;
        }
        else
        {
          if (this->print_output)
            (*this->Print.out) << "\n  " << this->Print.cbegin << "Do you really want to quit? Type \"yes\" to quit or \"no\" to continue." << this->Print.cend << "\n" << endl;
        }
      }
    }
    else
      this->ExecuteCommand(command);

    // end: execute commands
    // STEP 4 //////////////////////////////////////////////////////////////////////////////////////////////////////////

    usleep(100);

    if (this->online_mode && (exec_command == Commands.size()))
    {
      (*this->Print.out) << endl << "eof" << endl;

      usleep(1000);
      output_file.close();
      output_file_open = false;
      usleep(7500);
      cout << "Script executed successfully! The resulting output was saved in the file: " <<this->output_filename<< "." << endl;
      return true;
    }

    if (output_file_aux.is_open() && output_file_aux.good())
    {

      output_file_aux.close();

      this->Print.out = &output_file;
    }

    if (output_file_open && !this->keep_output_to_file)
    {
      (*this->Print.out) << flush;

      usleep(1000);
      output_file.close();
      output_file_open = false;

      this->Print.out = &cout;
      usleep(500);
    }

    if (this->exit)
    {
      if (output_file_open)
      {
        (*this->Print.out) << flush;

        usleep(1000);
        output_file.close();
        output_file_open = false;
      }

      cout << "\n  " << this->Print.cbegin << "End." << this->Print.cend << "\n" << endl;

      this->pre_exit = false;
      this->exit     = false;
      return true;
    }
  }
  return false;
}



/* ########################################################################################
######   ExecuteCommand(string command)                                              ######
######                                                                               ######
######   Version: 10.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) command   : a string containing a command to be executed                 ######
######   output:                                                                     ######
######   return value : command executed?                                            ######
###########################################################################################
######   description:                                                                ######
######   Interprets "command" and executes it. If "command" is not a system command  ######
######   call "ExecuteOrbifoldCommand(command)".                                     ######
######################################################################################## */
bool CPrompt::ExecuteCommand(string command)
{
  string parameter_string1 = "";
  string parameter_string2 = "";

  // ignore comments
  // updated on 17.09.2010
  if (this->FindCommandType1(command, "//", parameter_string1))
    return true;

  if (this->FindCommandType1(command, "@endif", parameter_string1))
    return true;

  if (this->FindCommandType1(command, "@else", parameter_string1))
    return true;

  // print enter
  // updated on 20.04.2011
  if (this->FindCommandType1(command, "@print enter", parameter_string1))
  {
    (*this->Print.out) << endl;

    if (this->print_output)
      this->MessageParameterNotKnown(parameter_string1);
    return true;
  }

  // print line
  // updated on 20.04.2011
  if (this->FindCommandType2(command, "@print(", parameter_string1, parameter_string2))
  {
    if (this->FindParameterType1(parameter_string2, "unformatted"))
      (*this->Print.out) << parameter_string1;
    else
      (*this->Print.out) << this->Print.cbegin << parameter_string1 << this->Print.cend << endl;

    if (this->print_output)
      this->MessageParameterNotKnown(parameter_string2);
    return true;
  }

  bool restore_old_output_type = false;
  OutputType old_output_type = Tstandard;

  if (this->FindParameterType1(command, "@mathematica"))
  {
    restore_old_output_type = true;
    old_output_type = this->Print.GetOutputType();
    this->Print.SetOutputType(Tmathematica);
  }
  else
  if (this->FindParameterType1(command, "@latex"))
  {
    restore_old_output_type = true;
    old_output_type = this->Print.GetOutputType();
    this->Print.SetOutputType(Tlatex);
  }
  else
  if (this->FindParameterType1(command, "@standard"))
  {
    restore_old_output_type = true;
    old_output_type = this->Print.GetOutputType();
    this->Print.SetOutputType(Tstandard);
  }

  bool     cend_printed = false;
  unsigned DurationWait = 5;

  size_t t1 = 0;
  size_t t2 = 0;
  size_t t3 = 0;
  unsigned i = 0;
  unsigned j = 0;
  unsigned k = 0;

  string tmp_string1 = "";
  string tmp_string2 = "";

  // step 1: child processes
  // step 2: prompt commands

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // wait until all child processes are finished
  // updated on 20.04.2011
  if (this->FindCommandType2(command, "wait(", parameter_string1, parameter_string2) || this->FindCommandType1(command, "wait", parameter_string2))
  {
    const unsigned waiting_intervall = 5;
    const unsigned MAX_waiting_intervall = 600;

    DurationWait = waiting_intervall;

    if (parameter_string1.size() != 0)
    {
      if (parameter_string1.find_first_not_of("0123456789") != string::npos)
      {
        if (this->print_output)
          (*this->Print.out) << "\n  " << this->Print.cbegin << "Duration of waiting interval ill defined: \"" << parameter_string1 << "\". Set to " << waiting_intervall << " sec." << this->Print.cend << "\n" << endl;

        DurationWait = waiting_intervall;
      }
      else
      {
        DurationWait = (unsigned)atoi(parameter_string1.c_str());
        if (DurationWait > MAX_waiting_intervall)
        {
          if (this->print_output)
            (*this->Print.out) << "\n  " << this->Print.cbegin << "Duration of waiting interval too long. Set to " << MAX_waiting_intervall << " sec." << this->Print.cend << "\n" << endl;

          DurationWait = MAX_waiting_intervall;
        }
      }
    }

    if (this->print_output)
      (*this->Print.out) << "\n  " << this->Print.cbegin << "Waiting for child processes to finish..." << flush;

    sleep(1);
    this->wait_for_processes = true;
  }

  CAnalyseModel Analyse;

  // STEP 1 //////////////////////////////////////////////////////////////////////////////////////////////////////////
  // begin: child processes
  t1 = this->Orbifolds.size();

  bool NoProcessRunning = true;
  do {
    NoProcessRunning = true;

    // begin: run through all orbifold directories of the prompt
    for (i = 0; i < t1; ++i)
    {
      COrbifold       &Orbifold   = this->Orbifolds[i];
      vector<SConfig> &VEVConfigs = this->AllVEVConfigs[i];

      // begin: run through all vev-configurations of the current orbifold directory
      t2 = VEVConfigs.size();
      for (j = 0; j < t2; ++j)
      {
        SConfig &VEVConfig = VEVConfigs[j];

        PID &PID_Data = VEVConfig.pid;
        const size_t p1 = PID_Data.PIDs.size();

        // begin: check status of child processes
        for (k = 0; k < p1; ++k)
        {
          if (!PID_Data.PID_Done[k])
          {
            pid_t currentPID = PID_Data.PIDs[k];

            int status = 0;
            if (waitpid(currentPID, &status, WNOHANG) > 0)
            {
              PID_Data.PID_Done[k] = true;

              std::ostringstream osID;
              osID << currentPID;
              const  string   Filename = PID_Data.PID_Filenames[k];
              const  unsigned JobIndex = PID_Data.PID_JobIndices[k];
              string          Command  = PID_Data.PID_Commands[k];

              const double diff = difftime (time(NULL), PID_Data.PID_StartingTimes[k]);
              // begin: process finished without problems
              if (WIFEXITED(status))
              {
                if (this->wait_for_processes)
                {
                  if (!cend_printed)
                    (*this->Print.out) << this->Print.cend << endl;
                  cend_printed = true;
                }

                if (i == this->OrbifoldIndex)
                  (*this->Print.out) << "\n  " << this->Print.cbegin << "PID " << currentPID << " done: \"" << Command << "\" TIME: " << setfill ('0') << setw (2) << (unsigned)(diff/3600.0) << ":" << setw (2) << (unsigned)(fmod((diff/60.0),60.0)) << ":" << setw (2) << (unsigned)fmod(diff, 60.0) << setfill (' ') << this->Print.cend << "\n" << endl;
                else
                  (*this->Print.out) << "\n  " << this->Print.cbegin << "For orbifold \"" << Orbifold.OrbifoldGroup.Label << "\" PID " << currentPID << " done: \"" << Command << "\" TIME: " << setfill ('0') << setw (2) << (unsigned)(diff/3600.0) << ":" << setw (2) << (unsigned)(fmod((diff/60.0),60.0)) << ":" << setw (2) << (unsigned)fmod(diff, 60.0) << setfill (' ') << this->Print.cend << "\n" << endl;

                vector<CField> &Fields = VEVConfig.Fields;

                std::ifstream in;
                if (JobIndex != 3)
                {
                  in.open(Filename.data());
                  if((!in.is_open()) || (!in.good()))
                  {
                    (*this->Print.out) << "\n  " << this->Print.cbegin << "File \"" << Filename << "\" not found." << this->Print.cend << "\n" << endl;
                    return true;
                  }
                }

                if (JobIndex != 3)
                  in.close();
              }
              // end: process finished without problems
              else
                (*this->Print.out) << "\n  " << this->Print.cbegin << "PID " << currentPID << " stopped: \"" << PID_Data.PID_Commands[k] << "\" TIME: " << setfill ('0') << setw (2) << (unsigned)(diff/3600.0) << ":" << setw (2) << (unsigned)(fmod((diff/60.0),60.0)) << ":" << setw (2) << (unsigned)fmod(diff, 60.0) << setfill (' ') << this->Print.cend << "\n" << endl;

              // updated on 24.08.2011
              if (JobIndex == 3) // create random orbifold
              {
                if (this->FindParameterType1(Command, "load when done"))
                  this->LoadOrbifolds(Filename, false, 0);

                if (!this->FindParameterType2(Command, "to(", tmp_string1))
                {
                  string rm_command = "rm " + Filename;
                  if (system(rm_command.data()) != 0)
                    cout << "  System(" << rm_command << ") failed." << endl;
                }
                else
                  (*this->Print.out) << "  " << this->Print.cbegin << "New orbifolds saved to file \"" << Filename << "\"." << this->Print.cend << "\n" << endl;
              }
            }
          }
        }
        // end: check status of child processes

        if (this->wait_for_processes && (find(PID_Data.PID_Done.begin(), PID_Data.PID_Done.end(), false) != PID_Data.PID_Done.end()))
          NoProcessRunning = false;
      }
      // end: run through all vev-configurations of the current orbifold directory
    }
    // end: run through all orbifold directories of the prompt

    if (this->wait_for_processes)
    {
      if (NoProcessRunning)
      {
        if (this->print_output)
        {
          if (cend_printed)
            (*this->Print.out) << "  " << this->Print.cbegin;
          (*this->Print.out) << "waiting done." << this->Print.cend << endl;
        }
        this->wait_for_processes = false;
        sleep(1);

        this->MessageParameterNotKnown(parameter_string2);
        return true;
      }
      else
      {
        if (this->print_output)
          (*this->Print.out) << "." << flush;

        sleep(DurationWait);
      }
    }
  } while(this->wait_for_processes);
  // end: child processes
  // STEP 1 //////////////////////////////////////////////////////////////////////////////////////////////////////////


  // STEP 2 //////////////////////////////////////////////////////////////////////////////////////////////////////////
  // begin: prompt commands

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // kill a processes
  // updated on 14.07.2011
  if (this->FindCommandType2(command, "kill(", parameter_string1, parameter_string2))
  {
    const bool KillAll = (parameter_string1 == "all");

    if (!KillAll && (parameter_string1.find_first_not_of("0123456789") != string::npos))
    {
      if (this->print_output)
      {
        (*this->Print.out) << "\n  " << this->Print.cbegin << "PID \"" << parameter_string1 << "\" ill defined." << this->Print.cend << "\n" << endl;
        this->MessageParameterNotKnown(parameter_string2);
      }
      return true;
    }


    pid_t killPID = 0;
    if (!KillAll)
      killPID = (unsigned)atoi(parameter_string1.c_str());

    bool PID_found = false;
    if (KillAll)
      PID_found = true;

    vector<unsigned> IndicesOfProcessesToKill;

    t1 = this->AllVEVConfigs.size();
    for (i = 0; i < t1; ++i)
    {
      vector<SConfig> &VEVConfigs = this->AllVEVConfigs[i];

      t2 = VEVConfigs.size();
      for (j = 0; j < t2; ++j)
      {
        PID &PID_Data = VEVConfigs[j].pid;

        IndicesOfProcessesToKill.clear();
        if (KillAll)
        {
          t3 = PID_Data.PIDs.size();
          for (k = 0; k < t3; ++k)
            IndicesOfProcessesToKill.push_back(k);
        }
        else
        {
          vector<pid_t>::iterator pos = find(PID_Data.PIDs.begin(), PID_Data.PIDs.end(), killPID);
          if (pos != PID_Data.PIDs.end())
          {
            PID_found = true;
            IndicesOfProcessesToKill.push_back(distance(PID_Data.PIDs.begin(), pos));
          }
        }

        t3 = IndicesOfProcessesToKill.size();
        if ((t3 != 0) && this->print_output)
          (*this->Print.out) << "\n";

        for (k = 0; k < t3; ++k)
        {
          if (!PID_Data.PID_Done[IndicesOfProcessesToKill[k]])
          {
            string command = "kill ";
            std::ostringstream os;
            os << PID_Data.PIDs[IndicesOfProcessesToKill[k]];
            command += os.str();

            if (system(command.data()) != 0)
              cout << "  System(" << command << ") failed." << endl;

            if (this->print_output)
              (*this->Print.out) << "  " << this->Print.cbegin << "Executing command \"" << command << "\"." << this->Print.cend << "\n";
          }
        }
      }
    }

    if (!PID_found && this->print_output)
    {
      (*this->Print.out) << "\n  " << this->Print.cbegin << "PID " << killPID << " unknown." << this->Print.cend << "\n" << endl;
      this->MessageParameterNotKnown(parameter_string2);
    }

    if (this->print_output)
    {
      (*this->Print.out) << endl;
      this->MessageParameterNotKnown(parameter_string2);
    }
    return true;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // show active processes
  // updated on 11.05.2011
  if (this->FindCommandType1(command, "ps", parameter_string1))
  {
    double diff = 0.0;
    bool printed = false;

    (*this->Print.out) << "\n";
    t1 = this->Orbifolds.size();
    for (i = 0; i < t1; ++i)
    {
      const COrbifold       &Orbifold   = this->Orbifolds[i];
      const vector<SConfig> &VEVConfigs = this->AllVEVConfigs[i];

      t2 = VEVConfigs.size();
      for (j = 0; j < t2; ++j)
      {
        const SConfig &VEVConfig = VEVConfigs[j];
        const PID     &PID_Data  = VEVConfig.pid;

        if (find(PID_Data.PID_Done.begin(), PID_Data.PID_Done.end(), false) != PID_Data.PID_Done.end())
        {
          printed = true;
          (*this->Print.out) << "  " << this->Print.cbegin << "Processes of model \"" << Orbifold.OrbifoldGroup.Label << "\" in vev-configuration \"" << VEVConfig.ConfigLabel << VEVConfig.ConfigNumber <<  "\":" << this->Print.cend << "\n";
          (*this->Print.out) << "  " << this->Print.cbegin << "  PID       TIME CMD" << this->Print.cend << endl;
          t3 = PID_Data.PIDs.size();
          for (k = 0; k < t3; ++k)
          {
            if (!PID_Data.PID_Done[k])
            {
              diff = difftime (time(NULL), PID_Data.PID_StartingTimes[k]);
              (*this->Print.out) << "  " << this->Print.cbegin << setw(5) << PID_Data.PIDs[k] << "   ";
              (*this->Print.out) << setfill ('0') << setw(2) << (unsigned)(diff/3600.0) << ":" << setw(2) << (unsigned)(fmod((diff/60.0),60.0)) << ":" << setw(2) << (unsigned)fmod(diff, 60.0);
              (*this->Print.out) << setfill (' ') << " " << PID_Data.PID_Commands[k] << this->Print.cend << "\n";
            }
          }
          (*this->Print.out) << endl;
        }
      }
    }
    if (!printed)
      (*this->Print.out) << "  " << this->Print.cbegin << "No processes are running." << this->Print.cend << "\n" << endl;

    this->MessageParameterNotKnown(parameter_string1);
    return true;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // exit orbifolder
  // updated on 12.10.2011
  if (this->FindCommandType1(command, "exit orbifolder", parameter_string1))
  {
    this->exit = true;
    this->MessageParameterNotKnown(parameter_string1);
    return true;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // exit
  // updated on 12.10.2011
  // !this->online_mode &&
  if ( this->FindCommandType1(command, "exit", parameter_string1))
  {
    this->pre_exit = true;

    t1 = this->AllVEVConfigs.size();
    for (i = 0; this->pre_exit && (i < t1); ++i)
    {
      const vector<SConfig> &VEVConfigs = this->AllVEVConfigs[i];

      t2 = VEVConfigs.size();
      for (j = 0; this->pre_exit && (j < t2); ++j)
      {
        const PID &PID_Data = VEVConfigs[j].pid;

        if (find(PID_Data.PID_Done.begin(), PID_Data.PID_Done.end(), false) != PID_Data.PID_Done.end())
        {
          if (this->print_output)
            (*this->Print.out) << "\n  " << this->Print.cbegin << "Cannot exit - some processes are stil running." << this->Print.cend << "\n" << endl;

          this->pre_exit = false;
        }
      }
    }
    if (this->pre_exit)
    {
      if (this->print_output)
        (*this->Print.out) << "\n  " << this->Print.cbegin << "Do you really want to quit? Type \"yes\" to quit or \"no\" to continue." << this->Print.cend << "\n" << endl;
    }
    this->MessageParameterNotKnown(parameter_string1);
    return true;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // system commands
  // updated on 20.04.2011
  if (this->FindCommandType1(command, "@", parameter_string1))
  {
    if (!this->online_mode && this->FindParameterType1(parameter_string1, "end print to file"))
    {
      this->output_filename     = "";
      this->keep_output_to_file = false;
      this->Print.out = &cout;
    }
    else
    if (this->FindParameterType1(parameter_string1, "status"))
    {
      (*this->Print.out) << "\n  status:\n";
      //if (this->online_mode)
      //  (*this->Print.out) << "    output is written to the web interface.\n";
      //else
      //{
      //  (*this->Print.out) << "    output is written to: ";
      //  if (this->keep_output_to_file)
      //    (*this->Print.out) << "file \"" << this->output_filename << "\".\n";
      //  else
      //    (*this->Print.out) << "screen.\n";
      //}

      (*this->Print.out) << "    typesetting: ";
      if (this->Print.GetOutputType() == Tstandard)
        (*this->Print.out) << "standard\n";
      if (this->Print.GetOutputType() == Tmathematica)
        (*this->Print.out) << "mathematica\n";
      if (this->Print.GetOutputType() == Tlatex)
        (*this->Print.out) << "latex\n";
      (*this->Print.out) << endl;
    }
    else
    if (this->FindParameterType2(parameter_string1, "typesetting(", tmp_string1))
    {
      if (tmp_string1 == "mathematica")
      {
        this->Print.SetOutputType(Tmathematica);

        if (this->print_output)
          (*this->Print.out) << "\n  " << this->Print.cbegin << "Typesetting set to mathematica style." << this->Print.cend << "\n" << endl;
      }
      else
      if (tmp_string1 == "latex")
      {
        this->Print.SetOutputType(Tlatex);

        if (this->print_output)
          (*this->Print.out) << "\n  " << this->Print.cbegin << "Typesetting set to latex style." << this->Print.cend << "\n" << endl;
      }
      else
      if (tmp_string1 == "standard")
      {
        this->Print.SetOutputType(Tstandard);

        if (this->print_output)
          (*this->Print.out) << "\n  " << this->Print.cbegin << "Typesetting set to standard style." << this->Print.cend << "\n" << endl;
      }
    }
    this->MessageParameterNotKnown(parameter_string1);
    return true;
  }
  // end: prompt commands
  // STEP 2 //////////////////////////////////////////////////////////////////////////////////////////////////////////

  this->ExecuteOrbifoldCommand(command);

  if (restore_old_output_type)
  {
    restore_old_output_type = false;
    this->Print.SetOutputType(old_output_type);
  }
  return true;
}



/* ########################################################################################
######   ExecuteOrbifoldCommand(string command)                                      ######
######                                                                               ######
######   Version: 12.01.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) command   : a string containing a command to be executed                 ######
######   output:                                                                     ######
######   return value : command executed?                                            ######
###########################################################################################
######   description:                                                                ######
######   Interprets "command" and executes it.                                       ######
######################################################################################## */
bool CPrompt::ExecuteOrbifoldCommand(string command)
{
  string tmp_string1 = "";
  string tmp_string2 = "";
  string tmp_string3 = "";
  string parameter_string1 = "";
  string parameter_string2 = "";
  string parameter_string3 = "";
  size_t t1 = 0;
  size_t t2 = 0;
  unsigned i = 0;
  unsigned j = 0;
  unsigned k = 0;

  bool case1 = false;
  bool case2 = false;

  CAnalyseModel Analyse;

  // step 6: child processes
  // step 4: commands available in all directories
  // step 5: commands available in specific folders


  // if the current directory is />
  if (this->current_folder[0] < 0)
  {

    if (command.substr(0,3) == "man")
    {
      if (command.length() >=4)
      {

        string path_doc = " ./doc/main/";
        // Insert path_doc en la position 4
        command.insert(4, path_doc);
        // Add ".man" to end
        command += ".man";


        // Redirect error output to a temporary file
        string temp_file = "temp_error.txt";
        string full_command = command + " 2>" + temp_file;

        // Execute the command
        int result = system(full_command.c_str());


        ifstream error_file(temp_file);
        stringstream error_stream;
        error_stream << error_file.rdbuf();
        string error_message = error_stream.str();

        // Delete temporary file
        error_file.close();
        remove(temp_file.c_str());

        // Verifity the result

        if (result !=0|| containsErrorKeywords(error_message,error_keywords)){
            cout << "Error: A command name is expected after the instruction man.\n"
                 << "Options are:\n* cd\n* create\n* delete\n* load\n* rename\n* save\n";
            return false;
        }


        return true;
      }
      else
        {
            cout << "Error: A command name is expected after the instruction man.\n"
                 << "Options are:\n* cd\n* create\n* delete\n* load\n* rename\n* save\n";
            return false;
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // change directory
    // updated on 10.10.2011
    if (this->FindCommandType1(command, "cd ", parameter_string1))
    {
      unsigned index = 0;
      if (this->MessageOrbifoldNotKnown(parameter_string1, index))
        return true;

      this->OrbifoldIndex = index;
      this->current_folder[0] = index;

      const COrbifold &Orbifold = this->Orbifolds[this->OrbifoldIndex];
      if (Orbifold.GetCheckStatus() == CheckedAndGood)
      {
        const SelfDualLattice Lattice = Orbifold.OrbifoldGroup.GetLattice();

        if (Lattice == E8xE8)
          SELFDUALLATTICE = 1;
        return true;
      }

      if (Orbifold.OrbifoldGroup.GetSpaceGroup_CheckStatus() == CheckedAndGood)
        this->MessageHelpCreateNewOrbifold(3);
      else
        this->MessageHelpCreateNewOrbifold();
      this->current_folder[1] = 1;
      return true;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // rename orbifold
    if (this->FindCommandType2(command, "rename orbifold(", parameter_string1, parameter_string2))
    {
      // begin: find parameters
      if (!this->FindParameterType2(parameter_string2, "to(", tmp_string1))
      {
        if (this->print_output)
          (*this->Print.out) << "\n  " << this->Print.cbegin << "New label of orbifold not specified." << this->Print.cend << "\n" << endl;

        this->MessageParameterNotKnown(parameter_string2);
        return true;
      }
      // end: find parameters

      // begin: from orbifold
      unsigned index = 0;
      if (this->MessageOrbifoldNotKnown(parameter_string1, index))
        return true;
      // end: from orbifold

      // begin: to orbifold
      if (this->MessageLabelError(tmp_string1))
      {
        this->MessageParameterNotKnown(parameter_string2);
        return true;
      }

      if (this->MessageOrbifoldAlreadyExists(tmp_string1))
      {
        if (!this->MessageParameterNotKnown(parameter_string2))
          (*this->Print.out) << endl;
        return true;
      }
      // end: to orbifold

      this->Orbifolds[index].OrbifoldGroup.Label = tmp_string1;

      if (this->print_output)
        (*this->Print.out) << "\n  " << this->Print.cbegin << "Orbifold \""  << parameter_string1 << "\" renamed to \"" << tmp_string1 << "\"." << this->Print.cend << "\n" << endl;

      this->MessageParameterNotKnown(parameter_string2);
      return true;
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // load orbifolds
    if (this->FindCommandType2(command, "load orbifolds(", parameter_string1, parameter_string2)
     || this->FindCommandType2(command, "load orbifold(", parameter_string1, parameter_string2))
    {
      const bool inequivalent      = this->FindParameterType1(parameter_string2, "inequivalent");
      (*this->Print.out) << "\n";
      this->LoadOrbifolds(parameter_string1, inequivalent, false);
      this->MessageParameterNotKnown(parameter_string2);
      return true;
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // save orbifolds //sept16
    if (this->FindCommandType2(command, "save orbifolds(", parameter_string1, parameter_string2)
     || this->FindCommandType2(command, "save orbifold(", parameter_string1, parameter_string2))
    {
      t1 = this->Orbifolds.size();
      if (t1 == 0)
      {
        if (this->print_output)
          (*this->Print.out) << "\n  " << this->Print.cbegin << "No orbifold to save." << this->Print.cend << "\n" << endl;
        this->MessageParameterNotKnown(parameter_string2);
        return true;
      }

      ofstream out(parameter_string1.data());
      if((!out.is_open()) || (!out.good()))
      {
        (*this->Print.out) << "\n  " << this->Print.cbegin << "File \"" << parameter_string1 << "\" not found." << this->Print.cend << "\n" << endl;

        this->MessageParameterNotKnown(parameter_string2);
        return true;
      }

      unsigned counter = 0;
      for (i = 0; i < t1; ++i)
      {
        if (this->Orbifolds[i].GetCheckStatus() == CheckedAndGood)
        {
          this->Orbifolds[i].OrbifoldGroup.PrintToFile(out);
          out << endl;
          ++counter;
        }
      }
      out.close();

      if (this->print_output)
      {
        if (counter == 0)
          (*this->Print.out) << "\n  " << this->Print.cbegin << "No orbifolds saved to file \"" << parameter_string1 << "\"." << this->Print.cend << "\n" << endl;
        else
        if (counter == 1)
          (*this->Print.out) << "\n  " << this->Print.cbegin << "One orbifold saved to file \"" << parameter_string1 << "\"." << this->Print.cend << "\n" << endl;
        else
          (*this->Print.out) << "\n  " << this->Print.cbegin << counter << " orbifolds saved to file \"" << parameter_string1 << "\"." << this->Print.cend << "\n" << endl;
      }

      this->MessageParameterNotKnown(parameter_string2);
      return true;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // create orbifold
    // updated on 21.10.2011
    if (this->FindCommandType2(command, "create orbifold(", parameter_string1, parameter_string2))
    {
      // begin: find parameters
      const bool b1 = this->FindParameterType2(parameter_string2, "from(", parameter_string3);
      const bool b2 = this->FindParameterType2(parameter_string2, "with point group(", parameter_string3);

      if ((!b1 && !b2) || (b1 && b2))
      {
        if (this->print_output)
          (*this->Print.out) << "\n  " << this->Print.cbegin << "One parameter needed: \"from(AnotherOrbifoldLabel)\" or \"with point group(M,N)\"." << this->Print.cend << "\n" << endl;

        this->MessageParameterNotKnown(parameter_string2);
        return true;
      }
      // end: find parameters

      // begin: new orbifold label
      if (this->MessageLabelError(parameter_string1))
      {
        this->MessageParameterNotKnown(parameter_string2);
        return true;
      }

      if (this->MessageOrbifoldAlreadyExists(parameter_string1))
      {
        (*this->Print.out) << endl;
        this->MessageParameterNotKnown(parameter_string2);
        return true;
      }
      // end: new orbifold label

      // case 1: "from(AnotherOrbifoldLabel)"
      if (b1)
      {
        // begin: from orbifold
        unsigned index = 0;
        if (this->MessageOrbifoldNotKnown(parameter_string3, index))
          return true;
        // end: from orbifold

        this->Orbifolds.push_back(this->Orbifolds[index]);
        this->Orbifolds[this->Orbifolds.size()-1].OrbifoldGroup.Label = parameter_string1;

        vector<SConfig> NewVEVConfigs = this->AllVEVConfigs[index];
        t1 = NewVEVConfigs.size();
        for (i = 0; i < t1; ++i)
        {
          PID &pid = NewVEVConfigs[i].pid;
          pid.PIDs.clear();
          pid.PID_Commands.clear();
          pid.PID_JobIndices.clear();
          pid.PID_StartingTimes.clear();
          pid.PID_Done.clear();
          pid.PID_Filenames.clear();
        }
        this->AllVEVConfigs.push_back(NewVEVConfigs);
        this->AllVEVConfigsIndex.push_back(this->AllVEVConfigsIndex[index]);

        if (this->print_output)
          (*this->Print.out) << "\n  " << this->Print.cbegin << "Orbifold \"" << parameter_string1 << "\" created from \"" << parameter_string3 << "\"." << this->Print.cend << "\n" << endl;

        this->MessageParameterNotKnown(parameter_string2);
        return true;
      }

      // case 2: "with point group(M,N)"
      vector<int> Orders;

      if (!convert_string_to_vector_of_int(parameter_string3, Orders) || ((Orders.size() != 1) && (Orders.size() != 2)) || (find(Orders.begin(), Orders.end(), 0) != Orders.end()))
      {
        if (this->print_output)
          (*this->Print.out) << "\n  " << this->Print.cbegin << "The point group of the orbifold is ill-defined." << this->Print.cend << "\n" << endl;

        this->MessageParameterNotKnown(parameter_string2);
        return true;
      }

       int ZM = 2;
       int ZN = 1;
       int ZK = 1;
      if (Orders.size() == 1)
       {
       ZN = Orders[0];
       }
      else
      {
        ZN = Orders[0];
        ZK = Orders[1];
      }

      COrbifold NewOrbifold;
      CSpaceGroup &SpaceGroup = NewOrbifold.OrbifoldGroup.AccessSpaceGroup();
      SpaceGroup.Clear();
      SpaceGroup.SetOrder(ZM, ZN, ZK);


      NewOrbifold.OrbifoldGroup.Label = parameter_string1;
      this->Orbifolds.push_back(NewOrbifold);

      SConfig TestConfig = NewOrbifold.StandardConfig;
      TestConfig.ConfigLabel = "TestConfig";
      TestConfig.ConfigNumber = 1;

      vector<SConfig> Configs;
      Configs.push_back(NewOrbifold.StandardConfig);
      Configs.push_back(TestConfig);
      this->AllVEVConfigs.push_back(Configs);
      this->AllVEVConfigsIndex.push_back(1);


     if (this->print_output)
      {
        (*this->Print.out) << "\n  " << this->Print.cbegin << "Orbifold with point group ";
        if (SpaceGroup.IsZMxZN())
          (*this->Print.out) << "Z" << SpaceGroup.GetN();
        if (SpaceGroup.IsZMxZNxZK())
          (*this->Print.out) << "xZ" << SpaceGroup.GetK();

        (*this->Print.out) << " created and stored in directory \"" << parameter_string1 << "\"." << this->Print.cend << "\n" << endl;
        (*this->Print.out) << "  " << this->Print.cbegin << "Use the command \"cd " << parameter_string1 << "\" to change the directory to the new orbifold." << this->Print.cend << "\n" << endl;
      }


      this->MessageParameterNotKnown(parameter_string2);
      return true;
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // create random orbifold
    // updated on 22.02.2012
    if (this->FindCommandType2(command, "create random orbifold from(", parameter_string1, parameter_string2))
    {
 // begin: find parameters

      bool save_if_new = false;
      bool save_all    = true;
      bool save_SM     = false;
      bool save_PS     = false;
      bool save_SU5    = false;

      unsigned number_of_generations = 3;
      bool generations_specified = false;


      bool Use_Filename = false;
      string Filename = "";

      const bool check_anomalies   = !this->FindParameterType1(parameter_string2, "do not check anomalies");
      const bool load_when_done    =  this->FindParameterType1(parameter_string2, "load when done");
      const bool print_info        =  this->FindParameterType1(parameter_string2, "print info");


      if (this->FindParameterType2(parameter_string2, "if(", parameter_string3))
      {
        save_SM     = this->FindParameterType1(parameter_string3, "SM");
        save_PS     = this->FindParameterType1(parameter_string3, "PS");
        save_SU5    = this->FindParameterType1(parameter_string3, "SU5");

        generations_specified = this->FindParameterType4(parameter_string3, "generations", number_of_generations);

        if (generations_specified && (!save_SM && !save_PS && !save_SU5))
        {
          if (this->print_output)
            (*this->Print.out) << "\n  " << this->Print.cbegin << "Cannot look for models with " << number_of_generations << " generations, because gauge group not specified." << this->Print.cend << flush;
          generations_specified = false;
        }

        if (save_SM || save_PS || save_SU5)
          save_all = false;

        save_if_new = this->FindParameterType1(parameter_string3, "inequivalent");
      }

      if (this->FindParameterType2(parameter_string2, "save to(", parameter_string3))
      {
        Use_Filename = true;
        Filename = parameter_string3;
      }

      const bool create_Orbifold = (!save_all || save_if_new || check_anomalies || print_info);

// begin: use original shifts and Wilson lines
      CRandomModel RandomModel(E8xE8);
      vector<bool> UseOrigShiftsAndWilsonLines(9, false);
      UseOrigShiftsAndWilsonLines[0]=true; //use Witten's shift


      if (this->FindParameterType2(parameter_string2, "use(", parameter_string3))
      {
        vector<string>   tmp_strings;
        vector<unsigned> tmp_unsigneds;
        global_DecomposeString(parameter_string3, tmp_strings, ",");
        global_ConvertVectorOfString2VectorOfUnsigned(tmp_strings, tmp_unsigneds);

        bool vector_ok = true;
        if (tmp_unsigneds.size() != 8)
          vector_ok = false;
        for (i = 0; vector_ok && (i < 8); ++i)
        {
          if ((tmp_unsigneds[i] != 0) && (tmp_unsigneds[i] != 1))
            vector_ok = false;
        }
        bool all_entries_one = true;
        for (i = 0; all_entries_one && (i < 8); ++i)
        {
          if (tmp_unsigneds[i] != 1)
            all_entries_one = false;
        }
        if (all_entries_one)
          vector_ok = false;

        if (!vector_ok)
        {
          if (this->print_output)
            (*this->Print.out) << "\n  " << this->Print.cbegin << "Parameter \"" << parameter_string3 << "\" of \"use(...)\" ill-defined." << this->Print.cend << "\n" << endl;
          return true;
        }

        for (i = 0; i < 8; ++i)
          UseOrigShiftsAndWilsonLines[i+1] = (tmp_unsigneds[i] == 1);
      }
// end: use original shifts and Wilson lines


      unsigned max_models = 1;

      if (this->FindParameterType2(parameter_string2, "#models(", parameter_string3))
      {
        if (parameter_string3 == "all")
          max_models = 250000000;
        else
        {
          if (parameter_string3.find_first_not_of("0123456789") != string::npos)
          {
            if (this->print_output)
              (*this->Print.out) << "\n  " << this->Print.cbegin << "Parameter \"" << parameter_string3 << "\" of \"#models(X)\" ill-defined." << this->Print.cend << "\n" << endl;
            return true;
          }
          max_models = (unsigned)atoi(parameter_string3.c_str());
        }
      }

      // end: find parameters

      unsigned OriginalOrbifoldIndex = 0;
      vector<unsigned> OrbifoldIndices;

      bool random_origin = false;

      if (parameter_string1 == "*")
      {
        random_origin = true;

        const size_t o1 = this->Orbifolds.size();
        for (i = 0; i < o1; ++i)
        {
          if (this->Orbifolds[i].GetCheckStatus() == CheckedAndGood)
            OrbifoldIndices.push_back(i);
        }
        if (OrbifoldIndices.size() == 0)
        {
          this->MessageParameterNotKnown(parameter_string2);
          if (this->print_output)
            (*this->Print.out) << "\n  " << this->Print.cbegin << "Cannot create orbifolds randomly from \"*\" because no valid orbifold model available." << this->Print.cend << "\n" << endl;
          return true;
        }
        OriginalOrbifoldIndex = OrbifoldIndices[0];
      }
      else
      {
        if (this->MessageOrbifoldNotKnown(parameter_string1, OriginalOrbifoldIndex))
        {
          this->MessageParameterNotKnown(parameter_string2);
          return true;
        }
        if (this->Orbifolds[OriginalOrbifoldIndex].GetCheckStatus() != CheckedAndGood)
        {
          this->MessageParameterNotKnown(parameter_string2);
          if (this->print_output)
            (*this->Print.out) << "\n  " << this->Print.cbegin << "Cannot create orbifolds randomly from \"" << parameter_string1 << "\" because model is not fully defined." << this->Print.cend << "\n" << endl;
          return true;
        }
      }

      if (!Use_Filename && !load_when_done && !print_info)
      {
        if (this->print_output)
          (*this->Print.out) << "\n  " << this->Print.cbegin << "Models neither saved nor displayed. Parameter \"load when done\" turned on." << this->Print.cend << flush;
        command += " load when done";
      }

      vector<SUSYMultiplet> Multiplets(2);
	  Multiplets[0]=Scalar;
	  Multiplets[1]=LeftFermi;

      int newPID = -1;
      newPID = fork();

      if (newPID < 0)  /* error occurred */
      {
        (*this->Print.out) << "\n  " << this->Print.cbegin << "Fork failed!" << this->Print.cend << "\n" << endl;
        return true;
      }

      if (newPID == 0) /* child process */
      {

        if (!Use_Filename)
        {
          Filename = "tmp_fileID";
          std::ostringstream osID;
          osID << getpid();
          Filename += osID.str();
          Filename += ".txt";
        }

        const double Rmax = RAND_MAX+0.000001;

        std::ofstream tmp_out(Filename.data());


        CInequivalentModels InequivModels;

        COrbifoldGroup NewOrbifoldGroup;
        vector<CVector> UnbrokenRoots;

        vector<CRandomModel>     RandomModels;
        vector<CSector>          CoreSpectrum;
        vector<vector<CSector> > CoreSpectra;

        unsigned SpecialIndex = 0;

        const size_t o1 = OrbifoldIndices.size();
        if (random_origin)
        {
          for (i = 0; i < o1; ++i)
          {
            NewOrbifoldGroup = this->Orbifolds[OrbifoldIndices[i]].OrbifoldGroup;
            NewOrbifoldGroup.LoadedU1Generators.clear();

            RandomModel.Initiate(NewOrbifoldGroup, UseOrigShiftsAndWilsonLines, UnbrokenRoots);
            RandomModels.push_back(RandomModel);

            CoreSpectrum.clear();
            COrbifoldCore OrbifoldCore(NewOrbifoldGroup, CoreSpectrum);
            CoreSpectra.push_back(CoreSpectrum);
          }
        }
        else
        {
          NewOrbifoldGroup = this->Orbifolds[OriginalOrbifoldIndex].OrbifoldGroup;
          NewOrbifoldGroup.LoadedU1Generators.clear();

           RandomModel.Initiate(NewOrbifoldGroup, UseOrigShiftsAndWilsonLines, UnbrokenRoots);
           RandomModels.push_back(RandomModel);

          CoreSpectrum.clear();
          COrbifoldCore OrbifoldCore(NewOrbifoldGroup, CoreSpectrum);
          CoreSpectra.push_back(CoreSpectrum);

          SpecialIndex = 0;
        }

        vector<SConfig> AllVEVConfigs;

        CSpectrum PrintSpectrum;

        bool save_current_model = false;
        string model_string = "";

        bool model_with_problem = false;

        unsigned problem_counter = 0;

        i = 1;
        while (i <= max_models)
        {
          if (random_origin)
          {
            OriginalOrbifoldIndex = (unsigned)(rand() * (double)o1/Rmax);
            if (OriginalOrbifoldIndex >= o1)
            {
              (*this->Print.out) << "\n  " << this->Print.cbegin << "Warning: orbifold-index out of range. Stop process." << this->Print.cend << endl;
              tmp_out.close();
              std::exit(0); // terminate child process
            }
            NewOrbifoldGroup = this->Orbifolds[OrbifoldIndices[OriginalOrbifoldIndex]].OrbifoldGroup;
            NewOrbifoldGroup.LoadedU1Generators.clear();

            SpecialIndex = OriginalOrbifoldIndex;
          }
          // create random shifts and Wilson lines
          if (NewOrbifoldGroup.CreateRandom(RandomModels[SpecialIndex], false))
          {

            model_with_problem = false;
            save_current_model = true;
            model_string = "Random";

            if (NewOrbifoldGroup.GetModularInvariance_CheckStatus() != CheckedAndGood)
            {
              (*this->Print.out) << "\n  " << this->Print.cbegin << "Warning: problems with modular invariance." << this->Print.cend << endl;
              model_with_problem = true;
              model_string += "_ProblemModInv";
            }

          if (create_Orbifold && !model_with_problem)
            {
              COrbifold NewOrbifold(NewOrbifoldGroup, CoreSpectra[SpecialIndex]);

              if (NewOrbifold.GetCheckStatus() != CheckedAndGood)
              {
                (*this->Print.out) << "\n  " << this->Print.cbegin << "Warning: problems with constructing the orbifold spectrum." << this->Print.cend << endl;
                model_with_problem = true;
                model_string += "_ProblemOrbi";
              }


             if (NewOrbifold.TachyonicStandardConfig.Fields.size() == 0)
             {
              if (!save_all && !model_with_problem)
              {
                AllVEVConfigs.clear();

                bool copy_save_SM  = save_SM;
                bool copy_save_PS  = save_PS;
                bool copy_save_SU5 = save_SU5;
                save_current_model = Analyse.AnalyseModel(NewOrbifold, NewOrbifold.StandardConfig, copy_save_SM, copy_save_PS, copy_save_SU5, AllVEVConfigs, this->Print, number_of_generations, false);
                if (save_current_model)
                {
                  model_string = "Model";
                  if (copy_save_SM)
                    model_string += "_SM";
                  if (copy_save_PS)
                    model_string += "_PS";
                  if (copy_save_SU5)
                    model_string += "_SU5_";
                }
              }
             }
             else
             {
			  save_current_model = false;
			 }

              if (save_current_model && (save_if_new || print_info) && !model_with_problem)
              {
                CSpectrum Spectrum(NewOrbifold.StandardConfig, Multiplets);

                if (save_if_new)
                {
                  save_current_model = InequivModels.IsSpectrumUnknown(Spectrum, true);
                }
                 if (save_current_model && print_info)
                 PrintSpectrum = Spectrum;
              }

              if (check_anomalies && !model_with_problem && (!save_if_new || (save_if_new && save_current_model)))
              {
                if ((!NewOrbifold.CheckAnomaly(NewOrbifold.StandardConfig, this->GaugeIndices, this->Print, false) ))
                {
                  //(*this->Print.out) << "\n  " << this->Print.cbegin << "Warning: problems with discrete and/or gauge anomalies." << this->Print.cend << endl;
                  //model_with_problem = true;
                  //model_string += "_ProblemAnomaly";
                  save_current_model = false;
                }
              }
            }

            // begin: save the current model
            if (save_current_model || model_with_problem)
            {
              std::ostringstream os;
              os << i;
              NewOrbifoldGroup.Label = model_string + os.str();
              NewOrbifoldGroup.PrintToFile(tmp_out);
              ++i;

              if (!model_with_problem && print_info)
              {
                (*this->Print.out) << "\n  " << this->Print.cbegin;
                if (random_origin)
                  NewOrbifoldGroup.GetSpaceGroup().PrintPointGroup(*this->Print.out);
                else
                  (*this->Print.out) << "Orbifold";

                (*this->Print.out) << " model \"" << NewOrbifoldGroup.Label << "\" from \"" << parameter_string1 << "\":" << this->Print.cend << "\n";

                this->Print.PrintSpectrum(PrintSpectrum);
                (*this->Print.out) << flush;
              }
            }

            if (model_with_problem)
            {
              (*this->Print.out) << "  " << this->Print.cbegin << "Orbifold label: " << NewOrbifoldGroup.Label << this->Print.cend << endl;
              ++problem_counter;
            }
            // end: save the current model
          }
        }

       tmp_out.close();

        if (problem_counter == 0)
          (*this->Print.out) << "\n  " << this->Print.cbegin << "Models created without problems." << this->Print.cend << flush;
        else
        {
          if (problem_counter == 1)
            (*this->Print.out) << "\n  " << this->Print.cbegin << "Warning: 1 model created with problems." << this->Print.cend << flush;
          else
            (*this->Print.out) << "\n  " << this->Print.cbegin << "Warning: " << problem_counter << " models created with problems." << this->Print.cend << flush;
        }

        // terminate child process
        std::exit(0);
      }
      /* parent process */
      usleep(100);

      if (!Use_Filename)
      {
        Filename = "tmp_fileID";
        std::ostringstream osID;
        osID << newPID;
        Filename += osID.str();
        Filename += ".txt";
      }

      (*this->Print.out) << "\n  " << this->Print.cbegin << "New child process \"PID " << newPID << "\" from command \"" << command << "\"." << this->Print.cend << "\n" << endl;

       if (print_info)
       (*this->Print.out) << "  " << this->Print.cbegin << "Note that details of newly created orbifold models can only be seen after process \"PID " << newPID << "\" has finished." << this->Print.cend << "\n" << endl;

// save info about the process to the orbifold that id the origin of the random models
      PID &PID_Data = this->AllVEVConfigs[OriginalOrbifoldIndex][this->AllVEVConfigsIndex[OriginalOrbifoldIndex]].pid;
      PID_Data.PIDs.push_back(newPID);
      PID_Data.PID_Commands.push_back(command);
      PID_Data.PID_JobIndices.push_back(3); // create random orbifold
      PID_Data.PID_StartingTimes.push_back(time (NULL));
      PID_Data.PID_Done.push_back(false);
      PID_Data.PID_Filenames.push_back(Filename);

      this->MessageParameterNotKnown(parameter_string2);
      return true;
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // delete orbifold
    if (this->FindCommandType1(command, "delete orbifold", parameter_string1))
    {
      bool process_running = false;
      if (this->FindParameterType2(parameter_string1, "(", tmp_string1))
      {
        unsigned index = 0;
        if (!this->MessageOrbifoldNotKnown(tmp_string1, index))
        {
          const vector<SConfig> &VEVConfigs = this->AllVEVConfigs[index];
          t1 = VEVConfigs.size();
          for (i = 0; !process_running && (i < t1); ++i)
          {
            const vector<bool> &PID_Done = VEVConfigs[i].pid.PID_Done;
            if (find(PID_Done.begin(), PID_Done.end(), false) != PID_Done.end())
              process_running = true;
          }
          if (process_running)
          {
            if (this->print_output)
              (*this->Print.out) << "\n  " << this->Print.cbegin << "Processes of orbifold \"" << tmp_string1 << "\" are still running. Cannot be deleted." << this->Print.cend << "\n" << endl;
          }
          else
          {
            this->Orbifolds.erase(this->Orbifolds.begin() + index);
            this->AllVEVConfigs.erase(this->AllVEVConfigs.begin() + index);
            this->AllVEVConfigsIndex.erase(this->AllVEVConfigsIndex.begin() + index);
            this->OrbifoldIndex = -1;
            if (this->print_output)
              (*this->Print.out) << "\n  " << this->Print.cbegin << "Orbifold \"" << tmp_string1 << "\" deleted." << this->Print.cend << "\n" << endl;
          }
        }
      }
      else
      if (this->FindParameterType1(parameter_string1, "s"))
      {
        t1 = this->Orbifolds.size();
        if (t1 == 0)
        {
          if (this->print_output)
            (*this->Print.out) << "\n  " << this->Print.cbegin << "No orbifold to delete." << this->Print.cend << "\n" << endl;
          this->MessageParameterNotKnown(parameter_string1);
          return true;
        }

        unsigned counter = 0;
        for (int i = t1-1; i >= 0; --i)
        {
          const vector<SConfig> &VEVConfigs = this->AllVEVConfigs[i];
          t2 = VEVConfigs.size();

          process_running = false;
          for (j = 0; !process_running && (j < t2); ++j)
          {
            const vector<bool> &PID_Done = VEVConfigs[j].pid.PID_Done;
            if (find(PID_Done.begin(), PID_Done.end(), false) != PID_Done.end())
              process_running = true;
          }
          if (process_running)
          {
            if (this->print_output)
              (*this->Print.out) << "\n  " << this->Print.cbegin << "Processes of orbifold \"" << this->Orbifolds[i].OrbifoldGroup.Label << "\" are still running. Cannot be deleted." << this->Print.cend;
          }
          else
          {
            ++counter;
            this->Orbifolds.erase(this->Orbifolds.begin() + i);
            this->AllVEVConfigs.erase(this->AllVEVConfigs.begin() + i);
            this->AllVEVConfigsIndex.erase(this->AllVEVConfigsIndex.begin() + i);
          }
        }

        this->OrbifoldIndex = -1;

        if (this->print_output)
        {
          if (t1 == counter)
            (*this->Print.out) << "\n  " << this->Print.cbegin << "All orbifolds deleted." << this->Print.cend << "\n" << endl;
          else
          {
            if (counter == 0)
              (*this->Print.out) << "\n  " << this->Print.cbegin << "No orbifold could be deleted." << this->Print.cend << "\n" << endl;
            else
            if (counter == 1)
              (*this->Print.out) << "\n  " << this->Print.cbegin << "Only one orbifold could be deleted." << this->Print.cend << "\n" << endl;
            else
              (*this->Print.out) << "\n  " << this->Print.cbegin << "Only " << counter << " orbifolds could be deleted." << this->Print.cend << "\n" << endl;
          }
        }
      }
      this->MessageParameterNotKnown(parameter_string1);
      return true;
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // show directories
    if (this->FindCommandType1(command, "dir", parameter_string1) || this->FindCommandType1(command, "help", parameter_string1) || this->FindCommandType1(command, "ll", parameter_string1))
    {
      if (this->FindParameterType1(parameter_string1, "processes"))
       this->PrintCommandsProcesses();
      else
      if (this->FindParameterType1(parameter_string1, "create random"))
      {
        (*this->Print.out) << "\n  create random orbifold from(OrbifoldLabel)\n";
        (*this->Print.out) << "  parameters:\n";
        (*this->Print.out) << "    \"save to(Filename)\"                     save result to model-file \"Filename\"\n";
        (*this->Print.out) << "    \"if(...)\"                               specify desired properties of the models to search for:\n";
        (*this->Print.out) << "                                                \"SM\", \"PS\", \"SU5\",\n";
        (*this->Print.out) << "                                                \"inequivalent\"\n";
        (*this->Print.out) << "    \"use(1,1,0,1,...)\"                      eight digits for two shifts and six Wilson lines:\n";
        (*this->Print.out) << "                                                \"1\" : take shift/WL from original model\n";
        (*this->Print.out) << "                                                \"0\" : create shift/WL randomly\n";
        (*this->Print.out) << "    \"#models(X)\"                            create \"X\" models with specified properties;\n";
        (*this->Print.out) << "                                              use \"X\" = \"all\" to create as many as possible\n";
        (*this->Print.out) << "    \"print info\"                            print summary of spectrum\n";
        (*this->Print.out) << "    \"load when done\"                        after process finished, load new models\n";
        (*this->Print.out) << "    \"do not check anomalies\"                speeds up the process\n";
        (*this->Print.out) << "  examples:\n";
        (*this->Print.out) << "    create random orbifold from(A) if(inequivalent SM) save to(File1.txt) #models(all)\n";
        (*this->Print.out) << "    create random orbifold from(B) if(inequivalent) save to(File2.txt) #models(100) use(1,1,0,0,0,0,0,0)\n\n" << flush;
      }
      else
      if (this->FindParameterType1(parameter_string1, "system commands"))
      {
        (*this->Print.out) << "\n  system commands:\n";
        (*this->Print.out) << "    @status\n";
        //if (!this->online_mode)
        //{
          (*this->Print.out) << "    @print enter                              print new line\n";
          (*this->Print.out) << "    @print(string)                            optional: unformatted\n";
          (*this->Print.out) << "    @begin print to file(Filename)\n";
          (*this->Print.out) << "    @end print to file\n";
        //}
        (*this->Print.out) << "    @typesetting(Type)                       \"Type\" can be \"mathematica\", \"latex\" or \"standard\" \n\n" << flush;
      }
      else
      {
        const bool PrintSubDir = !this->FindParameterType1(parameter_string1, "no subdirectories");

        (*this->Print.out) << "\n  commands of this directory:\n";
        (*this->Print.out) << "    load program(Filename)                    load commands from \"Filename\"\n\n";
        (*this->Print.out) << "    load orbifolds(Filename)                  optional: \"inequivalent\"\n";
        (*this->Print.out) << "    save orbifolds(Filename)\n\n";
        (*this->Print.out) << "    create orbifold(OrbifoldLabel) with point group(M,N)\n";
        (*this->Print.out) << "    create orbifold(OrbifoldLabel) from(AnotherOrbifoldLabel)\n";
        (*this->Print.out) << "    create random orbifold from(OrbifoldLabel)\n";
        (*this->Print.out) << "                                              various parameters, see \"help create random\"\n";
        (*this->Print.out) << "    rename orbifold(OldOrbifoldLabel) to(NewOrbifoldLabel)\n";
        (*this->Print.out) << "    delete orbifold(OrbifoldLabel)\n";
        (*this->Print.out) << "    delete orbifolds\n\n";

        if (PrintSubDir)
        {
          t1 = this->Orbifolds.size();
          if (t1 != 0)
          {
            (*this->Print.out) << "  change directory:\n";
            for (i = 0; i < t1; ++i)
              (*this->Print.out) << "    cd " << this->Orbifolds[i].OrbifoldGroup.Label << "\n";
            (*this->Print.out) << "\n";
          }
        }

        (*this->Print.out) << "  general commands:\n";
        (*this->Print.out) << "    dir                                       show commands; optional: \"no subdirectories\"\n";
        (*this->Print.out) << "    help                                      optional: \"create random\", \"system commands\", \"processes\"\n";
        if (!this->online_mode)
          (*this->Print.out) << "    exit                                      exit program\n";
        (*this->Print.out) << "\n" << flush;;
      }

      this->MessageParameterNotKnown(parameter_string1);
      return true;
    }
  }
  else
  {
    // STEP 4 //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // begin: commands available in all orbifold directories
    COrbifold          &Orbifold       = this->Orbifolds[this->OrbifoldIndex];
    vector<SConfig>    &VEVConfigs     = this->AllVEVConfigs[this->OrbifoldIndex];
    unsigned           &VEVConfigIndex = this->AllVEVConfigsIndex[this->OrbifoldIndex];
    SConfig            &VEVConfig      = VEVConfigs[VEVConfigIndex];
    vector<CField>     &Fields         = VEVConfig.Fields;

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // change directory
    // updated on 29.06.2011
    if (this->FindCommandType1(command, "cd ~", parameter_string1))
    {
      this->current_folder[0] = -1;
      this->current_folder[1] = 0;
      this->current_folder[2] = 0;

      this->MessageParameterNotKnown(parameter_string1);
      return true;
    }

    if (this->FindCommandType0(command, "m"))
    {
      if (Orbifold.GetCheckStatus() != CheckedAndGood)
      {
        if (Orbifold.OrbifoldGroup.GetSpaceGroup_CheckStatus() == CheckedAndGood)
          this->MessageHelpCreateNewOrbifold(3);
        else
          this->MessageHelpCreateNewOrbifold();

        this->current_folder[1] = 1;
        return true;
      }
      this->current_folder[1] = 1;
      return true;
    }

    if (this->FindCommandType0(command, "gg"))
    {
      if (Orbifold.GetCheckStatus() != CheckedAndGood)
      {
        if (Orbifold.OrbifoldGroup.GetSpaceGroup_CheckStatus() == CheckedAndGood)
          this->MessageHelpCreateNewOrbifold(3);
        else
          this->MessageHelpCreateNewOrbifold();

        this->current_folder[1] = 1;
        return true;
      }
      this->current_folder[1] = 2;
      return true;
    }

    if (this->FindCommandType0(command, "s"))
    {
      if (Orbifold.GetCheckStatus() != CheckedAndGood)
      {
        if (Orbifold.OrbifoldGroup.GetSpaceGroup_CheckStatus() == CheckedAndGood)
          this->MessageHelpCreateNewOrbifold(3);
        else
          this->MessageHelpCreateNewOrbifold();

        this->current_folder[1] = 1;
        return true;
      }
      this->current_folder[1] = 3;
      return true;
    }

    if (this->FindCommandType0(command, "v"))
    {
      if (Orbifold.GetCheckStatus() != CheckedAndGood)
      {
        if (Orbifold.OrbifoldGroup.GetSpaceGroup_CheckStatus() == CheckedAndGood)
          this->MessageHelpCreateNewOrbifold(3);
        else
          this->MessageHelpCreateNewOrbifold();

        this->current_folder[1] = 1;
        return true;
      }
      this->current_folder[1] = 5;
      this->current_folder[2] = 0;
      return true;
    }

    if (this->FindCommandType0(command, "l"))
    {
      if (Orbifold.GetCheckStatus() != CheckedAndGood)
      {
        if (Orbifold.OrbifoldGroup.GetSpaceGroup_CheckStatus() == CheckedAndGood)
          this->MessageHelpCreateNewOrbifold(3);
        else
          this->MessageHelpCreateNewOrbifold();

        this->current_folder[1] = 1;
        return true;
      }
      this->current_folder[1] = 5;
      this->current_folder[2] = 1;
      return true;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // create set or monomial
    // updated on 26.01.2011
    if ((case1 = this->FindCommandType1(command, "create set", parameter_string1)))
    {
      unsigned index = 0;

      // find the label of the new set or new monomial and save it into the string "tmp_string1"
      if (this->FindParameterType2(parameter_string1, "(", tmp_string1))
      {
        if (this->MessageLabelError(tmp_string1))
          return true;

        vector<string> &NamesOfSetsOfFields = VEVConfig.NamesOfSetsOfFields;
        vector<vector<unsigned> > &SetsOfFields   = VEVConfig.SetsOfFields;

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // create an empty set or create a set from a monomial
        // updated on 19.04.2011
        if (case1)
        {
          if (this->MessageXAlreadyExists(NamesOfSetsOfFields, "Set", tmp_string1))
            return true;

          vector<unsigned> SetOfIndices;
          vector<unsigned> tmp_SetOfFields;

          t1 = tmp_SetOfFields.size();
          if (t1 == 0)
          {
            if (this->print_output)
              (*this->Print.out) << "\n  " << this->Print.cbegin << "Empty set \"" << tmp_string1 << "\"";
          }
          else
          {
            stable_sort(tmp_SetOfFields.begin(), tmp_SetOfFields.end());
            if (this->print_output)
            {
              (*this->Print.out) << "  " << this->Print.cbegin << "Set \"" << tmp_string1 << "\" with " << t1 << " field";
              if (t1 != 1) (*this->Print.out) << "s";
            }
          }
          if (this->print_output)
            (*this->Print.out) << " created." << this->Print.cend << "\n" << endl;

          NamesOfSetsOfFields.push_back(tmp_string1);
          SetsOfFields.push_back(tmp_SetOfFields);
        }
      }
      this->MessageParameterNotKnown(parameter_string1);
      return true;
    }



    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // delete set or monomial
    // updated on 19.04.2011
    if ((case1 = this->FindCommandType1(command, "delete set", parameter_string1)))
    {
      unsigned index = 0;
      // delete a set or all sets
      if (case1)
      {
        vector<string>            &NamesOfSetsOfFields = VEVConfig.NamesOfSetsOfFields;
        vector<vector<unsigned> > &SetsOfFields        = VEVConfig.SetsOfFields;

        // delete set
        // updated on 19.04.2011
        if (this->FindParameterType2(parameter_string1, "(", tmp_string1))
        {
          if (this->MessageXNotKnown(NamesOfSetsOfFields, "Set", tmp_string1, index))
            return true;

          NamesOfSetsOfFields.erase(NamesOfSetsOfFields.begin() + index);
          SetsOfFields.erase(SetsOfFields.begin() + index);

          if (this->print_output)
            (*this->Print.out) << "\n  " << this->Print.cbegin << "Set \"" << tmp_string1 << "\" deleted." << this->Print.cend << "\n" << endl;
        }
        else
        // delete sets
        // updated on 19.04.2011
        if (this->FindParameterType1(parameter_string1, "s"))
        {
          NamesOfSetsOfFields.clear();
          SetsOfFields.clear();

          if (this->print_output)
            (*this->Print.out) << "\n  " << this->Print.cbegin << "All sets deleted." << this->Print.cend << "\n" << endl;
        }
      }
      this->MessageParameterNotKnown(parameter_string1);
      return true;
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // insert/remove fields into/from a set (monomial)
    // updated on 19.04.2011
    if ((case1 = this->FindCommandType2(command, "insert(", parameter_string1, parameter_string2))
     || (case2 = this->FindCommandType2(command, "remove(", parameter_string1, parameter_string2)))
    {
      unsigned index = 0;
      // begin: set
      if (this->FindParameterType2(parameter_string2, "set(", tmp_string1) &&
          ((case1 && this->FindParameterType1(parameter_string2, "into")) || (case2 && this->FindParameterType1(parameter_string2, "from"))))
      {
        vector<string>            &NamesOfSetsOfFields = VEVConfig.NamesOfSetsOfFields;
        vector<vector<unsigned> > &SetsOfFields        = VEVConfig.SetsOfFields;

        if (this->MessageXNotKnown(NamesOfSetsOfFields, "Set", tmp_string1, index))
          return true;

        vector<unsigned> &CurrentSet = SetsOfFields[index];

        vector<SUSYMultiplet> Multiplets(2);
        Multiplets[0]=Scalar;
	    Multiplets[1]=LeftFermi;

        vector<string> FieldLabels;
        ExtractLabels(Multiplets, parameter_string1, FieldLabels);
        vector<unsigned> FieldIndices = GetIndices(FieldLabels);

        // find conditions and apply them
        if (!this->FindConditionsAndFilterFieldIndices(parameter_string2, FieldIndices))
        {
          cout << "Warning in bool CPrompt::ExecuteCommand(...): Could not apply the conditions. Return false." << endl;
          return false;
        }

        //insert or remove fields
        unsigned field_counter = 0;

        t1 = FieldIndices.size();
        for (i = 0; i < t1; ++i)
        {
          index = FieldIndices[i];
          vector<unsigned>::iterator pos = find(CurrentSet.begin(), CurrentSet.end(), index);

          if (case1 == (pos == CurrentSet.end()))
          {
            ++field_counter;
            if (case1)
              CurrentSet.push_back(index);
            else
            if (case2)
              CurrentSet.erase(pos);
          }
        }

        if (this->print_output)
        {
          (*this->Print.out) << "\n  " << this->Print.cbegin << field_counter << " field";
          if (field_counter != 1) (*this->Print.out) << "s";
        }

        if (case1)
        {
          stable_sort(CurrentSet.begin(), CurrentSet.end());
          if (this->print_output)
            (*this->Print.out) << " inserted into";
        }
        else
        if (case2)
        {
          if (this->print_output)
            (*this->Print.out) << " removed from";
        }

        if (this->print_output)
          (*this->Print.out) << " set \"" << tmp_string1 << "\"." << this->Print.cend << "\n" << endl;
      }
      // end: set
      this->MessageParameterNotKnown(parameter_string2);
      return true;
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // print set or monomial
    // updated on 11.05.2011
    if ((case1 = this->FindCommandType1(command, "print set", parameter_string1)))
    {
      string           empty_string = "";
      unsigned         max_length   = 0;
      vector<unsigned> SetOfIndicesToPrint;
      unsigned         index = 0;

      // print a set or all sets
      if (case1)
      {
        vector<string>            &NamesOfSetsOfFields = VEVConfig.NamesOfSetsOfFields;
        vector<vector<unsigned> > &SetsOfFields        = VEVConfig.SetsOfFields;

        const unsigned break_line = 20;

        // collect all sets
        // updated on 20.04.2011
        if ((parameter_string1.find("s", 0) != string::npos) && (parameter_string1.find("(", 0) == string::npos))
        {
          this->FindParameterType1(parameter_string1, "s");

          const bool if_not_empty = this->FindParameterType1(parameter_string1, "if not empty");

          t1 = NamesOfSetsOfFields.size();
          if (if_not_empty)
          {
            for (i = 0; i < t1; ++i)
            {
              if (SetsOfFields[i].size() != 0)
                SetOfIndicesToPrint.push_back(i);
            }
          }
          else
          {
            for (i = 0; i < t1; ++i)
              SetOfIndicesToPrint.push_back(i);
          }
        }
        else
        // find the set with label "tmp_string1"
        // updated on 25.01.2011
        if (this->FindParameterType2(parameter_string1, "(", tmp_string1))
        {
          if (this->MessageXNotKnown(NamesOfSetsOfFields, "Set", tmp_string1, index))
            return true;

          SetOfIndicesToPrint.push_back(index);
        }

        // print the set(s)
        t1 = SetOfIndicesToPrint.size();
        if (t1 == 0)
        {
          (*this->Print.out) << "\n  " << this->Print.cbegin << "No sets available." << this->Print.cend << "\n" << endl;
          return true;
        }

        for (i = 0; i < t1; ++i)
        {
          index = SetOfIndicesToPrint[i];
          if (NamesOfSetsOfFields[index].size() > max_length)
            max_length = NamesOfSetsOfFields[index].size();
        }

        (*this->Print.out) << "\n";

        empty_string = "";
        for (i = 0; i < t1; ++i)
        {
          index = SetOfIndicesToPrint[i];

          empty_string.resize(max_length - NamesOfSetsOfFields[index].size(), ' ');
          (*this->Print.out) << "  " << empty_string << NamesOfSetsOfFields[index] << " = {";
          const vector<unsigned> &CurrentSet = SetsOfFields[index];

          t2 = CurrentSet.size();
          for (j = 0; j < t2; ++j)
          {
            this->Print.PrintLabel(Fields[CurrentSet[j]], VEVConfig.use_Labels);
            if (j+1 < t2)
              (*this->Print.out) << this->Print.separator;

            if (((j+1) % break_line) == 0)
            {
              empty_string.resize(max_length + 3, ' ');
              (*this->Print.out) << "\n  " << empty_string;
            }
          }
          (*this->Print.out) << "}" << this->Print.endofset << "\n" << endl;
        }
      }
      this->MessageParameterNotKnown(parameter_string1);
      return true;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // count the number of fields in a set
    // updated on 19.05.2011
    if (this->FindCommandType2(command, "#fields in set(", parameter_string1, parameter_string2))
    {
      vector<string>            &NamesOfSetsOfFields = VEVConfig.NamesOfSetsOfFields;
      vector<vector<unsigned> > &SetsOfFields        = VEVConfig.SetsOfFields;

      unsigned index = 0;
      if (this->MessageXNotKnown(NamesOfSetsOfFields, "Set", parameter_string1, index))
        return true;

      (*this->Print.out) << "\n  " << this->Print.cbegin << "#fields in set \"" << parameter_string1 << "\" = " << SetsOfFields[index].size() << this->Print.cend << "\n" << endl;

      this->MessageParameterNotKnown(parameter_string2);
      return true;
    }

    // end: commands available in all orbifold directories
    // STEP 4 //////////////////////////////////////////////////////////////////////////////////////////////////////////

    const bool UsingStandardConfig = ((VEVConfig.ConfigLabel == "StandardConfig") && (VEVConfig.ConfigNumber == 1));

    // STEP 5 //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // begin: commands available in sub folders
    switch (this->current_folder[1])
    { //
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // ..
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////
      case 0:
      { //

        if (command.substr(0,3) == "man")
        {

        if (command.length() >= 4)
        {

        string path_doc = " ./doc/orbidir/";
        command.insert(4, path_doc);
        command += ".man";


        string temp_file = "temp_error.txt";
        string full_command = command + " 2>" + temp_file;

        int result = system(full_command.c_str());

        ifstream error_file(temp_file);
        stringstream error_stream;
        error_stream << error_file.rdbuf();
        string error_message = error_stream.str();

        error_file.close();
        remove(temp_file.c_str());

        if (result !=0|| containsErrorKeywords(error_message,error_keywords)){
            cout << "Error: A command name is expected after the instruction man.\n"
                 << "Options are:\n* cd\n";
            return false;
        }

            return true;
        }
        else
          {
            cout << "Error: A command name is expected after the instruction man.\n"
                 << "Options are:\n* cd\n";
            return false;
          }
        }

        // updated on 29.06.2011
        if (this->FindCommandType1(command, "cd ", parameter_string1))
        {
          if (this->FindParameterType1(parameter_string1, ".."))
          {
            this->current_folder[0] = -1;
            return true;
          }

          if (Orbifold.GetCheckStatus() != CheckedAndGood)
          {
            if (Orbifold.OrbifoldGroup.GetSpaceGroup_CheckStatus() == CheckedAndGood)
              this->MessageHelpCreateNewOrbifold(3);
            else
              this->MessageHelpCreateNewOrbifold();

            this->current_folder[1] = 1;
            return true;
          }

          if (this->FindParameterType1(parameter_string1, "model"))
          {
            this->current_folder[1] = 1;
            return true;
          }

          if (this->FindParameterType1(parameter_string1, "gauge group"))
          {
            this->current_folder[1] = 2;
            return true;
          }

          if (this->FindParameterType1(parameter_string1, "spectrum"))
          {
            this->current_folder[1] = 3;
            return true;
          }

          if (this->FindParameterType1(parameter_string1, "vev-config/labels"))
          {
            this->current_folder[1] = 5;
            this->current_folder[2] = 1;
            return true;
          }

          if (this->FindParameterType1(parameter_string1, "vev-config"))
          {
            this->current_folder[1] = 5;
            this->current_folder[2] = 0;
            return true;
          }
          this->MessageParameterNotKnown(parameter_string1);
          return true;
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // updated on 09.09.2011
        if (this->FindCommandType1(command, "dir", parameter_string1) || this->FindCommandType1(command, "help", parameter_string1) || this->FindCommandType1(command, "ll", parameter_string1))
        {
          if (this->FindParameterType1(parameter_string1, "processes"))
          this->PrintCommandsProcesses();
          else
          if (this->FindParameterType1(parameter_string1, "short cuts"))
          {
            (*this->Print.out) << "\n  short cuts:\n";
            (*this->Print.out) << "    m   change directory to /model>\n";
            (*this->Print.out) << "    gg  change directory to /gauge group>\n";
            (*this->Print.out) << "    s   change directory to /spectrum>\n";
            (*this->Print.out) << "    v   change directory to /vev-config>\n";
            (*this->Print.out) << "    l   change directory to /vev-config/labels>\n\n" << flush;
          }
          else
          if (this->FindParameterType1(parameter_string1, "conditions"))
            this->PrintCommandsConditions();
          else
          if (this->FindParameterType1(parameter_string1, "sets"))
            this->PrintCommandsSets();
          else
          {
            (*this->Print.out) << "\n  special commands of this directory:\n";
            (*this->Print.out) << "  change directory:\n";
            (*this->Print.out) << "    cd model                                  change directory to /model>\n";
            (*this->Print.out) << "    cd gauge group                            change directory to /gauge group>\n";
            (*this->Print.out) << "    cd spectrum                               change directory to /spectrum>\n";
            (*this->Print.out) << "    cd vev-config                             change directory to /vev-config>\n";
            (*this->Print.out) << "    cd vev-config/labels                      change directory to /labels>\n\n";
            (*this->Print.out) << "  general commands:\n";
            (*this->Print.out) << "    dir                                       show commands\n";
            (*this->Print.out) << "    help                                      optional: \"conditions\", \"processes\", \"sets\", \"short cuts\"\n";
            (*this->Print.out) << "    cd ..                                     leave this directory\n";
            if (!this->online_mode)
              (*this->Print.out) << "    exit                                      exit program\n";
            (*this->Print.out) << "\n" << flush;;
          }
          this->MessageParameterNotKnown(parameter_string1);
          return true;
        }
        break;
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // model
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////
      case 1:
      {
        if (command.substr(0,3) == "man")
        {

        if (command.length() >= 4)
        {

        string path_doc = " ./doc/model/";
        command.insert(4, path_doc);
        command += ".man";



        string temp_file = "temp_error.txt";
        string full_command = command + " 2>" + temp_file;

        int result = system(full_command.c_str());

        ifstream error_file(temp_file);
        stringstream error_stream;
        error_stream << error_file.rdbuf();
        string error_message = error_stream.str();

        error_file.close();
        remove(temp_file.c_str());


        if (result !=0|| containsErrorKeywords(error_message,error_keywords)){
            cout << "Error: A command name is expected after the instruction man.\n"
                 << "Options are:\n* cd\n* print\n* set\n* use\n";
            return false;
        }

            return true;
        }
        else
          {
            cout << "Error: A command name is expected after the instruction man.\n"
                 << "Options are:\n* cd\n* print\n* set\n* use\n";
            return false;
          }
        }





        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // updated on 24.01.2011
        if (this->FindCommandType1(command, "cd ..", parameter_string1))
        {
          this->current_folder[1] = 0;

          this->MessageParameterNotKnown(parameter_string1);
          return true;
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // print heterotic string type
        // updated on 07.09.2011
        if (this->FindCommandType1(command, "print heterotic string type", parameter_string1))
        {
          if (this->print_output)
          {
            const SelfDualLattice Lattice = Orbifold.OrbifoldGroup.GetLattice();
            if (Lattice == E8xE8)
              (*this->Print.out) << "\n  " << this->Print.cbegin << "Using the SO(16)xSO(16) heterotic string." << this->Print.cend << "\n" << endl;
           // else
           //   if (Lattice == Spin32)
           //     (*this->Print.out) << "\n  " << this->Print.cbegin << "Using the Spin(32)/Z_2 heterotic string." << this->Print.cend << "\n" << endl;
            else
              (*this->Print.out) << "\n  " << this->Print.cbegin << "Heterotic string not defined." << this->Print.cend << "\n" << endl;
          }

          this->MessageParameterNotKnown(parameter_string1);
          return true;
        }

        // print available space groups
        // updated on 16.09.2011
        if (this->FindCommandType1(command, "print available space groups", parameter_string1))
        {

          const unsigned N = Orbifold.OrbifoldGroup.GetOrderZN();
          const unsigned K = Orbifold.OrbifoldGroup.GetOrderZK();
          this->FindSpaceGroupsInDirectory(N, K, "Geometry/");


          const bool PrintFilenames = true;

          // begin: print the space groups
          (*this->Print.out) << "\n  " << this->Print.cbegin << "available ";
          if (K != 1)
            (*this->Print.out) << "Z_" << N << " x Z_" << K;
          else
            (*this->Print.out) << "Z_" << N;
          (*this->Print.out) << " space groups: ";

          const size_t s1 = this->PV.AvailableLatticesFilenames.size();

          // no space group
          if (s1 == 0)
          {
            (*this->Print.out) << "none" << this->Print.cend << "\n" << endl;
            this->MessageParameterNotKnown(parameter_string1);
            return true;
          }

          (*this->Print.out) << this->Print.cend << "\n";

          const string space = "                                             ";
          string emptyspace = "";

          int max_length1 = 15;
          int max_length2 = 18;
          for (unsigned i = 0; PrintFilenames && (i < s1); ++i)
          {
            if (this->PV.AvailableLatticesLabels[i].size() > max_length1)
              max_length1 = this->PV.AvailableLatticesLabels[i].size();

            if (this->PV.AvailableAdditionalLabels[i].size() > max_length2)
              max_length2 = this->PV.AvailableAdditionalLabels[i].size();
          }

          (*this->Print.out) << this->Print.cbegin << "     # | lattice   ";

          emptyspace = space;
          emptyspace.resize(max_length1-10);
          (*this->Print.out) << emptyspace << " | additional label ";
          emptyspace = space;
          emptyspace.resize(max_length2-18);
          (*this->Print.out) << emptyspace << "  | geometry file" << this->Print.cend << "\n";
          (*this->Print.out) << "  " << this->Print.cbegin << "  ----------------------------------------------------------------------------------------------------- " << this->Print.cend << "\n";

          for (unsigned i = 0; i < s1; ++i)
          {
            (*this->Print.out) << "    " << this->Print.cbegin << setw(2) << i+1 << " | " << this->PV.AvailableLatticesLabels[i];
            if (this->PV.AvailableLatticesLabels[i].size() < max_length1)
            {
              emptyspace = space;
              emptyspace.resize(max_length1 - this->PV.AvailableLatticesLabels[i].size());
              (*this->Print.out) << emptyspace;
            }

            (*this->Print.out) << " | " << this->PV.AvailableAdditionalLabels[i];
            if (this->PV.AvailableAdditionalLabels[i].size() < max_length2)
            {
              emptyspace = space;
              emptyspace.resize(max_length2 - this->PV.AvailableAdditionalLabels[i].size());
              (*this->Print.out) << emptyspace;
            }

            if (PrintFilenames)
              (*this->Print.out) << " | \"" << this->PV.AvailableLatticesFilenames[i] << "\"";

            (*this->Print.out) << this->Print.cend << "\n";
          }

          // if only one space group exists and none was set before, set it automatically
          if ((s1 == 1) && (Orbifold.OrbifoldGroup.GetSpaceGroup_CheckStatus() != CheckedAndGood))
          {
            (*this->Print.out) << "\n  " << this->Print.cbegin << "Only one space group, hence chosen automatically." << this->Print.cend << "\n";
            this->ExecuteOrbifoldCommand("use space group(1)");
          }
          else
            (*this->Print.out) << "\n";

          (*this->Print.out) << flush;
          // end: print the space groups

          this->MessageParameterNotKnown(parameter_string1);
          return true;
        }


        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // print ...
        // updated on 24.08.2011
        if (this->FindCommandType1(command, "print", parameter_string1))
        {
          const COrbifoldGroup &OrbifoldGroup = Orbifold.OrbifoldGroup;
          const CSpaceGroup    &SpaceGroup    = OrbifoldGroup.GetSpaceGroup();


          if (this->FindParameterType1(parameter_string1, "orbifold label"))
          {
            (*this->Print.out) << "\n  " << this->Print.cbegin << "Orbifold \"" << OrbifoldGroup.Label << "\"." << this->Print.cend << "\n" << endl;
          }
          else
          // updated on 10.08.2011
          if (this->FindParameterType1(parameter_string1, "point group"))
          {
            (*this->Print.out) << "\n  " << this->Print.cbegin << "Point group is ";
            if (SpaceGroup.IsZMxZNxZK())
              (*this->Print.out) << "Z_" << SpaceGroup.GetN() << " x Z_" << SpaceGroup.GetK();
            else
              (*this->Print.out) << "Z_" << SpaceGroup.GetN();

            if (SpaceGroup.additional_label != "")
              (*this->Print.out) << " - " << SpaceGroup.additional_label;

            (*this->Print.out) << "." << this->Print.cend << "\n" << endl;
          }
          else
          // updated on 10.08.2011
          if (this->FindParameterType1(parameter_string1, "space group"))
          {
            if (SpaceGroup.lattice_label == "")
              (*this->Print.out) << "\n  " << this->Print.cbegin << "Space group not chosen yet. Choose it using \"use space group\". See \"help\"." << this->Print.cend << "\n" << endl;
            else
            {
              (*this->Print.out) << "\n  " << this->Print.cbegin << "Space group based on ";
              if (SpaceGroup.IsZMxZNxZK())
                (*this->Print.out) << "Z_" << SpaceGroup.GetN() << " x Z_" << SpaceGroup.GetK();
              else
                (*this->Print.out) << "Z_" << SpaceGroup.GetN();

              if (SpaceGroup.additional_label != "")
                (*this->Print.out) << " - " << SpaceGroup.additional_label;

              (*this->Print.out) << " point group and root-lattice of " << SpaceGroup.lattice_label << "." << this->Print.cend << "\n";
              (*this->Print.out) << "  " << this->Print.cbegin << "Generators are:" << this->Print.cend << "\n" << endl;

              const vector<CSpaceGroupElement> &SG_Generators_Twist = SpaceGroup.GetSG_Generators_Twist();
              const vector<CSpaceGroupElement> &SG_Generators_Shift = SpaceGroup.GetSG_Generators_Shift();

              t1 = SG_Generators_Twist.size();
              for (i = 0; i < t1; ++i)
              {
                if (i == 0)
                  (*this->Print.out) << "  " << this->Print.set_open;
                else
                  (*this->Print.out) << this->Print.separator << "\n  ";

                this->Print.PrintSGElement(SG_Generators_Twist[i]);
              }
              t1 = SG_Generators_Shift.size();
              for (i = 0; i < t1; ++i)
              {
                (*this->Print.out) << this->Print.separator << "\n  ";
                this->Print.PrintSGElement(SG_Generators_Shift[i]);
              }
              (*this->Print.out) << this->Print.set_close << this->Print.endofset << "\n" << endl;
            }
          }
          else
          // updated on 10.08.2011
          if (this->FindParameterType1(parameter_string1, "twists") || this->FindParameterType1(parameter_string1, "twist"))
          {
            if (SpaceGroup.IsZMxZNxZK())
            {
              (*this->Print.out) << "\n  v" << this->Print.underscore << "1 = ";
              this->Print.PrintRational(SpaceGroup.GetTwist(1), SO8);
              (*this->Print.out) << this->Print.endofset << "\n  v" << this->Print.underscore << "2 = ";
              this->Print.PrintRational(SpaceGroup.GetTwist(2), SO8);
              (*this->Print.out) << this->Print.endofset;
            }
            else
            {
              (*this->Print.out) << "\n  v" << this->Print.underscore << "1 = ";
              this->Print.PrintRational(SpaceGroup.GetTwist(1), SO8);
              (*this->Print.out) << this->Print.endofset;
            }
            (*this->Print.out) << "\n" << endl;
          }
          else
          // updated on 04.07.2011
          if (this->FindParameterType1(parameter_string1, "#SUSY"))
          {
            (*this->Print.out) << "\n  " << this->Print.cbegin << "N = " << Orbifold.OrbifoldGroup.GetNumberOfSupersymmetry() << " SUSY in 4d." << this->Print.cend << "\n" << endl;
          }
          else
         //
            if (this->FindParameterType1(parameter_string1, "shifts") || this->FindParameterType1(parameter_string1, "shift"))
          {
            const SelfDualLattice Lattice = OrbifoldGroup.GetLattice();

            if (SpaceGroup.IsZMxZNxZK())
            {
              (*this->Print.out) << "\n  V" << this->Print.underscore << "1 = ";
              this->Print.PrintRational(OrbifoldGroup.GetShift(1), Lattice);
              (*this->Print.out) << this->Print.endofset << "\n  V" << this->Print.underscore << "2 = ";
              this->Print.PrintRational(OrbifoldGroup.GetShift(2), Lattice);
              (*this->Print.out) << this->Print.endofset;
            }
            else
            {
              (*this->Print.out) << "\n  V" << this->Print.underscore << "1 = ";
              this->Print.PrintRational(OrbifoldGroup.GetShift(1), Lattice);
              (*this->Print.out) << this->Print.endofset;
            }
            (*this->Print.out) << "\n" << endl;
          }
          else
          // updated on 10.08.2011
            if (this->FindParameterType1(parameter_string1, "Wilson lines") || this->FindParameterType1(parameter_string1, "WLs"))
          {
            const SelfDualLattice   Lattice          = OrbifoldGroup.GetLattice();
            const vector<unsigned> &WL_AllowedOrders = SpaceGroup.GetWL_AllowedOrders();

            (*this->Print.out) << "\n";
            this->Print.PrintIdentifiedWilsonLines(SpaceGroup.GetWL_Relations());
            (*this->Print.out) << "\n  " << this->Print.cbegin << "Allowed orders of the Wilson lines: ";
            for (i = 0; i < LatticeDim; ++i)
              (*this->Print.out) << WL_AllowedOrders[i] << " ";
            (*this->Print.out) << this->Print.cend << "\n\n";

            this->Print.PrintRational(OrbifoldGroup.GetWilsonLines());
            (*this->Print.out) << flush;
            if (OrbifoldGroup.UseFreelyActingWL)
            {
              (*this->Print.out) << "  " << this->Print.cbegin << "freely acting Wilson line:" << this->Print.cend << "\n  ";
              this->Print.PrintRational(OrbifoldGroup.FreelyActingWilsonLine, Lattice);
              (*this->Print.out) << "\n";
            }
            (*this->Print.out) << endl;
          }

          this->MessageParameterNotKnown(parameter_string1);
          return true;
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // use space group
        // updated on 16.09.2011
        if (this->FindCommandType2(command, "use space group(", parameter_string1, parameter_string2))
        {
          (*this->Print.out) << "\n";
          if (parameter_string1.find_first_not_of("0123456789") != string::npos)
          {
            if (this->print_output)
              (*this->Print.out) << "  " << this->Print.cbegin << "Space group #" << parameter_string1 << " not known." << this->Print.cend << "\n" << endl;
            return true;
          }
          const size_t a1 = this->PV.AvailableLatticesFilenames.size();
          const unsigned lattice_number = (unsigned)atoi(parameter_string1.c_str());
          if ((lattice_number == 0) || (lattice_number > a1))
          {
            if (this->print_output)
              (*this->Print.out) << "  " << this->Print.cbegin << "Space group #" << parameter_string1 << " not known." << this->Print.cend << "\n" << endl;
            return true;
          }

          if (Orbifold.OrbifoldGroup.GetSpaceGroup_CheckStatus() == CheckedAndGood)
          {
            const CSpaceGroup &OldSpaceGroup = Orbifold.OrbifoldGroup.GetSpaceGroup();

            if ((OldSpaceGroup.lattice_label == PV.AvailableLatticesLabels[lattice_number-1]) && (OldSpaceGroup.additional_label == PV.AvailableAdditionalLabels[lattice_number-1]))
            {
              if (this->print_output)
                (*this->Print.out) << "  " << this->Print.cbegin << "Space group #" << parameter_string1 << " is already in use." << this->Print.cend << "\n" << endl;
              return true;
            }
          }

          CSpaceGroup NewSpaceGroup;
          if (!NewSpaceGroup.LoadSpaceGroup(this->PV.AvailableLatticesFilenames[lattice_number-1]))
          {
            if (this->print_output)
              (*this->Print.out) << "  " << this->Print.cbegin << "New space group could not be loaded. Thus, the model is not changed." << this->Print.cend << "\n" << endl;
            return true;
          }
          NewSpaceGroup.Check();

          if (this->print_output)
            (*this->Print.out) << "  " << this->Print.cbegin << "Now using space group #" << parameter_string1 << "." << this->Print.cend << endl;

          // begin: clear old model
          Orbifold.OrbifoldGroup.AccessSpaceGroup() = NewSpaceGroup;
          Orbifold.Reset(true, true, true, true);

          VEVConfigs.clear();
          VEVConfigIndex = 0;

          if (this->print_output)
            (*this->Print.out) << "  " << this->Print.cbegin << "Orbifold model \"" << Orbifold.OrbifoldGroup.Label << "\" cleared." << this->Print.cend << "\n";
          // end: clear old model

          // begin: print
          if (this->print_output)
          {
            (*this->Print.out) << "  " << this->Print.cbegin << "T^6 from root-lattice of " << NewSpaceGroup.lattice_label << " and ";
            if (NewSpaceGroup.IsZMxZN())
            {
              (*this->Print.out) << this->Print.cend << "\n  " << this->Print.cbegin << "twist vector" << this->Print.cend << " v1 = ";
              this->Print.PrintRational(NewSpaceGroup.GetTwist(1), SO8);
              (*this->Print.out) << this->Print.endofset;
            }
            else
            if (NewSpaceGroup.IsZMxZNxZK())
            {

              (*this->Print.out) << this->Print.cend << "\n  " << this->Print.cbegin << "twist vector" << this->Print.cend << " v1 = ";
              this->Print.PrintRational(NewSpaceGroup.GetTwist(1), SO8);
              (*this->Print.out) << this->Print.endofset;
              (*this->Print.out) << this->Print.cend << "\n  " << this->Print.cbegin << "twist vector" << this->Print.cend << " v2 = ";
              this->Print.PrintRational(NewSpaceGroup.GetTwist(2), SO8);
              (*this->Print.out) << this->Print.endofset;
            }
            else
            {
              (*this->Print.out) << "twist vector" << this->Print.cend << " v = ";
              this->Print.PrintRational(NewSpaceGroup.GetTwist(0), SO8);
              (*this->Print.out) << this->Print.endofset;
            }
            (*this->Print.out) << "\n\n";
            this->Print.PrintIdentifiedWilsonLines(NewSpaceGroup.GetWL_Relations());
            (*this->Print.out) << "\n";
          }

          this->MessageHelpCreateNewOrbifold(2);
          // end: print

          this->MessageParameterNotKnown(parameter_string2);
          return true;
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // set shift
        // updated on 20.09.2011
        if (this->FindCommandType1(command, "set shift", parameter_string1))
        {
          COrbifoldGroup     NewOrbifoldGroup = Orbifold.OrbifoldGroup;
          const CSpaceGroup &SpaceGroup       = NewOrbifoldGroup.GetSpaceGroup();

          const SelfDualLattice Lattice = NewOrbifoldGroup.GetLattice();
          const bool            ZMxZN   = SpaceGroup.IsZMxZN();
          const bool            ZMxZNxZK   = SpaceGroup.IsZMxZNxZK();

          if (SpaceGroup.GetNumberOfSectors() == 0)
          {
            if (this->print_output)
              (*this->Print.out) << "\n  " << this->Print.cbegin << "Please choose a space group first, using \"use space group(i)\". See \"help\"." << this->Print.cend << "\n" << endl;
            return true;
          }

          // begin: find parameters
          const bool b1 = this->FindParameterType1(parameter_string1, "standard embedding");
          const bool b2 = this->FindParameterType1(parameter_string1, "V");

          if ((!b1 && !b2) || (b1 && b2))
          {
            if (this->print_output)
            {
              if (ZMxZN)
                (*this->Print.out) << "\n  " << this->Print.cbegin << "One parameter needed: \"standard embedding\" or \"V(i) = <16D vector>\"." << this->Print.cend << "\n" << endl;
              else
                (*this->Print.out) << "\n  " << this->Print.cbegin << "One parameter needed: \"standard embedding\" or \"V = <16D vector>\"." << this->Print.cend << "\n" << endl;
            }

            this->MessageParameterNotKnown(parameter_string1);
            return true;
          }
          // end: find parameters

          NewOrbifoldGroup.AccessWilsonLines().SetToZero(Lattice);

         if (b1)
          {
            CShiftVector newShift(Lattice);
            newShift[0] = 0;
            newShift[1] = 0;
            newShift[2] = 0;
            newShift[3] = 1;
            newShift[4] = 0;
            newShift[5] = 0;
            newShift[6] = 0;
            newShift[7] = 0;
            newShift[8] = 0;
            newShift[9] = 0;
            newShift[10] = 0;
            newShift[11] = 1;
            newShift[12] = 0;
            newShift[13] = 0;
            newShift[14] = 0;
            newShift[15] = 0;

            NewOrbifoldGroup.AccessShift(0) = newShift;

           const CTwistVector &ZN_Twist = SpaceGroup.GetTwist(1);
              newShift[0] = ZN_Twist[1];
              newShift[1] = ZN_Twist[2];
              newShift[2] = ZN_Twist[3];
              for (i = 3; i < 16; ++i)
                newShift[i] = 0;

              NewOrbifoldGroup.AccessShift(1) = newShift;

             if (ZMxZNxZK)
            {
              const CTwistVector &ZK_Twist = SpaceGroup.GetTwist(2);
              newShift[0] = ZK_Twist[1];
              newShift[1] = ZK_Twist[2];
              newShift[2] = ZK_Twist[3];
              for (i = 3; i < 16; ++i)
                newShift[i] = 0;

              NewOrbifoldGroup.AccessShift(2) = newShift;
            }
          }

          // set shift to non-standard embedding
          if (b2)
          {
            CShiftVector newShift(Lattice);
            newShift[0] = 0;
            newShift[1] = 0;
            newShift[2] = 0;
            newShift[3] = 1;
            newShift[4] = 0;
            newShift[5] = 0;
            newShift[6] = 0;
            newShift[7] = 0;
            newShift[8] = 0;
            newShift[9] = 0;
            newShift[10] = 0;
            newShift[11] = 1;
            newShift[12] = 0;
            newShift[13] = 0;
            newShift[14] = 0;
            newShift[15] = 0;

            NewOrbifoldGroup.AccessShift(0) = newShift;

            int shift_number = -1;

            if (ZMxZNxZK)
            {
              if (this->FindParameterType2(parameter_string1, "(", tmp_string1))
              {
                if (tmp_string1 == "1")
                  shift_number = 1;
                else
                  if (tmp_string1 == "2")
                    shift_number = 2;
              }
            }
            else
              shift_number = 1;

            rationalVector RationalVector;

            string::size_type loc1 = 0;
            bool command_ok = (shift_number != -1);
            if (command_ok)
            {
              loc1 = parameter_string1.find("=", 0);
              if (loc1 == string::npos)
                command_ok = false;
              else
                command_ok = convert_string_to_vector_of_rational(parameter_string1.substr(loc1+1, string::npos), RationalVector);
            }

            if (RationalVector.size() != 16)
              command_ok = false;

            if (!command_ok)
            {
              if (this->print_output)
                (*this->Print.out) << "\n  " << this->Print.cbegin << "Shift \"" << parameter_string1.substr(loc1+1, string::npos) << "\" failed. Nothing changed." << this->Print.cend << "\n" << endl;
              return true;
            }

            NewOrbifoldGroup.AccessShift(shift_number) = RationalVector;

            if (ZMxZNxZK && (shift_number == 1))
            {
              CShiftVector ZeroShift(Lattice);
              NewOrbifoldGroup.AccessShift(2) = ZeroShift;

              Orbifold.OrbifoldGroup = NewOrbifoldGroup;
              Orbifold.Reset(false, false, false, false);

              if (this->print_output)
                (*this->Print.out) << "\n  " << this->Print.cbegin << "Next, select the second shift V_2 using \"set shift V(2) = <16D vector>\"." << this->Print.cend << "\n" << endl;
              return true;
            }
          }

          // begin: create and check model dependent part of the orbifold group
          NewOrbifoldGroup.CreateModelDependentPart(this->Print, this->print_output);
          if (NewOrbifoldGroup.GetModularInvariance_CheckStatus() != CheckedAndGood)
          {
            if (this->print_output)
              (*this->Print.out) << "\n  " << this->Print.cbegin << "Modular invariance failed. Using old shift." << this->Print.cend << "\n" << endl;

            return true;
          }
          // end: create and check model dependent part of the orbifold group

          if (this->print_output)
            (*this->Print.out) << "\n  " << this->Print.cbegin << "Modular invariance ok." << this->Print.cend << "\n";

          // begin: clear old model
          VEVConfigs.clear();
          VEVConfigIndex = 0;

          if (this->print_output)
            (*this->Print.out) << "  " << this->Print.cbegin << "Orbifold model \"" << NewOrbifoldGroup.Label << "\" cleared. Computing new spectrum ." << flush;
          // end: clear old model

          // begin: create new model
          COrbifold NewOrbifold(NewOrbifoldGroup);
          Orbifold = NewOrbifold;

          if (this->print_output)
            (*this->Print.out) << "." << flush;

          Orbifold.CheckAnomaly(Orbifold.StandardConfig, this->GaugeIndices, this->Print, false);

          if (this->print_output)
            (*this->Print.out) << ". done." << this->Print.cend << "\n" << endl;


          SConfig TestConfig = Orbifold.StandardConfig;
          TestConfig.ConfigLabel = "TestConfig";
          TestConfig.ConfigNumber = 1;

          VEVConfigs.push_back(Orbifold.StandardConfig);
          VEVConfigs.push_back(TestConfig);
          VEVConfigIndex = 1;
          // end: create new model

          return true;
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // set Wilson line
        // updated on 20.09.2011
        if (this->FindCommandType2(command, "set WL W(", parameter_string1, parameter_string2))
        {
          COrbifoldGroup NewOrbifoldGroup = Orbifold.OrbifoldGroup;

          const CSpaceGroup               &SpaceGroup       = NewOrbifoldGroup.GetSpaceGroup();
          const SelfDualLattice            Lattice          = NewOrbifoldGroup.GetLattice();
          const vector<vector<unsigned> > &WL_Relations     = SpaceGroup.GetWL_Relations();
          const vector<unsigned>          &WL_AllowedOrders = SpaceGroup.GetWL_AllowedOrders();

          int WL_number = -1;

          if (parameter_string1 == "1")
            WL_number = 0;
          else
          if (parameter_string1 == "2")
            WL_number = 1;
          else
          if (parameter_string1 == "3")
            WL_number = 2;
          else
          if (parameter_string1 == "4")
            WL_number = 3;
          else
          if (parameter_string1 == "5")
            WL_number = 4;
          else
          if (parameter_string1 == "6")
            WL_number = 5;

          rationalVector RationalVector;

          string::size_type loc1 = 0;
          bool command_ok = (WL_number != -1);
          if (command_ok)
          {
            loc1 = parameter_string2.find("=", 0);
            if (loc1 == string::npos)
              command_ok = false;
            else
              command_ok = convert_string_to_vector_of_rational(parameter_string2.substr(loc1+1, string::npos), RationalVector);
          }

          if (RationalVector.size() != 16)
            command_ok = false;

          if (!command_ok)
          {
            if (this->print_output)
              (*this->Print.out) << "\n  " << this->Print.cbegin << "Wilson line \"" << parameter_string2.substr(loc1+1, string::npos) << "\" failed. Nothing changed." << this->Print.cend << "\n" << endl;
            return true;
          }

          CWilsonLine newWL(Lattice);
          newWL = RationalVector;

          const unsigned WLOrder = WL_AllowedOrders[WL_number];
          // check only if the allowed order of the wilson line is known
          if ((WLOrder != 0) && ((WLOrder % newWL.GetOrder()) != 0))
          {
            if (this->print_output)
              (*this->Print.out) << "\n  " << this->Print.cbegin << "Order of the new Wilson line failed. Using old Wilson lines." << this->Print.cend << "\n" << endl;

            return true;
          }

          NewOrbifoldGroup.AccessWilsonLines().SetWilsonLine(WL_number, newWL);

          if (this->print_output)
            (*this->Print.out) << "\n";
          // begin: change all Wilson lines that are identified on the orbifold
          t1 = WL_Relations.size();
          for (i = 0; i < t1; ++i)
          {
            const vector<unsigned> &WL_Relation = WL_Relations[i];
            if (find(WL_Relation.begin(), WL_Relation.end(), WL_number) != WL_Relation.end())
            {
              if (this->print_output)
                (*this->Print.out) << "  " << this->Print.cbegin << "Identified Wilson lines: W_" << WL_number+1;

              t2 = WL_Relation.size();
              for (j = 0; j < t2; ++j)
              {
                NewOrbifoldGroup.AccessWilsonLines().SetWilsonLine(WL_Relation[j], newWL);
                if (this->print_output && (WL_Relation[j] != WL_number))
                  (*this->Print.out) << " = W_" << WL_Relation[j] + 1;
              }
              if (this->print_output)
                (*this->Print.out) << this->Print.cend << "\n";
            }
          }
          // end: change all Wilson lines that are identified on the orbifold

          NewOrbifoldGroup.AccessWilsonLines().Check(WL_Relations, WL_AllowedOrders);

          // begin: create and check model dependent part of the orbifold group
          NewOrbifoldGroup.CreateModelDependentPart(this->Print, this->print_output);
          if (NewOrbifoldGroup.GetModularInvariance_CheckStatus() != CheckedAndGood)
          {
            if (this->print_output)
              (*this->Print.out) << "  " << this->Print.cbegin << "Modular invariance failed. Using old Wilson lines." << this->Print.cend << "\n" << endl;

            return true;
          }
          // end: create and check model dependent part of the orbifold group

          if (this->print_output)
            (*this->Print.out) << "  " << this->Print.cbegin << "Modular invariance ok." << this->Print.cend << "\n";

          // begin: clear old model
          VEVConfigs.clear();
          VEVConfigIndex = 0;

          if (this->print_output)
            (*this->Print.out) << "  " << this->Print.cbegin << "Orbifold model \"" << NewOrbifoldGroup.Label << "\" cleared. Computing new spectrum ." << flush;
          // end: clear old model

          // begin: create new model
          COrbifold NewOrbifold(NewOrbifoldGroup);
          Orbifold = NewOrbifold;

          if (this->print_output)
            (*this->Print.out) << "." << flush;

          Orbifold.CheckAnomaly(Orbifold.StandardConfig, this->GaugeIndices, this->Print, false);

          if (this->print_output)
            (*this->Print.out) << ". done." << this->Print.cend << "\n" << endl;

          SConfig TestConfig = Orbifold.StandardConfig;
          TestConfig.ConfigLabel = "TestConfig";
          TestConfig.ConfigNumber = 1;

          VEVConfigs.push_back(Orbifold.StandardConfig);
          VEVConfigs.push_back(TestConfig);
          VEVConfigIndex = 1;
          // end: create new model

          return true;
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // updated on 23.09.2011
        if (this->FindCommandType1(command, "dir", parameter_string1) || this->FindCommandType1(command, "help", parameter_string1) || this->FindCommandType1(command, "ll", parameter_string1))
        {
          if (this->FindParameterType1(parameter_string1, "processes"))
          this->PrintCommandsProcesses();
          else
          if (this->FindParameterType1(parameter_string1, "print"))
          {
            (*this->Print.out) << "\n  print commands:\n";
            (*this->Print.out) << "    print orbifold label\n";
            (*this->Print.out) << "    print heterotic string type\n";
            (*this->Print.out) << "    print available space groups\n";
            (*this->Print.out) << "    print point group\n";
            (*this->Print.out) << "    print space group\n";
            (*this->Print.out) << "    print twist\n";
            (*this->Print.out) << "    print #SUSY\n";
            (*this->Print.out) << "    print shift\n";
            (*this->Print.out) << "    print Wilson lines\n\n" << flush;
          }
          else
          if (this->FindParameterType1(parameter_string1, "short cuts"))
          {
            (*this->Print.out) << "\n  short cuts:\n";
            (*this->Print.out) << "    gg  change directory to /gauge group>\n";
            (*this->Print.out) << "    s   change directory to /spectrum>\n";
            (*this->Print.out) << "    v   change directory to /vev-config>\n";
            (*this->Print.out) << "    l   change directory to /vev-config/labels>\n\n" << flush;
          }
          else
          if (this->FindParameterType1(parameter_string1, "conditions"))
            this->PrintCommandsConditions();
          else
          if (this->FindParameterType1(parameter_string1, "sets"))
            this->PrintCommandsSets();
          else
          {
            const CSpaceGroup &SpaceGroup = Orbifold.OrbifoldGroup.GetSpaceGroup();
            const bool         ZMxZNxZK      = SpaceGroup.IsZMxZNxZK();

            (*this->Print.out) << "\n  special commands of this directory:\n";
            (*this->Print.out) << "    print ...                                 various parameters, see \"help print\"\n\n";
            (*this->Print.out) << "    use space group(i)                        ";

            this->FindSpaceGroupsInDirectory(SpaceGroup.GetN(), SpaceGroup.GetK(), "Geometry/");
            if (this->PV.AvailableLatticesFilenames.size() == 1)
            {
              (*this->Print.out) << "with i = 1 for " << this->PV.AvailableLatticesLabels[0];
              if (this->PV.AvailableAdditionalLabels[0] != "")
                (*this->Print.out) << " - " << this->PV.AvailableAdditionalLabels[0];
              (*this->Print.out) << "\n";
            }
            else
              this->PrintFor(this->PV.AvailableLatticesFilenames.size(), "space group", "i");
            (*this->Print.out) << "\n";

              if (ZMxZNxZK)
              (*this->Print.out) << "    set shift V(i) = <16D vector>             for i = 1,2\n";
            else
              (*this->Print.out) << "    set shift V = <16D vector>\n";

            (*this->Print.out) << "    set shift standard embedding\n";
            (*this->Print.out) << "    set WL W(i) = <16D vector>                for i = 1,..," << LatticeDim << "\n\n";


            (*this->Print.out) << "  general commands:\n";
            (*this->Print.out) << "    dir                                       show commands\n";
            (*this->Print.out) << "    help                                      optional: \"conditions\", \"sets\", \"print\", \"short cuts\"\n";
            (*this->Print.out) << "    cd ..                                     leave this directory\n";
            if (!this->online_mode)
              (*this->Print.out) << "    exit                                      exit program\n";
            (*this->Print.out) << "\n" << flush;;
          }
          this->MessageParameterNotKnown(parameter_string1);
          return true;
        }
        break;
      }

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // gauge group
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////
      case 2:
      {

        if (command.substr(0,3) == "man")
        {

        if(command.length()>= 4)
        {

        string path_doc = " ./doc/gauge_group/";
        // Insertar path_doc en la posición 4
        command.insert(4, path_doc);
        // Concatenar ".man" al final
        command += ".man";


        // Redirigir la salida de error estándar a un archivo temporal
        string temp_file = "temp_error.txt";
        string full_command = command + " 2>" + temp_file;

        // Ejecutar el comando utilizando system
        int result = system(full_command.c_str());


        // Leer el archivo temporal para verificar el mensaje de error
        ifstream error_file(temp_file);
        stringstream error_stream;
        error_stream << error_file.rdbuf();
        string error_message = error_stream.str();

        // Eliminar el archivo temporal
        error_file.close();
        remove(temp_file.c_str());

        // Verificar el resultado de la ejercución y el mensaje de error

        if (result !=0|| containsErrorKeywords(error_message,error_keywords)){
            cout << "Error: A command name is expected after the instruction man.\n"
                 << "Options are:\n* cd\n* print\n* set\n";
            return false;
        }


        return true;
      }
      else
        {
            cout << "Error: A command name is expected after the instruction man.\n"
                 << "Options are:\n* cd\n* print\n* set\n";
            return false;
        }
        }




        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // updated on 18.02.2011
        if (this->FindCommandType1(command, "cd ..", parameter_string1))
        {
          this->current_folder[1] = 0;

          this->MessageParameterNotKnown(parameter_string1);
          return true;
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // updated on 29.02.2012
        // begin print
        if (this->FindCommandType1(command, "print", parameter_string1))
        {
	      // commands:
          if (this->FindParameterType1(parameter_string1, "gg") || this->FindParameterType1(parameter_string1, "gauge group"))
          {
            (*this->Print.out) << "\n";
            this->Print.PrintGaugeGroup(VEVConfig, true);
          }
          else
          if (this->FindParameterType1(parameter_string1, "beta coefficients"))
          {
            (*this->Print.out) << "\n  " << this->Print.cbegin << "Beta function coefficients in vev-configuration \"" << VEVConfig.ConfigLabel << VEVConfig.ConfigNumber << "\" (N = " << VEVConfig.InvariantSupercharges.size() << " SUSY):" << this->Print.cend << "\n";

            t1 = VEVConfig.SymmetryGroup.observable_sector_GGs.size();
            for (i = 0; i < t1; ++i)
            {
              t2 = VEVConfig.SymmetryGroup.observable_sector_GGs[i];

              (*this->Print.out) << "    b_{";
              this->Print.PrintGaugeGroupFactor(VEVConfig, t2);
              (*this->Print.out) << "} = " << Analyse.ComputeBetaFunctionCoefficient(this->GaugeIndices, VEVConfig, t2) << endl;
            }
            (*this->Print.out) << endl;
          }
          else
          if (this->FindParameterType1(parameter_string1, "simple roots"))
          {
            (*this->Print.out) << "\n";
            this->Print.PrintSimpleRoots(VEVConfig.SymmetryGroup.GaugeGroup, -1);
            (*this->Print.out) << endl;
          }
          else
          if (this->FindParameterType2(parameter_string1, "simple root(", tmp_string1))
          {
            (*this->Print.out) << "\n";
            if (tmp_string1.find_first_not_of("0123456789") != string::npos)
            {
              (*this->Print.out) << "  " << this->Print.cbegin << "Simple root #" << tmp_string1 << " not known." << this->Print.cend << "\n" << endl;
              return true;
            }
            int sr_number = (unsigned)atoi(tmp_string1.c_str()) - 1;

            this->Print.PrintSimpleRoots(VEVConfig.SymmetryGroup.GaugeGroup, (unsigned)sr_number);
            (*this->Print.out) << endl;
          }
          else
          if (this->FindParameterType1(parameter_string1, "FI term"))
          {
            (*this->Print.out) << "\n  ";
            if (VEVConfig.SymmetryGroup.IsFirstU1Anomalous)
            {
              switch(this->Print.GetOutputType())
              {
                case Tstandard:
                {
                  (*this->Print.out) << "tr Q_anom = " << setw(5) << setprecision(2) << VEVConfig.SymmetryGroup.D0_FI_term << "\n" << endl;
                  break;
                }
                case Tmathematica:
                {
                  (*this->Print.out) << "TrQanom = " << setw(5) << setprecision(2) << VEVConfig.SymmetryGroup.D0_FI_term << ";\n" << endl;
                  break;
                }
                case Tlatex:
                {
                  (*this->Print.out) << "$\\text{tr} Q_\\text{anom} = " << setw(5) << setprecision(2) << VEVConfig.SymmetryGroup.D0_FI_term << "$\n" << endl;
                  break;
                }
              }
            }
            else
              (*this->Print.out) << this->Print.cbegin << "No anomalous U(1)." << this->Print.cend << "\n" << endl;
          }
          else
          if (this->FindParameterType1(parameter_string1, "anomaly info"))
          {
            if (this->Print.GetOutputType() == Tmathematica)
            (*this->Print.out) << "(* bash info: option \"for mathematica\" not possible for this command. *)" << endl;
            Orbifold.CheckAnomaly(VEVConfig, this->GaugeIndices, this->Print, true);
          }
          else
          if (this->FindParameterType1(parameter_string1, "B-L generator"))
          {
            const CVector Null(16);
            if ((VEVConfig.SymmetryGroup.BmLGenerator.GetSize() == 16) && (VEVConfig.SymmetryGroup.BmLGenerator != Null))
            {
              (*this->Print.out) << "\n  " << this->Print.cbegin << "U(1)_B-L generator:" << this->Print.cend << "\n    ";
              this->Print.PrintRational(VEVConfig.SymmetryGroup.BmLGenerator, Orbifold.OrbifoldGroup.GetLattice());
              (*this->Print.out) << "\n" << endl;
            }
            else
            {
              if (this->print_output)
                (*this->Print.out) << "  " << this->Print.cbegin << "U(1)_B-L generator not defined." << this->Print.cend << "\n" << endl;
            }
          }
          else
          if (this->FindParameterType1(parameter_string1, "U1 generators"))
          {
            (*this->Print.out) << "\n";
            this->Print.PrintU1Directions(VEVConfig.SymmetryGroup, -1);
            (*this->Print.out) << endl;
          }
          else
          if (this->FindParameterType2(parameter_string1, "U1 generator(", tmp_string1))
          {
            (*this->Print.out) << "\n";
            if (tmp_string1.find_first_not_of("0123456789") != string::npos)
            {
              (*this->Print.out) << "  " << this->Print.cbegin << "U(1) generator #" << tmp_string1 << " not known." << this->Print.cend << "\n" << endl;
              return true;
            }
            int U1_number = (unsigned)atoi(tmp_string1.c_str()) - 1;

            this->Print.PrintU1Directions(VEVConfig.SymmetryGroup, (unsigned)U1_number);
            (*this->Print.out) << endl;
          }
          this->MessageParameterNotKnown(parameter_string1);
          return true;
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // set generator of some U(1) of or U(1)_B-L and compute new charges
        if ((case1 = this->FindCommandType2(command, "set U1(", parameter_string1, parameter_string2)) || (case2 = this->FindCommandType1(command, "set B-L", parameter_string2)))
        {
          if (UsingStandardConfig)
          {
            if (this->print_output)
              (*this->Print.out) << "\n  " << this->Print.cbegin << "Vev-configuration \"StandardConfig1\" cannot be changed." << this->Print.cend << "\n" << endl;
            return true;
          }

          bool NeedsToBeNonAnomalous = true;
          if (case2)
          {
            if (this->FindParameterType1(parameter_string2, "allow for anomalous B-L"))
              NeedsToBeNonAnomalous = false;
          }

          int U1_number = -1;

          if (case1)
          {
            if (parameter_string1.find_first_not_of("0123456789") != string::npos)
            {
              if (this->print_output)
                (*this->Print.out) << "\n  " << this->Print.cbegin << "Index \"" << parameter_string1 << "\" of U(1) generator ill-defined." << this->Print.cend << "\n" << endl;
              return true;
            }

            U1_number = (unsigned)atoi(parameter_string1.c_str()) - 1;
            const size_t u1 = VEVConfig.SymmetryGroup.GaugeGroup.u1directions.size();

            if (U1_number >= u1)
            {
              if (this->print_output)
                (*this->Print.out) << "\n  " << this->Print.cbegin << "U(1) generator #" << U1_number+1 << " does not exist. Number of U(1)s: " << u1 << this->Print.cend << "\n";
              return true;
            }

            if ((U1_number == 0) && VEVConfig.SymmetryGroup.IsFirstU1Anomalous)
            {
              if (this->print_output)
                (*this->Print.out) << "\n  " << this->Print.cbegin << "Cannot change the anomalous U(1) generator." << this->Print.cend << "\n" << endl;
              return true;
            }
          }
          if (case2)
            U1_number = 23;

          rationalVector RationalVector;

          bool command_ok = (U1_number != -1);
          if (command_ok)
          {
            string::size_type loc1 = parameter_string2.find("=", 0);
            if (loc1 == string::npos)
              command_ok = false;
            else
              command_ok = convert_string_to_vector_of_rational(parameter_string2.substr(loc1+1, string::npos), RationalVector);
          }

          if (RationalVector.size() != 16)
            command_ok = false;

          if (!command_ok)
          {
            if (this->print_output)
              (*this->Print.out) << "\n  " << this->Print.cbegin << "U(1) generator failed. Nothing changed." << this->Print.cend << "\n" << endl;
            return true;
          }

          CVector Generator;
          Generator = RationalVector;
          if (case1)
          {
            if (Orbifold.Config_SetU1Direction(VEVConfig, Generator, (unsigned)U1_number))
            {
              if (this->print_output)
                (*this->Print.out) << "\n  " << this->Print.cbegin << U1_number+1 << "-th U(1) generator changed." << this->Print.cend << "\n" << endl;
            }
            else
            {
              if (this->print_output)
                (*this->Print.out) << "\n  " << this->Print.cbegin << "U(1) generator failed. Nothing changed." << this->Print.cend << "\n" << endl;
            }
          }
          if (case2)
          {
            if (Analyse.SetBmLGenerator(Orbifold, VEVConfig, Generator, NeedsToBeNonAnomalous))
            {
              if (this->print_output)
                (*this->Print.out) << "\n  " << this->Print.cbegin << "U(1)_B-L generator changed." << this->Print.cend << "\n" << endl;
            }
            else
            {
              if (this->print_output)
                (*this->Print.out) << "\n  " << this->Print.cbegin << "U(1)_B-L generator failed. Nothing changed." << this->Print.cend << "\n" << endl;
            }
          }
          return true;
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // updated on 29.02.2012
        if (this->FindCommandType1(command, "dir", parameter_string1) || this->FindCommandType1(command, "help", parameter_string1) || this->FindCommandType1(command, "ll", parameter_string1))
        {
          const size_t s1 = VEVConfig.SymmetryGroup.GaugeGroup.u1directions.size();
          size_t s2 = 0;
          for (unsigned i = 0; i < VEVConfig.SymmetryGroup.GaugeGroup.factor.size(); ++i)
            s2 += VEVConfig.SymmetryGroup.GaugeGroup.factor[i].simpleroots.size();

          if (this->FindParameterType1(parameter_string1, "processes"))
            this->PrintCommandsProcesses();
          else
          if (this->FindParameterType1(parameter_string1, "print"))
          {
            (*this->Print.out) << "\n  print commands:\n";
            (*this->Print.out) << "    print gauge group\n";
            (*this->Print.out) << "    print beta coefficients\n";
            if (s2 != 0)
            {
              (*this->Print.out) << "    print simple root(i)                      ";
              this->PrintFor(s2, "simple roots", "i");
              (*this->Print.out) << "    print simple roots\n";
            }
            if (s1 != 0)
            {
              (*this->Print.out) << "    print anomaly info\n";
            }
            if (s1 != 0)
            {
              (*this->Print.out) << "    print B-L generator\n";
              (*this->Print.out) << "    print U1 generator(i)                     ";
              this->PrintFor(s1, "U(1)s", "i");
              (*this->Print.out) << "    print U1 generators\n";
            }
            (*this->Print.out) << "\n" << flush;
          }
          else
          if (this->FindParameterType1(parameter_string1, "short cuts"))
          {
            (*this->Print.out) << "\n  short cuts:\n";
            (*this->Print.out) << "    m   change directory to /model>\n";
            (*this->Print.out) << "    s   change directory to /spectrum>\n";
            (*this->Print.out) << "    v   change directory to /vev-config>\n";
            (*this->Print.out) << "    l   change directory to /vev-config/labels>\n\n" << flush;
          }
          else
          if (this->FindParameterType1(parameter_string1, "conditions"))
            this->PrintCommandsConditions();
          else
           if (this->FindParameterType1(parameter_string1, "sets"))
            this->PrintCommandsSets();
          else
          {
            (*this->Print.out) << "\n  special commands of this directory:\n";
            (*this->Print.out) << "    print ...                                 various parameters, see \"help print\"\n\n";

            if (!UsingStandardConfig)
            {
              if (s1 != 0)
              {
                (*this->Print.out) << "    set U1(i) = <16D vector>                  ";
                this->PrintFor(s1, "U(1)s", "i");
                (*this->Print.out) << "\n";
                (*this->Print.out) << "    set B-L = <16D vector>                    optional: \"allow for anomalous B-L\"\n\n";
              }
            }
            (*this->Print.out) << "  general commands:\n";
            (*this->Print.out) << "    dir                                       show commands\n";
            (*this->Print.out) << "    help                                      optional: \"conditions\", \"print\", \"processes\", \n";
            (*this->Print.out) << "                                                        \"sets\", \"short cuts\"\n";
            (*this->Print.out) << "    cd ..                                     leave this directory\n";
            if (!this->online_mode)
              (*this->Print.out) << "    exit                                      exit program\n";
            (*this->Print.out) << "\n" << flush;;
          }
          this->MessageParameterNotKnown(parameter_string1);
          return true;
        }
        break;
      }

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // spectrum
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////
      case 3:
      {



        if (command.substr(0,3) == "man")
        {

        if (command.length()>=4)
        {
        string path_doc = " ./doc/spectrum/";
        command.insert(4, path_doc);
        command += ".man";


        string temp_file = "temp_error.txt";
        string full_command = command + " 2>" + temp_file;

        int result = system(full_command.c_str());


        ifstream error_file(temp_file);
        stringstream error_stream;
        error_stream << error_file.rdbuf();
        string error_message = error_stream.str();

        error_file.close();
        remove(temp_file.c_str());


        if (result !=0|| containsErrorKeywords(error_message,error_keywords)){
            cout << "Error: A command name is expected after the instruction man.\n"
                 << "Options are:\n* cd\n* print\n* textable\n* if\n* sets\n";
            return false;
        }


        return true;
        }

        else
        {
            cout << "Error: A command name is expected after the instruction man.\n"
                 << "Options are:\n* cd\n* print\n* textable\n* if\n* sets\n";
            return false;
        }

        }




        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // updated on 18.02.2011
        if (this->FindCommandType1(command, "cd ..", parameter_string1))
        {
          this->current_folder[1] = 0;

          this->MessageParameterNotKnown(parameter_string1);
          return true;
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // print details for some fields
        // updated on 27.09.2011
        if (this->FindCommandType2(command, "print(", parameter_string1, parameter_string2))
        {

          vector<SUSYMultiplet> Multiplets(2); //Particle types to be printed
          Multiplets[0]=Scalar;				//and to be given for equivalence check
	      Multiplets[1]=LeftFermi;

          vector<string> FieldLabels;
          ExtractLabels(Multiplets, parameter_string1, FieldLabels);
          vector<unsigned> FieldIndices = GetIndices(FieldLabels);

          // find conditions and apply them
          if (!this->FindConditionsAndFilterFieldIndices(parameter_string2, FieldIndices))
          {
            (*this->Print.out) << "\n  " << this->Print.cbegin << "Conditions failed." << this->Print.cend << "\n" << endl;
            this->MessageParameterNotKnown(parameter_string2);
            return true;
          }
          if (FieldIndices.size() == 0)
          {
            (*this->Print.out) << "\n  " << this->Print.cbegin << "No state to print." << this->Print.cend << "\n" << endl;
            this->MessageParameterNotKnown(parameter_string2);
            return true;
          }

          const bool PrintInternalInformation = this->FindParameterType1(parameter_string2, "with internal information");

          this->Print.PrintStates(Orbifold, VEVConfig, FieldIndices, PrintInternalInformation);

          this->MessageParameterNotKnown(parameter_string2);
          return true;
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // print summary ...
        // updated on 12.09.2011
        if (this->FindCommandType1(command, "print summary", parameter_string1))
        {
          const bool print_NoU1Charges = this->FindParameterType1(parameter_string1, "no U1s");
          const bool print_labels      = this->FindParameterType1(parameter_string1, "with labels");

          vector<SUSYMultiplet> Multiplets(2);			//Particle types to be printed
	      Multiplets[0]=Scalar;							//and to be given for equivalence check
	      Multiplets[1]=LeftFermi;

          const vector<unsigned> Orig_observable_sector_U1s = VEVConfig.SymmetryGroup.observable_sector_U1s;
          if (print_NoU1Charges)
            VEVConfig.SymmetryGroup.observable_sector_U1s.clear();

         if (this->FindParameterType1(parameter_string1, " of sectors"))
          {
            (*this->Print.out) << "\n";
            this->Print.PrintSummaryOfSectors(Orbifold, VEVConfig, Multiplets, print_labels);
          }
          else
          if (this->FindParameterType1(parameter_string1, " of fixed points"))
          {
            (*this->Print.out) << "\n";
            this->Print.PrintSummaryOfFixedBranes(Orbifold, VEVConfig, Multiplets, print_labels);
          }
          else
          // print summary of twisted sector T(k,l) (use T(0,0) for untwisted sector)
          if (this->FindParameterType2(parameter_string1, "of sector T(", tmp_string1))
          {
            (*this->Print.out) << "\n";

            vector<int> unsigned_Sector;
            convert_string_to_vector_of_int(tmp_string1, unsigned_Sector);

            if (unsigned_Sector.size() == 3)
            {
              t1 = Orbifold.GetNumberOfSectors();
              for (i = 0; i < t1; ++i)
              {
                const CSector &current_Sector = Orbifold.GetSector(i);
                if ((current_Sector.Get_m() == unsigned_Sector[0]) && (current_Sector.Get_n() == unsigned_Sector[1]) && (current_Sector.Get_k() == unsigned_Sector[2]))
                {
                  for (int j=0; j<2; j++)
                  {
                  this->Print.PrintSummary(current_Sector, VEVConfig, Multiplets[j], print_labels);
                  (*this->Print.out) << endl;
			      }
                  return true;
                }
              }
            }
            (*this->Print.out) << "  " << this->Print.cbegin << "Sector \"" << tmp_string1 << "\" not known." << this->Print.cend << "\n" << endl;
          }
          else
          // print summary of fixed point
          if (this->FindParameterType2(parameter_string1, "of fixed point(", tmp_string1))
          {
            bool FixedBraneFound = false;
            const CFixedBrane &FixedBrane = this->AccessFixedBrane(tmp_string1, Orbifold, FixedBraneFound);

            (*this->Print.out) << "\n";
            if (FixedBraneFound)
              for (int j=0; j<2; j++)
              {
              this->Print.PrintSummary(FixedBrane, Orbifold.OrbifoldGroup, VEVConfig, Multiplets[j], print_labels);
		      }
            else
              (*this->Print.out) << "  " << this->Print.cbegin << "Fixed point \"" << tmp_string1 << "\" not known." << this->Print.cend << "\n";
            (*this->Print.out) << endl;

          }
          else
          {
            (*this->Print.out) << "\n";
            this->Print.PrintSummaryOfVEVConfig(VEVConfig, Multiplets, print_labels);
            (*this->Print.out) << "\n";
          }
          if (print_NoU1Charges)
            VEVConfig.SymmetryGroup.observable_sector_U1s = Orig_observable_sector_U1s;

          this->MessageParameterNotKnown(parameter_string1);
          return true;
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // print all states   //sept23
        // updated on 11.05.2011
        if (this->FindCommandType1(command, "print all states", parameter_string1))
        {
          this->Print.PrintStates(Orbifold, VEVConfig);

          this->MessageParameterNotKnown(parameter_string1);
          return true;
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // tex details for some fields
        // updated on 21.09.2011
        if (this->FindCommandType2(command, "tex table(", parameter_string1, parameter_string2))
        {

          vector<SUSYMultiplet> Multiplets(2);
          Multiplets[0]=Scalar;
	      Multiplets[1]=LeftFermi;

          vector<string> FieldLabels;
          ExtractLabels(Multiplets, parameter_string1, FieldLabels);
          vector<unsigned> FieldIndices = GetIndices(FieldLabels);

          // find conditions and apply them
          if (!this->FindConditionsAndFilterFieldIndices(parameter_string2, FieldIndices))
          {
            cout << "Warning in bool CPrompt::ExecuteCommand(...): Could not apply the conditions. Return false." << endl;
            return false;
          }

          // begin: choose the details that will be listed
          this->FindParameterType2(parameter_string2, "print labels(", tmp_string1);

          vector<int> int_PrintLabels;
          convert_string_to_vector_of_int(tmp_string1, int_PrintLabels);

          unsigned number_of_Labels = 0;
          const size_t f1 = VEVConfig.Fields.size();
          if (f1 != 0)
            number_of_Labels = VEVConfig.Fields[0].Labels.size();

          t1 = int_PrintLabels.size();
          vector<unsigned> PrintLabels;
          for (i = 0; i < t1; ++i)
          {
            if ((int_PrintLabels[i] > 0) && (int_PrintLabels[i] <= number_of_Labels))
              PrintLabels.push_back((unsigned)(int_PrintLabels[i]-1));
            else
              (*this->Print.out) << "\n  " << this->Print.cbegin << "Label #" << int_PrintLabels[i] << " not known." << this->Print.cend << endl;
          }
          // end: choose the details that will be listed

          if (FieldIndices.size() == 0)
          {
            (*this->Print.out) << "\n  " << this->Print.cbegin << "No state to print." << this->Print.cend << "\n" << endl;
            return true;
          }

          (*this->Print.out) << "\\documentclass[a4paper,12pt,twoside]{article}\n";
          (*this->Print.out) << "\\usepackage{longtable}\n";
          (*this->Print.out) << "\\usepackage{amsmath}\n\n";
          (*this->Print.out) << "\\newcommand{\\rep}[1]{\\ensuremath\\boldsymbol{#1}}\n";
          (*this->Print.out) << "\\newcommand{\\crep}[1]{\\ensuremath\\overline{\\boldsymbol{#1}}}\n\n";
          (*this->Print.out) << "\\begin{document}\n";
          (*this->Print.out) << "\\addtolength{\\hoffset}{-1.5cm}\n";

          this->Print.TexSpectrum(Orbifold, VEVConfig, FieldIndices, PrintLabels);
          (*this->Print.out) << "\n\\end{document}\n";

          this->MessageParameterNotKnown(parameter_string2);
          return true;
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // print list of charges
        // updated on 21.06.2011
        if (this->FindCommandType2(command, "print list of charges(", parameter_string1, parameter_string2))
        {

          vector<SUSYMultiplet> Multiplets(2);
          Multiplets[0]=Scalar;
	      Multiplets[1]=LeftFermi;

          vector<string> FieldLabels;
          ExtractLabels(Multiplets, parameter_string1, FieldLabels);
          vector<unsigned> FieldIndices = GetIndices(FieldLabels);

          // find conditions and apply them
          if (!this->FindConditionsAndFilterFieldIndices(parameter_string2, FieldIndices))
          {
            cout << "Warning in bool CPrompt::ExecuteCommand(...): Could not apply the conditions. Return false." << endl;
            return false;
          }

          tmp_string1 = "";
          this->FindParameterType2(parameter_string2, "label of list(", tmp_string1);

          if (FieldIndices.size() == 0)
          {
            (*this->Print.out) << "\n  " << this->Print.cbegin << "No state to print." << this->Print.cend << "\n" << endl;

            this->MessageParameterNotKnown(parameter_string2);
            return true;
          }
          this->Print.PrintListOfCharges(Orbifold, VEVConfig, FieldIndices, tmp_string1);
          (*this->Print.out) << endl;

          this->MessageParameterNotKnown(parameter_string2);
          return true;
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // updated on 16.09.2011
        if (this->FindCommandType1(command, "dir", parameter_string1) || this->FindCommandType1(command, "help", parameter_string1) || this->FindCommandType1(command, "ll", parameter_string1))
        {
          if (this->FindParameterType1(parameter_string1, "processes"))
            this->PrintCommandsProcesses();
          else
          if (this->FindParameterType1(parameter_string1, "print summary"))
          {
            (*this->Print.out) << "\n  print summary\n";
            (*this->Print.out) << "    print a summary table of massless fields\n\n";
            (*this->Print.out) << "  parameters:\n";
            (*this->Print.out) << "    \"of sectors\"                              group table by sectors\n";
            (*this->Print.out) << "    \"of fixed points\"                         group table by fixed points\n";
            (*this->Print.out)  << "    \"of sector T(k,m,n)\"                      print sector T(k,m,n) only\n";
            (*this->Print.out)  << "    \"of fixed point(label)\"                   print fixed point \"label\" only\n";
            (*this->Print.out) << "    \"of fixed point(k,m,n,n1,n2,n3,n4,n5,n6)\" print specified fixed point only\n";
            (*this->Print.out) << "    \"no U1s\"                                  omit the U(1) charges\n";
            (*this->Print.out) << "    \"with labels\"                             print the field labels\n\n" << flush;
          }
          else
          if (this->FindParameterType1(parameter_string1, "short cuts"))
          {
            (*this->Print.out) << "\n  short cuts:\n";
            (*this->Print.out) << "    m   change directory to /model>\n";
            (*this->Print.out) << "    gg  change directory to /gauge group>\n";
            (*this->Print.out) << "    v   change directory to /vev-config>\n";
            (*this->Print.out) << "    l   change directory to /vev-config/labels>\n\n" << flush;
          }
          else
          if (this->FindParameterType1(parameter_string1, "conditions"))
            this->PrintCommandsConditions();
          else
          if (this->FindParameterType1(parameter_string1, "sets"))
            this->PrintCommandsSets();
          else
          {
            (*this->Print.out) << "\n  special commands of this directory:\n";
            (*this->Print.out) << "    print(fields)                             optional: \"with internal information\"\n";
            (*this->Print.out) << "    print all states\n";
            (*this->Print.out) << "    print summary                             various parameters, see \"help print summary\"\n";
            (*this->Print.out) << "    print list of charges(fields)             optional: \"label of list(Label)\"\n";
            (*this->Print.out) << "    tex table(fields)                         optional: print labels(i,j,..)\n\n";
            (*this->Print.out) << "  optional for many commands of this directory:\n";
            (*this->Print.out) << "    if(condition)                             only if \"condition\" is fulfiled \n";
            (*this->Print.out) << "  general commands:\n";
            (*this->Print.out) << "    dir (or help)                             show commands\n";
            (*this->Print.out) << "                                              optional: \"conditions\", \"sets\", \"short cuts\" , \"print summary\"\n";
            (*this->Print.out) << "    cd ..                                     leave this directory\n";
            if (!this->online_mode)
              (*this->Print.out) << "    exit                                      exit program\n";
            (*this->Print.out) << "\n" << flush;;
          }
          this->MessageParameterNotKnown(parameter_string1);
          return true;
        }
        break;
      }

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // vev-config
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////
      case 5:
      {
        switch (this->current_folder[2])
        {
          ////////////////////////////////////////////////////////////////////////////////////////////////////////////
          // ..
          ////////////////////////////////////////////////////////////////////////////////////////////////////////////
          case 0:
          {  //l
            bool print_configs = false;
            bool print_help    = false;



            if (command.substr(0,3) == "man")
            {

            if (command.length()>=4)
            {
            string path_doc = " ./doc/vev-config/";
            command.insert(4, path_doc);
            command += ".man";


            string temp_file = "temp_error.txt";
            string full_command = command + " 2>" + temp_file;

            int result = system(full_command.c_str());


            ifstream error_file(temp_file);
            stringstream error_stream;
            error_stream << error_file.rdbuf();
            string error_message = error_stream.str();

            error_file.close();
            remove(temp_file.c_str());


            if (result !=0|| containsErrorKeywords(error_message,error_keywords)){
            cout << "Error: A command name is expected after the instruction man.\n"
                 << "Options are:\n* cd\n* analyze\n* create\n* delete\n* print\n* rename\n* select\n* use\n";
            return false;
            }


            return true;
            }
            else
            {
            cout << "Error: A command name is expected after the instruction man.\n"
                 << "Options are:\n* cd\n* analyze\n* create\n* delete\n* print\n* rename\n* select\n* use\n";
            return false;
            }
            }




            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // updated on 04.09.2011
            if (this->FindCommandType1(command, "dir", parameter_string1) || this->FindCommandType1(command, "help", parameter_string1) || this->FindCommandType1(command, "ll", parameter_string1))
            {
              if (parameter_string1 == "")
                print_configs = true;

              print_help  = true;
            }

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // updated on 18.02.2011
            if (this->FindCommandType1(command, "cd labels", parameter_string1))
            {
              this->current_folder[2] = 1;

              this->MessageParameterNotKnown(parameter_string1);
              return true;
            }

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // updated on 18.02.2011
            if (this->FindCommandType1(command, "cd ..", parameter_string1))
            {
              this->current_folder[1] = 0;

              this->MessageParameterNotKnown(parameter_string1);
              return true;
            }

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // updated on 04.09.2011
            if (print_configs || this->FindCommandType1(command, "print configurations", parameter_string1) || this->FindCommandType1(command, "print configs", parameter_string1))
            {
              int max_length = 0;
              int length     = 0;

              t1 = VEVConfigs.size();
              if (t1 > 0)
              {
                std::ostringstream os;
                os << VEVConfigs[0].ConfigNumber;
                tmp_string1 = os.str();
                max_length = VEVConfigs[0].ConfigLabel.size() + tmp_string1.size();
              }

              for (i = 1; i < t1; ++i)
              {
                std::ostringstream os;
                os << VEVConfigs[i].ConfigNumber;
                tmp_string1 = os.str();
                length = VEVConfigs[i].ConfigLabel.size() + tmp_string1.size();

                if (length > max_length)
                  max_length = length;
              }

              (*this->Print.out) << "\n  " << this->Print.cbegin << "list of vev-configurations: " << this->Print.cend << "\n";
              (*this->Print.out) << "  " << this->Print.cbegin << "   label ";
              for (i = 3; i < max_length; ++i)
                (*this->Print.out) << " ";
              (*this->Print.out) << "| field label # |" << this->Print.cend << "\n";
              (*this->Print.out) << "  " << this->Print.cbegin << "  -------------------------------------- " << this->Print.cend << "\n";

              unsigned counter_printed_VEVs = 0;
              for (i = 0; i < t1; ++i)
              {
                const SConfig &VEVConfig = VEVConfigs[i];
                if (i == VEVConfigIndex)
                  (*this->Print.out) << "  -> ";
                else
                  (*this->Print.out) << "     ";

                (*this->Print.out) << "\"" << VEVConfig.ConfigLabel << VEVConfig.ConfigNumber << "\"";

                std::ostringstream os;
                os << VEVConfig.ConfigNumber;
                tmp_string1 = os.str();
                length = VEVConfig.ConfigLabel.size() + tmp_string1.size();
                for (j = length; j < max_length; ++j)
                  (*this->Print.out) << " ";

                const size_t f1 = VEVConfig.Fields.size();
                if (f1 == 0)
                {
                  cout << "Warning in bool CPrompt::ExecuteCommand(...): \"Fields\" is empty. Return true." << endl;
                  return true;
                }
                const unsigned number_of_labels = VEVConfig.Fields[0].Labels.size();
                (*this->Print.out) << " |      " << setw(3) << VEVConfig.use_Labels + 1 << " /" << setw(3) << number_of_labels << " | ";
                (*this->Print.out) << "\n";
              }
              if (print_help)
                (*this->Print.out) << flush;
              else
              {
                (*this->Print.out) << endl;
                this->MessageParameterNotKnown(parameter_string1);
                return true;
              }
            }

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // updated on 26.08.2011
            if (this->FindCommandType2(command, "use config(", parameter_string1, parameter_string2))
            {
              string   useVEVConfigLabel  = parameter_string1;
              unsigned useVEVConfigNumber = 1;
              this->SplitVEVConfigLabel(useVEVConfigLabel, useVEVConfigNumber);

              unsigned index = 0;
              if (this->MessageVEVConfigNotKnown(VEVConfigs, useVEVConfigLabel, useVEVConfigNumber, index))
              {
                this->MessageParameterNotKnown(parameter_string2);
                return true;
              }
              if (VEVConfigIndex == index)
              {
                if (this->print_output)
                  (*this->Print.out) << "\n  " << this->Print.cbegin << "Vev-configuration \"" << useVEVConfigLabel << useVEVConfigNumber << "\" is already in use." << this->Print.cend << "\n" << endl;
              }
              else
              {
                VEVConfigIndex = index;

                if (this->print_output)
                  (*this->Print.out) << "\n  " << this->Print.cbegin << "Now using vev-configuration \"" << useVEVConfigLabel << useVEVConfigNumber << "\"." << this->Print.cend << "\n" << endl;
              }

              this->MessageParameterNotKnown(parameter_string2);
              return true;
            }

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // updated on 15.07.2011
            if (this->FindCommandType2(command, "create config(", parameter_string1, parameter_string2))
            {
              // begin: find parameters
              string   fromVEVConfigLabel  = "StandardConfig";
              unsigned fromVEVConfigNumber = 1;

              // if no other origin is specified use the standard config to create the new config
              if (this->FindParameterType2(parameter_string2, "from(", fromVEVConfigLabel))
                this->SplitVEVConfigLabel(fromVEVConfigLabel, fromVEVConfigNumber);
              // end: find parameters

              // begin: new config label
              if (this->MessageLabelError(parameter_string1))
              {
                this->MessageParameterNotKnown(parameter_string2);
                return true;
              }

              string   newVEVConfigLabel  = parameter_string1;
              unsigned newVEVConfigNumber = 1;
              this->SplitVEVConfigLabel(newVEVConfigLabel, newVEVConfigNumber);

              while (this->MessageVEVConfigAlreadyExists(VEVConfigs, newVEVConfigLabel, newVEVConfigNumber))
                ++newVEVConfigNumber;
              // end: new config label

              unsigned index = 0;
              if (this->MessageVEVConfigNotKnown(VEVConfigs, fromVEVConfigLabel, fromVEVConfigNumber, index))
              {
                this->MessageParameterNotKnown(parameter_string2);
                return true;
              }

              SConfig NewVEVConfig = VEVConfigs[index];
              NewVEVConfig.ConfigLabel  = newVEVConfigLabel;
              NewVEVConfig.ConfigNumber = newVEVConfigNumber;

              VEVConfigIndex = VEVConfigs.size();
              VEVConfigs.push_back(NewVEVConfig);

              if (this->print_output)
              {
                (*this->Print.out) << "\n  " << this->Print.cbegin << "Vev-configuration \"" << newVEVConfigLabel << newVEVConfigNumber << "\" created from \"" << fromVEVConfigLabel << fromVEVConfigNumber << "\"." << this->Print.cend;
                (*this->Print.out) << "\n  " << this->Print.cbegin << "Now using new configuration \"" << newVEVConfigLabel << newVEVConfigNumber << "\"." << this->Print.cend << "\n" << endl;
              }

              this->MessageParameterNotKnown(parameter_string2);
              return true;
            }

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // updated on 15.07.2011
            if (this->FindCommandType2(command, "rename config(", parameter_string1, parameter_string2))
            {
              // begin: find parameters
              if (!this->FindParameterType2(parameter_string2, "to(", tmp_string1))
              {
                if (this->print_output)
                  (*this->Print.out) << "\n  " << this->Print.cbegin << "New label of vev-configuration not specified." << this->Print.cend << "\n" << endl;

                this->MessageParameterNotKnown(parameter_string2);
                return true;
              }
              // end: find parameters

              // begin: from config
              string   fromVEVConfigLabel  = parameter_string1;
              unsigned fromVEVConfigNumber = 1;
              this->SplitVEVConfigLabel(fromVEVConfigLabel, fromVEVConfigNumber);

              if ((fromVEVConfigLabel == "StandardConfig") && (fromVEVConfigNumber == 1))
              {
                if (this->print_output)
                  (*this->Print.out) << "\n  " << this->Print.cbegin << "Vev-configuration \"StandardConfig1\" cannot be changed." << this->Print.cend << "\n" << endl;

                this->MessageParameterNotKnown(parameter_string2);
                return true;
              }

              unsigned index = 0;
              if (this->MessageVEVConfigNotKnown(VEVConfigs, fromVEVConfigLabel, fromVEVConfigNumber, index))
              {
                this->MessageParameterNotKnown(parameter_string2);
                return true;
              }
              // end: from config

              // begin: to config
              if (this->MessageLabelError(tmp_string1))
              {
                this->MessageParameterNotKnown(tmp_string1);
                return true;
              }

              string   newVEVConfigLabel  = tmp_string1;
              unsigned newVEVConfigNumber = 1;
              this->SplitVEVConfigLabel(newVEVConfigLabel, newVEVConfigNumber);

              while (this->MessageVEVConfigAlreadyExists(VEVConfigs, newVEVConfigLabel, newVEVConfigNumber))
                ++newVEVConfigNumber;
              // end: to config

              SConfig &RenamedVEVConfig = VEVConfigs[index];
              RenamedVEVConfig.ConfigLabel  = newVEVConfigLabel;
              RenamedVEVConfig.ConfigNumber = newVEVConfigNumber;

              if (this->print_output)
                (*this->Print.out) << "\n  " << this->Print.cbegin << "Vev-configuration \""  << fromVEVConfigLabel << fromVEVConfigNumber << "\" renamed to \"" << newVEVConfigLabel << newVEVConfigNumber << "\"." << this->Print.cend << "\n" << endl;

              this->MessageParameterNotKnown(parameter_string2);
              return true;
            }

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // updated on 24.10.2011
            if (this->FindCommandType2(command, "delete config(", parameter_string1, parameter_string2))
            {
              string   deleteVEVConfigLabel  = parameter_string1;
              unsigned deleteVEVConfigNumber = 1;
              this->SplitVEVConfigLabel(deleteVEVConfigLabel, deleteVEVConfigNumber);

              if ((deleteVEVConfigLabel == "StandardConfig") && (deleteVEVConfigNumber == 1))
              {
                if (this->print_output)
                  (*this->Print.out) << "\n  " << this->Print.cbegin << "Vev-configuration \"StandardConfig1\" cannot be changed." << this->Print.cend << "\n" << endl;

                this->MessageParameterNotKnown(parameter_string2);
                return true;
              }

              unsigned index = 0;
              if (this->MessageVEVConfigNotKnown(VEVConfigs, deleteVEVConfigLabel, deleteVEVConfigNumber, index))
              {
                this->MessageParameterNotKnown(parameter_string2);
                return true;
              }

              const vector<bool> &PID_Done = VEVConfigs[index].pid.PID_Done;
              if (find(PID_Done.begin(), PID_Done.end(), false) != PID_Done.end())
              {
                if (this->print_output)
                  (*this->Print.out) << "\n  " << this->Print.cbegin << "Processes of vev-configuration \"" << VEVConfigs[index].ConfigLabel << VEVConfigs[index].ConfigNumber << "\" are still running. Cannot be deleted." << this->Print.cend << "\n" << endl;
              }
              else
              {
                VEVConfigs.erase(VEVConfigs.begin() + index);
                VEVConfigIndex = index-1;

                if (this->print_output)
                {
                  (*this->Print.out) << "\n  " << this->Print.cbegin << "Vev-configuration \"" << deleteVEVConfigLabel << deleteVEVConfigNumber << "\" deleted." << this->Print.cend;
                  (*this->Print.out) << "\n  " << this->Print.cbegin << "Now using vev-configuration \"" << VEVConfigs[VEVConfigIndex].ConfigLabel << VEVConfigs[VEVConfigIndex].ConfigNumber << "\"." << this->Print.cend << "\n" << endl;
                }
              }

              this->MessageParameterNotKnown(parameter_string2);
              return true;
            }


            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // updated on 15.07.2011
            if (this->FindCommandType1(command, "print gauge group", parameter_string1))
            {
              (*this->Print.out) << "\n";
              this->Print.PrintGaugeGroup(VEVConfig, true);

              if ((VEVConfig.SymmetryGroup.observable_sector_GGs.size() != VEVConfig.SymmetryGroup.GaugeGroup.factor.size())
               || (VEVConfig.SymmetryGroup.observable_sector_U1s.size() != VEVConfig.SymmetryGroup.GaugeGroup.u1directions.size()))
                (*this->Print.out) << "  " << this->Print.cbegin << "(factors in brackets, e.g. [SU(2)], belong to the hidden sector in this vev-configuration)" << this->Print.cend << "\n\n";

              (*this->Print.out) << flush;

              this->MessageParameterNotKnown(parameter_string1);
              return true;
            }


            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // updated on 27.01.2012
            if (this->FindCommandType1(command, "select observable sector:", parameter_string1))
            {
              if (UsingStandardConfig)
              {
                if (this->print_output)
                  (*this->Print.out) << "\n  " << this->Print.cbegin << "Vev-configuration \"StandardConfig1\" cannot be changed." << this->Print.cend << "\n" << endl;
                return true;
              }

              bool ok = true;

              // select gauge group factors
              if (this->FindParameterType1(parameter_string1, "full gauge group"))
              {
                t1 = VEVConfig.SymmetryGroup.GaugeGroup.factor.size();
                VEVConfig.SymmetryGroup.observable_sector_GGs.clear();
                for (i = 0; i < t1; ++i)
                  VEVConfig.SymmetryGroup.observable_sector_GGs.push_back(i);

                if (this->print_output)
                {
                  (*this->Print.out) << "\n  " << this->Print.cbegin << "Observable sector changed." << this->Print.cend << "\n";
                  this->Print.PrintGaugeGroup(VEVConfig, true);
                }
              }
              else
              if (this->FindParameterType1(parameter_string1, "no gauge groups"))
              {
                VEVConfig.SymmetryGroup.observable_sector_GGs.clear();

                if (this->print_output)
                {
                  (*this->Print.out) << "\n  " << this->Print.cbegin << "Observable sector changed." << this->Print.cend << "\n";
                  this->Print.PrintGaugeGroup(VEVConfig, true);
                }
              }
              else
              if (this->FindParameterType2(parameter_string1, "gauge group(", parameter_string2))
              {
                t1 = VEVConfig.SymmetryGroup.GaugeGroup.factor.size();
                VEVConfig.SymmetryGroup.observable_sector_GGs.clear();

                vector<int> ggs;
                t2 = 0;

                ok = true;
                if ((parameter_string2 != "") && (parameter_string2.find_first_not_of(",0123456789") == string::npos))
                {
                  convert_string_to_vector_of_int(parameter_string2, ggs);
                  t2 = ggs.size();

                  for (i = 0; ok && (i < t2); ++i)
                  {
                    if ((ggs[i] <= 0) || (ggs[i] > t1))
                      ok = false;
                  }
                }

                if (ok)
                {
                  for (i = 0; i < t2; ++i)
                    VEVConfig.SymmetryGroup.observable_sector_GGs.push_back(ggs[i]-1);

                  if (this->print_output)
                  {
                    (*this->Print.out) << "\n  " << this->Print.cbegin << "Observable sector changed:" << this->Print.cend;
                    this->Print.PrintGaugeGroup(VEVConfig, true);
                  }
                }
                else
                {
                  if (this->print_output)
                    (*this->Print.out) << "\n  " << this->Print.cbegin << "Observable sector not valid." << this->Print.cend << endl;
                }
              }

              // select all U(1) factors
              if (this->FindParameterType1(parameter_string1, "all U1s"))
              {
                t1 = VEVConfig.SymmetryGroup.GaugeGroup.u1directions.size();
                VEVConfig.SymmetryGroup.observable_sector_U1s.clear();
                for (i = 0; i < t1; ++i)
                  VEVConfig.SymmetryGroup.observable_sector_U1s.push_back(i);

                if (this->print_output)
                {
                  (*this->Print.out) << "\n  " << this->Print.cbegin << "Observable sector changed." << this->Print.cend;
                  this->Print.PrintGaugeGroup(VEVConfig, true);
                }
              }
              else
              // select U(1) factors
              if (this->FindParameterType1(parameter_string1, " no U1s"))
              {
                VEVConfig.SymmetryGroup.observable_sector_U1s.clear();

                if (this->print_output)
                {
                  (*this->Print.out) << "\n  " << this->Print.cbegin << "Observable sector changed." << this->Print.cend;
                  this->Print.PrintGaugeGroup(VEVConfig, true);
                }
              }
              else
              if (this->FindParameterType2(parameter_string1, " U1s(", parameter_string2))
              {
                t1 = VEVConfig.SymmetryGroup.GaugeGroup.u1directions.size();
                VEVConfig.SymmetryGroup.observable_sector_U1s.clear();

                vector<int> u1s;
                t2 = 0;

                ok = true;
                if ((parameter_string2 != "") && (parameter_string2.find_first_not_of(",0123456789") == string::npos))
                {
                  convert_string_to_vector_of_int(parameter_string2, u1s);
                  t2 = u1s.size();

                  for (i = 0; ok && (i < t2); ++i)
                  {
                    if ((u1s[i] <= 0) || (u1s[i] > t1))
                      ok = false;
                  }
                }

                if (ok)
                {
                  for (i = 0; i < t2; ++i)
                    VEVConfig.SymmetryGroup.observable_sector_U1s.push_back(u1s[i]-1);

                  if (this->print_output)
                  {
                    (*this->Print.out) << "\n  " << this->Print.cbegin << "Observable sector changed." << this->Print.cend;
                    this->Print.PrintGaugeGroup(VEVConfig, true);
                  }
                }
                else
                {
                  if (this->print_output)
                    (*this->Print.out) << "\n  " << this->Print.cbegin << "Observable sector not valid." << this->Print.cend << endl;
                }
              }

              this->MessageParameterNotKnown(parameter_string1);
              return true;
            }

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // analyze the current config
            // updated on 15.07.2011
            if (this->FindCommandType1(command, "analyze config", parameter_string1) || this->FindCommandType1(command, "analyse config", parameter_string1))
            {
              // begin: find parameters
              const bool SM_PrintSU5SimpleRoots = this->FindParameterType1(parameter_string1, "print SU(5) simple roots");

              unsigned NumberOfGenerations = 3;
              this->FindParameterType4(parameter_string1, "generations", NumberOfGenerations);
              // end: find parameters

              bool SM = true;
              bool PS = true;
              bool SU5 = true;

              if (this->print_output)
                (*this->Print.out) << "\n  " << this->Print.cbegin << "Analyzing vev-configuration \"" << VEVConfig.ConfigLabel << VEVConfig.ConfigNumber << "\" ..." << this->Print.cend << "\n";

              vector<SConfig> NewVEVConfigs;
              if (Analyse.AnalyseModel(Orbifold, VEVConfig, SM, PS, SU5, NewVEVConfigs, this->Print, NumberOfGenerations, SM_PrintSU5SimpleRoots))
              {
                t1 = NewVEVConfigs.size();
                for (i = 0; i < t1; ++i)
                {
                  SConfig &NewVEVConfig = NewVEVConfigs[i];
                  while (this->MessageVEVConfigAlreadyExists(VEVConfigs, NewVEVConfig.ConfigLabel, NewVEVConfig.ConfigNumber))
                    ++NewVEVConfig.ConfigNumber;

                  VEVConfigs.push_back(NewVEVConfig);
                }
                // select last config
                VEVConfigIndex = VEVConfigs.size()-1;

                if (this->print_output)
                {
                  if (t1 == 1)
                    (*this->Print.out) << "  " << this->Print.cbegin << "One vev-configuration identified, labeled \"" << VEVConfigs[VEVConfigIndex].ConfigLabel << VEVConfigs[VEVConfigIndex].ConfigNumber << "\" and selected." << this->Print.cend;
                  else
                  {
                    (*this->Print.out) << "\n  " << this->Print.cbegin << t1 << " vev-configurations identified. Print all configurations:" << this->Print.cend << endl;
                    this->ExecuteOrbifoldCommand("print configs");
                    (*this->Print.out) << "  " << this->Print.cbegin << "Now using vev-configuration \"" << VEVConfigs[VEVConfigIndex].ConfigLabel << VEVConfigs[VEVConfigIndex].ConfigNumber << "\"." << this->Print.cend;
                  }
                  this->current_folder[2] = 1;
                  this->ExecuteOrbifoldCommand("print labels");
                  this->current_folder[2] = 0;
                }
              }
              else
              {
                if (this->print_output)
                  (*this->Print.out) << "  " << this->Print.cbegin << "No vev-configuration identified." << this->Print.cend << "\n" << endl;
              }

              this->MessageParameterNotKnown(parameter_string1);
              return true;
            }

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // updated on 21.09.2011
            if (print_help)
            {
              if (this->FindParameterType1(parameter_string1, "short cuts"))
              {
                (*this->Print.out) << "\n  short cuts:\n";
                (*this->Print.out) << "    m   change directory to /model>\n";
                (*this->Print.out) << "    gg  change directory to /gauge group>\n";
                (*this->Print.out) << "    s   change directory to /spectrum>\n";
                (*this->Print.out) << "    l   change directory to /vev-config/labels>\n\n" << flush;
              }
              else
              if (this->FindParameterType1(parameter_string1, "conditions"))
                this->PrintCommandsConditions();
              else
              if (this->FindParameterType1(parameter_string1, "sets"))
                this->PrintCommandsSets();
              else
              {
                (*this->Print.out) << "\n  special commands of this directory:\n";
                (*this->Print.out) << "    use config(ConfigLabel)                   change to configuration \"ConfigLabel\"\n";
                (*this->Print.out) << "    create config(ConfigLabel)                optional: \"from(AnotherConfigLabel)\"\n";
                (*this->Print.out) << "    rename config(OldConfigLabel) to(NewConfigLabel)\n";
                (*this->Print.out) << "    delete config(ConfigLabel)\n\n";
                (*this->Print.out) << "    print configs\n";
                (*this->Print.out) << "    print gauge group\n";

                if (!UsingStandardConfig)
                {
                  t1 = VEVConfig.SymmetryGroup.GaugeGroup.factor.size();
                  t2 = VEVConfig.SymmetryGroup.GaugeGroup.u1directions.size();

                  if ((t1 != 0) || (t2 != 0))
                  {
                    (*this->Print.out) << "\n";
                    if (t1 != 0)
                    {
                      if (t1 == 1)
                      {
                        (*this->Print.out) << "    select observable sector: ...             parameters: \"gauge group(i)\" ";
                        this->PrintFor(t1, "gauge groups", "i");
                      }
                      else
                      if (t1 > 1)
                      {
                        (*this->Print.out) << "    select observable sector: ...             parameters: \"gauge group(i,j,...)\" ";
                        this->PrintFor(t1, "gauge groups", "i,j");
                      }
                      (*this->Print.out) << "                                                          \"full gauge group\"\n";
                      (*this->Print.out) << "                                                          \"no gauge groups\"\n";
                    }
                    if (t2 != 0)
                    {
                      if (t2 == 1)
                      {
                        (*this->Print.out) << "    select observable sector: ...             parameters: \"U1s(i)\" ";
                        this->PrintFor(t2, "U(1)s", "i");
                      }
                      else
                      if (t2 > 1)
                      {
                        (*this->Print.out) << "    select observable sector: ...             parameters: \"U1s(i,j,...)\" ";
                        this->PrintFor(t2, "U(1)s", "i,j");
                      }
                      (*this->Print.out) << "                                                          \"all U1s\"\n";
                      (*this->Print.out) << "                                                          \"no U1s\"\n";
                    }
                  }
                  (*this->Print.out) << "\n";
                }
                (*this->Print.out) << "    analyze config                            optional: \"print SU(5) simple roots\"\n";
                (*this->Print.out) << "  change directory:\n";
                (*this->Print.out) << "    cd labels                                 change directory to /labels>\n\n";
                (*this->Print.out) << "  general commands:\n";
                (*this->Print.out) << "    dir                                       show commands\n";
                (*this->Print.out) << "    help                                      optional: \"conditions\", \"processes\", \"sets\", \"short cuts\"\n";
                (*this->Print.out) << "    cd ..                                     leave this directory\n";
                if (!this->online_mode)
                  (*this->Print.out) << "    exit                                      exit program\n";
                (*this->Print.out) << "\n" << flush;;
              }

              this->MessageParameterNotKnown(parameter_string1);
              return true;
            }
            break;
          }

          ////////////////////////////////////////////////////////////////////////////////////////////////////////////
          // labels
          ////////////////////////////////////////////////////////////////////////////////////////////////////////////
          case 1:
          {

            if (command.substr(0,3) == "man")
            {

            if (command.length()>=4)
            {
            string path_doc = " ./doc/labels/";
            // Insertar path_doc en la posición 4
            command.insert(4, path_doc);
            // Concatenar ".man" al final
            command += ".man";


            // Redirigir la salida de error estándar a un archivo temporal
            string temp_file = "temp_error.txt";
            string full_command = command + " 2>" + temp_file;

            // Ejecutar el comando utilizando system
            int result = system(full_command.c_str());


            // Leer el archivo temporal para verificar el mensaje de error
            ifstream error_file(temp_file);
            stringstream error_stream;
            error_stream << error_file.rdbuf();
            string error_message = error_stream.str();

            // Eliminar el archivo temporal
            error_file.close();
            remove(temp_file.c_str());

            // Verificar el resultado de la ejercución y el mensaje de error

            if (result !=0|| containsErrorKeywords(error_message,error_keywords)){
            cout << "Error: A command name is expected after the instruction man.\n"
                 << "Options are:\n* cd\n* assign\n* change\n* create\n* load\n* print\n* save\n* use\n";
            return false;
            }


            return true;
            }
            else
            {
            cout << "Error: A command name is expected after the instruction man.\n"
                 << "Options are:\n* cd\n* assign\n* change\n* create\n* load\n* print\n* save\n* use\n";
            return false;
            }

            }


            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // updated on 18.02.2011
            if (this->FindCommandType1(command, "cd ..", parameter_string1))
            {
              this->current_folder[2] = 0;

              this->MessageParameterNotKnown(parameter_string1);
              return true;
            }


            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // updated on 26.09.2011
            if (this->FindCommandType2(command, "use label(", parameter_string1, parameter_string2))
            {
              if (parameter_string1.find_first_not_of("0123456789") != string::npos)
              {
                if (this->print_output)
                  (*this->Print.out) << "\n  " << this->Print.cbegin << "Label #" << parameter_string1 << " not known." << this->Print.cend << "\n" << endl;
                return true;
              }

              if (UsingStandardConfig)
              {
                if (this->print_output)
                  (*this->Print.out) << "\n  " << this->Print.cbegin << "Vev-configuration \"StandardConfig1\" cannot be changed." << this->Print.cend << "\n" << endl;
                return true;
              }

              unsigned tmp = (unsigned)atoi(parameter_string1.c_str());
              if ((tmp == 0) || (tmp > Fields[0].Labels.size()))
              {
                if (this->print_output)
                  (*this->Print.out) << "\n  " << this->Print.cbegin << "Label #" << parameter_string1 << " not known." << this->Print.cend << "\n" << endl;

                return true;
              }

              if (VEVConfig.use_Labels == tmp - 1)
              {
                if (this->print_output)
                  (*this->Print.out) << "\n  " << this->Print.cbegin << "Label #" << parameter_string1 << " is already in use." << this->Print.cend << "\n" << endl;

                return true;
              }

              VEVConfig.use_Labels = tmp - 1;
              if (this->print_output)
                (*this->Print.out) << "\n  " << this->Print.cbegin << "Now using label #" << parameter_string1 << "." << this->Print.cend << "\n" << endl;

              this->MessageParameterNotKnown(parameter_string2);
              return true;
            }


            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // change the label of a single field
            // updated on 26.09.2011
            if (this->FindCommandType2(command, "change label(", parameter_string1, parameter_string2))
            {
              if (UsingStandardConfig)
              {
                if (this->print_output)
                  (*this->Print.out) << "\n  " << this->Print.cbegin << "Vev-configuration \"StandardConfig1\" cannot be changed." << this->Print.cend << "\n" << endl;
                return true;
              }

             vector<SUSYMultiplet> Multiplets(2);
	         Multiplets[0]=Scalar;
	         Multiplets[1]=LeftFermi;


              if (this->print_output)
                (*this->Print.out) << "\n";

              bool stop = false;

              if (!this->FindParameterType2(parameter_string2, "to(", tmp_string1))
                stop = true;
              else
              {
                vector<string> FieldLabels;
                ExtractLabels(Multiplets, parameter_string1, FieldLabels);
                vector<unsigned> FieldIndices = GetIndicesOnlyFieldWithNumber(FieldLabels);

                if (FieldIndices.size() != 1)
                  stop = true;

                if (tmp_string1.size() == 0)
                  stop = true;

                string::size_type loc1 = tmp_string1.find("_");
                if (loc1 == string::npos)
                  stop = true;

                unsigned NewIndex = 0;
                if (!stop)
                {
                  tmp_string2 = tmp_string1.substr(0, loc1);
                  tmp_string3 = tmp_string1.substr(loc1+1, string::npos);

                  if (this->MessageLabelError(tmp_string2) || (tmp_string3.size() == 0) || (tmp_string3.find_first_not_of("0123456789") != string::npos))
                    stop = true;
                  else
                  {
                    NewIndex = (unsigned)atoi(tmp_string3.c_str());

                    if (NewIndex == 0)
                    {
                      stop = true;
                      (*this->Print.out) << "  " << this->Print.cbegin << "Label with index 0 not allowed." << this->Print.cend << endl;
                    }
                  }
                }

                const size_t f1 = Fields.size();
                for (i = 0; !stop && (i < f1); ++i)
                {
                  if (i != FieldIndices[0])
                  {
                    const CField &Field = Fields[i];
                    if ((Field.Labels[VEVConfig.use_Labels] == tmp_string2) && (Field.Numbers[VEVConfig.use_Labels] == NewIndex))
                    {
                      if (this->print_output)
                        (*this->Print.out) << "  " << this->Print.cbegin << "Label used by another field." << this->Print.cend << endl;

                      stop = true;
                    }
                  }
                }
                if (!stop)
                {
                  CField &Field = Fields[FieldIndices[0]];
                  Field.Labels[VEVConfig.use_Labels]  = tmp_string2;
                  Field.Numbers[VEVConfig.use_Labels] = NewIndex;

                  if (this->print_output)
                    (*this->Print.out) << "  " << this->Print.cbegin << "Label changed from \"" << parameter_string1 << "\" to \"" << tmp_string2 << "_" << NewIndex << "\"." << this->Print.cend << "\n" << endl;
                }
              }
              if (stop)
              {
                if (this->print_output)
                  (*this->Print.out) << "  " << this->Print.cbegin << "Label not changed." << this->Print.cend << "\n" << endl;
              }

              this->MessageParameterNotKnown(parameter_string2);
              return true;
            }


            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // create labels for fields, sorting the fields with respect to the i-th U(1) charge
            // updated on 26.09.2011
            if (this->FindCommandType1(command, "create labels", parameter_string1))
            {

              if (UsingStandardConfig)
              {
                if (this->print_output)
                  (*this->Print.out) << "\n  " << this->Print.cbegin << "Vev-configuration \"StandardConfig1\" cannot be changed." << this->Print.cend << "\n" << endl;
                return true;
              }

              vector<SUSYMultiplet> Multiplets(2);
	          Multiplets[0]=Scalar;
	          Multiplets[1]=LeftFermi;


              (*this->Print.out) << "\n";

               if (!Analyse.Labels_Create(std::cin, VEVConfig, this->Print, Multiplets, true))
              {
                if (this->print_output)
                  (*this->Print.out) << "\n  " << this->Print.cbegin << "Creating new labels failed." << this->Print.cend << "\n" << endl;
                return true;
              }
              VEVConfig.use_Labels = Fields[0].Labels.size()-1;

              if (this->print_output)
                (*this->Print.out) << "\n  " << this->Print.cbegin << "Now using new labels (i.e. #" << VEVConfig.use_Labels+1 << ") of the fields." << this->Print.cend << "\n" << endl;

              this->MessageParameterNotKnown(parameter_string1);
              return true;
            }

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // assign a label to a fixed point
            // updated on 26.09.2011
            if (this->FindCommandType2(command, "assign label(", parameter_string1, parameter_string2))
            {
              if (UsingStandardConfig)
              {
                if (this->print_output)
                  (*this->Print.out) << "\n  " << this->Print.cbegin << "Vev-configuration \"StandardConfig1\" cannot be changed." << this->Print.cend << "\n" << endl;
                return true;
              }

              if (this->FindParameterType2(parameter_string2, "to fixed point(", tmp_string1))
              {
                CSpaceGroupElement SGElement;
                if (this->GetLocalization(tmp_string1, SGElement))
                {

                  bool FixedBraneFound = false;
                  CFixedBrane &FixedBrane = Orbifold.AccessFixedBrane(SGElement, FixedBraneFound);

                  if (FixedBraneFound)
                  {

                    FixedBrane.AccessFixedBraneLabel() = parameter_string1;

                    if (this->print_output)
                    {
                      (*this->Print.out) << "\n  " << this->Print.cbegin << "Label \"" << FixedBrane.GetFixedBraneLabel() << "\" assigned to fixed point ";
                      this->Print.PrintSGElement(SGElement);
                      (*this->Print.out) << ")." << this->Print.cend << endl;
                    }
                  }
                }
              }
              this->MessageParameterNotKnown(parameter_string2);
              return true;
            }

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // updated on 12.10.2011
            if (this->FindCommandType1(command, "print labels", parameter_string1))
            {
              (*this->Print.out) << "\n  " << this->Print.cbegin << "Using label #" << VEVConfig.use_Labels+1 << " of the fields." << this->Print.cend << "\n" << endl;

              vector<SUSYMultiplet> Multiplets(2);
     	      Multiplets[0]=Scalar;
	          Multiplets[1]=LeftFermi;

              this->Print.PrintSummaryOfVEVConfig(VEVConfig, Multiplets, true);
              (*this->Print.out) << endl;

              this->MessageParameterNotKnown(parameter_string1);
              return true;
            }

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // updated on 26.09.2011
            if (this->FindCommandType2(command, "load labels(", parameter_string1, parameter_string2))
            {
              if (UsingStandardConfig)
              {
                if (this->print_output)
                  (*this->Print.out) << "\n  " << this->Print.cbegin << "Vev-configuration \"StandardConfig1\" cannot be changed." << this->Print.cend << "\n" << endl;
                return true;
              }

              vector<SUSYMultiplet> Multiplets(2);
	          Multiplets[0]=Scalar;
	          Multiplets[1]=LeftFermi;

              //if (this->online_mode)
              //{
              //  if (this->print_output)
              //    (*this->Print.out) << "\n  " << this->Print.cbegin << "Command disabled in the web interface." << this->Print.cend << "\n" << endl;
              //  return true;
              //}

              if (Analyse.Labels_Load(parameter_string1, VEVConfig))
              {
                if (this->print_output)
                {
                  (*this->Print.out) << "\n  " << this->Print.cbegin << "New labels loaded. Now using label #" << VEVConfig.use_Labels+1 << " of the fields." << this->Print.cend << "\n" << endl;
                  this->Print.PrintSummaryOfVEVConfig(VEVConfig, Multiplets, true);
                  (*this->Print.out) << endl;
                }
              }
              else
              {
                if (this->print_output)
                  (*this->Print.out) << "\n  " << this->Print.cbegin << "The file \"" << parameter_string1 << "\" does not exist or does not contain valid labels." << this->Print.cend << "\n" << endl;
              }

              this->MessageParameterNotKnown(parameter_string2);
              return true;
            }

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // updated on 30.08.2011
            if (this->FindCommandType2(command, "save labels(", parameter_string1, parameter_string2))
            {
              //if (this->online_mode)
              //{
              //  if (this->print_output)
              //    (*this->Print.out) << "\n  " << this->Print.cbegin << "Command disabled in the web interface." << this->Print.cend << "\n" << endl;
              //  return true;
              //}

              std::ofstream out(parameter_string1.data());
              if((!out.is_open()) || (!out.good()))
              {
                (*this->Print.out) << "\n  " << this->Print.cbegin << "File \"" << parameter_string1 << "\" not found." << this->Print.cend << "\n" << endl;

                this->MessageParameterNotKnown(parameter_string2);
                return true;
              }

              const size_t f1 = Fields.size();
              for (i = 0; i < f1; ++i)
              {
                if (i != 0) out << "\n";
                out << i << " " << Fields[i].Labels[VEVConfig.use_Labels] << "_" << Fields[i].Numbers[VEVConfig.use_Labels];
              }
              out.close();

              if (this->print_output)
                (*this->Print.out) << "  " << this->Print.cbegin << f1 << " labels saved to file \"" << parameter_string1 << "\"." << this->Print.cend << "\n" << endl;

              this->MessageParameterNotKnown(parameter_string2);
              return true;
            }

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // updated on 26.09.2011
            if (this->FindCommandType1(command, "dir", parameter_string1) || this->FindCommandType1(command, "help", parameter_string1) || this->FindCommandType1(command, "ll", parameter_string1))
            {
              if (this->FindParameterType1(parameter_string1, "processes"))
                this->PrintCommandsProcesses();
              else
              if (this->FindParameterType1(parameter_string1, "short cuts"))
              {
                (*this->Print.out) << "\n  short cuts:\n";
                (*this->Print.out) << "    m   change directory to /model>\n";
                (*this->Print.out) << "    gg  change directory to /gauge group>\n";
                (*this->Print.out) << "    s   change directory to /spectrum>\n";
                (*this->Print.out) << "    v   change directory to /vev-config>\n\n" << flush;
              }
              else
              if (this->FindParameterType1(parameter_string1, "conditions"))
                this->PrintCommandsConditions();
              else
              if (this->FindParameterType1(parameter_string1, "sets"))
                this->PrintCommandsSets();
              else
              {
                (*this->Print.out) << "\n  special commands of this directory:\n";
                if (!UsingStandardConfig)
                {
                  (*this->Print.out) << "    change label(A_i) to(B_j)\n";
                  (*this->Print.out) << "    create labels\n";
                  (*this->Print.out) << "    assign label(Label) to fixed point(k,m,n,n1,n2,n3,n4,n5,n6)\n";
                }
                (*this->Print.out) << "    print labels                              \n";
                if (!UsingStandardConfig)
                {
                  (*this->Print.out) << "    use label(i)                              ";

                  const size_t f1 = VEVConfig.Fields.size();
                  if (f1 == 0)
                  {
                    cout << "Warning in bool CPrompt::ExecuteCommand(...): \"Fields\" is empty. Return true." << endl;
                    return true;
                  }
                  this->PrintFor(VEVConfig.Fields[0].Labels.size(), "labels", "i");
                  (*this->Print.out) << "\n";
                }
                //if (!this->online_mode)
                //{
                  if (!UsingStandardConfig)
                    (*this->Print.out) << "    load labels(Filename)\n";
                  (*this->Print.out) << "    save labels(Filename)\n";
                //}
                (*this->Print.out) << "\n";
                (*this->Print.out) << "  general commands:\n";
                (*this->Print.out) << "    dir                                       show commands\n";
                (*this->Print.out) << "    help                                      optional: \"conditions\", \"processes\", \"sets\", \"short cuts\"\n";
                (*this->Print.out) << "    cd ..                                     leave this directory\n";
                if (!this->online_mode)
                  (*this->Print.out) << "    exit                                      exit program\n";
                (*this->Print.out) << "\n" << flush;;
              }
              this->MessageParameterNotKnown(parameter_string1);
              return true;
            }

            break;
          }
        }
        break;
      }

    }
    // end: commands available in sub folders
    // STEP 5 //////////////////////////////////////////////////////////////////////////////////////////////////////////
  }
  this->MessageParameterNotKnown(command);
  return true;
}


/* ########################################################################################
######   LoadProgram(const string &Filename, vector<string> &Commands)               ######
######                                                                               ######
######   Version: 03.05.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Filename  : filename of the script file                                  ######
######   output:                                                                     ######
######   2) Commands  : list of commanmds read from "Filename"                       ######
######   return value : finished successfully?                                       ######
###########################################################################################
######   description:                                                                ######
######   Reads commands from file "Filename" and stores them in "Commands".          ######
######################################################################################## */
bool CPrompt::LoadProgram(const string &Filename, vector<string> &Commands)
{
  struct stat stFileInfo;
  int intStat = stat(Filename.c_str(),&stFileInfo);
  if (intStat == 0)
  {
    usleep(200);
    std::ifstream input(Filename.data());

    string command = "";
    unsigned counter = 0;
    while (GetSaveLine(input, command))
    {
      ++counter;
      Commands.push_back(command);
    }
    input.close();
    (*this->Print.out) << "\n  " << this->Print.cbegin << counter << " commands loaded from file \"" << Filename << "\"." << this->Print.cend << "\n" << endl;
  }
  return true;
}



/* ########################################################################################
######   LoadOrbifolds(const string &Filename, bool inequivalent, ...)               ######
######                                                                               ######
######   Version: 19.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Filename                      : filename of the model file               ######
######   2) inequivalent                  : only load inequivalent orbifold models   ######
######   3) compare_couplings_up_to_order : refines the comparision method           ######
######   output:                                                                     ######
######   return value                     : finished successfully?                   ######
###########################################################################################
######   description:                                                                ######
######   Load orbifold models from file "Filename" into the prompt.                  ######
######################################################################################## */
bool CPrompt::LoadOrbifolds(const string &Filename, bool inequivalent, unsigned compare_couplings_up_to_order)
{
 if (this->print_output)
  {
    (*this->Print.out) << "  " << this->Print.cbegin << "Load ";
    if (inequivalent)
      (*this->Print.out) << "inequivalent ";
    (*this->Print.out) << "orbifolds from file \"" << Filename << "\"";
    (*this->Print.out) << "." << this->Print.cend << flush;
  }

  string ProgramFilename = "";
  std::ifstream in(Filename.data());
  if((!in.is_open()) || (!in.good()))
  {
    if (this->print_output)
      (*this->Print.out) << "\n  " << this->Print.cbegin << "File \"" << Filename << "\" not found." << this->Print.cend << "\n" << endl;
    return false;
  }

   vector<SUSYMultiplet> Multiplets(2);
   Multiplets[0]=Scalar;
   Multiplets[1]=LeftFermi;

  CInequivalentModels InequivModels;
  vector<string> FieldLabels;
  unsigned j = 0;

  const size_t o1 = this->Orbifolds.size();
  if (inequivalent)
  {
    for (unsigned i = 0; i < o1; ++i)
    {
      COrbifold &Orbifold = this->Orbifolds[i];

      CSpectrum Spectrum(Orbifold.StandardConfig, Multiplets);

      InequivModels.IsSpectrumUnknown(Spectrum, true);
    }
  }

  unsigned counter = 0;
  const unsigned MAXcounter = 2000;
  bool Orbifold_loaded = true;

  unsigned NewNumber = 1;
  string NewLabelPart1 = "";
  string NewLabel = "";

  COrbifoldGroup NewOrbifoldGroup;


  while ((counter < MAXcounter) && !in.eof())
  {
    if (NewOrbifoldGroup.LoadOrbifoldGroup(in, ProgramFilename))
    {

      if (!this->MessageOrbifoldAlreadyExists(NewOrbifoldGroup.Label))
      {
        COrbifold NewOrbifold(NewOrbifoldGroup);
        NewOrbifold.CheckAnomaly(NewOrbifold.StandardConfig, this->GaugeIndices, this->Print, false);

        Orbifold_loaded = true;

        if (inequivalent)
        {
          CSpectrum Spectrum(NewOrbifold.StandardConfig, Multiplets);

          Orbifold_loaded = InequivModels.IsSpectrumUnknown(Spectrum, true);
        }

        if (Orbifold_loaded)
        {
          this->Orbifolds.push_back(NewOrbifold);
          ++counter;

          SConfig TestConfig = NewOrbifold.StandardConfig;
          TestConfig.ConfigLabel = "TestConfig";
          TestConfig.ConfigNumber = 1;

          vector<SConfig> Configs;
          Configs.push_back(NewOrbifold.StandardConfig);
          Configs.push_back(TestConfig);
          this->AllVEVConfigs.push_back(Configs);
          this->AllVEVConfigsIndex.push_back(1);
        }
      }
    }
  }
  in.close();

  this->current_folder.assign(10,0);
  this->current_folder[0]  = -1;
  this->current_folder[1]  = 0;
  this->current_folder[2]  = 0;

  if (this->print_output)
  {
    if (this->Orbifolds.size() == MAXcounter)
      (*this->Print.out) << "\n  " << this->Print.cbegin << "Maximal limit of " << MAXcounter << " orbifold models reached." << this->Print.cend;
    else
    {
      if (counter == 0)
        (*this->Print.out) << "\n  " << this->Print.cbegin << "No orbifolds loaded." << this->Print.cend;
      else
        if (counter == 1)
          (*this->Print.out) << "\n  " << this->Print.cbegin << "Orbifold \"" << this->Orbifolds[this->Orbifolds.size() - 1].OrbifoldGroup.Label << "\" loaded." << this->Print.cend;
      else
        (*this->Print.out) << "\n  " << this->Print.cbegin << counter << " orbifolds loaded." << this->Print.cend;
    }
    (*this->Print.out) << "\n\n";
  }
  (*this->Print.out) << flush;
  return true;
}


/* ########################################################################################
######   FindCommandType0(const string &inputstring, const string &command) const    ######
######                                                                               ######
######   Version: 06.01.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) inputstring : contains the command given by the user                     ######
######   2) command     : search for this command                                    ######
######   output:                                                                     ######
######   return value   : execute "command"?                                         ######
###########################################################################################
######   description:                                                                ######
######   Compares "command" and the input "inputstring".                             ######
######################################################################################## */
bool CPrompt::FindCommandType0(const string &inputstring, const string &command) const
{
  if (inputstring == command)
    return true;

  return false;
}



/* ########################################################################################
######   FindCommandType1(const string &inputstring, ...) const                      ######
######                                                                               ######
######   Version: 06.01.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) inputstring : contains the command given by the user                     ######
######   2) command     : search for this command                                    ######
######   output:                                                                     ######
######   3) parameters  : the part of "inputstring" not being "command" is stored as ######
######                    "parameters"                                               ######
######   return value   : execute "command"?                                         ######
###########################################################################################
######   description:                                                                ######
######   Compares "command" and the input "inputstring".                             ######
######################################################################################## */
bool CPrompt::FindCommandType1(const string &inputstring, const string &command, string &parameters) const
{
  size_t c1 = command.size();

  if (inputstring.substr(0, c1) == command)
  {
    parameters = inputstring.substr(c1, string::npos);
    return true;
  }
  return false;
}



/* ########################################################################################
######   FindCommandType2(const string &inputstring, ...) const                      ######
######                                                                               ######
######   Version: 06.01.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) inputstring       : contains the command given by the user               ######
######   2) command           : search for this command                              ######
######   output:                                                                     ######
######   3) command_parameter : see example                                          ######
######   4) other_parameters  : see example                                          ######
######   return value         : execute "command"?                                   ######
###########################################################################################
######   description:                                                                ######
######   Example: inputstring = "print(something) with something else"               ######
######   then: command = "print("                                                    ######
######         command_parameter = "something"                                       ######
######         other_parameters = "with something else"                              ######
######################################################################################## */
bool CPrompt::FindCommandType2(const string &inputstring, const string &command, string &command_parameter, string &other_parameters) const
{
  string::size_type loc1 = inputstring.find_first_not_of(" ");
  if (loc1 == string::npos)
    return false;

  string::size_type loc2 = 0;

  if (inputstring.substr(loc1, command.size()) == command)
  {
    loc1 = inputstring.find("(", 0);
    if (loc1 != string::npos)
    {
      loc2 = inputstring.find( ")", loc1);
      if (loc2 != string::npos)
        command_parameter = inputstring.substr(loc1+1, loc2-loc1-1);
      else
        return false;
    }
    else
      return false;

    other_parameters = inputstring.substr(loc2+1, string::npos);

    return true;
  }
  return false;
}



/* ########################################################################################
######   FindParameterType1(string &inputstring, const string &parameter) const      ######
######                                                                               ######
######   Version: 06.01.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) inputstring       : contains the parameter given by the user             ######
######   2) parameter         : search for this parameter                            ######
######   output:                                                                     ######
######   return value         : "parameter" contained in "inputstring"?              ######
###########################################################################################
######   description:                                                                ######
######   If "parameter" is contained in "inputstring" remove it and return true.     ######
######################################################################################## */
bool CPrompt::FindParameterType1(string &inputstring, const string &parameter) const
{
  string::size_type loc1 = inputstring.find(parameter, 0);
  if (loc1 != string::npos)
  {
    string tmp = inputstring;
    inputstring = tmp.substr(0, loc1);
    inputstring += tmp.substr(loc1 + parameter.size(), string::npos);

    return true;
  }
  return false;
}



/* ########################################################################################
######   FindParameterType2(string &inputstring, const string &parameter, ...) const ######
######                                                                               ######
######   Version: 06.01.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) inputstring     : contains the parameter given by the user               ######
######   2) parameter       : search for this parameter                              ######
######   output:                                                                     ######
######   3) parameter_value : see example                                            ######
######   return value       : "parameter" contained in "inputstring"?                ######
###########################################################################################
######   description:                                                                ######
######   If "parameter" is contained in "inputstring" remove it and return true.     ######
######   Example: inputstring = "if(Q_4 == even) and more"                           ######
######   then: parameter = "if("                                                     ######
######         parameter_value = "Q_4 == even"                                       ######
######         inputstring = "and more"                                              ######
######################################################################################## */
bool CPrompt::FindParameterType2(string &inputstring, const string &parameter, string &parameter_value) const
{
  size_t c1 = parameter.size();

  string::size_type loc1 = inputstring.find(parameter, 0);
  if (loc1 != string::npos)
  {
    string::size_type loc2 = inputstring.find(")", loc1+c1);
    if (loc2 != string::npos)
    {
      parameter_value = inputstring.substr(loc1+c1, loc2-loc1-c1);
      string tmp = inputstring;
      inputstring = tmp.substr(0, loc1);
      inputstring += tmp.substr(loc2 + 1, string::npos);
      return true;
    }
  }
  return false;
}



/* ########################################################################################
######   FindParameterType3(string &inputstring, const string &parameter, ...) const ######
######                                                                               ######
######   Version: 06.01.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) inputstring : contains the parameter given by the user                   ######
######   2) parameter   : search for this parameter                                  ######
######   output:                                                                     ######
######   3) index       : see example                                                ######
######   return value   : "parameter" contained in "inputstring"?                    ######
###########################################################################################
######   description:                                                                ######
######   If "parameter" is contained in "inputstring" remove it and return true.     ######
######   Example: inputstring = "Q_4 == even"                                        ######
######   then: parameter = "Q_"                                                      ######
######         index = 4                                                             ######
######         inputstring = "== even"                                               ######
######################################################################################## */
bool CPrompt::FindParameterType3(string &inputstring, const string &parameter, unsigned &index) const
{
  size_t c1 = parameter.size();

  string::size_type loc1 = inputstring.find(parameter, 0);
  if (loc1 != string::npos)
  {
    string::size_type loc2 = inputstring.find(" ", loc1+c1);
    if (loc2 != string::npos)
    {
      index = atoi(inputstring.substr(loc1+c1, loc2-loc1-c1).c_str());
      string tmp = inputstring;
      inputstring = tmp.substr(0, loc1);
      inputstring += tmp.substr(loc2, string::npos);
      return true;
    }
  }
  return false;
}



/* ########################################################################################
######   FindParameterType4(string &inputstring, const string &parameter, ...) const ######
######                                                                               ######
######   Version: 06.01.2012                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) inputstring : contains the parameter given by the user                   ######
######   2) parameter   : search for this parameter                                  ######
######   output:                                                                     ######
######   3) number      : see example                                                ######
######   return value   : "parameter" contained in "inputstring"?                    ######
###########################################################################################
######   description:                                                                ######
######   If "parameter" is contained in "inputstring" remove it and return true.     ######
######   Example: inputstring = "3generations"                                       ######
######   then: parameter = "generations"                                             ######
######         number = 3                                                            ######
######         inputstring = ""                                                      ######
######################################################################################## */
bool CPrompt::FindParameterType4(string &inputstring, const string &parameter, unsigned &number) const
{
  string::size_type loc1 = inputstring.find(parameter, 0);
  if (loc1 != string::npos)
  {
    int loc2 = loc1-1;
    while (loc2 >= 0)
    {
      if ((inputstring.substr(loc2,1)).find_first_of("1234567890") == string::npos)
        break;

      --loc2;
    }
    string tmp = inputstring.substr(loc2+1,loc1-loc2-1);
    if ((tmp.size() == 0) || (tmp.find_first_not_of("1234567890") != string::npos))
      return false;

    number = atoi(tmp.c_str());

    tmp = inputstring;
    inputstring = tmp.substr(0, loc1);
    inputstring += tmp.substr(loc1 + parameter.size(), string::npos);

    return true;
  }
  return false;
}



/* ########################################################################################
######   MessageHelpCreateNewOrbifold(unsigned StartWithLine) const                  ######
######                                                                               ######
######   Version: 19.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) StartWithLine : write the message starting from the "StartWithLine"-th   ######
######                      line                                                     ######
######   output:                                                                     ######
######   -                                                                           ######
###########################################################################################
######   description:                                                                ######
######   Displays a help explaining how to create a new orbifold model.              ######
######################################################################################## */
void CPrompt::MessageHelpCreateNewOrbifold(unsigned StartWithLine) const
{
  if (!this->print_output)
    return;

  (*this->Print.out) << "\n  " << this->Print.cbegin << "Heterotic string with SO(16)xSO(16) lattice is assigned. Input data for orbifold model \"" << this->Orbifolds[this->OrbifoldIndex].OrbifoldGroup.Label << "\" is needed:" << this->Print.cend << "\n";

  vector<string> Text;

  Text.push_back(") print available space groups    : Print a list of possible space groups.");
  Text.push_back(") use space group(i)              : Choose the space group.");

  string SetShift = ") set shift standard embedding      or";
  if (this->Orbifolds[this->OrbifoldIndex].OrbifoldGroup.GetOrderZN() != 1)
    SetShift += this->Print.cend + "\n  " + this->Print.cbegin + "     set shift V(i) = <16D vector>   : Set the gauge shift.";
  else
    SetShift += this->Print.cend + "\n  " + this->Print.cbegin + "     set shift V = <16D vector>      : Set the gauge shift.";
  Text.push_back(SetShift);

  Text.push_back(") set WL W(i) = <16D vector>      : Set the i-th Wilson line.");
  const size_t t1 = Text.size();

  for (unsigned i = StartWithLine; i < t1; ++i)
    (*this->Print.out) << "  " << this->Print.cbegin << "  " << i+1-StartWithLine << Text[i] << this->Print.cend << "\n";

  (*this->Print.out) << endl;
}



/* ########################################################################################
######   MessageLabelError(const string &Label) const                                ######
######                                                                               ######
######   Version: 19.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Label     : label of, for example, a field or orbifold to be checked for ######
######                  correctness                                                  ######
######   output:                                                                     ######
######   return value : is the label correct?                                        ######
###########################################################################################
######   description:                                                                ######
######   Checks whether "Label" is correct and displays a message if not.            ######
######################################################################################## */
bool CPrompt::MessageLabelError(const string &Label) const
{
  if ((Label.size() == 0) || (Label.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789") != string::npos))
  {
    if (this->print_output)
      (*this->Print.out) << "  " << this->Print.cbegin << "Label Error: The label is only allowed to contain characters and numbers." << this->Print.cend << endl;
    return true;
  }
  return false;
}



/* ########################################################################################
######   MessageOrbifoldAlreadyExists(const string &OrbifoldLabel, ...) const        ######
######                                                                               ######
######   Version: 19.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) OrbifoldLabel : label of an orbifold be checked                          ######
######   output:                                                                     ######
######   return value     : does an orbifold with label "OrbifoldLabel" already      ######
######                      exist?                                                   ######
###########################################################################################
######   description:                                                                ######
######   Checks whether "OrbifoldLabel" already exists.                              ######
######################################################################################## */
bool CPrompt::MessageOrbifoldAlreadyExists(const string &OrbifoldLabel, bool PrintOutput) const
{
  size_t s1 = this->Orbifolds.size();
  for (unsigned i = 0; i < s1; ++i)
  {
    if (OrbifoldLabel == this->Orbifolds[i].OrbifoldGroup.Label)
    {
      if (this->print_output && PrintOutput)
        (*this->Print.out) << "\n  " << this->Print.cbegin << "Orbifold labeled \"" << OrbifoldLabel << "\" already exists." << this->Print.cend << endl;
      return true;
    }
  }
  return false;
}



/* ########################################################################################
######   MessageOrbifoldNotKnown(const string &OrbifoldLabel, unsigned &index) const ######
######                                                                               ######
######   Version: 19.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) OrbifoldLabel : label of an orbifold be checked                          ######
######   output:                                                                     ######
######   2) index         : the index of "OrbifoldLabel" in the member variable      ######
######                      "Orbifolds"                                              ######
######   return value     : is the orbifold with label "OrbifoldLabel" unknown?      ######
###########################################################################################
######   description:                                                                ######
######   Finds the index of the orbifold "OrbifoldLabel" in the member variable      ######
######   "Orbifolds". If orbifold not found, display a message.                      ######
######################################################################################## */
bool CPrompt::MessageOrbifoldNotKnown(const string &OrbifoldLabel, unsigned &index) const
{
  size_t s1 = this->Orbifolds.size();
  for (unsigned i = 0; i < s1; ++i)
  {
    if (OrbifoldLabel == this->Orbifolds[i].OrbifoldGroup.Label)
    {
      index = i;
      return false;
    }
  }
  if (this->print_output)
    (*this->Print.out) << "\n  " << this->Print.cbegin << "Orbifold \"" << OrbifoldLabel << "\" not known." << this->Print.cend << "\n" << endl;
  return true;
}



/* ########################################################################################
######   MessageParameterNotKnown(const string &parameter_string) const              ######
######                                                                               ######
######   Version: 19.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) parameter_string : a parameter given by the user                         ######
######   output:                                                                     ######
######   return value        : is the parameter not empty?                           ######
###########################################################################################
######   description:                                                                ######
######   Display a message if "parameter_string" is not empty.                       ######
######################################################################################## */
bool CPrompt::MessageParameterNotKnown(const string &parameter_string) const
{
  if (parameter_string.find_first_not_of(" ") != string::npos)
  {
    if (this->print_output)
      (*this->Print.out) << "  " << this->Print.cbegin << "bash: " << parameter_string << ": parameters not known." << this->Print.cend << "\n" << endl;
    return true;
  }
  return false;
}



/* ########################################################################################
######   MessageVEVConfigAlreadyExists(const vector<SConfig> &VEVConfigs, ...) const ######
######                                                                               ######
######   Version: 19.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) VEVConfigs      : search for "VEVConfigLabel" in here                    ######
######   2) VEVConfigLabel  : 1st part of the label of vev-config to be checked      ######
######   3) VEVConfigNumber : 2nd part of the label of vev-config to be checked      ######
######   output:                                                                     ######
######   return value       : does a vev-config with given label already exist?      ######
###########################################################################################
######   description:                                                                ######
######   Checks whether the given vev-config label already exists.                   ######
######################################################################################## */
bool CPrompt::MessageVEVConfigAlreadyExists(const vector<SConfig> &VEVConfigs, const string &VEVConfigLabel, const unsigned &VEVConfigNumber) const
{
  const size_t s1 = VEVConfigs.size();
  for (unsigned i = 0; i < s1; ++i)
  {
    if ((VEVConfigs[i].ConfigLabel == VEVConfigLabel) && (VEVConfigs[i].ConfigNumber == VEVConfigNumber))
    {
      if (this->print_output)
        (*this->Print.out) << "\n  " << this->Print.cbegin << "Vev-configuration \"" << VEVConfigLabel << VEVConfigNumber << "\" already exists." << this->Print.cend << endl;
      return true;
    }
  }
  return false;
}



/* ########################################################################################
######   MessageVEVConfigNotKnown(const vector<SConfig> &VEVConfigs, ...) const      ######
######                                                                               ######
######   Version: 19.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) VEVConfigs      : search for "VEVConfigLabel" in here                    ######
######   2) VEVConfigLabel  : 1st part of the label of a vev-config be checked       ######
######   3) VEVConfigNumber : 2nd part of the label of a vev-config be checked       ######
######   output:                                                                     ######
######   4) index           : the index of the vev-config label in "VEVConfigs"      ######
######   return value       : is the vev-config label unknown?                       ######
###########################################################################################
######   description:                                                                ######
######   Finds the index of the vev-config "VEVConfigLabel""VEVConfigNumber" in      ######
######   "VEVConfigs". If vev-config not found, display a message.                   ######
######################################################################################## */
bool CPrompt::MessageVEVConfigNotKnown(const vector<SConfig> &VEVConfigs, const string &VEVConfigLabel, const unsigned &VEVConfigNumber, unsigned &index) const
{
  string tmp = "";
  const size_t s1 = VEVConfigs.size();
  for (unsigned i = 0; i < s1; ++i)
  {
    if ((VEVConfigs[i].ConfigLabel == VEVConfigLabel) && (VEVConfigs[i].ConfigNumber == VEVConfigNumber))
    {
      index = i;
      return false;
    }
  }
  if (this->print_output)
    (*this->Print.out) << "\n  " << this->Print.cbegin << "Vev-configuration \"" << VEVConfigLabel << VEVConfigNumber << "\" not known." << this->Print.cend << "\n" << endl;
  return true;
}



/* ########################################################################################
######   MessageVEVConfigAlreadyExists(const vector<SConfig> &VEVConfigs, ...) const ######
######                                                                               ######
######   Version: 19.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) NamesOfSetsOfX  : search for "X_Label" in here                           ######
######   2) X               : the type of label, e.g. "Set" or "Monomial"            ######
######   3) X_Label         : label to be checked                                    ######
######   output:                                                                     ######
######   return value       : does NamesOfSetsOfX contain "X_Label"?                 ######
###########################################################################################
######   description:                                                                ######
######   Checks whether the given vev-config label already exists.                   ######
######################################################################################## */
bool CPrompt::MessageXAlreadyExists(const vector<string> &NamesOfSetsOfX, const string &X, const string &X_Label) const
{
  if (find(NamesOfSetsOfX.begin(), NamesOfSetsOfX.end(), X_Label) != NamesOfSetsOfX.end())
  {
    if (this->print_output)
      (*this->Print.out) << "\n  " << this->Print.cbegin << X << " \"" << X_Label << "\" already exists." << this->Print.cend << endl;
    return true;
  }
  return false;
}



/* ########################################################################################
######   MessageXNotKnown(const vector<string> &NamesOfSetsOfX, ...) const           ######
######                                                                               ######
######   Version: 19.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) NamesOfSetsOfX : search for "X_Label" in this set                        ######
######   2) X              : the type of label, e.g. "Set" or "Monomial"             ######
######   3) X_Label        : label to search for                                     ######
######   output:                                                                     ######
######   4) index          : the index of "X_Label" in "NamesOfSetsOfX"              ######
######   return value      : is the label unknown?                                   ######
###########################################################################################
######   description:                                                                ######
######   Finds the index of "X_Label" in "NamesOfSetsOfX". If not found, display a   ######
######   message.                                                                    ######
######################################################################################## */
bool CPrompt::MessageXNotKnown(const vector<string> &NamesOfSetsOfX, const string &X, const string &X_Label, unsigned &index) const
{
  vector<string>::const_iterator pos = find(NamesOfSetsOfX.begin(), NamesOfSetsOfX.end(), X_Label);
  if (pos == NamesOfSetsOfX.end())
  {
    if (this->print_output)
      (*this->Print.out) << "\n  " << this->Print.cbegin << X << " \"" << X_Label << "\" not known." << this->Print.cend << "\n" << endl;
    return true;
  }
  index = distance(NamesOfSetsOfX.begin(), pos);
  return false;
}



/* ########################################################################################
######   PrintCurrentDirectory(string &output) const                                 ######
######                                                                               ######
######   Version: 10.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   1) output    : write current directory string here                          ######
######   return value : is current directory valid?                                  ######
######################################################################################## */
bool CPrompt::PrintCurrentDirectory(string &output) const
{
  if (this->current_folder[0] == -1)
  {
    output = "> "; // Remove prompt symbol
    return true;
  }
  else
  {
    output = "/" + this->Orbifolds[this->OrbifoldIndex].OrbifoldGroup.Label;
    //output = "";
    switch (this->current_folder[1])
    {
      case 0:
      {
        output += "> "; // Fix double prompt
        return true;
      }
      case 1:
      {
        output += "/model> ";
        return true;
      }
      case 2:
      {
        output += "/gauge group> ";
        return true;
      }
      case 3:
      {
        output += "/spectrum> ";
        return true;
      }
      case 4:
      {
        output += "/couplings> ";
        return true;
      }
      case 5:
      {
        output += "/vev-config";
        switch (this->current_folder[2])
        {
          case 0:
          {
            output += "> ";
            return true;
          }
          case 1:
          {
            output += "/labels> ";
            return true;
          }
        }
      }
    }
  }
  return false;
}







/* ########################################################################################
######   PrintCommandsConditions() const                                             ######
######                                                                               ######
######   Version: 19.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
######################################################################################## */
void CPrompt::PrintCommandsConditions() const
{
  if (!this->print_output)
    return;

  (*this->Print.out) << "\n  if(condition):\n";
  (*this->Print.out) << "    A condition consists of three parts, e.g. if(length == even):\n\n";
  (*this->Print.out) << "      1) left hand side: the variable:\n";
  (*this->Print.out) << "        \"length\"  : the length-square of the left-moving momentum, (p_sh)^2\n";
  (*this->Print.out) << "        \"B-L\"     : the B-L charge\n";
  (*this->Print.out) << "        \"Q_i\"     : i-th U(1) charge, i = 1,2,3,...\n";
  (*this->Print.out) << "        \"p_sh_i\"  : i-th component of the left-moving momentum p_sh, i = 1,..,16\n";
  (*this->Print.out) << "        \"q_sh_i\"  : i-th component of the right-moving momentum q_sh, i = 1,..,4\n";
  (*this->Print.out) << "        \"#osci.\"  : number of oscillators acting on the left-mover\n";
  (*this->Print.out) << "        \"label\"   : will compare the field label\n\n";
  (*this->Print.out) << "      2) in between: the comparison operator\n";
  (*this->Print.out) << "        \"==\", \"!=\", \">\", \">=\", \"<\" and \"<=\" or \n";
  (*this->Print.out) << "        \"involves\", \"!involves\" for field labels\n";
  (*this->Print.out) << "      3) right hand side: the value\n";
  (*this->Print.out) << "        \"even\"    : only with comparison \"==\" and \"!=\"\n";
  (*this->Print.out) << "        \"odd\"     : only with comparison \"==\" and \"!=\"\n";
  (*this->Print.out) << "        \"even/odd\": only with comparison \"==\" and \"!=\"\n";
  (*this->Print.out) << "        string    : check whether the field labels involve \"string\" or not\n";
  (*this->Print.out) << "        rational  : a rational number, e.g. \"0\" and \"1/2\"\n\n";
  (*this->Print.out) << "    More examples:\n";
  (*this->Print.out) << "      1) if(length == 3/2)\n";
  (*this->Print.out) << "      2) if(B-L == even)\n";
  (*this->Print.out) << "      3) if(Q_1 > 1/3)\n";
  (*this->Print.out) << "      4) if(#osci. != 0)\n";
  (*this->Print.out) << "      5) if(label involves X)\n\n" << flush;
}


/* ########################################################################################
######   PrintCommandsProcesses() const                                              ######
######                                                                               ######
######   Version: 19.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
######################################################################################## */
void CPrompt::PrintCommandsProcesses() const
{
  if (!this->print_output)
    return;

  (*this->Print.out) << "\n  processes:\n";
  (*this->Print.out) << "    ps                                        list all active processes\n";
  (*this->Print.out) << "    kill(A)                                   terminate process with ID \"A\"; use \"kill(all)\" to terminate all\n";
  (*this->Print.out) << "    wait(X)                                   wait until all processes have been terminated (check every \"X\" seconds)\n\n" << flush;
}


/* ########################################################################################
######   PrintCommandsSets() const                                                   ######
######                                                                               ######
######   Version: 19.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   -                                                                           ######
######   output:                                                                     ######
######   -                                                                           ######
######################################################################################## */
void CPrompt::PrintCommandsSets() const
{
  if (!this->print_output)
    return;

  (*this->Print.out) << "\n  sets of fields:\n";
  (*this->Print.out) << "    create set(SetLabel)\n";
  (*this->Print.out) << "    delete set(SetLabel)\n";
  (*this->Print.out) << "    delete sets\n";
  (*this->Print.out) << "    insert(fields) into set(SetLabel)         optional: \"if(condition)\"\n";
  (*this->Print.out) << "    remove(fields) from set(SetLabel)         optional: \"if(condition)\"\n";
  (*this->Print.out) << "    print sets                                optional: \"if not empty\"\n";
  (*this->Print.out) << "    print set(SetLabel)\n";
  (*this->Print.out) << "    #fields in set(SetLabel)\n\n" << flush;
}


/* ########################################################################################
######   PrintFor(unsigned number_of_Type, const string &Type, ...) const            ######
######                                                                               ######
######   Version: 10.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) number_of_Type : max number for the index "Var"                          ######
######   2) Type           : what the index runs through, e.g. "space group"         ######
######   3) Var            : the index label, e.g. "i"                               ######
######   output:                                                                     ######
######   -                                                                           ######
######################################################################################## */
void CPrompt::PrintFor(unsigned number_of_Type, const string &Type, const string &Var) const
{
  switch (number_of_Type)
  {
    case 0:
      (*this->Print.out) << "Warning: no " << Type << " available!\n";
      break;
    case 1:
      (*this->Print.out) << "for " << Var << " = 1\n";
      break;
    case 2:
      (*this->Print.out) << "for " << Var << " = 1,2\n";
      break;
    case 3:
      (*this->Print.out) << "for " << Var << " = 1,2,3\n";
      break;
    default:
      (*this->Print.out) << "for " << Var << " = 1,...," << number_of_Type << "\n";
  }
}


/* ########################################################################################
######   FindSpaceGroupsInDirectory(const unsigned &M, const unsigned &N, ...)       ######
######                                                                               ######
######   Version: 29.06.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) M         : order of Z_M                                                 ######
######   2) N         : order of Z_N                                                 ######
######   3) directory : directory of local computer where the geometry files are     ######
######                  stored, default should be "Geometry/"                        ######
######   output:                                                                     ######
######   return value : finished succesfully?                                        ######
######################################################################################## */
bool CPrompt::FindSpaceGroupsInDirectory(const unsigned &N, const unsigned &K, const string &directory)
{
  this->PV.AvailableLatticesFilenames.clear();
  this->PV.AvailableLatticesLabels.clear();
  this->PV.AvailableAdditionalLabels.clear();

  std::ostringstream os1, os2;
  os1 << N;

  string tmp = "dir " + directory + "Geometry_Z" + os1.str();
  if (K != 1)
  {
    os2 << K;
    tmp += "xZ" + os2.str();
  }
  tmp += "*.txt";

  // begin: search in the directory "Geometry" for possible files
  vector<string> PossibleFilenames;

  string::size_type loc1 = 0;
  string::size_type loc2 = 0;
  char line[130];
  FILE *fptr = popen(tmp.c_str(), "r");
  while (fgets( line, sizeof line, fptr))
  {
    tmp = line;

    loc1 = 0;
    loc2 = 0;
    while (loc1 != string::npos)
    {
      loc1 = tmp.find("Geometry_", loc2);
      if (loc1 != string::npos)
      {
        tmp = tmp.substr(loc1, string::npos);
        loc2 = tmp.find(".txt", 0);
        if (loc2 != string::npos)
        {
          PossibleFilenames.push_back(directory + tmp.substr(0, loc2+4));
          loc1 = loc2 + 1;
        }
      }
    }
  }
  pclose(fptr);
  // end: search in the directory "Geometry" for possible files

  string tmp_additional_label = "";
  string tmp_lattice_label    = "";

  bool go_on = true;
  bool PointGroupFound = false;

  const size_t s0 = PossibleFilenames.size();
  for (unsigned i = 0; i < s0; ++i)
  {
    std::ifstream in;
    in.open(PossibleFilenames[i].data(), ifstream::in);
    if(!in.is_open() || !in.good())
    {
      cout << "\n  Warning in bool CPrompt::FindSpaceGroupsInDirectory(...) : Could not find the file \"" << PossibleFilenames[i] << "\". Return false." << endl;
      return false;
    }

    // begin: find if current file fits to the chosen orbifold and get the lattice_label
    tmp_additional_label = "";
    tmp_lattice_label    = "";

    go_on = true;
    PointGroupFound = false;

    while (go_on && (!PointGroupFound || (tmp_additional_label == "") || (tmp_lattice_label == "")) && GetSaveLine(in, tmp))
    {
      if (tmp.substr(0,11) == "point group")
      {
        unsigned tmp_M = 0;
        unsigned tmp_N = 0;
        unsigned tmp_K = 0;

        GetSaveLine(in, tmp);
        std::istringstream line1(tmp);
        line1 >> tmp_M;

        GetSaveLine(in, tmp);
        std::istringstream line2(tmp);
        line2 >> tmp_N;

        GetSaveLine(in, tmp);
        std::istringstream line3(tmp);
        line3 >> tmp_K;

        PointGroupFound = true;
        if ((tmp_N != N) || (tmp_K != K))
          go_on = false;
      }
      else
        if (tmp.substr(0,16) == "additional label")
          GetSaveLine(in, tmp_additional_label);
      else
        if (tmp.substr(0,13) == "lattice label")
          GetSaveLine(in, tmp_lattice_label);
    }
    in.close();
    // end: find if current file fits to the chosen orbifold and get the lattice_label

    if (!PointGroupFound || (tmp_lattice_label == ""))
      go_on = false;

    if (go_on)
    {
      this->PV.AvailableLatticesFilenames.push_back(PossibleFilenames[i]);
      this->PV.AvailableLatticesLabels.push_back(tmp_lattice_label);
      this->PV.AvailableAdditionalLabels.push_back(tmp_additional_label);
    }
  }
    return true;
}


/* ########################################################################################
######   SplitVEVConfigLabel(string &VEVConfigLabel, unsigned &VEVConfigNumber) const######
######                                                                               ######
######   Version: 29.06.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) VEVConfigLabel  : the full vev-config label, e.g. "Test1"                ######
######   output:                                                                     ######
######   1) VEVConfigLabel  : the label, e.g. "Test"                                 ######
######   2) VEVConfigNumber : the number, e.g. 1                                     ######
######   return value : finished succesfully?                                        ######
######################################################################################## */
bool CPrompt::SplitVEVConfigLabel(string &VEVConfigLabel, unsigned &VEVConfigNumber) const
{
  if (VEVConfigLabel.size() == 0)
  {
    if (this->print_output)
      (*this->Print.out) << "\n  " << this->Print.cbegin << "Empty string is not a valid configuration label. Correct labels look like \"ExampleVEVConfig12\"." << this->Print.cend << endl;
    VEVConfigLabel  = "ErrorLabel";
    VEVConfigNumber = 1;
    return false;
  }

  const size_t pos = VEVConfigLabel.find_first_of("0123456789", 0);
  if (pos == string::npos)
  {
    VEVConfigNumber = 1;
    return true;
  }

  if ((pos != string::npos) && (VEVConfigLabel.find_first_not_of("0123456789", pos) != string::npos))
  {
    if (this->print_output)
      (*this->Print.out) << "\n  " << this->Print.cbegin << "\"" << VEVConfigLabel << "\" is not a valid configuration label. Correct labels look like \"ExampleVEVConfig12\"." << this->Print.cend << endl;
    VEVConfigLabel  = "ErrorLabel";
    VEVConfigNumber = 1;
    return false;
  }

  VEVConfigNumber = atoi((VEVConfigLabel.substr(pos)).c_str());
  VEVConfigLabel  = VEVConfigLabel.substr(0, pos);
  return true;
}



/* ########################################################################################
######   GetIndices(const vector<string> &FieldLabels) const                         ######
######                                                                               ######
######   Version: 19.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) FieldLabels : some field labels, e.g. "N Q_1 Q_2"                        ######
######   output:                                                                     ######
######   return value   : the corresponding indices, e.g. the indices of all fields  ######
######                    "N" and the indices of "Q_1" and "Q_2"                     ######
######################################################################################## */
vector<unsigned> CPrompt::GetIndices(const vector<string> &FieldLabels) const
{
  const SConfig &VEVConfig = this->AllVEVConfigs[this->OrbifoldIndex][this->AllVEVConfigsIndex[this->OrbifoldIndex]];

  return global_GetIndices(FieldLabels, VEVConfig.use_Labels, VEVConfig.Fields);
}



/* ########################################################################################
######   GetIndicesOnlyFieldWithNumber(const vector<string> &FieldLabels) const      ######
######                                                                               ######
######   Version: 19.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) FieldLabels : some field labels, e.g. "Q_1 Q_2"                          ######
######   output:                                                                     ######
######   return value   : the corresponding indices, e.g. the indices of "Q_1" and   ######
######                    "Q_2"                                                      ######
######################################################################################## */
vector<unsigned> CPrompt::GetIndicesOnlyFieldWithNumber(const vector<string> &FieldLabels) const
{
  const SConfig &VEVConfig = this->AllVEVConfigs[this->OrbifoldIndex][this->AllVEVConfigsIndex[this->OrbifoldIndex]];

  return global_GetIndicesOnlyFieldWithNumber(FieldLabels, VEVConfig.use_Labels, VEVConfig.Fields);
}


/* ########################################################################################
######   GetLocalization(const string &Localization, CSpaceGroupElement &...) const  ######
######                                                                               ######
######   Version: 19.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Localization : string                                                    ######
######   output:                                                                     ######
######   2) result       : the space-group element of the fixed point specified by   ######
######                     the string "Localization"                                 ######
######   return value    : localization found?                                       ######
######################################################################################## */
bool CPrompt::GetLocalization(const string &Localization, CSpaceGroupElement &result) const
{
  if (Localization.substr(0,8) == "loc of ")
  {
    vector<string> FieldLabels;
    FieldLabels.push_back(Localization.substr(8, string::npos));
    vector<unsigned> FieldIndices = GetIndices(FieldLabels);

    if (FieldIndices.size() == 1)
    {
      result = this->AllVEVConfigs[this->OrbifoldIndex][this->AllVEVConfigsIndex[this->OrbifoldIndex]].Fields[FieldIndices[0]].SGElement;
      return true;
    }
  }
  else
  {
    vector<int> tmp_result;
    convert_string_to_vector_of_int(Localization, tmp_result);

    if (tmp_result.size() == 9)
    {
      result.Set_m(tmp_result[0]);
      result.Set_n(tmp_result[1]);
      result.Set_k(tmp_result[2]);
      result.Set_n_alpha(0, tmp_result[3]);
      result.Set_n_alpha(1, tmp_result[4]);
      result.Set_n_alpha(2, tmp_result[5]);
      result.Set_n_alpha(3, tmp_result[6]);
      result.Set_n_alpha(4, tmp_result[7]);
      result.Set_n_alpha(5, tmp_result[8]);
      return true;
    }
    else
    {
      if (tmp_result.size() == 3)
      {
        result.Set_m(tmp_result[0]);
        result.Set_n(tmp_result[1]);
        result.Set_k(tmp_result[2]);
        result.Set_n_alpha(0, 0);
        result.Set_n_alpha(1, 0);
        result.Set_n_alpha(2, 0);
        result.Set_n_alpha(3, 0);
        result.Set_n_alpha(4, 0);
        result.Set_n_alpha(5, 0);
        return true;
      }
    }
  }
  if (this->print_output)
    (*this->Print.out) << "\n  " << this->Print.cbegin << "Localization \"" << Localization << "\" not found." << this->Print.cend << endl;
  return false;
}


/* ########################################################################################
######   GetLocalization(const string &Localization, CSpaceGroupElement &...) const  ######
######                                                                               ######
######   Version: 19.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) input           : a string either containing a fixed point / brane label ######
######                        or a constructing element                              ######
######   2) Orbifold        : the orbifold to look for "input"                       ######
######   output:                                                                     ######
######   3) FixedBraneFound : fixed point / brane "input" found in "Orbifold"?       ######
######   return value       : access to the fixed point / brane "input"              ######
######################################################################################## */
CFixedBrane &CPrompt::AccessFixedBrane(const string &input, COrbifold &Orbifold, bool &FixedBraneFound) const
{
  unsigned j = 0;

  size_t t1 = Orbifold.GetNumberOfSectors();
  size_t t2 = 0;

  // print summary of fixed point using the fixed point label
  for (unsigned i = 0; i < t1; ++i)
  {
    CSector &Sector = Orbifold.AccessSector(i);
    t2 = Sector.GetNumberOfFixedBranes();

    for (j = 0; j < t2; ++j)
    {
      CFixedBrane &FixedBrane = Sector.AccessFixedBrane(j);
      if (FixedBrane.GetFixedBraneLabel() == input)
      {
        FixedBraneFound = true;
        return FixedBrane;
      }
    }
  }

  // print summary of fixed point using the constructing element
  CSpaceGroupElement SGElement;

  if (this->GetLocalization(input, SGElement))
  {
    FixedBraneFound = true;
    return Orbifold.AccessFixedBrane(SGElement, FixedBraneFound);
  }
  FixedBraneFound = false;
  return Orbifold.AccessSector(0).AccessFixedBrane(0);
}

/* ########################################################################################
######   ExtractLabels(const SUSYMultiplet &Multiplet, string input, ...)            ######
######                                                                               ######
######   Version: 28.03.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Multiplet             : LeftChiral, RightChiral, Vector, ...             ######
######   2) input                 : string                                           ######
######                                                                               ######
######   output:                                                                     ######
######   2) FieldLabels           : vector of field labels                           ######
######                                                                               ######
######################################################################################## */
void CPrompt::ExtractLabels(const vector<SUSYMultiplet>  &Multiplet, string input, vector<string> &FieldLabels)
{
  // begin: remove VEV symbols "<" and ">" from the input string
  string::size_type loc1 = 0;
  while (loc1 != string::npos)
  {
    loc1 = input.find("<", 0);
    if (loc1 != string::npos)
      input.erase(loc1,1);

    loc1 = input.find(">", 0);
    if (loc1 != string::npos)
      input.erase(loc1,1);
  }
  // end: remove VEV symbols "<" and ">" from the input string

  const SConfig &VEVConfig = this->AllVEVConfigs[this->OrbifoldIndex][this->AllVEVConfigsIndex[this->OrbifoldIndex]];

  unsigned i = 0;

  string tmp = "";
  string tmp_string1 = "";
  string tmp_string2 = "";
  loc1 = 0;
  string::size_type loc2 = 0;
  string::size_type loc3 = 0;
  size_t t1 = 0;

  vector<string>::const_iterator pos;

  // run through the string
  while (loc1 != string::npos)
  {
    loc1 = input.find(" ", loc2);

    // extract label
    tmp_string1 = input.substr(loc2, loc1 - loc2);
    loc2 = loc1 + 1;
    //cout << "Label->" << tmp_string1 << "<-" << endl;

    // begin: find operators
    loc3  = tmp_string1.find_last_of("+-\\");
    if (loc3 != string::npos)
    {
      vector<string> FieldLabelsA;
      vector<string> FieldLabelsB;
      ExtractLabels(Multiplet, tmp_string1.substr(0, loc3),              FieldLabelsA);
      ExtractLabels(Multiplet, tmp_string1.substr(loc3+1, string::npos), FieldLabelsB);

      string op = tmp_string1.substr(loc3, 1);
      // begin: find difference A\B or A-B
      if ((op.substr(0,1) == "-") || (op.substr(0,1) == "\\"))
      {
        //cout << "minus" << endl;
        t1 = FieldLabelsA.size();
        for (i = 0; i < t1; ++i)
        {
          tmp = FieldLabelsA[i];
          // add label only if not contained before and if label not contained in set(B)
          if ((find(FieldLabelsB.begin(), FieldLabelsB.end(), tmp) == FieldLabelsB.end()) && (find(FieldLabels.begin(), FieldLabels.end(), tmp) == FieldLabels.end()))
            FieldLabels.push_back(tmp);
        }
      }
      // end: find difference A\B or A-B
      else
      // begin: find union A+B
      if (op.substr(0,1) == "+")
      {
        //cout << "plus" << endl;
        t1 = FieldLabelsA.size();
        for (i = 0; i < t1; ++i)
        {
          tmp = FieldLabelsA[i];
          // add label only if not contained before
          if (find(FieldLabels.begin(), FieldLabels.end(), tmp) == FieldLabels.end())
            FieldLabels.push_back(tmp);
        }
        t1 = FieldLabelsB.size();
        for (i = 0; i < t1; ++i)
        {
          tmp = FieldLabelsB[i];
          // add label only if not contained before
          if (find(FieldLabels.begin(), FieldLabels.end(), tmp) == FieldLabels.end())
            FieldLabels.push_back(tmp);
        }
      }
      // end: find union A+B
    }
    // end: find operators
    else
    {
      //cout << "no operator" << endl;
      // begin: insert all fields
      if (tmp_string1 == "*")
      {
      t1 = VEVConfig.Fields.size();

        for (i = 0; i < t1; ++i)
        {
		 for (int j=0; j<Multiplet.size(); j++)
         {
          const CField &Field = VEVConfig.Fields[i];

          if (Field.Multiplet == Multiplet[j])
          {
            tmp_string2 = Field.Labels[VEVConfig.use_Labels];
            std::ostringstream os;
            os << Field.Numbers[VEVConfig.use_Labels];
            tmp_string2 += "_";
            tmp_string2 += os.str();

            // add label only if not contained before
            if (find(FieldLabels.begin(), FieldLabels.end(), tmp_string2) == FieldLabels.end())
              FieldLabels.push_back(tmp_string2);
          }
         }
        }
        return;
      }
      // end: insert all fields

      pos = find(VEVConfig.NamesOfSetsOfFields.begin(), VEVConfig.NamesOfSetsOfFields.end(), tmp_string1);
      // begin: "tmp_string1" labels a field or several fields
      if (pos == VEVConfig.NamesOfSetsOfFields.end())
      {
        // a single field
        if (tmp_string1.find("_") != string::npos)
        {
          t1 = VEVConfig.Fields.size();
          for (i = 0; i < t1; ++i)
          {
            for (int j=0; j<Multiplet.size(); j++)
            {
             const CField &Field = VEVConfig.Fields[i];

             if (Field.Multiplet == Multiplet[j])
             {
              tmp_string2 = Field.Labels[VEVConfig.use_Labels];
              std::ostringstream os;
              os << Field.Numbers[VEVConfig.use_Labels];
              tmp_string2 += "_";
              tmp_string2 += os.str();

              // add label only if not contained before
              if ((tmp_string2 == tmp_string1) && (find(FieldLabels.begin(), FieldLabels.end(), tmp_string2) == FieldLabels.end()))
                FieldLabels.push_back(tmp_string2);
             }
           }
          }
        }
        // several fields
        else
        {
          t1 = VEVConfig.Fields.size();
          for (i = 0; i < t1; ++i)
          {
            for (int j=0; j<Multiplet.size(); j++)
            {
             const CField &Field = VEVConfig.Fields[i];
             if (Field.Multiplet == Multiplet[j])
             {
              tmp_string2 = Field.Labels[VEVConfig.use_Labels];
              std::ostringstream os;
              os << Field.Numbers[VEVConfig.use_Labels];
              tmp_string2 += "_";
              tmp_string2 += os.str();

              // add label only if not contained before
              if ((Field.Labels[VEVConfig.use_Labels] == tmp_string1) && (find(FieldLabels.begin(), FieldLabels.end(), tmp_string2) == FieldLabels.end()))
                FieldLabels.push_back(tmp_string2);
             }
           }
          }
        }
      }
      // end: "tmp_string1" labels a field or several fields
      // begin: "tmp_string1" labels a set
      else
      {
        const vector<unsigned> &CurrentSet = VEVConfig.SetsOfFields[distance(VEVConfig.NamesOfSetsOfFields.begin(), pos)];

        t1 = CurrentSet.size();
        for (i = 0; i < t1; ++i)
        {
        for (int j=0; j<Multiplet.size(); j++)
         {
          const CField &Field = VEVConfig.Fields[CurrentSet[i]];
          if (Field.Multiplet == Multiplet[j])
          {
            tmp_string2 = Field.Labels[VEVConfig.use_Labels];
            std::ostringstream os;
            os << Field.Numbers[VEVConfig.use_Labels];
            tmp_string2 += "_";
            tmp_string2 += os.str();

            // add label only if not contained before
            if (find(FieldLabels.begin(), FieldLabels.end(), tmp_string2) == FieldLabels.end())
              FieldLabels.push_back(tmp_string2);
          }
         }
        }
      }
      // end: "tmp_string1" labels a set
    }
  }
}


/* ########################################################################################
######   FindSUSYType(string &input_string, ...) const                               ######
######                                                                               ######
######   Version: 04.07.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) input_string          : input command string                             ######
######   2) NumberOfSupersymmetry : number of SUSY                                   ######
######   output:                                                                     ######
######   3) Multiplet             : the SUSY multiplet type, e.g. "LeftChiral"       ######
######   return value             : did a problem occur?                             ######
###########################################################################################
######   description:                                                                ######
######   Finds SUSY type in "input_string".                                          ######
######################################################################################## */
bool CPrompt::FindSUSYType(string &input_string, int NumberOfSupersymmetry, SUSYMultiplet &Multiplet) const
{
  switch (NumberOfSupersymmetry)
  {
    case 1:
    {
      Multiplet = LeftChiral;
      break;
    }
    case 2:
    {
      Multiplet = AnyKind;
      break;
    }
  }

  if (this->FindParameterType1(input_string, "left-chiral"))
  {
    Multiplet = LeftChiral;
    return true;
  }
  if (this->FindParameterType1(input_string, "right-chiral"))
  {
    Multiplet = RightChiral;
    return true;
  }
  if (this->FindParameterType1(input_string, "anykind"))
  {
    Multiplet = AnyKind;
    return true;
  }
  if (this->FindParameterType1(input_string, "vectorcc"))
  {
    Multiplet = VectorCC;
    return true;
  }
  if (this->FindParameterType1(input_string, "vector"))
  {
    Multiplet = Vector;
    return true;
  }
  if (this->FindParameterType1(input_string, "halfhyper"))
  {
    Multiplet = Halfhyper;
    return true;
  }
  if (this->FindParameterType1(input_string, "hyper"))
  {
    Multiplet = Hyper;
    return true;
  }
  if (this->FindParameterType1(input_string, "gravitycc"))
  {
    Multiplet = GravityCC;
    return true;
  }
  if (this->FindParameterType1(input_string, "gravity"))
  {
    Multiplet = Gravity;
    return true;
  }
  if (this->FindParameterType1(input_string, "modulus"))
  {
    Multiplet = LCModulus;
    return true;
  }
  if (this->FindParameterType1(input_string, "moduluscc"))
  {
    Multiplet = RCModulus;
    return true;
  }
  return false;
}



/* ########################################################################################
######   FindConditionsAndFilterFieldIndices(...) const                              ######
######                                                                               ######
######   Version: 26.01.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) input_string : if-conditions are searched in and removed from this       ######
######                     string                                                    ######
######   2) FieldIndices : vector of field indices; the ones that do not fulfill all ######
######                     all conditions will be removed                            ######
######   output:                                                                     ######
######   return value    : did a problem occur?                                      ######
###########################################################################################
######   description:                                                                ######
######   First call "FindConditions(...)" then "ApplyConditions(...)".               ######
######################################################################################## */
bool CPrompt::FindConditionsAndFilterFieldIndices(string &input_string, vector<unsigned> &FieldIndices) const
{
  vector<SCondition> Conditions;
  if (this->FindConditions(input_string, Conditions))
  {
    if (!this->ApplyConditions(Conditions, FieldIndices))
    {
      cout << "Warning in bool CPrompt::FindConditionsAndFilterFieldIndices(...) const: Could not apply the conditions. Return false." << endl;
      return false;
    }
  }
  return true;
}



/* ########################################################################################
######   FindConditions(string &input_string, vector<SCondition> &Conditions) const  ######
######                                                                               ######
######   Version: 16.05.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) input_string : if-conditions are searched in and removed from this       ######
######                     string                                                    ######
######   2) Conditions   : set of conditions extracted from the input string         ######
######   output:                                                                     ######
######   return value    : are conditions found in the input string?                 ######
###########################################################################################
######   description:                                                                ######
######   Interpret if-conditions in "input_string" as "SCondition"s and store them   ######
######   in "Conditions".                                                            ######
######################################################################################## */
bool CPrompt::FindConditions(string &input_string, vector<SCondition> &Conditions) const
{
  if (Conditions.size() != 0)
  {
    cout << "Warning in bool CPrompt::FindConditions(...) const: Set of conditions is not empty. Now cleared." << endl;
    Conditions.clear();
  }

  const COrbifold      &Orbifold   = this->Orbifolds[this->OrbifoldIndex];
  const SSymmetryGroup &GaugeGroup = Orbifold.StandardConfig.SymmetryGroup;

  unsigned index = 0;
  bool condition_ok = true;
  string tmp_string = "";

  // begin: find conditions
  while (this->FindParameterType2(input_string, "if(", tmp_string))
  {
    condition_ok = true;
    SCondition NewCondition;
    NewCondition.ith_entry   = 0;
    NewCondition.value       = 0;
    NewCondition.field_label = "";

    // begin: left hand side of the condition
    // (p_sh)^2 of the left-moving shifted momentum p_sh of the field
    if (this->FindParameterType1(tmp_string, "length"))
      NewCondition.ConditionTypeId = 0;
    else
    // vev of the field
    if (this->FindParameterType1(tmp_string, "vev"))
      NewCondition.ConditionTypeId = 1;
    else
    // B-L charge of the field
    if (this->FindParameterType1(tmp_string, "B-L"))
      NewCondition.ConditionTypeId = 2;
    else
    // i = 0,1,2,.. i-th accidental U(1) charge of the field
    if (this->FindParameterType3(tmp_string, "acc. Q_", index))
    {
      NewCondition.ConditionTypeId = 4;
      NewCondition.ith_entry = index-1;
    }
    else
    // i = 1,2,... i-th U(1) charge of the field
    if (this->FindParameterType3(tmp_string, "Q_", index))
    {
      if ((index < 1) || (index > GaugeGroup.GaugeGroup.u1directions.size()))
      {
        condition_ok = false;
        (*this->Print.out) << "  " << this->Print.cbegin << "i-th U(1), i = " << index << " out of range." << this->Print.cend << endl;
      }

      NewCondition.ConditionTypeId = 3;
      NewCondition.ith_entry = index-1;
    }
    else
    // dimension of the irrep with respect to the i-th gauge group factor
    if (this->FindParameterType3(tmp_string, "rep_", index))
    {
      NewCondition.ConditionTypeId = 5;
      NewCondition.ith_entry = index-1;
    }
    else
    // i = 1, ..., 16 i-th entry of the left-moving shifted momentum p_sh of the field (first p_sh if non-Abelian representation)
    if (this->FindParameterType3(tmp_string, "p_sh_", index))
    {
      if ((index < 1) || (index > 16))
      {
        condition_ok = false;
        (*this->Print.out) << "  " << this->Print.cbegin << "p_sh_i, i = " << index << " out of range." << this->Print.cend << endl;
      }

      NewCondition.ConditionTypeId = 6;
      NewCondition.ith_entry = index-1;
    }
    else
    // i = 1,2,3,4 i-th entry of the fields q_sh-charge
    if (this->FindParameterType3(tmp_string, "q_sh_", index))
    {
      if ((index < 1) || (index > 4))
      {
        condition_ok = false;
        (*this->Print.out) << "  " << this->Print.cbegin << "q_sh_i, i = " << index << " out of range." << this->Print.cend << endl;
      }

      NewCondition.ConditionTypeId = 7;
      NewCondition.ith_entry = index-1;
    }
    else
    // number of string oscillators acing on the fields ground state
    if (this->FindParameterType1(tmp_string, "#osci."))
      NewCondition.ConditionTypeId = 8;
    else
    // field label
    if (this->FindParameterType1(tmp_string, "label"))
      NewCondition.ConditionTypeId = 9;
    else
    {
      condition_ok = false;
      (*this->Print.out) << "  " << this->Print.cbegin << "Condition \"" << tmp_string << "\": left hand side not known." << this->Print.cend << endl;
    }
    // end: left hand side of the condition

    // begin: equality sign of the condition
    if (this->FindParameterType1(tmp_string, " == "))
      NewCondition.ComparisonTypeId = 0;
    else
    if (this->FindParameterType1(tmp_string, " != "))
      NewCondition.ComparisonTypeId = 1;
    else
    if (this->FindParameterType1(tmp_string, " > "))
      NewCondition.ComparisonTypeId = 2;
    else
    if (this->FindParameterType1(tmp_string, " >= "))
      NewCondition.ComparisonTypeId = 3;
    else
    if (this->FindParameterType1(tmp_string, " < "))
      NewCondition.ComparisonTypeId = 4;
    else
    if (this->FindParameterType1(tmp_string, " <= "))
      NewCondition.ComparisonTypeId = 5;
    else
    if (this->FindParameterType1(tmp_string, " involves "))
      NewCondition.ComparisonTypeId = 6;
    else
    if (this->FindParameterType1(tmp_string, " !involves "))
      NewCondition.ComparisonTypeId = 7;
    else
    {
      condition_ok = false;
      (*this->Print.out) << "  " << this->Print.cbegin << "Condition \"" << tmp_string << "\": equality sign not known." << this->Print.cend << endl;
    }
    // end: equality sign of the condition

    // begin: right hand side of the condition
    if (this->FindParameterType1(tmp_string, "even/odd"))
    {
      if ((NewCondition.ComparisonTypeId != 0) && (NewCondition.ComparisonTypeId != 1))
      {
        (*this->Print.out) << "  " << this->Print.cbegin << "Condition \"" << tmp_string << "\": even/odd only works with == or != ." << this->Print.cend << endl;
        condition_ok = false;
      }
      NewCondition.ValueTypeId = 2;
    }
    else
    if (this->FindParameterType1(tmp_string, "even"))
    {
      if ((NewCondition.ComparisonTypeId != 0) && (NewCondition.ComparisonTypeId != 1))
      {
        (*this->Print.out) << "  " << this->Print.cbegin << "Condition \"" << tmp_string << "\": even only works with == or != ." << this->Print.cend << endl;
        condition_ok = false;
      }
      NewCondition.ValueTypeId = 0;
    }
    else
    if (this->FindParameterType1(tmp_string, "odd"))
    {
      if ((NewCondition.ComparisonTypeId != 0) && (NewCondition.ComparisonTypeId != 1))
      {
        (*this->Print.out) << "  " << this->Print.cbegin << "Condition \"" << tmp_string << "\": odd only works with == or != ." << this->Print.cend << endl;
        condition_ok = false;
      }
      NewCondition.ValueTypeId = 1;
    }
    else
    {
      if (tmp_string.find_first_not_of("+-0123456789/") != string::npos)
      {
        if (NewCondition.ConditionTypeId != 9)
        {
          (*this->Print.out) << "  " << this->Print.cbegin << "Condition \"" << tmp_string << "\": rational number expected." << this->Print.cend << endl;
          condition_ok = false;
        }
        else
        {
          NewCondition.field_label = tmp_string;
          NewCondition.ValueTypeId = 4;
        }
      }
      else
      {
        // if value is integer (i.e. not a rational number)
        if (tmp_string.find("/") == string::npos)
          NewCondition.value = rational<int>(atoi(tmp_string.c_str()), 1);
        // if value is a rational number
        else
        {
          istringstream myStream(tmp_string);

          if (!(myStream >> NewCondition.value))
          {
            (*this->Print.out) << "  " << this->Print.cbegin << "Condition \"" << tmp_string << "\": rational number expected." << this->Print.cend << endl;
            condition_ok = false;
          }
        }
        NewCondition.ValueTypeId = 3;
      }
    }
    // end: right hand side of the condition

    if ((NewCondition.ConditionTypeId == 9)
    && (((NewCondition.ComparisonTypeId != 0) && (NewCondition.ComparisonTypeId != 1) && (NewCondition.ComparisonTypeId != 6) && (NewCondition.ComparisonTypeId != 7)) || (NewCondition.ValueTypeId != 4)))
    {
      (*this->Print.out) << "  " << this->Print.cbegin << "Condition for field labels failed." << this->Print.cend << endl;
      condition_ok = false;
    }
    if ((NewCondition.ConditionTypeId != 9)
    && ((NewCondition.ComparisonTypeId == 6) || (NewCondition.ComparisonTypeId == 7) || (NewCondition.ValueTypeId == 4)))
    {
      (*this->Print.out) << "  " << this->Print.cbegin << "Condition only works for field labels." << this->Print.cend << endl;
      condition_ok = false;
    }

    if (condition_ok)
    {
      Conditions.push_back(NewCondition);
      /*cout << "new condition" << endl;
      cout << "ConditionTypeId:  " << NewCondition.ConditionTypeId << endl;
      cout << "ith_entry:        " << NewCondition.ith_entry << endl;
      cout << "ComparisonTypeId: " << NewCondition.ComparisonTypeId << endl;
      cout << "ValueTypeId:      " << NewCondition.ValueTypeId << endl;
      cout << "value:            " << NewCondition.value << endl;
      cout << "field_label:      " << NewCondition.field_label << endl;*/
    }
  }
  if (Conditions.size() != 0)
    return true;
  else
    return false;
}



/* ########################################################################################
######   ApplyConditions(vector<string> &Commands, unsigned &exec_command)           ######
######                                                                               ######
######   Version: 10.10.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Commands     : list of prompt commands                                   ######
######   2) exec_command : index of the command that shall be executed next          ######
######   output:                                                                     ######
######   return value    : did a problem occur?                                      ######
###########################################################################################
######   description:                                                                ######
######   Evaluate "@if(..)", "@else", "@endif".                                      ######
######################################################################################## */
bool CPrompt::ApplyConditions(vector<string> &Commands, unsigned &exec_command)
{
  if (exec_command >= Commands.size())
    return true;

  string current_command = Commands[exec_command];
  string condition = "";
  if (!this->FindCommandType1(current_command, "@if(", condition))
    return true;

  condition = "";
  vector<string> PartsOfCondition;

  const SConfig &VEVConfig = this->AllVEVConfigs[this->OrbifoldIndex][this->AllVEVConfigsIndex[this->OrbifoldIndex]];

  bool condition_fulfilled = true;
  while (condition_fulfilled && this->FindParameterType2(current_command, "if(", condition))
  {
    condition_fulfilled = false;

    PartsOfCondition.clear();
    global_DecomposeString(condition, PartsOfCondition, " ");

    if (PartsOfCondition[0] == "set")
    {
      if (PartsOfCondition.size() == 4)
      {
        unsigned index = 0;
        if (this->MessageXNotKnown(VEVConfig.NamesOfSetsOfFields, "Set", PartsOfCondition[1], index))
          return false;

        if (PartsOfCondition[2] == "is")
        {
          if (PartsOfCondition[3] == "empty")
            condition_fulfilled = (VEVConfig.SetsOfFields[index].size() == 0);
          else
          if (PartsOfCondition[3] == "!empty")
            condition_fulfilled = (VEVConfig.SetsOfFields[index].size() != 0);
        }
      }
    }
    else
    if ((PartsOfCondition[0] == "hidden") && (PartsOfCondition[1] == "gauge") && (PartsOfCondition[2] == "group"))
    {
      string hidden_GG  = "";
      string comparison = "";
      unsigned index_from = 0;
      unsigned index_to   = VEVConfig.SymmetryGroup.GaugeGroup.factor.size();

      // for example: @if(hidden gauge group involves A1)
      if (PartsOfCondition.size() == 5)
      {
        comparison = PartsOfCondition[3];
        hidden_GG  = PartsOfCondition[4];
      }
      else
      // for example: @if(hidden gauge group from second E8 !involves E)
      if ((PartsOfCondition.size() == 8) && (PartsOfCondition[3] == "from") && (PartsOfCondition[5] == "E8"))
      {
        comparison = PartsOfCondition[6];
        hidden_GG  = PartsOfCondition[7];

        bool error = false;
        if ((PartsOfCondition[4] == "first"))
          index_to = VEVConfig.SymmetryGroup.Position_of_and_in_GaugeGroup;
        else
        if ((PartsOfCondition[4] == "second"))
          index_from = VEVConfig.SymmetryGroup.Position_of_and_in_GaugeGroup;
        else
          error = true;

        if (error || (this->Orbifolds[this->OrbifoldIndex].OrbifoldGroup.GetLattice() != E8xE8))
        {
          (*this->Print.out) << "  " << this->Print.cbegin << "Warning in @if condition: hidden gauge group cannot originate from E8." << this->Print.cend << endl;
          return false;
        }
      }
      string   algebra = "";
      unsigned rank    = 0;

      size_t pos = hidden_GG.find_first_of("0123456789", 0);

      const bool with_rank = (pos != string::npos);
      if (with_rank)
      {
        algebra = hidden_GG.substr(0,pos);
        rank = (unsigned)atoi((hidden_GG.substr(pos, string::npos)).c_str()) - 1;
        if (algebra == "E")
          ++rank;
      }
      else
        algebra = hidden_GG;

      if ((algebra != "SU") && (algebra != "SO") && (algebra != "E"))
      {
        (*this->Print.out) << "  " << this->Print.cbegin << "Warning in @if condition: Gauge group \"" << hidden_GG << "\" not of ADE type." << this->Print.cend << endl;
        return false;
      }

      if ((comparison != "involves") && (comparison != "!involves"))
      {
        (*this->Print.out) << "  " << this->Print.cbegin << "Warning in @if condition: Comparision must be \"involves\" or \"!involves\"." << this->Print.cend << endl;
        return false;
      }

      if (comparison == "!involves")
      {
        condition_fulfilled = true;
        for (unsigned i = index_from; condition_fulfilled && (i < index_to); ++i)
        {
          // analyse hidden gauge group factors only
          if (find(VEVConfig.SymmetryGroup.observable_sector_GGs.begin(), VEVConfig.SymmetryGroup.observable_sector_GGs.end(), i) == VEVConfig.SymmetryGroup.observable_sector_GGs.end())
          {
            const gaugeGroupFactor<double> &ggf = VEVConfig.SymmetryGroup.GaugeGroup.factor[i];

            if (((ggf.algebra[0] == 'A') && (algebra == "SU")) ||
                ((ggf.algebra[0] == 'D') && (algebra == "SO")) ||
                ((ggf.algebra[0] == 'E') && (algebra == "E")))
            {
              if (with_rank)
              {
                if (ggf.rank == rank)
                  condition_fulfilled = false;
              }
              else
                condition_fulfilled = false;
            }
          }
        }
      }
      else
      if (comparison == "involves")
      {
        condition_fulfilled = false;
        for (unsigned i = index_from; !condition_fulfilled && (i < index_to); ++i)
        {
          // analyse hidden gauge group factors only
          if (find(VEVConfig.SymmetryGroup.observable_sector_GGs.begin(), VEVConfig.SymmetryGroup.observable_sector_GGs.end(), i) == VEVConfig.SymmetryGroup.observable_sector_GGs.end())
          {
            const gaugeGroupFactor<double> &ggf = VEVConfig.SymmetryGroup.GaugeGroup.factor[i];

            if (((ggf.algebra[0] == 'A') && (algebra == "SU")) ||
                ((ggf.algebra[0] == 'D') && (algebra == "SO")) ||
                ((ggf.algebra[0] == 'E') && (algebra == "E")))
            {
              if (with_rank)
              {
                if (ggf.rank == rank)
                  condition_fulfilled = true;
              }
              else
                condition_fulfilled = true;
            }
          }
        }
      }
    }
    else
    if (PartsOfCondition[0] == "#fields")
    {
      if ((PartsOfCondition.size() == 6) && (PartsOfCondition[1] == "in") && (PartsOfCondition[2] == "set"))
      {
        unsigned index = 0;
        if (this->MessageXNotKnown(VEVConfig.NamesOfSetsOfFields, "Set", PartsOfCondition[3], index))
          return false;

        if (PartsOfCondition[5].find_first_not_of(" 0123456789") != string::npos)
        {
          if (this->print_output)
            (*this->Print.out) << "\n  " << this->Print.cbegin << PartsOfCondition[5] << " not a number." << this->Print.cend << "\n" << endl;

          return false;
        }

        const unsigned number_of_fields = (unsigned)atoi(PartsOfCondition[5].c_str());

        if (PartsOfCondition[4] == ">")
          condition_fulfilled = (VEVConfig.SetsOfFields[index].size() > number_of_fields);
        else
        if (PartsOfCondition[4] == "<")
          condition_fulfilled = (VEVConfig.SetsOfFields[index].size() < number_of_fields);
        else
        if (PartsOfCondition[4] == ">=")
          condition_fulfilled = (VEVConfig.SetsOfFields[index].size() >= number_of_fields);
        else
        if (PartsOfCondition[4] == "<=")
          condition_fulfilled = (VEVConfig.SetsOfFields[index].size() <= number_of_fields);
        else
        if (PartsOfCondition[4] == "==")
          condition_fulfilled = (VEVConfig.SetsOfFields[index].size() == number_of_fields);
        else
        if (PartsOfCondition[4] == "!=")
          condition_fulfilled = (VEVConfig.SetsOfFields[index].size() != number_of_fields);
      }
    }
  }
  ++exec_command;

  if (condition_fulfilled)
  {
    unsigned number_of_open_if = 1;

    for (unsigned i = exec_command; i < Commands.size(); ++i)
    {
      if (this->FindCommandType1(Commands[i], "@if(", condition))
        ++number_of_open_if;
      else
      if ((number_of_open_if == 1) && this->FindCommandType0(Commands[i], "@else"))
      {
        unsigned j = 0;
        for (j = i; (number_of_open_if != 0) && (j < Commands.size()); ++j)
        {
          if (this->FindCommandType1(Commands[j], "@if(", condition))
            ++number_of_open_if;

          if (this->FindCommandType0(Commands[j], "@endif"))
            --number_of_open_if;
        }
        Commands.erase(Commands.begin()+i, Commands.begin()+j);

        return true;
      }
      else
      if (this->FindCommandType0(Commands[i], "@endif"))
      {
        --number_of_open_if;

        if (number_of_open_if == 0)
        {
          Commands.erase(Commands.begin()+i);
          return true;
        }
      }
    }
  }
  else
  {
    unsigned number_of_open_if = 1;

    for (unsigned i = exec_command; i < Commands.size(); ++i)
    {
      if ((number_of_open_if == 1) && this->FindCommandType0(Commands[i], "@else"))
      {
        Commands.erase(Commands.begin()+exec_command, Commands.begin()+i+1);

        // search and delete corresponding "@endif"
        unsigned j = 0;
        for (j = exec_command; (number_of_open_if != 0) && (j < Commands.size()); ++j)
        {
          if (this->FindCommandType1(Commands[j], "@if(", condition))
            ++number_of_open_if;

          if (this->FindCommandType0(Commands[j], "@endif"))
            --number_of_open_if;
        }
        Commands.erase(Commands.begin()+j);

        return true;
      }
      else
      if (this->FindCommandType1(Commands[i], "@if(", condition))
        ++number_of_open_if;
      else
      if (this->FindCommandType0(Commands[i], "@endif"))
      {
        if (number_of_open_if == 1)
        {
          Commands.erase(Commands.begin()+exec_command, Commands.begin()+i+1);
          return true;
        }
        else
          --number_of_open_if;
      }
    }
  }

  return false;
}



/* ########################################################################################
######   ApplyConditions(const vector<SCondition> &Conditions, ...) const            ######
######                                                                               ######
######   Version: 19.08.2011                                                         ######
######   Check-Level: 1                                                              ######
######                                                                               ######
###########################################################################################
######   input:                                                                      ######
######   1) Conditions   : set of conditions to be imposed on the fields             ######
######   2) FieldIndices : vector of field indices; the ones that do not fulfill all ######
######                     conditions will be removed                                ######
######   output:                                                                     ######
######   return value    : did a problem occur?                                      ######
###########################################################################################
######   description:                                                                ######
######   Remove those field indices from "FieldIndices" that do not fulfill all      ######
######   conditions "Conditions".                                                    ######
######################################################################################## */
bool CPrompt::ApplyConditions(const vector<SCondition> &Conditions, vector<unsigned> &FieldIndices) const
{
  // Set the precision
  const double prec = 0.0001;

  const COrbifold       &Orbifold  = this->Orbifolds[this->OrbifoldIndex];
  const vector<CSector> &Sectors   = Orbifold.GetSectors();
  const SConfig         &VEVConfig = this->AllVEVConfigs[this->OrbifoldIndex][this->AllVEVConfigsIndex[this->OrbifoldIndex]];
  const vector<CField>  &Fields    = VEVConfig.Fields;

  unsigned j = 0;

  const size_t c1 = Conditions.size();
  if (c1 == 0)
    return true;

  bool equal = true;
  bool field_ok = true;
  double TestObject = 0.0;

  vector<unsigned> NewFieldIndices;

  const size_t s1 = FieldIndices.size();
  for (unsigned i = 0; i < s1; ++i)
  {
    const CField &Field = Fields[FieldIndices[i]];

    field_ok = true;
    for (j = 0; field_ok && (j < c1); ++j)
    {
      const SCondition &Condition = Conditions[j];

      // set "TestObject", specified by "ConditionTypeId"
      if (Condition.ConditionTypeId == 9)
      {
        const string &label = Field.Labels[VEVConfig.use_Labels];
        switch (Condition.ComparisonTypeId)
        {
          case 0: // equal
          {
            if (label != Condition.field_label)
              field_ok = false;
            break;
          }
          case 1: // not equal
          {
            if (label == Condition.field_label)
              field_ok = false;
            break;
          }
          case 6: // involves
          {
            if (label.find(Condition.field_label, 0) == string::npos)
              field_ok = false;
            break;
          }
          case 7: // not involves
          {
            if (label.find(Condition.field_label, 0) != string::npos)
              field_ok = false;
            break;
          }
        }
      }
      else
      {
        TestObject = 0.0;

        switch (Condition.ConditionTypeId)
        {
          case 0: // length
          {
            TestObject = Field.GetLMWeight(0, Sectors).GetSqrTo(16);
            break;
          }
          case 1: // vev
          {
            TestObject = Field.VEVs.GetLength();
            break;
          }
          case 2: // B-L
          {
            TestObject = Field.BmLCharge;
            break;
          }
          case 3: // Q_i
          {
            TestObject = Field.U1Charges[Condition.ith_entry];
            break;
          }
          case 4: // accidental Q_i
          {
            const rational<CHugeInt> &AccU1Charge = Field.AccU1Charges[Condition.ith_entry];
            TestObject = ((long double)AccU1Charge.numerator().ToLongLongInt())/((long double)AccU1Charge.denominator().ToLongLongInt());
            break;
          }
          case 5: // rep_i
          {
            TestObject = Field.Dimensions[Condition.ith_entry].Dimension;
            break;
          }
          case 6: // left-moving shifted momentum p_sh_i
          {
            TestObject = Field.GetLMWeight(0, Sectors)[Condition.ith_entry];
            break;
          }
          case 7: // R-charge R_i
          {
            TestObject = Field.GetRMWeight(0, Sectors)[Condition.ith_entry];
            break;
          }
          case 8: // number of oscillators acting on the state
          {
            TestObject =  Field.GetNumberOfOscillators(Sectors);
            break;
          }
          default:
          {
            cout << "Warning in bool CPrompt::ApplyConditions(...) const: ConditionTypeId not known. Return false." << endl;
            return false;
          }
        }

        // == and !=
        if ((equal = (Condition.ComparisonTypeId == 0)) || (Condition.ComparisonTypeId == 1))
        {
          const rational<int> RationalTestObject = D2Rat(TestObject);

          switch (Condition.ValueTypeId)
          {
            case 0: // equal or not-equal to even
            {
              // field is not ok if "TestObject" is not an integer or not even
              if (equal && ((RationalTestObject.denominator() != 1) || (RationalTestObject.numerator() % 2 != 0)))
                field_ok = false;
              else
              // field is not ok if "TestObject" is even
              if (!equal && (RationalTestObject.denominator() == 1) && (RationalTestObject.numerator() % 2 == 0))
                field_ok = false;

              break;
            }
            case 1: // equal or not-equal to odd
            {
              // field is not ok if "TestObject" is not an integer or even
              if (equal && ((RationalTestObject.denominator() != 1) || (RationalTestObject.numerator() % 2 == 0)))
                field_ok = false;
              else
              // field is not ok if "TestObject" is odd
              if (!equal && (RationalTestObject.denominator() == 1) && (RationalTestObject.numerator() % 2 != 0))
                field_ok = false;

              break;
            }
            case 2: // equal or not-equal to (even/odd)
            {
              // field is not ok if "TestObject" is the numerator is odd or the denominator is even
              if (equal && ((RationalTestObject.numerator() % 2 != 0) || (RationalTestObject.denominator() % 2 == 0)))
                field_ok = false;
              else
              // field is not ok if "TestObject" is the numerator is even and the denominator is odd
              if (!equal && (RationalTestObject.numerator() % 2 == 0) && (RationalTestObject.denominator() % 2 != 0))
                field_ok = false;

              break;
            }
            case 3: // equal or not-equal to rational number
            {
              rational<int> diff = RationalTestObject - Condition.value;
              long double tmp_diff = fabs(((long double)diff.numerator())/((long double)diff.denominator()));

              // field is not ok if "TestObject" is not equal to "value"
              if (equal && (tmp_diff > prec))
                field_ok = false;
              else
              // field is not ok if "TestObject" is equal to "value"
              if (!equal && (tmp_diff < prec))
                field_ok = false;

              break;
            }
            default:
            {
              cout << "Warning in bool CPrompt::ApplyConditions(...) const: ValueTypeId not known. Return false." << endl;
              return false;
            }
          }
        }
        // the inequalities only work with a rational number to compare
        else
        {
          if (Condition.ValueTypeId != 3)
          {
            cout << "Warning in bool CPrompt::ApplyConditions(...) const: Inequalities only work with rational numbers. Return false." << endl;
            return false;
          }
          long double value = ((long double)Condition.value.numerator())/((long double)Condition.value.denominator());

          // >
          if (Condition.ComparisonTypeId == 2)
          {
            if ((fabs(TestObject - value) < prec) || (TestObject + prec < value))
              field_ok = false;
          }
          else
          // >=
          if (Condition.ComparisonTypeId == 3)
          {
            if (TestObject + prec < value)
              field_ok = false;
          }
          else
          // <
          if (Condition.ComparisonTypeId == 4)
          {
            if ((fabs(TestObject - value) < prec) || (TestObject > value + prec))
              field_ok = false;
          }
          else
          // <=
          if (Condition.ComparisonTypeId == 5)
          {
            if (TestObject > value + prec)
              field_ok = false;
          }
        }
      }
    }
    if (field_ok)
      NewFieldIndices.push_back(FieldIndices[i]);
  }
  FieldIndices = NewFieldIndices;
  return true;
}
