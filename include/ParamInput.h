/*
  Name: ParamInput.h
  Copyright: DDDC
  Author: Xiaofeng Liu
  Date: 09-01-08 10:39
  Description: Some defination for the global input parameters used in Cyndi.
*/
#ifndef INPUT_H
#define INPUT_H
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cassert>
#include "utility.h"
using namespace std;

enum FFTYPE
{
	FF_TAFF = 1,
	FF_MMFF94 = 2
};

// MOGA global parameters
class MOGAParam
{
public:
	// MOGA parameter;
	// If generate conformers on the fly, here define the maximum number of conformers used
	// for following matching
	int NumObjects_;
	int MaxNumConformer_;
    int MaxNumGen_;
    int PopSize_;
	int MaxNumRun_;
    float PXover_;
    float PMutation_;
    float n_distribution_c;
    float n_distribution_m;
	std::vector<float> EPSILON_;
	// toogle if reserve the input conformer
	bool KeepInputConformer_;
	// toogle if optimize the input conformer befor MOGA conformer generation
	bool OptimizeInputConformer_;
    // Force Field related parameter
    bool OptimizeConformer_;
	// force filed type: TAFF or MMFF94
	int FFType_;
	float EnergyCutoff_;
	float rmsdScaleFactor_;
    int MaxNumIteration_;
    float MaxGrd_;
	string InputFile_;
	string OutputFile_;
	float BasicSeed_;
	bool UseInputRandomSeed_;
	// v3 (2026): index of the first molecule to process in the input file
	// (0-based). Lets a batch driver split one multi-molecule file across
	// several independent Cyndi processes (one chunk each) and derive each
	// molecule's RNG seed from its GLOBAL index, so the union of chunk
	// outputs is identical to a single run over the whole file.
	int StartIndex_;
	int MaxMolecules_;   // 0 = all molecules from StartIndex_ to EOF
    MOGAParam():
	InputFile_("input.mol2"),OutputFile_("output.mol2"),NumObjects_(4),MaxNumConformer_(500),MaxNumGen_(500),MaxNumRun_(2),PopSize_(500),PXover_(0.8),PMutation_(0.08),n_distribution_c(15),rmsdScaleFactor_(0.3),n_distribution_m(20),EPSILON_(4,0.),UseInputRandomSeed_(false),KeepInputConformer_(false), OptimizeInputConformer_(false), EnergyCutoff_(20.0),OptimizeConformer_(false),MaxNumIteration_(100),MaxGrd_(0.01),BasicSeed_(0.34), FFType_(FF_MMFF94), StartIndex_(0), MaxMolecules_(0){}
};
extern MOGAParam MOGAParam_;

inline void CyndiUsage()
{
	cout<<"Cyndi -input input.mol2 -output output.mol2 [-parm Cyndi.parm]"<<endl;
}

bool ReadParameter(string name);

void CyndiOptionParser(int argc, char* argv[]);



#endif
       
