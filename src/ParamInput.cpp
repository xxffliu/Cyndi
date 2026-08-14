#include "../include/ParamInput.h"

MOGAParam MOGAParam_;

bool ReadParameter(string name)
{
	ifstream infile;
	infile.open(name.c_str(), ifstream::in);
	if(!infile)
	{
		cout<<"Error! Cannot access the parameter file with the name of "<<name<<endl;
		return false;
	}
	while(!infile.eof())
	{
		string line;
		while(getline(infile, line))
		{
			if(line.find("#") == 0 || line == "")
				continue;
			else
				cout<<line<<endl;
			vector<string> fields;
			stringstream sin;
			sin.str(line);
			string token;
			while(sin>>token)
				fields.push_back(token);
			if(fields.size() < 2)
				continue;
			//MOGA parameters start from here
			if(fields[0] == "MOGA_Max_Conformers")
				MOGAParam_.MaxNumConformer_ = str2int(fields[1]);
			else if(fields[0] == "MOGA_Num_Objectives")
				MOGAParam_.NumObjects_ = str2int(fields[1]);
			else if(fields[0] == "MOGA_Random_Seed")
			{
				MOGAParam_.BasicSeed_ = str2double(fields[1]);
				MOGAParam_.UseInputRandomSeed_ = true;
			}

			else if (fields[0] == "MOGA_Max_Generation")

				MOGAParam_.MaxNumGen_ = str2int(fields[1]);

			else if (fields[0] == "MOGA_Population_Size")

				MOGAParam_.PopSize_ = str2int(fields[1]);
			else if(fields[0] == "MOGA_Max_Run")
				MOGAParam_.MaxNumRun_ = str2int(fields[1]);

			else if (fields[0] == "MOGA_Crossover_Probability")

				MOGAParam_.PXover_ = str2double(fields[1]);

			else if (fields[0] == "MOGA_Mutation_Probability")

				MOGAParam_.PMutation_ = str2double(fields[1]);

			else if (fields[0] == "MOGA_SBX")
			{
				MOGAParam_.n_distribution_c = str2double(fields[1]);
				MOGAParam_.n_distribution_m = str2double(fields[1]);
			}

			else if(fields[0] == "MOGA_VDW_Energy_Epsilon")
				MOGAParam_.EPSILON_[0] = str2double(fields[1]);
			else if(fields[0] == "MOGA_Torsion_Energy_Epsilon")
				MOGAParam_.EPSILON_[1] = str2double(fields[1]);
			else if(fields[0] == "MOGA_RMSD_Epsilon")
				MOGAParam_.EPSILON_[2] = str2double(fields[1]);
			else if(fields[0] == "MOGA_Gyration_Radius_Epsilon")
				MOGAParam_.EPSILON_[3] = str2double(fields[1]);

			else if(fields[0] == "MOGA_Keep_Input_Conformer[Y/N]")

				if(fields[1] == "Y" ||fields[1] == "y" || fields[1] == "yes")
					MOGAParam_.KeepInputConformer_ = true;
				else
					MOGAParam_.KeepInputConformer_ = false;

			else if(fields[0] == "MOGA_ForceField_Type[TAFF/MMFF94]")
				if(fields[1] == "MMFF94")
					MOGAParam_.FFType_ = FF_MMFF94;
				else
					MOGAParam_.FFType_ = FF_TAFF;

			else if(fields[0] == "MOGA_Optimize_Input_Conformer[Y]/N]")

				if(fields[1] == "Y" ||fields[1] == "y" || fields[1] == "yes")
					MOGAParam_.OptimizeInputConformer_ = true;
				else
					MOGAParam_.OptimizeInputConformer_ = false;



			else if (fields[0] == "MOGA_Optimize_Conformer[Y/N]")

				if(fields[1] == "Y" ||fields[1] == "y" || fields[1] == "yes")
					MOGAParam_.OptimizeConformer_ = true;
				else
					MOGAParam_.OptimizeConformer_ = false;

			else if (fields[0] == "MOGA_Max_Opt_Iteration")
				MOGAParam_.MaxNumIteration_ = str2int(fields[1]);

			else if (fields[0] == "MOGA_Max_Opt_Gradient")
				if(MOGAParam_.OptimizeConformer_)
					MOGAParam_.MaxGrd_ = str2double(fields[1]);
				else
					continue;
			else if(fields[0] == "MOGA_Energy_Cutoff")
				MOGAParam_.EnergyCutoff_ = str2double(fields[1]);

			else if(fields[0] == "MOGA_RMSD_Scale_Factor")
				MOGAParam_.rmsdScaleFactor_ = str2double(fields[1]);

			else
				continue;
		}
	}
	return true;
}

void CyndiOptionParser(int argc, char* argv[])
{
	if(argc == 1)
	{
		CyndiUsage();
		exit(1);
	}
	else
	{
		for(int i = 1; i < argc;)
		{
			// input molecule name
			if(string(argv[i]) == "-input")
			{
				string molfile(argv[i+1]);
				if(molfile.find(".mol2") == string::npos)
				{
					cout<<"Inappropriate MOL2 file specified, program exiting..."<<endl;
					exit(1);
				}
				MOGAParam_.InputFile_ = molfile;
				i += 2;
			}
			// output molecule name
			else if(string(argv[i]) == "-output")
			{
				string molfile(argv[i+1]);
				if(molfile.find(".mol2") == string::npos)
				{
					cout<<"Output file must be in MOL2 format, program exiting..."<<endl;
					exit(1);
				}
				MOGAParam_.OutputFile_ = molfile;
				i += 2;
			}
			// max confs
			else if(string(argv[i]) == "-maxconfs")
			{
				MOGAParam_.MaxNumConformer_ = str2int(argv[i+1]);
				i += 2;
			}
			// MOGA generation
			else if(string(argv[i]) == "-gagen")
			{
				MOGAParam_.MaxNumGen_ = str2int(argv[i+1]);
				i += 2;
			}
			// MOGA population size
			else if(string(argv[i]) == "-gapop")
			{
				MOGAParam_.PopSize_ = str2int(argv[i+1]);
				i += 2;
			}
			// MOGA run times
			else if(string(argv[i]) == "-garun")
			{
				MOGAParam_.MaxNumRun_ = str2int(argv[i+1]);
				i += 2;
			}
			// MOGA crossover rate
			else if(string(argv[i]) == "-gaxover")
			{
				MOGAParam_.PXover_ = str2double(argv[i+1]);
				i += 2;
			}
			// MOGA mutation rate
			else if(string(argv[i]) == "-gamut")
			{
				MOGAParam_.PMutation_ = str2double(argv[i+1]);
				i += 2;
			}
			// objectives
			else if(string(argv[i]) == "-objs")
			{
				MOGAParam_.NumObjects_ = str2int(argv[i+1]);
				i += 2;
			}
			else if(string(argv[i]) == "-vdw-epsilon")
			{
				MOGAParam_.EPSILON_[0] = str2double(argv[i+1]);
				i += 2;
			}
			else if(string(argv[i]) == "-torsion-epsilon")
			{
				MOGAParam_.EPSILON_[1] = str2double(argv[i+1]);
				i += 2;
			}
			else if(string(argv[i]) == "-rmsd-epsilon")
			{
				MOGAParam_.EPSILON_[2] = str2double(argv[i+1]);
				i += 2;
			}
			else if(string(argv[i]) == "-gr-epsilon")
			{
				MOGAParam_.EPSILON_[3] = str2double(argv[i+1]);
				i += 2;
			}
			// forcefield
			else if(string(argv[i]) == "-ff")
			{
				string fftype = string(argv[i+1]);
				if(fftype == "MMFF94")
					MOGAParam_.FFType_ = FF_MMFF94;
				else if(fftype == "TAFF")
					MOGAParam_.FFType_ = FF_TAFF;
				else
				{
					cout<<"Unknown forcefield type "<<fftype<<". Only MMFF94 or TAFF supported."<<endl;
					exit(1);
				}
				i += 2;
			}
			// keep input conformer
			else if(string(argv[i]) == "-keepinput")
			{
				MOGAParam_.KeepInputConformer_ = true;
				i += 1;
			}
			// optimize input conformer
			else if(string(argv[i]) == "-optinput")
			{
				MOGAParam_.OptimizeInputConformer_ = true;
				i += 1;
			}
			// optimize output conformer
			else if(string(argv[i]) == "-optoutput")
			{
				MOGAParam_.OptimizeConformer_ = true;
				i += 1;
			}
			// energy cutoff
			else if(string(argv[i]) == "-energycutoff")
			{
				MOGAParam_.EnergyCutoff_ = str2double(argv[i+1]);
				i += 2;
			}
			// RMSD cutoff
			else if(string(argv[i]) == "-rmsd")
			{
				MOGAParam_.rmsdScaleFactor_ = str2double(argv[i+1]);
				i += 2;
			}
			// parameter file
			else if(string(argv[i]) == "-parm")
			{
				string para_name = string(argv[i+1]);
				if(ReadParameter(para_name) == false)
				{
					cout<<"Fail to read parameter file, program exiting..."<<endl;
					exit(1);
				}
				i += 2;
			}
			else
			{
				cout<<"Unknown argument: "<<argv[i]<<endl;
				CyndiUsage();
				exit(1);
			}
		}
	}
	return;
}