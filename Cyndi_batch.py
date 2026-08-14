#!/usr/bin/env python

import os, sys, getopt
from subprocess import call


# global variables
InputFile_ = None
OutputPrefix_ = "Cyndi_"
NumConf_ = 100
NumGen_ = 200
NumRun_ = 1
PopSize_ = 100
CrossProb_ = 0.85
MutProb_ = 0.1
EnergyCutOff_ = 20.0
rmsdScaleFactor_ = 0.2
KeepInput_ = 'N'
MinConf_ = 'N'
SplitOut_ = 'N'
InputMode_ = "mmol2"
RandomSeed_ = 1000
# ---------------------------------------------------------------------------------------------------------------
def usage():
    print """
    Generate conformers for the listed ligands in batch mode.
    Usage: Cyndi_batch.py [options] input-mol2.list
    
    options:
    -h,     --help=                 Print help message
    -o,     --output-prefix=        Prefix for the resultant mol2 files
    -c,     --num-conf=             Maximum number of conformers for each ligands to be output
    -g,     --num-gen=              Maximum number of generations in MOGA
    -r,     --num-run=              Maximum number of MOGA
    -p,     --pop-size=             Population size in MOGA
    --cross-prob=                   Crossover probability in MOGA
    --mut-prob=                     Mutation probability in MOGA
    -k,     --keep-input=           Keep input conformer or not
    -m,     --min-conf=             Minimize generated conformers or not
    --energy-cutoff=                The energy cutoff for the output conformers
    --rmsd-factor=                  The rmsd scale factor to remove redundant conformers
    -i,     --input-mode=           Mol2 file containing multiple ligands (mmol2) or list file containing names of different mol2 files (list)
    -s,     --split-outfile=        Split the conformers ensembles for defferent ligands into seperate files.
    --random-seed=                  The seed (0-1)of random number generator, if unspecified, a seed will be generated from current time.
    """
def generate_parameters():
    paramfile = open("CyndiParam.in", 'w')
    paramfile.write("MOGA_Max_Conformers\t\t%i\n\
MOGA_Num_Objectives\t\t3\n\
MOGA_Max_Generation\t\t%i\n\
MOGA_Population_Size\t\t%i\n\
MOGA_Max_Run\t\t%i\n\
MOGA_Crossover_Probability\t\t%f\n\
MOGA_Mutation_Probability\t\t%f\n\
MOGA_SBX\t\t15 20\n\
MOGA_Epsilon_Quaternion\t\t 3 0.3 0.1 2\n\
MOGA_Keep_Input_Conformer[Y/N]\t\t%s\n\
MOGA_Optimize_Conformer[Y/N]\t\t%s\n\
MOGA_Energy_Cutoff\t\t%i\n\
MOGA_RMSD_Scale_Factor\t\t%f\n\
MOGA_Max_Opt_Iteration\t\t500\n\
MOGA_Split_Output[Y/N]\t\t%s\n"\
%(NumConf_, NumGen_, PopSize_, NumRun_, CrossProb_, MutProb_, KeepInput_, MinConf_, EnergyCutOff_, rmsdScaleFactor_, SplitOut_))
    paramfile.close()
    if(RandomSeed_ != 1000):
        paramfile = open("CyndiParam.in", 'a')
        paramfile.write("MOGA_Random_Seed[0-1]\t\t%f\n" %(RandomSeed_))
    
if __name__ == "__main__":
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'ho:c:g:p:kmi:sr:',
                                  ['help', 'output-prefix=', 'num-conf=',
                                   'num-gen=', 'pop-size=', 'cross-prob=',
                                   'mut-prob=', 'keep-input','min-conf', 'input-mode=',
                                   'split-outfile', 'num-run=', 'random-seed=',
                                   'energy-cutoff=', 'rmsd-factor='])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
# ---------------------------------------------------------------------------------------------------------
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-o", "--output-prefix"):
            OutputPrefix_ = a + '_'
        elif o in ("-c","--num-conf"):
            NumConf_ = int(a)
        elif o in ("-g", "--num-gen"):
            NumGen_ = int(a)
        elif o in ("-r", "--num-run"):
            NumRun_ = int(a)
        elif o in ("-p", "--pop-size"):
            PopSize_ = int(a)
        elif o == "-cross-prob":
            CrossProb_ = float(a)
        elif o == "--mut-prob":
            MutProb_ = float(a)
        elif o in ("-k", "--keep-input"):
            KeepInput_ = 'Y'
        elif o in ("-m", "--min-conf"):
            MinConf_ = 'Y'
        elif o == "--energy-cutoff":
            EnergyCutOff_ = float(a)
        elif o == "--rmsd-factor":
            rmsdScaleFactor_ = float(a)
        elif o in ("-i", "--input-mode"):
            if a != "mmol2" and a != "list":
                print "Please provide appropriate input-mode argument value: mmol2 or list"
                sys.exit()
            else:
                InputMode_ = a
        elif o in ("-s", "--split-output"):
            SplitOut_ = 'Y'
        elif o in ("-seed", "--random-seed"):
            RandomSeed_ = float(a)
    
    if len(args) != 1:
        usage()
        sys.exit()
    else:
        InputFile_ = args[0]
    
    generate_parameters()
    
    if InputMode_ == "mmol2":
        os.system("..\\Release\\Cyndi.exe -i %s -o %s -p CyndiParam.in" %(InputFile_, OutputPrefix_))
    elif InputMode_ == "list":
        for name in open(InputFile_, 'r').readlines():
            name = name.strip()
            if name.find(".mol2") == -1:
                print "Error! No appropriate mol2 file names in the input list file"
                sys.exit()
            os.system("..\\Release\\Cyndi.exe -i %s -o %s -p CyndiParam.in" %(name, OutputPrefix_))