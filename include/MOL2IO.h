// This class reads/writes information of molecule from/into mol2(Tripos) file and parse 
// the information to pack into/from a MOL object

#ifndef MOL2IO_H
#define MOL2IO_H

#include <fstream>
#include <iostream>
#include <vector>

#include "Mol.h"
#include "utility.h"
#include "Vector3.h"

class MOL2IO{
      public:
             // default constructor
             MOL2IO();
             // constructor
             MOL2IO(const string&, const string&);
             // destructor
             virtual ~MOL2IO();
             void clear();
			 // open the specified mol2 file
			 bool open(const string& name, const string& mode);
             //void bond_to_mol(MOL::MOL& mol);
             // read a molecule from the mol2 file;
             bool read(MOL* mol);
			 // read a serial of conformers of one molecule from t he mol2 file, the name must have the format of "mol_name_N", where N = 1, 2, 3, 4...
			 bool readConfs(deque<MOL>& conf);
             //write a molecule into the mol2 file
             bool write(const MOL& mol);
             // load the first( or next )molecule from the mol2 file
             //extern MOL* read();
             
      private:
                void read_atom_section_();
                void read_bond_section_();
                void read_molecule_section_();
                void read_set_section_();
                void read_substructure_section_();
                bool build_all_(MOL* mol);
                void stream_clear_();
                struct AtomInfo{
                       int id;
                       string name;
                       vector3 position;
                       string symbolic_type;
                       int substructure_id;
                       int substructure_name;
                       double charge;
                       };
                struct BondInfo{
                       int id;
                       int atom1;
                       int atom2;
                       string type;
                       };
                struct MolInfo{
                       string name;
                       int num_of_atoms;
                       int num_of_bonds;
                       int num_of_substructures;
                       int num_of_features;
                       int num_of_sets;
                       string type;
                       string charge_type;
                       string annotate;
                       };
                struct SetInfo{
                       string name;
                       string type;
                       string sub_type;
                       string annotate;
                       int num_of_members;
                       std::vector<int> members;
                       };
                struct SubstrInfo{
                       string name;
                       int root_atom;
                       string type;
                       int dict_type;
                       string chain;
                       string sub_type;
                       int inter_bonds;
                       string comment;
                       };          
                std::vector<AtomInfo> atoms_;
                std::vector<BondInfo> bonds_;
                std::vector<SetInfo> sets_;
                std::vector<SubstrInfo> substructures_;
                MolInfo molecule_;
                
                int num_of_lines_;
                //string line_;
                //char buffer_[4096];
                
                ifstream input_;
                ofstream output_;
                // define a stream position mark, so we can read the molecules
                // sequentialliy from the multiple mol2 file.
                ifstream::pos_type mark_;
                // if we have read the first molecule, we set this flag value to true
                bool first_mol_flag_;
                // 
                bool has_sets_info_;
                bool has_substr_info_;
				// maps for fast access the radius and atomic weight for each type of element
				map<string, float> AtomRadius_;
				map<string, float> AtomWeight_;
};
#endif
