#ifndef FFPARAMETER_H
#define FFPARAMETER_H

#include <vector>
#include <map>
#include <string>


using namespace std;




class FFParameter{
      public:
             struct AtomTypes{
             AtomTypes();
             AtomTypes(const AtomTypes&);
             const AtomTypes& operator=(const AtomTypes&);
             virtual ~AtomTypes();
             void clear();
             // numeric representation for each atom type
             map<string, int> type_map_;
             // symbolic names of each atom type
             vector<string> type_;
             };
             FFParameter();
             FFParameter(const string file_name);
             FFParameter(const FFParameter& to_copy);
             virtual ~FFParameter();
             const FFParameter& operator=(const FFParameter&);
             virtual void clear();
             inline bool is_initialized(){ return is_initialized_; }
             void read_parameter(const string& file_name);
             bool initialize_atom_type();
             AtomTypes& get_atomtypes();
             int get_type(const string& name) const;
             string get_type_name(const int type) const;
             int get_num_types() const;
             bool is_valid();
             // a map holding each line in parameter file
             map<string, vector<vector<string> > >params_in_each_section; 
      protected:
                AtomTypes atom_types_;
                bool is_initialized_;
                bool is_valid_;
                
      };
#endif             
