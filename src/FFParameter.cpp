#include "../include/FFParameter.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include "../include/utility.h"
FFParameter::AtomTypes::AtomTypes():
                       type_map_(),
                       type_(){}
                       
FFParameter::AtomTypes::AtomTypes(const AtomTypes& to_copy):
                           type_map_(to_copy.type_map_),
                           type_(to_copy.type_){}
                           
const FFParameter::AtomTypes& FFParameter::AtomTypes::operator=(const AtomTypes& to_assign){
      if(!(this == &to_assign)){
           type_map_ = to_assign.type_map_;
           type_ = to_assign.type_;
           return *this;
           }
	  return *this;
      }
      
FFParameter::AtomTypes::~AtomTypes(){}
void FFParameter::AtomTypes::clear(){
        type_map_.clear();
        type_.clear();
        }
FFParameter::FFParameter():
                           atom_types_(),
                           is_initialized_(false),
                           is_valid_(false){}

FFParameter::FFParameter(const string file_name):
                               is_initialized_(false),
                               is_valid_(false){
                               read_parameter(file_name);
                               }

FFParameter::FFParameter(const FFParameter& to_copy):
                               atom_types_(to_copy.atom_types_),
                               //atom_types_.type_(to_copy.atom_types_.type_),
                               is_initialized_(to_copy.is_initialized_){}

FFParameter::~FFParameter(){
                            clear();
                            }

void FFParameter::clear(){
                     atom_types_.clear();
                     is_initialized_ = false;
                     }

const FFParameter& FFParameter::operator=(const FFParameter& to_assign){
      if (!(this == &to_assign)){
         atom_types_ = to_assign.atom_types_;
      //atom_types_.type_ = to_assign.atom_types_.type;
         is_initialized_ = to_assign.is_initialized_;
      return *this;
      }
	  return *this;
}
      
FFParameter::AtomTypes& FFParameter::get_atomtypes(){
           return atom_types_;
           }
// get the atom numeric type according to given symbolic type name
// if failed returned -1           
int FFParameter::get_type(const string& symbol) const{
    map<string, int>::const_iterator iter = atom_types_.type_map_.find(symbol);
    //debug
    //cout<<symbol<<" "<<iter->second<<endl;
    if (iter != atom_types_.type_map_.end())
        return iter->second;
    else
        return -1;
}

string FFParameter::get_type_name(const int type) const{
       static const string empty;
       if (type<0 || type > int(atom_types_.type_.size()))
           return empty;
       else 
           return atom_types_.type_[type-1];
}

int FFParameter::get_num_types() const{
               return int(atom_types_.type_.size());
               }

bool FFParameter::is_valid(){
     return is_valid_;
     }

// read atom type information from param_in_each_section and store them in AtomTypes
bool FFParameter::initialize_atom_type(){
     // check if the param_map is empty or key value "atom" exists
     if (params_in_each_section.empty() || params_in_each_section.count("atom") == 0){
         cout<<"Error: FFParameter::initialize_atom_type(): No parameters are read or missing Atom type information"<<endl;
         is_valid_ = false;
         return false;
         }
     else{
          // pack the symbolic and numeric atom types defined in the force field into the atom_types_
         atom_types_.clear();
         vector<vector<string> >::iterator iter = params_in_each_section["atom"].begin();
         for ( ; iter!=params_in_each_section["atom"].end(); ++iter){
             atom_types_.type_.push_back((*iter)[1]);
             atom_types_.type_map_.insert(make_pair((*iter)[1],str2int((*iter)[0])));
             }
         is_valid_ = true;
         return true;
         }
     }
         
// parse each line of parameters in the parameter file into the map
// using the component name as the key
void FFParameter::read_parameter(const string& file_name){
     params_in_each_section.clear();
     ifstream fin(file_name.c_str());
     if(!fin){
        cout<<"FFParameter::read_parameter(): Parameter opening error! File may not exist!"<<endl;
        exit(1);
        }
     // read the file line by line
     string line;
     string word, key, old_key;
     bool first_time_reading = true;
     vector<string> inner_vec;
     vector<vector<string> > outer_vec;
     while (getline(fin, line)){
           if (line[0] == '#')
               continue;
           inner_vec.clear();
           istringstream stream(line);
           old_key = key;
           stream >> key;
           while (stream>>word)
               inner_vec.push_back(word);
           if (old_key == key || first_time_reading){
               outer_vec.push_back(inner_vec);
               first_time_reading = false;
               }
           else{
                params_in_each_section.insert(make_pair(old_key, outer_vec));
                outer_vec.clear();
                outer_vec.push_back(inner_vec);
                first_time_reading = true;
                }
     }
     if (initialize_atom_type())
         is_initialized_ = true;
     //cout<<params_in_each_section["bond"].size()<<" "<<params_in_each_section["angle"].size()
        // <<" "<<params_in_each_section["torsion"].size()<<" "<<params_in_each_section["vdw"].size()<<endl;
     }

 
                  
     
