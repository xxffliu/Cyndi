#include "../include/TAFFStretch.h"
#include "../include/TAFF.h"
#include "../include/utility.h"
#include <cmath>
using namespace std;

// constructors
TAFFStretch::TAFFStretch():
                           FFComponent(),
                           fast_access(),
                           num_of_atom_types_(0),
                           stretch_data_holder_(){
                                                  // set component name
                                                  set_name("TAFF Stretch");
                                                  }
TAFFStretch::TAFFStretch(ForceField& ff):
                                     FFComponent(ff),
                                     fast_access(),
                                     num_of_atom_types_(0),
                                     stretch_data_holder_(){
                                                            set_name("TAFF Stretch");
                                                            }
// copy constructor
TAFFStretch::TAFFStretch(const TAFFStretch& to_copy):
                               FFComponent(to_copy){
                                                    fast_access = to_copy.fast_access;
                                                    num_of_atom_types_ = to_copy.num_of_atom_types_;
                                                    stretch_data_holder_ = to_copy.stretch_data_holder_;
                                                    }
// destructor
TAFFStretch::~TAFFStretch(){
                            fast_access.clear();
                            stretch_data_holder_.clear();
                            }

// extract bond stretch parameters from FFParameter object bonded to the force field
// and establish a hash table for fast access
bool TAFFStretch::extract_BS_parameters(FFParameter& ffp){
     if (!ffp.is_valid())
         return false;
     // build a two dim array of atom types and loop variable
     //FFParameter::AtomTypes& atom_types = ffp.get_atomtypes();
     int num_types = ffp.get_num_types();
     num_of_atom_types_ = num_types;
     for (int i = 1; i <=5; i++){
         //fast_access[i].num_of_atom_types = num_types;
         // allocate two onedimensional fields for the two parameters
         // and a two dimensional field of boolean variables.
         // we might have to delete old stuff lying around
         /*if (fast_access[i].k)
            delete [] fast_acess[i].k;
         if (fast_access[i].r0)
            delete [] fast_access[i].r0;
         if (fast_access[i].is_defined)
            delete [] fast_access[i].is_defined;
         fast_access[i].k  = new float[num_types * num_types];
         fast_access[i].r0 = new float[num_types * num_types];
         fast_access[i].is_defined = new bool[num_types * num_types];*/
         fast_access[i].k.clear();
         fast_access[i].k.resize(num_types * num_types);
         fast_access[i].r0.clear();
         fast_access[i].r0.resize(num_types * num_types);
         fast_access[i].is_defined.clear();
         fast_access[i].is_defined.resize(num_types * num_types);
         for (int j = 0; j < num_types * num_types; j++) 
			      fast_access[i].is_defined[j] = false;
     }
     int type1, type2, index, bond_order;
     string name_type1, name_type2;   
     for (vector<vector<string> >::iterator it = force_field_->get_parameters().params_in_each_section["bond"].begin();
         it != force_field_->get_parameters().params_in_each_section["bond"].end(); ++it){
               name_type1 = (*it)[0];
               name_type2 = (*it)[1];
               bond_order = get_bond_order((*it)[2]);
               type1 = ffp.get_type(name_type1);
               type2 = ffp.get_type(name_type2);
               if ((type1 == -1) || (type2 == -1)){
                   cout<<"error! no numeric atom type defined for atom "<<name_type1<<" or "<<name_type2<<endl;
                   return false;
                   }
               index = ((type1-1) * num_types + (type2-1));
               fast_access[bond_order].is_defined[index] = true;
               fast_access[bond_order].r0[index] = str2double((*it)[3]);
               fast_access[bond_order].k[index] = str2double((*it)[4]);
               //build the symmetric index
               index = ((type1-1) + num_types * (type2-1));
               fast_access[bond_order].is_defined[index] = true;
               fast_access[bond_order].r0[index] = str2double((*it)[3]);
               fast_access[bond_order].k[index] = str2double((*it)[4]);
     }
     return true;
}

// query a set of parameters has defined for a given combination of atoms;
bool TAFFStretch::has_params(int i, int j, int bond_order){
     if ((i < 0 || i > num_of_atom_types_) || (bond_order < 1 || bond_order > 5))
         return false;
     if ((j < 0 || j > num_of_atom_types_) || (bond_order < 1 || bond_order > 5))
         return false;
     return fast_access[bond_order].is_defined[(i-1)* num_of_atom_types_ + (j-1)];
     }
     
// return the parameters for a given atom type combination
TAFFStretch::ForceValues TAFFStretch::get_params(int i, int j, int bond_order){
                         TAFFStretch::ForceValues param;
                         if (assign_params(param, i, j, bond_order))
                             return param;
                         else{
                             cout<<"TAFFStretch::get_params():Warning: no bond stretch parameters are assigned for atom "<<
                             i<<" and "<<j<<endl;
                             param.k = 0.0;
                             param.r0 = 0.0;
                             return param;
                             } 
                         }

// assign the parameters to a given atom type combination
// if no parameters are define for this 2 combination, return false and nothing is changed
bool TAFFStretch::assign_params(TAFFStretch::ForceValues& param, int i, int j, int bond_order){
     if (has_params(i, j, bond_order)){
                       param.k = fast_access[bond_order].k[(i-1)*num_of_atom_types_+(j-1)];
                       param.r0 = fast_access[bond_order].r0[(i-1)*num_of_atom_types_+(j-1)];
                       return true;
                       }
     else
         return false;
}

// set up method
bool TAFFStretch::setup(){
     if (force_field_ == NULL){
                          cout<<"TAFFStretch::setup(): force field bound component can not be found"<<endl;
                          return false;
                          }
     // clear the parameter holder
     stretch_data_holder_.clear();
     // tempararily set this component enabled;
     setenabled(true);
     
     TAFF* taff_force_field = dynamic_cast<TAFF*>(force_field_);
     if ((taff_force_field == NULL) || taff_force_field->has_initialized_param()){
                          bool result = extract_BS_parameters(force_field_->get_parameters());
                          if (!result){
                                       cout << "TAFFStretch::setup(): can not found bond stretch section"<<endl;
                                       return false;
                                       }
     }                    
                          
     TAFFStretch::ForceValues value;
     // retrieve all stretch parameters
     for (BONDVec::iterator biter = force_field_->get_bonds().begin();
          biter != force_field_->get_bonds().end(); ++biter){
                int bond_order = get_bond_order((*biter)->get_type());
                int atom_type_A = (*biter)->get_first_atom()->get_type();
                int atom_type_B = (*biter)->get_second_atom()->get_type();
                stretch_data_holder_.push_back(TAFFStretch::StretchData());
                stretch_data_holder_.back().atom1 = (*biter)->get_first_atom();
                stretch_data_holder_.back().atom2 = (*biter)->get_second_atom();
                stretch_data_holder_.back().bond_order = bond_order;
                if(!assign_params(value, atom_type_A, atom_type_B, bond_order))
					if(!assign_params(value, atom_type_B, atom_type_A, bond_order))
						if(!assign_params(value, ATOM::ANY_TYPE, atom_type_B, bond_order))
						{
#ifdef DEBUG
							cout<<"TAFFStretch::setup(): cannot find bond stretch parameters for "<<
							force_field_->get_parameters().get_type_name(atom_type_A)<<"-"<<
							force_field_->get_parameters().get_type_name(atom_type_B)<<" with bond type: "
							<<(*biter)->get_type()<<" and default parameters are assigned"<<endl;
#endif
							// if we cannot assign proper type for the bond, default parameters are assigned (form tripos)
							value.k = 600.0;
							value.r0 = 1.50;
							force_field_->add_unassigned_atom((*biter)->get_first_atom());
							force_field_->add_unassigned_atom((*biter)->get_second_atom());
							//continue;
                     }
                stretch_data_holder_.back().value = value;
     }
// everything goes well
     return true;
}
// update methods
double TAFFStretch::update_energy(){
       //energy initializion
       energy_ = 0.0;
       //debug
       //cout<<stretch_data_holder_.size()<<endl;
       // iterate all bonds and summarize the stretch energies
       for (vector<StretchData>::size_type i = 0; i<stretch_data_holder_.size(); ++i){
           double distance = stretch_data_holder_[i].atom1->get_position().dist(stretch_data_holder_[i].atom2->get_position());
           //debug
           //cout<<stretch_data_holder_[i].atom1->get_symbol_type()<<" "<<stretch_data_holder_[i].atom2->get_symbol_type()<<" "<<distance<<endl;
           energy_ += 0.5*stretch_data_holder_[i].value.k*pow((distance - stretch_data_holder_[i].value.r0),2);
       }
       return energy_;
}
// calculate current forces imposed by stretch and add them to the force field;
void TAFFStretch::update_forces(){
     if (get_force_field() == 0){
                           cout<<"TAFFStretch::update_force(): error! this component doesn't bond to any force field"<<endl;
                           return;
                           }    
     // iterate all bonds and update forces
     for (vector<StretchData>::size_type i = 0; i<stretch_data_holder_.size(); ++i){
         vector3 force_1(stretch_data_holder_[i].atom1->get_force()), force_2(stretch_data_holder_[i].atom2->get_force());
         vector3 direction(stretch_data_holder_[i].atom1->get_position()-stretch_data_holder_[i].atom2->get_position());
         double distance = direction.length();
		 direction.normalize();
         if (distance != 0)
		 {
                      // convert the units from kcal/mol A to N
                      // kcal -> J: 1e3 * 4.2
                      // A -> m: 1e-10
                      // J/mol -> J: Avogadro
			 double dE = stretch_data_holder_[i].value.k * (distance - stretch_data_holder_[i].value.r0) * FORCE_FACTOR;
			 direction *= dE;
			 force_1 -= direction;
			 force_2 += direction;
		 }
         stretch_data_holder_[i].atom1->set_force(force_1);
         stretch_data_holder_[i].atom2->set_force(force_2);
         }
}


                           



                                                            
                           
                           
