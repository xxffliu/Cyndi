#include "../include/TAFFBend.h"
#include "../include/TAFF.h"
#include "../include/Atom.h"
#include "../include/Bond.h"
//#include <cmath>
#include "../include/utility.h"

/*#ifndef AVOGADRO
#define AVOGADRO  6.0221367E+23L;
#endif
*/

//using namespace std;
//clear function of BendHashData
void TAFFBend::BendHashData::clear(){
     k.clear();
     theta0.clear();
     is_defined.clear();
     }
     
// consturctors
TAFFBend::TAFFBend():
                     FFComponent(),
                     fast_access(),
                     num_of_atom_types_(0),
                     bend_data_holder_(){
                                         // set component name
                                         set_name("TAFF Bend");
                                         }
TAFFBend::TAFFBend(ForceField& ff):
                               FFComponent(ff),
                               fast_access(),
                               num_of_atom_types_(0),
                               bend_data_holder_(){
                                                   set_name("TAFF Bend");
                                                   }
// copy constructor
TAFFBend::TAFFBend(const TAFFBend& to_copy):
                         FFComponent(to_copy){
                                              fast_access = to_copy.fast_access;
                                              bend_data_holder_ = to_copy.bend_data_holder_;
                                              num_of_atom_types_ = to_copy.num_of_atom_types_;
                                              }
// destructor
TAFFBend::~TAFFBend(){
                      fast_access.clear();
                      bend_data_holder_.clear();
                      }

// extract angle bend parameters from the bonded FFParamter object and build a hashtable for fast access
bool TAFFBend::extract_AB_parameters(FFParameter& ffp){
     if(!ffp.is_valid()){
         //debug
         cout<<"force field parameters is not valid"<<endl;
         return false;
         }
     // build a 2 dim array of atom types and loop variables
     //FFParameter::AtomTypes atom_types = ffp.get_atomtypes();
     int num_types = ffp.get_num_types();
     int num_entries = num_types * num_types * num_types;
     num_of_atom_types_ = num_types;
     // clear internal data structures
     /*if (fast_access.k)
            delete [] fast_access.k;
     if (fast_access.theta0)
            delete [] fast_access.theta0;
     if (fast_access.is_defined)
            delete [] fast_access.is_defined;
     fast_access.k  = new double[num_entries];
     fast_access.theta0 = new double[num_entries];
     fast_access.is_defined = new bool[num_entries];
     */
     fast_access.clear();
     fast_access.k.resize(num_entries);
     fast_access.theta0.resize(num_entries);
     fast_access.is_defined.resize(num_entries);
     
     for (int i = 0; i < num_entries; i++) 
         fast_access.is_defined[i] = false;
     //debug
     //cout<<num_types<<endl;
     // start pack the parameters into the vector fast_access
     int type1, type2, type3,index;
     string name_type1, name_type2, name_type3;
     for (vector<vector<string> >::iterator it = force_field_->get_parameters().params_in_each_section["angle"].begin();
         it != force_field_->get_parameters().params_in_each_section["angle"].end(); ++it){
               name_type1 = (*it)[0];
               name_type2 = (*it)[1];
               name_type3 = (*it)[2];
               type1 = ffp.get_type(name_type1);
               type2 = ffp.get_type(name_type2);
               type3 = ffp.get_type(name_type3);
               //debug
               //cout<<type1<<" "<<type2<<" "<<type3<<" "<<endl;
               if ((type1 == -1) || (type2 == -1) || (type3 == -1)){
                   cout<<"error! no numeric atom type defined for atom "<<name_type1<<" or "<<name_type2<<" or "<<name_type3<<endl;
                   return false;
                   }
               index = ((type1-1) + (type2-1) * num_types + (type3-1)*num_types*num_types);
               fast_access.is_defined[index] = true;
               fast_access.theta0[index] = str2double((*it)[3]);
               fast_access.k[index] = str2double((*it)[4]);
               //build the symmetric index
               index = ((type3-1) + (type2-1) * num_types + (type1-1)*num_types*num_types);
               fast_access.is_defined[index] = true;
               fast_access.theta0[index] = str2double((*it)[3]);
               fast_access.k[index] = str2double((*it)[4]);
         }
     return true;
}

// query a set of angle bend parameters has defined for a given combination of atom types
bool TAFFBend::has_params(int i, int j, int k) const{
     if ((i > 0 && i <= num_of_atom_types_) && (j>0 && j <= num_of_atom_types_) && (k>0 && k <= num_of_atom_types_))
         return fast_access.is_defined[(i-1) + (j-1) * num_of_atom_types_ + (k-1) * num_of_atom_types_ * num_of_atom_types_];
     else
         return false;
     }
// return the parameters for a given combination of atom types
TAFFBend::ForceValues TAFFBend::get_params(int i, int j, int k) const{
                      TAFFBend::ForceValues value;
                      if (assign_params(value, i, j, k))
                         return value;
                      else{
                          cout<<"TAFFBend::get_params():Warning: no angle bend parameters are assigned for atom "<<
                          i<<" and "<<j<<" and "<<k<<endl;
                          value.k = 0.0;
                          value.theta0 = 0.0;
                          return value;
                          }
}
// assign the parameters for a given combination of atom types;
bool TAFFBend::assign_params(TAFFBend::ForceValues& param, int i, int j, int k) const{
     if (has_params(i, j, k)){
                       //debug
                       //cout<<i<<" "<<j<<" "<<k<<" is defined"<<endl;
                       param.k = fast_access.k[(i-1) + (j-1)*num_of_atom_types_ + (k-1)*num_of_atom_types_*num_of_atom_types_];
                       param.theta0 = fast_access.theta0[(i-1) + (j-1)*num_of_atom_types_ + (k-1)*num_of_atom_types_*num_of_atom_types_];
                       return true;
                       }
     else
         return false;
}

// set up method
bool TAFFBend::setup(){
     if (get_force_field() == 0){
                          cout<<"TAFFBend::setup(): force field bound component can not be found"<<endl;
                          return false;
                          }
     // clear the parameter container
     bend_data_holder_.clear();
     // tempararily set this component enabled;
     setenabled(true);
     
     TAFF* taff_force_field = dynamic_cast<TAFF*>(force_field_);
     if (taff_force_field == 0 || taff_force_field->has_initialized_param()){
                          bool result = extract_AB_parameters(get_force_field()->get_parameters());
                          if (!result){
                                       cout << "TAFFBend::setup(): can not found angle bend section"<<endl;
                                       return false;
                                       }
                          }
     TAFFBend::ForceValues value;
     // retrieve all bend parameters
     BONDVec::iterator biter, biter2;
     for (ATOMVec::iterator aiter = force_field_->get_atoms().begin();
          aiter != force_field_->get_atoms().end(); ++aiter){
                if((*aiter)->get_symbol_type() == "H")
                    continue;
                //cout<<(*aiter)->get_id()<<" "<<*aiter<<endl;
                for (biter = (*aiter)->get_bond_list().begin(); biter != (*aiter)->get_bond_list().end(); ++biter){
                    for (biter2 = biter, ++biter2; biter2 != (*aiter)->get_bond_list().end(); ++biter2){
                        //cout<<(*biter)->get_first_atom()<<" "<<(*biter)->get_second_atom()<<endl;
                        int atom_type_A = (*biter)->get_partner(*aiter)->get_type();
                        int atom_type_B = (*aiter)->get_type();
                        int atom_type_C = (*biter2)->get_partner(*aiter)->get_type();
                        bend_data_holder_.push_back(TAFFBend::BendData());
                        bend_data_holder_.back().atom1 = (*biter)->get_partner(*aiter);
                        bend_data_holder_.back().atom2 = *aiter;
                        bend_data_holder_.back().atom3 = (*biter2)->get_partner(*aiter);
                        if(!assign_params(value, atom_type_A, atom_type_B, atom_type_C))
							if(!assign_params(value, ATOM::ANY_TYPE, atom_type_B, atom_type_C))
								if(!assign_params(value, atom_type_A, atom_type_B, ATOM::ANY_TYPE))
									if(!assign_params(value, ATOM::ANY_TYPE, atom_type_B, ATOM::ANY_TYPE))
									{
                            //cout<<ATOM::ANY_TYPE<<" "<<atom_type_B<<" "<<ATOM::ANY_TYPE<<" is not defined"<<endl;
#ifdef DEBUG
										cout<<"TAFFBend::setup(): cannot find angle bend parameters for "<<
										force_field_->get_parameters().get_type_name(atom_type_A)<<"-"<<
										force_field_->get_parameters().get_type_name(atom_type_B)<<"-"<<
										force_field_->get_parameters().get_type_name(atom_type_C)<<" and default parameterss are assigned"<<endl;
#endif
										value.k = 0.024;
										value.theta0 = 109.5;
										get_force_field()->add_unassigned_atom((*biter)->get_partner(*aiter));
										get_force_field()->add_unassigned_atom((*biter2)->get_partner(*aiter));
										get_force_field()->add_unassigned_atom(*aiter);
										//continue;
									}
                        bend_data_holder_.back().value = value;
                    }
                }
          }
// everything goes well
     return true;
}

// update methods
double TAFFBend::update_energy(){
       //energy initializion
       energy_ = 0.0;
       if (bend_data_holder_.size() == 0)
           return 0.0;
       // iterate all angles and summarize the stretch energies
       for (vector<BendData>::size_type i = 0; i<bend_data_holder_.size(); ++i){
           //cout<<bend_data_holder_.size()<<endl;
           vector3 v1, v2;
           v1 = bend_data_holder_[i].atom1->get_position() - bend_data_holder_[i].atom2->get_position();
           v2 = bend_data_holder_[i].atom3->get_position() - bend_data_holder_[i].atom2->get_position();
           
           v1 = v1.normalize();
           v2 = v2.normalize();
           if((v1.length()<1e-3) && (v2.length()<1e-3)){
                               cout<<v1.length()<<" "<<v2.length()<<endl;
                               exit(1);
                               }
                    double costheta = v1*v2;
				double theta;
				if (costheta > 1.0) 
					costheta = 1.0;
				else if (costheta < -1.0) 
					costheta = -1.0;
         theta = acos(costheta)*180.0/PI;
           double angle = vectorAngle(v1, v2);
           energy_ += 0.5 * bend_data_holder_[i].value.k * pow((theta - bend_data_holder_[i].value.theta0), 2);
       }
       return energy_;
}

// calculate current forces imposed by bend and add them to the force field;
void TAFFBend::update_forces(){
     if (get_force_field() == 0){
                           cout<<"TAFFBend::update_force(): error! this component doesn't bond to any force field"<<endl;
                           return;
                           }
     // iterate all angles and update forces
     for (vector<BendData>::size_type i = 0; i < bend_data_holder_.size(); ++i){
         vector3 force_1 = bend_data_holder_[i].atom1->get_force(), force_2 = bend_data_holder_[i].atom2->get_force(), force_3 = bend_data_holder_[i].atom3->get_force();
         // calculate the vectors between atom1 and atom2, atom3 and atom2 then normalize them
         vector3 v1 = bend_data_holder_[i].atom1->get_position() - bend_data_holder_[i].atom2->get_position();
         vector3 v2 = bend_data_holder_[i].atom3->get_position() - bend_data_holder_[i].atom2->get_position();
         double length = v1.length();
         if (length == 0.0) continue;
         double inverse_v1 = 1/length;
         v1 = v1.normalize();
         length = v2.length();
         if (length == 0.0) continue;
         double inverse_v2 = 1/length;
         v2 = v2.normalize();
         double angle = vectorAngle(v1, v2);
         double costheta = v1 * v2;
				double theta;
				if (costheta > 1.0) 
					costheta = 1.0;
				else if (costheta < -1.0) 
					costheta = -1.0;
         theta = acos(costheta)*180.0/PI;
         // calculate the cross product of v1 and v2 and then normalize it
         //vector3 cross_product = cross(v1, v2);
         //cross_product = cross_product.normalize();
         // convert the units from kcal/mol A to N
         // kcal -> J: 1e3 * 4.2
         // A -> m: 1e-10
         // J/mol -> J: Avogadro
         double factor = 4.184*1e13 / AVOGADRO * 180.0/PI * bend_data_holder_[i].value.k * (theta - bend_data_holder_[i].value.theta0);
         vector3 n1 = factor*(-1.0/(sin(acos(costheta))+0.00001)*inverse_v1*(v2-v1*costheta));
         vector3 n2 = factor*(-1.0/(sin(acos(costheta))+0.00001)*inverse_v2*(v1-v2*costheta));
         //vector3 n1 = factor*(v1%cross_product);
         //vector3 n2 = factor*(v2%cross_product);
         force_1 -= n1;
         force_2 += n1;
         force_2 += n2;
         force_3 -= n2;
         bend_data_holder_[i].atom1->set_force(force_1);
         bend_data_holder_[i].atom2->set_force(force_2);
         bend_data_holder_[i].atom3->set_force(force_3);
         }
}
                      

     
     
