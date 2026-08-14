#include "../include/TAFFVDW.h"
#include "../include/TAFF.h"
#include "../include/utility.h"

using namespace std;

// constructors
TAFFVDW::TAFFVDW():
                           FFComponent(),
                           fast_access(),
                           num_of_atom_types_(0),
                           vdw_data_holder_(),
                           vdw_energy_(0.0){
                                            //set component name
                                            set_name("TAFF VDW");
                                            }
TAFFVDW::TAFFVDW(ForceField& ff):
                           FFComponent(ff),
                           fast_access(),
                           num_of_atom_types_(0),
                           vdw_data_holder_(),
                           vdw_energy_(0.0){
                                            //set component name
                                            set_name("TAFF VDW");
                                            }
//copy constructor
TAFFVDW::TAFFVDW(const TAFFVDW& to_copy):
                           FFComponent(to_copy){
                                                fast_access = to_copy.fast_access;
                                                num_of_atom_types_ = to_copy.num_of_atom_types_;
                                                vdw_data_holder_ = to_copy.vdw_data_holder_;
                                                vdw_energy_ = to_copy.vdw_energy_;
                                                }
// destructor
TAFFVDW::~TAFFVDW(){
                            //fast_access.clear();
                            vdw_data_holder_.clear();
                            }
// extract vdw parameters from the non bonded FFParamter object and build a hashtable for fast access
bool TAFFVDW::extract_VDW_parameters(FFParameter& ffp){
     if(!ffp.is_valid())
         return false;
     // 
     // build a 2 dim array of atom types and loop variables
     //FFParameter::AtomTypes atom_types = ffp.get_atomtypes();
     int num_types = ffp.get_num_types();
     int num_entries = num_types * num_types;
     num_of_atom_types_ = num_types;
     // clear internal data structures
     fast_access.Aij.clear();
     fast_access.Aij.resize(num_entries);
     fast_access.Bij.clear();
     fast_access.Bij.resize(num_entries);
     fast_access.is_defined.clear();
     fast_access.is_defined.resize(num_types);
     for (int i = 0; i < num_types; i++) 
         fast_access.is_defined[i] = false;
     // start pack the parameters into the vector fast_access;
     int type;
     string name_type;
     vector<double> epsilon(num_types);
     vector<double> r(num_types);
     //debug
     //cout<<get_force_field()->get_parameters().params_in_each_section["vdw"].size()<<endl;
     for (vector<vector<string> >::iterator it = force_field_->get_parameters().params_in_each_section["vdw"].begin();
         it != force_field_->get_parameters().params_in_each_section["vdw"].end(); ++it){
               name_type = (*it)[0];
               type = ffp.get_type(name_type);
               //debug
               //cout<<"type is "<<type<<endl;
               if (type == -1){
                   cout<<"TAFFVDW::extract_VDW_parameters(): error! no numeric atom type defined for atom "<<name_type<<endl;
                   return false;
                   }
               fast_access.is_defined[type-1] = true;
               epsilon[type-1] = str2double((*it)[2]);
               r[type-1] = str2double((*it)[1]);
         }

     // now assemble all vdw parameters for all known atom type pairs
     for (int i = 0; i< num_types; i++)
         for (int j = 0; j<num_types; j++){
               // calculate index for Aij and Bij fields
               int index = i +j * num_types;
               int sym_index = i * num_types + j;
               //debug
               //cout<<fast_access.is_defined[i]<<" "<<fast_access.is_defined[i]<<endl;
               if (fast_access.is_defined[i] && fast_access.is_defined[j]){
                                            // calculate Aij and Bij from epsilon and r values
                                            // average method is arithmatic average
                                            double R = r[i] + r[j];
                                            double R6 = pow(R,6);
                                            double R12 = pow(R,12);
                                            double Epsilon = sqrt(epsilon[i]*epsilon[j]);
                                            fast_access.Aij[index] = Epsilon*R12;
                                            fast_access.Bij[index] = 2.0*Epsilon*R6;
                                            fast_access.Aij[sym_index] = fast_access.Aij[index];
                                            fast_access.Bij[sym_index] = fast_access.Bij[index];
               }
               // ignore atom type *
               else if((i+1) == ATOM::ANY_TYPE || (j+1) == ATOM::ANY_TYPE){
                    fast_access.Aij[index] = 0.0;
                    fast_access.Bij[index] = 0.0;
                    fast_access.Aij[sym_index] = 0.0;
                    fast_access.Bij[sym_index] = 0.0;
                    }
               else{
                    cout<<"Warning:TAFFVDW::extract_VDW_parameters(): no Aij and Bij parameters has been defined for: "<<
                          ffp.get_type_name(i+1)<<" and "<<ffp.get_type_name(j+1)<<endl;
                    fast_access.Aij[index] = 0.0;
                    fast_access.Bij[index] = 0.0;
                    fast_access.Aij[sym_index] = 0.0;
                    fast_access.Bij[sym_index] = 0.0;
               }
         }
     return true;
}

// query a set of vdw parameters has defined for a given combination of atom types
bool TAFFVDW::has_params(int i, int j) const{
     if ((i > 0 && i <= num_of_atom_types_) && (j>0 && j <= num_of_atom_types_))
         return (fast_access.is_defined[i-1] && fast_access.is_defined[j-1]);
     else
         return false;
     }
// return the parameters for a given combination of atom types
TAFFVDW::VDWForceValues TAFFVDW::get_params(int i, int j) const{
                      TAFFVDW::VDWForceValues value;
                      if (assign_params(value, i, j))
                         return value;
                      else{
                          cout<<"TAFFVDW::get_params():Warning: no vdw parameters are assigned for atom "<<
                          i<<" and "<<j<<endl;
                          value.A = 0.0;
                          value.B = 0.0;
                          return value;
                          }
}
// assign the parameters for a given combination of atom types;
bool TAFFVDW::assign_params(TAFFVDW::VDWForceValues& param, int i, int j) const{

     if (has_params(i, j)){
                       param.A = fast_access.Aij[(i-1) + (j-1)*num_of_atom_types_];
                       param.B = fast_access.Bij[(i-1) + (j-1)*num_of_atom_types_];
                       return true;
                       }
     else
         return false;
}

// set up method
bool TAFFVDW::setup(){
     if (force_field_ == NULL){
                          cout<<"TAFFVDW::setup(): force field bound component can not be found"<<endl;
                          return false;
                          }
     // clear the vdw parameter container
     vdw_data_holder_.clear();
     // tempararily set this component enabled;
     setenabled(true);
     // extract the L-J vdw parameters
     TAFF* taff_force_field = dynamic_cast<TAFF*>(force_field_);
     if (taff_force_field == NULL || taff_force_field->has_initialized_param()){
                          bool result = extract_VDW_parameters(force_field_->get_parameters());
                          if (!result){
                                       cout << "TAFFVDW::setup(): can not access L-J vdw section"<<endl;
                                       return false;
                                       }
                          }
     // now iterate all atom pairs and remove 1-2, 1-3 pairs
     VDWForceValues value;
     bool is_14 = false;
     for (ATOMVec::iterator aiter1 = force_field_->get_atoms().begin();
          aiter1 != force_field_->get_atoms().end(); ++aiter1){
                 for (ATOMVec::iterator aiter2 = aiter1;
                     aiter2 != force_field_->get_atoms().end(); ++aiter2){
                            // first check if is the same atom
                            if (aiter1 == aiter2) continue;
                            // second check if is bonded
                            else if ((*aiter2)->is_bonded_to(*aiter1)) continue;
                            // then check if is geminal
                            else if ((*aiter2)->is_geminal_to(*aiter1)) continue;
                            else{
                                 // check if is vincinal and switch the mark;
                                 if ((*aiter2)->is_vicinal_to(*aiter1))
                                    is_14 = true;
                                 int atom_type_A = (*aiter1)->get_type();
                                 int atom_type_B = (*aiter2)->get_type();
                                 //debug
                                 //cout<<(*aiter1)->get_id()<<" "<<(*aiter2)->get_id()<<endl;
                                 vdw_data_holder_.push_back(TAFFVDW::VDWData());
                                 vdw_data_holder_.back().atom1 = *aiter1;
                                 vdw_data_holder_.back().atom2 = *aiter2;
                                 vdw_data_holder_.back().is_14_interaction = is_14;
                                 // check if the given atom pair has predefined vdw parameters
                                 if(!assign_params(value, atom_type_A, atom_type_B))
                                 {
                                      cout<<"TAFFVDW::setup(): cannot find L-J vdw parameters for "<<
                                      force_field_->get_parameters().get_type_name(atom_type_A)<<"-"<<
                                      force_field_->get_parameters().get_type_name(atom_type_B)<<endl;
                                      value.A = 0.0;
                                      value.B = 0.0;
                                      force_field_->add_unassigned_atom(*aiter1);
                                      force_field_->add_unassigned_atom(*aiter2);
                                      }
                                 vdw_data_holder_.back().value = value;
                            }
                     }
          }
          return true;
     }

double TAFFVDW::update_energy(){
       energy_ = 0.0;
       vdw_energy_ = 0.0;
       // temperarily use distance independent dielectric constant 
       double vdw_1_4_scale_factor = 1.0;
       // loop all non bond atom pairs and calculate vdw and ele energies
       for (vector<VDWData>::size_type i = 0; i != vdw_data_holder_.size();
           ++i){
                   double distance = vdw_data_holder_[i].atom1->get_position().dist(vdw_data_holder_[i].atom2->get_position());

                   double inver_dis = 1/distance;
                   double inver_dis_2 = inver_dis*inver_dis;
                   double inver_dis_3 = inver_dis_2*inver_dis;
                   double inver_dis_6 = inver_dis_3*inver_dis_3;
                   if (!vdw_data_holder_[i].is_14_interaction){
                      vdw_energy_ += inver_dis_6*(inver_dis_6*vdw_data_holder_[i].value.A - vdw_data_holder_[i].value.B);
                      }
                   else{
                   // ATTENTION!! the format of 1-4 special non bond interaction
                       
                       vdw_energy_ += vdw_1_4_scale_factor*inver_dis_6*(inver_dis_6*vdw_data_holder_[i].value.A - vdw_data_holder_[i].value.B);
                       }
                   
           }

       energy_ = vdw_energy_;
       //energy_ = ele_energy_;
       return energy_;            
       }
// get vdw energy
double TAFFVDW::get_vdw_energy()const{
       return vdw_energy_;
       }

// calculate the forces imposed on each atoms by vdw and ele interaction
void TAFFVDW::update_forces(){
     if (force_field_ == NULL){
                           cout<<"TAFFVDW::update_force(): error! this component doesn't bond to any force field"<<endl;
                           return;
                           }
     double vdw_1_4_scale_factor = 1.0;
     // iterate all non bond pairs and update forces
     for (vector<VDWData>::size_type i = 0; i != vdw_data_holder_.size();++i){
         vector3 force_1 = vdw_data_holder_[i].atom1->get_force(), force_2 = vdw_data_holder_[i].atom2->get_force();
         vector3 direction(vdw_data_holder_[i].atom1->get_position()-vdw_data_holder_[i].atom2->get_position());
         double distance = direction.length();
         double distance_2 = direction.length_2();
         //direction = direction.normalize();
         if (distance != 0.0){
                      double inverse_dist = 1/distance;
                      double inverse_dist_2 = inverse_dist*inverse_dist;
                      //define a convert factor of vdw interaction
                      // convert the units from kcal/mol A to N
                      // kcal -> J: 1e3 * 4.2
                      // A -> m: 1e-10
                      // J/mol -> J: Avogadro
                      double inverse_dist_6 = inverse_dist_2*inverse_dist_2*inverse_dist_2;
                      double inverse_dist_12 = inverse_dist_6*inverse_dist_6;
                      double vdw_factor;
                      if (!vdw_data_holder_[i].is_14_interaction)
                         vdw_factor = 4.184*(1e13/AVOGADRO)*inverse_dist_6*inverse_dist_2 * vdw_1_4_scale_factor * (12*vdw_data_holder_[i].value.A * inverse_dist_6 - 6*vdw_data_holder_[i].value.B);
                      else
                          vdw_factor = 4.184*(1e13/AVOGADRO)*inverse_dist_6 * inverse_dist_2 *( 12*vdw_data_holder_[i].value.A * inverse_dist_6 - 6*vdw_data_holder_[i].value.B);
                      //vdw_factor *= 4.2*1e13/AVOGADRO;
                      //direction.normalize();
                      force_1 += direction*vdw_factor;
                      force_2 -= direction*vdw_factor;
                      vdw_data_holder_[i].atom1->set_force(force_1);
                      vdw_data_holder_[i].atom2->set_force(force_2);
         }
     }
}


                      
void TAFFVDW::update(){;}                      


       
                                                                       
                                 
                                 
                            
                
