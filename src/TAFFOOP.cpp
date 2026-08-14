#include "../include/TAFF.h"
#include "../include/TAFFOOP.h"
#include "../include/Vector3.h"
#include "../include/utility.h"
using namespace std;
class vector3;
// constructors
TAFFOOP::TAFFOOP():
                   FFComponent(),
                   fast_access(),
                   oop_data_holder_(NUM_OOP_TYPES){
                                      set_name("TAFF OOP");
                                      }
TAFFOOP::TAFFOOP(ForceField& ff):
                             FFComponent(ff),
                             fast_access(),
                             oop_data_holder_(NUM_OOP_TYPES){
                                                set_name("TAFF OOP");
                                                }

// copy constructor
TAFFOOP::TAFFOOP(const TAFFOOP& to_copy):
                       FFComponent(to_copy){
                                            fast_access = to_copy.fast_access;
                                            oop_data_holder_ = to_copy.oop_data_holder_;
                                            }
// destructor
TAFFOOP::~TAFFOOP(){
                    fast_access.k0.clear();
                    fast_access.is_defined.clear();
                    oop_data_holder_.clear();
                    }

// extract angle oop ff parameters from the TAFF parameter file
// and build some data structure for fast access these data 
bool TAFFOOP::extract_OOP_parameters(FFParameter& ffp){
     // check the force field is valid
     if (!ffp.is_valid()){
         //debug
         cout<<"TAFFOOP::extract_OOP_parameters(): ffp is not valid"<<endl;
         return false;
         }
     // build a one dim array of atom types
     //FFParameter::AtomTypes atom_types = ffp.get_atomtypes();
     int num_types = ffp.get_num_types();
     num_of_atom_types_ = num_types;
     // clear internal data structures
     /*if (fast_access.k)
            delete [] fast_acess.k;
     if (fast_access.is_defined)
            delete [] fast_access.is_defined;
     fast_access.k  = new double[num_types];
     fast_access.is_defined = new bool[num_types];*/
     fast_access.k0.clear();
     fast_access.k0.resize(num_types);
     fast_access.is_defined.clear();
     fast_access.is_defined.resize(num_types);
     for (int i = 0; i < num_types; i++) 
         fast_access.is_defined[i] = false;
     // start to pack the parameters into the vector of fast access
     int root_type;
     string name_type;
     for (vector<vector<string> >::iterator it = force_field_->get_parameters().params_in_each_section["oopbend"].begin();
         it != force_field_->get_parameters().params_in_each_section["oopbend"].end(); ++it){
               name_type = (*it)[0];
               root_type = ffp.get_type(name_type);
               if ((root_type == -1)){
                   cout<<"error! no numeric atom type defined for atom "<<name_type<<endl;
                   return false;
                   }
               int index = root_type;
               fast_access.is_defined[index] = true;
               fast_access.k0[index] = str2double((*it)[1]);

         }
     return true;
}

// query a root atom parameters has defined for a given root atom type
bool TAFFOOP::has_params(int i) const{
     if (i < 0 || i > num_of_atom_types_)
         return false;
     else
         return fast_access.is_defined[i];
     }
// return the parameters for a given combination of atom types
double TAFFOOP::get_params(int i) const{
                      double value;
                      if (assign_params(value, i))
                         return value;
                      else{
                          cout<<"TAFFOOP::get_params():Warning: no OOP parameters are assigned for atom "<<
                             i<<endl;
                          return 0.0;
                          }
}
// assign the parameters for a given combination of atom types;
bool TAFFOOP::assign_params(double& param, int i) const{
     if (has_params(i)){
                       param = fast_access.k0[i];
                       return true;
                       }
     else
         return false;
}
// set up method
bool TAFFOOP::setup(){
     if (force_field_ == NULL){
                          cout<<"TAFFOOP::setup(): force field bound component can not be found"<<endl;
                          return false;
                          }
     // clear the parameter holder
     oop_data_holder_.clear();
     // tempararily set this component enabled;
     setenabled(true);
     
     TAFF* taff_force_field = dynamic_cast<TAFF*>(force_field_);
     if (taff_force_field == NULL || taff_force_field->has_initialized_param()){
                          bool result = extract_OOP_parameters(get_force_field()->get_parameters());

                          if (!result){
                                       cout << "TAFFOOP::setup(): can not found OOP section"<<endl;
                                       return false;
                                       }
                          }
                          
     double value;
     // retrieve all stretch parameters
     for (ATOMVec::iterator aiter = force_field_->get_atoms().begin();
          aiter != force_field_->get_atoms().end(); ++aiter){
                string type = (*aiter)->get_symbol_type();
                int atom_type_root = (*aiter)->get_type();
                //check if the root atom is defined in the parameters list
                if ((type != "C.2") && (type != "C.cat") && (type != "N.2") && (type != "N.am") && (type != "N.ar") && (type != "N.pl3")&&(type != "C.ar")){
                   continue;
                   }
                // check if the root atom has 3 neighbor atoms

                if ((*aiter)->get_num_neighbor_bond() != 3){

                   continue;
                   }
                oop_data_holder_.push_back(TAFFOOP::OOPData());
                BONDVec::iterator biter = (*aiter)->get_bond_list().begin();
                int atom_type_A = (*biter)->get_partner(*aiter)->get_type();
                oop_data_holder_.back().atom1 = (*biter)->get_partner(*aiter);
                biter++;
                int atom_type_B = (*biter)->get_partner(*aiter)->get_type();
                oop_data_holder_.back().atom2 = (*biter)->get_partner(*aiter);
                biter++;
                int atom_type_C = (*biter)->get_partner(*aiter)->get_type();
                oop_data_holder_.back().atom3 = (*biter)->get_partner(*aiter);
                oop_data_holder_.back().root = *aiter;
                if (!assign_params(value, atom_type_root))
                {
#ifdef DEBUG
                     cout<<"TAFFOOP::setup(): cannot find oop parameters for "<<
                     force_field_->get_parameters().get_type_name(atom_type_root)<<endl;
#endif
                     value = 0.0;
                     force_field_->add_unassigned_atom(*aiter);
                     }
                oop_data_holder_.back().k = value;
     }
// everything goes well
     return true;
}
// update methods
double TAFFOOP::update_energy(){
       //energy initializion
       energy_ = 0.0;
       // iterate all OOP and summarize the oop energies

       for (vector<OOPData>::size_type i = 0; i<oop_data_holder_.size(); ++i){
           double distance = Point2Plane(oop_data_holder_[i].root->get_position(),oop_data_holder_[i].atom1->get_position(),
                                                  oop_data_holder_[i].atom2->get_position(),oop_data_holder_[i].atom3->get_position());

           energy_ += oop_data_holder_[i].k*distance*distance/2;
       }
       return energy_;
}
// calculate current forces imposed by stretch and add them to the force field;
void TAFFOOP::update_forces(){
     if (get_force_field() == 0){
                           cout<<"TAFFOOP::update_force(): error! this component doesn't bond to any force field"<<endl;
                           return;
                           }
     // iterate all OOP and update forces
     for (vector<OOPData>::size_type i = 0; i<oop_data_holder_.size(); ++i){
         vector3 force_1(oop_data_holder_[i].root->get_force()), force_2(oop_data_holder_[i].atom1->get_force()), force_3(oop_data_holder_[i].atom2->get_force()), force_4(oop_data_holder_[i].atom3->get_force());
         // root is the central atom, define 3 vectors from 3 bonded atoms
         /*vector3 ad = oop_data_holder_[i].atom1->get_position() - oop_data_holder_[i].root->get_position();
         vector3 bd = oop_data_holder_[i].atom2->get_position() - oop_data_holder_[i].root->get_position();
         vector3 cd = oop_data_holder_[i].atom3->get_position() - oop_data_holder_[i].root->get_position();
         */
         vector3 v01 = oop_data_holder_[i].root->get_position() - oop_data_holder_[i].atom1->get_position();
         vector3 v21 = oop_data_holder_[i].atom2->get_position() - oop_data_holder_[i].atom1->get_position();
         vector3 v31 = oop_data_holder_[i].atom3->get_position() - oop_data_holder_[i].atom1->get_position();
         /*double length = v1.length();
         if (length == 0.0) continue;
         double inverse_v1 = 1/length;
         v1 = v1.normalize();
         length = v2.length();
         if (length == 0.0) continue;
         double inverse_v2 = 1/length;
         v2 = v2.normalize();
         length = v3.length();
         if (length == 0.0) continue;
         double inverse_v3 = 1/length;
         v3 = v3.normalize();*/
         
         vector3 t = cross(v21,v31);
         double length_t2 = t.length_2();
         double v_dot = dot(t, v01);
         if (v_dot > 0)
             t = -t;
         if (length_t2 == 0)
             continue;
         t = t.normalize();
         // calculate the distance between root atom and plane defined by other 3 atoms
         double distance = Point2Plane(oop_data_holder_[i].root->get_position(),oop_data_holder_[i].atom1->get_position(),
                                                  oop_data_holder_[i].atom2->get_position(),oop_data_holder_[i].atom3->get_position());
         /*double costheta_1 = distance*inverse_v1;
         double costheta_2 = distance*inverse_v2;
         double costheta_3 = distance*inverse_v3;*/
         

         // convert the units from kcal/(mol*A) to N
         // kcal -> J: 1e3 * 4.184
         // A -> m: 1e-10
         // J/mol -> J: Avogadro
         double factor = 4.184*1e13/AVOGADRO*oop_data_holder_[i].k*distance;
         //double factor = 4.184*1e13/AVOGADRO*oop_data_holder_[i].k*ad.length()*v_dot;
         //debug
         vector3 v23 = oop_data_holder_[i].atom2->get_position() - oop_data_holder_[i].atom3->get_position();
         vector3 v12 = oop_data_holder_[i].atom1->get_position() - oop_data_holder_[i].atom2->get_position();
         //vector3 td = ad - t*(v_dot/length_t2);
         vector3& v0 = oop_data_holder_[i].root->get_position();
         vector3& v1 = oop_data_holder_[i].atom1->get_position();
         vector3& v2 = oop_data_holder_[i].atom2->get_position();
         vector3& v3 = oop_data_holder_[i].atom3->get_position();
         double t1, t2, t3;
         t1 = factor*MinDist_lines(v0, t, v2, v23)/MinDist_lines(v1, -t, v2, v23);
         t2 = factor*MinDist_lines(v0, t, v3, v31)/MinDist_lines(v2, -t, v3, v31);
         t3 = factor*MinDist_lines(v0, t, v1, v12)/MinDist_lines(v3, -t, v1, v12);

         force_1 += factor*t;
         force_2 -= t1*t;
         force_3 -= t2*t;
         force_4 -= t3*t;

         oop_data_holder_[i].root->set_force(force_1);
         oop_data_holder_[i].atom1->set_force(force_2);
         oop_data_holder_[i].atom2->set_force(force_3);
         oop_data_holder_[i].atom3->set_force(force_4);
         }
}
