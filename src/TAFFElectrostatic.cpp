#include "../include/TAFFElectrostatic.h"
#include "../include/TAFF.h"
#include "../include/utility.h"

// constructors
TAFFEle::TAFFEle():
                           FFComponent(),
						   ele_pair_list_(),
                           ele_energy_(0.0){
                                            //set component name
                                            set_name("TAFF ELE");
                                            }
TAFFEle::TAFFEle(ForceField& ff):
                           FFComponent(ff),
							   ele_pair_list_(),
                           ele_energy_(0.0){
                                            //set component name
                                            set_name("TAFF ELE");
                                            }
//copy constructor
TAFFEle::TAFFEle(const TAFFEle& to_copy):
                           FFComponent(to_copy),
							   ele_pair_list_(to_copy.ele_pair_list_),
						ele_energy_(to_copy.ele_energy_)                                             
						   {}
// destructor
TAFFEle::~TAFFEle(){
ele_energy_ = 0.0;
ele_pair_list_.clear();
}

// set up method
bool TAFFEle::setup()
{
     if (force_field_ == NULL)
	 {
                          cout<<"TAFFEle::setup(): force field bound component can not be found"<<endl;
                          return false;
	 }
	 ele_pair_list_.clear();
     // tempararily set this component enabled;
     setenabled(true);
     // now iterate all atom pairs and remove 1-2, 1-3 pairs
     bool is_14 = false;
     for (ATOMVec::iterator aiter1 = force_field_->get_atoms().begin();
          aiter1 != force_field_->get_atoms().end(); ++aiter1)
	 {
                 for (ATOMVec::iterator aiter2 = aiter1; aiter2 != force_field_->get_atoms().end(); ++aiter2)
				 {
                            // first check if is the same atom
                            if (aiter1 == aiter2) continue;
                            // second check if is bonded
                            else if ((*aiter2)->is_bonded_to(*aiter1)) continue;
                            // then check if is geminal
                            else if ((*aiter2)->is_geminal_to(*aiter1)) continue;
                            else
							{
                                 // check if is vincinal and switch the mark;
                                 if ((*aiter2)->is_vicinal_to(*aiter1))
                                    is_14 = true;
                                 //debug
								 ele_pair_list_.push_back(TAFFEle::ElePair());
                                 ele_pair_list_.back().atom1 = *aiter1;
                                 ele_pair_list_.back().atom2 = *aiter2;
                                 ele_pair_list_.back().is_14_interaction = is_14;
                            }
                     }
          }
          return true;
     }

double TAFFEle::update_energy(){
       energy_ = 0.0;
       ele_energy_ = 0.0;
       // temperarily use distance independent dielectric constant 
       double ele_1_4_scale_factor = 1.0;
       // loop all non bond atom pairs and calculate vdw and ele energies
       for (vector<ElePair>::size_type i = 0; i != ele_pair_list_.size();
           ++i){
                   double distance = ele_pair_list_[i].atom1->get_position().dist(ele_pair_list_[i].atom2->get_position());
                   double inver_dis = 1/distance;
				   // Caution!!!!!!!!
                   // temerarilly using constant dielectric constant;
                   if (!ele_pair_list_[i].is_14_interaction){
                      ele_energy_ += 332.17*inver_dis*(ele_pair_list_[i].atom1->get_charge()*ele_pair_list_[i].atom2->get_charge());
                      }
                   else{
                   // ATTENTION!! the format of 1-4 special non bond interaction
                       
                       ele_energy_ += 332.17*ele_1_4_scale_factor*inver_dis*(ele_pair_list_[i].atom1->get_charge()*ele_pair_list_[i].atom2->get_charge());
              
                       }
                   
           }

       energy_ = ele_energy_;
       return energy_;            
       }

// get electrostatistic energy
double TAFFEle::get_ele_energy()const{
       return ele_energy_;
       }

// calculate the forces imposed on each atoms by vdw and ele interaction
void TAFFEle::update_forces(){
     if (force_field_ == NULL){
                           cout<<"TAFFEle::update_force(): error! this component doesn't bond to any force field"<<endl;
                           return;
                           }
     double ele_1_4_scale_factor = 0.75;
     // iterate all non bond pairs and update forces
     for (vector<ElePair>::size_type i = 0; i != ele_pair_list_.size();++i){
         vector3 force_1 = ele_pair_list_[i].atom1->get_force(), force_2 = ele_pair_list_[i].atom2->get_force();
         vector3 direction(ele_pair_list_[i].atom1->get_position()-ele_pair_list_[i].atom2->get_position());
         double distance = direction.length();
         //double distance_2 = direction.length_2();
         direction.normalize();
         if (distance != 0.0){
                      double inverse_dist = 1/distance;
                      double inverse_dist_2 = inverse_dist*inverse_dist;
                      // first define a convert factor of electrostatistic interaction
                      // convert the units from kcal/mol A to N
                      // kcal -> J: 1e3 * 4.2
                      // A -> m: 1e-10
                      // J/mol -> J: Avogadro
                      // temperarily use distance independent dielectric constant
                      double ele_factor;
                      if (!ele_pair_list_[i].is_14_interaction)
                         ele_factor = 332.17 * FORCE_FACTOR; 
                      else
                         ele_factor = 332.17 * FORCE_FACTOR * ele_1_4_scale_factor;
                      ele_factor *= (ele_pair_list_[i].atom1->get_charge()*ele_pair_list_[i].atom2->get_charge()*inverse_dist_2);
                      force_1 += direction * ele_factor;
                      force_2 -= direction * ele_factor;

                      ele_pair_list_[i].atom1->set_force(force_1);
                      ele_pair_list_[i].atom2->set_force(force_2);
         }
     }
}
void TAFFEle::update(){;} 