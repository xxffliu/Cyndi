#include "../include/Atom.h"
//#include "Bond.h"

ATOM::ATOM():id(0),name("X"),symbol_type_("DU"),element_("X"),type(-1),subst_id(0),subst_name("UNK"),charge(0.0),
			 radius_(0.), weight_(0.), num_of_bond_(0),num_of_rotor_bond_(0),_p(), AromaticRingAtom_(false),mmff94_type(0),mmff94_symbol_type_("DU"), mmff94_partial_charge(0.), orig_position_(),_f(),is_ring(false),ring_id(),is_root_node(false),rotation_mark_(false){}

ATOM::ATOM(const ATOM& atom0):
           id(atom0.id),
           name(atom0.name),
           symbol_type_(atom0.symbol_type_),
           type(atom0.type),
		   mmff94_type(atom0.mmff94_type),
		   mmff94_symbol_type_(atom0.mmff94_symbol_type_),
		   mmff94_partial_charge(atom0.mmff94_partial_charge),
           element_(atom0.element_),
		   radius_(atom0.radius_),
		   weight_(atom0.weight_),
           subst_id(atom0.subst_id),
           subst_name(atom0.subst_name),
           charge(atom0.charge),
           num_of_bond_(atom0.num_of_bond_),
		   num_of_rotor_bond_(atom0.num_of_rotor_bond_),
		   num_of_bonded_atom_(atom0.num_of_bonded_atom_),
		   num_of_bonded_heavy_atom_(atom0.num_of_bonded_heavy_atom_),
           _p(atom0._p),
           _f(atom0._f),
           orig_position_(atom0.orig_position_),
           bond_(atom0.bond_),
           is_ring(atom0.is_ring),
           ring_id(atom0.ring_id),
           is_root_node(atom0.is_root_node),
		   AromaticRingAtom_(atom0.AromaticRingAtom_),
           rotation_mark_(atom0.rotation_mark_)
           {}

ATOM& ATOM::operator=(const ATOM& atom0){
      if (this == &atom0)
          return *this;
      id = atom0.id;
      name = atom0.name;
      symbol_type_ = atom0.symbol_type_;
      type = atom0.type;
	  mmff94_type = atom0.mmff94_type;
	  mmff94_symbol_type_ = atom0.mmff94_symbol_type_;
      element_ = atom0.element_;
	  radius_ = atom0.radius_;
	  weight_ = atom0.weight_;
      subst_id = atom0.subst_id;
      subst_name = atom0.subst_name;
      charge = atom0.charge;
	  mmff94_partial_charge = atom0.mmff94_partial_charge;
      num_of_bond_ = (atom0.num_of_bond_);
	  num_of_rotor_bond_ = (atom0.num_of_rotor_bond_);
	  num_of_bonded_atom_ = (atom0.num_of_bonded_atom_);
	  num_of_bonded_heavy_atom_ = (atom0.num_of_bonded_heavy_atom_);
      bond_ = (atom0.bond_);
      _p = atom0._p;
      _f = atom0._f;
      orig_position_ = atom0.orig_position_;
      is_ring = atom0.is_ring;
      ring_id = atom0.ring_id;
      is_root_node = atom0.is_root_node;
      rotation_mark_ = atom0.rotation_mark_;
	  AromaticRingAtom_=atom0.AromaticRingAtom_;
      return *this;
      }
      
void ATOM::clear(){
              id = 0;
              name = "X";
              symbol_type_ = "DU";
              type = -1;
			  mmff94_type = 0;
			  mmff94_symbol_type_ = "DU";
              element_ = "X";
			  radius_ = 0.;
			  weight_  = 0.;
              subst_id = 0;
              subst_name = "UNK";
              charge = 0.0;
			  mmff94_partial_charge = 0.0;
              num_of_bond_ = 0;
			  num_of_rotor_bond_ = 0;
			  num_of_bonded_atom_ = 0;
			  num_of_bonded_heavy_atom_ = 0;
              bond_.clear();
              _p = VZero;
              _f = VZero;
              orig_position_ = VZero;
              is_ring = false;
              ring_id.clear();
              is_root_node = false;
              rotation_mark_ = false;
			  AromaticRingAtom_=false;
              }

ATOM::~ATOM(){
              clear();
              }

int ATOM::get_atomic_num() const
{
	if(element_=="H")
		return 1;
	else if(element_=="C")
		return 6;
	else if(element_=="O")
		return 8;
	else if(element_=="N")
		return 7;
	else if(element_=="P")
		return 15;
	else if(element_=="S")
		return 16;
	else if(element_=="F")
		return 9;
	else if(element_=="Cl")
		return 17;
	else if(element_=="Br")
		return  35;
	else if(element_=="I")
		return 53;
	else if(element_ == "Si")
		return 14;
	else return 0;
	
}
int ATOM::get_num_PI_electrons()
{
	if(is_carbon())
	{
		if(symbol_type_ == "C.ar" || symbol_type_ == "C.2")
			return 1;
	}
	else if(is_nitrogen())
	{
		if(symbol_type_ == "N.ar" || symbol_type_ == "N.2")
		// fixed by xfliu. 20090227
		{
			if(num_of_bonded_heavy_atom_ == 2)
				return 1;
			else if(num_of_bonded_heavy_atom_ == 3)
			{
				//fixed by xfliu, 20090227
				return 1;
			}
		}
		else if(symbol_type_ == "N.pl3")
			return 2;
		else if(symbol_type_ == "N.3" && num_of_bond_ == 3)
			return 2;
	}
	else if(element_ == "O" || element_ == "S")
		return 2;
	return 0;
}
void ATOM::add_neighbor_bond_list(BOND* pbond){
     bond_.push_back(pbond);
     if (bond_.size() > MAX_NUM_BOND){
         cout<<"ATOM::add_neighbor_bond_list(): the number of bonds exceeds the maximum connecting number of this atom type"<<endl;
         return ;
         }
     }

void ATOM::add_neighbor_atom_list(ATOM* patom){
     atom_.push_back(patom);
     if (atom_.size() > MAX_NUM_BOND){
         cout<<"Warning: ATOM::add_neighbor_atom_list(): the number of bonded atoms exceeds the maximum connecting number of this atom type"<<endl;
         return ;
         }
     }
void ATOM::clear_neighbor_bond_list(){
     bond_.clear();
     }
     
void ATOM::clear_neighbor_atom_list(){
     atom_.clear();
     }
     
void ATOM::set_num_neighbor_bond(int i){
     num_of_bond_ = i;
     }
void ATOM::set_num_neighbor_atom(int i){
     num_of_bonded_atom_ = i;
     }
void ATOM::set_num_neighbor_heavy_atom(int i){
     num_of_bonded_heavy_atom_ = i;
     }
bool ATOM::is_bonded_to(ATOM* atom){
     if (bond_.empty()){
        //debug
        //cout<<"bond is empty"<<endl;
        return false;
        }
     if (this == atom)
        return false;
     for (BONDVec::iterator biter = bond_.begin(); biter != bond_.end(); ++biter){
         if ((*biter)->has_atom(atom->get_id()))
            return true;
     }
     return false;
}
     
bool ATOM::is_geminal_to(ATOM* atom){
     if(bond_.empty())
        return false;
     if(this == atom)
        return false;
     //ATOM& atom1 = const_cast<ATOM&>(atom);
     for (vector<BOND*>::iterator biter = atom->get_bond_list().begin(); biter != atom->get_bond_list().end(); ++biter){
         ATOM* atom_tmp = (*biter)->get_partner(atom);
         if (this->is_bonded_to(atom_tmp))
            return true;
     }
     return false;
}

bool ATOM::is_vicinal_to(ATOM* atom){
     if(bond_.empty())
        return false;
     if(this == atom)
        return false;
     //ATOM& atom1 = const_cast<ATOM&>(atom);
     for (vector<BOND*>::iterator biter = atom->get_bond_list().begin(); biter != atom->get_bond_list().end(); ++biter){
         ATOM* atom_tmp = (*biter)->get_partner(atom);
         if (this->is_geminal_to(atom_tmp))
            return true;
     }
     return false;
}
     
void ATOM::swap(ATOM& atom1, ATOM& atom2){
     ATOM temp_atom = atom1;
     atom1 = atom2;
     atom2 = temp_atom;
     }
     
bool ATOM::is_polar_hydrogen(){
     if (element_ == "H"){
                  for(ATOMVec::iterator aiter = atom_.begin(); aiter != atom_.end(); ++aiter){
                                        if((*aiter)->is_hetero())
                                            return true;
                  }
     }
     return false;
}

bool ATOM::is_nonpolar_hydrogen(){
     return !is_polar_hydrogen();
}

bool ATOM::bond_to_hydrogen(){
     for(ATOMVec::iterator aiter = atom_.begin(); aiter != atom_.end(); ++aiter)
         {
              if((*aiter)->is_hydrogen())
                  return true;
         }
	 return false;
}

bool ATOM::is_carboxyl(){
     if (symbol_type_ == "C.2" && num_of_bonded_atom_ == 3)
     {
         int oxygen_counter = 0;
         for(ATOMVec::iterator aiter = atom_.begin(); aiter != atom_.end(); ++aiter)
         {
              if((*aiter)->is_oxygen() && (*aiter)->get_num_bonded_heavy_atom() == 1)
              {
                  oxygen_counter += 1;
              }   
         }
         if(oxygen_counter == 2)
             return true;
         else
             return false;
     }
     return false;
}

bool ATOM::is_carboxyl_oxygen(){
     if (symbol_type_ == "O.co2" && num_of_bonded_atom_ == 1)
         return true;
     else if (is_oxygen() && num_of_bonded_heavy_atom_ == 1){
          for (ATOMVec::iterator aiter = atom_.begin(); aiter!=atom_.end(); ++aiter){
              if ((*aiter)->is_carboxyl())
                  return true;
              else
                  return false;
          }
     }
     return false;
}

bool ATOM::is_carbonyl()
{
     if (symbol_type_ == "C.2" && num_of_bonded_atom_ == 3)
     {
         int oxygen_counter = 0;
         for(ATOMVec::iterator aiter = atom_.begin(); aiter != atom_.end(); ++aiter)
         {
              if((*aiter)->is_oxygen() && (*aiter)->get_num_bonded_heavy_atom() == 1)
              {
                  oxygen_counter += 1;
              }
			  else if((*aiter)->is_sulfur() && (*aiter)->get_num_bonded_heavy_atom() == 1)
			  {
				  oxygen_counter += 1;
			  }
         }
         if(oxygen_counter == 1)
             return true;
         else
             return false;
     }
     return false;
}

bool ATOM::is_phosphate()
{
     if (is_phosphorus())
     {
         int oxygen_counter = 0;
         for(ATOMVec::iterator aiter = atom_.begin(); aiter != atom_.end(); ++aiter)
         {
              if((*aiter)->is_oxygen() && (*aiter)->get_num_bonded_heavy_atom() == 1)
                  oxygen_counter += 1;
         }
         if (oxygen_counter == 2 || oxygen_counter == 3 || oxygen_counter == 4)
             return true;
         else 
             return false;
     }
     return false;
}

bool ATOM::is_phosphate_oxygen(){
     if (is_oxygen() && num_of_bonded_heavy_atom_ == 1){
          for (ATOMVec::iterator aiter = atom_.begin(); aiter!=atom_.end(); ++aiter){
              if ((*aiter)->is_phosphate())
                  return true;
          }
     }
     return false;
}

bool ATOM::is_sulphate(){
     if (symbol_type_ == "S.O" || symbol_type_ == "S.O2")
     {
         int oxygen_counter = 0;
         for(ATOMVec::iterator aiter = atom_.begin(); aiter != atom_.end(); ++aiter)
         {
              if((*aiter)->is_oxygen() && (*aiter)->get_num_bonded_heavy_atom() == 1)
                  oxygen_counter += 1;
         }
         if(oxygen_counter == 1 || oxygen_counter ==2)
             return true;
         else
             return false;
     }
     return false;
}

bool ATOM::is_sulphate_oxygen(){
     if (is_oxygen() && num_of_bonded_heavy_atom_ == 1){
          for (ATOMVec::iterator aiter = atom_.begin(); aiter!=atom_.end(); ++aiter){
              if ((*aiter)->is_sulphate())
                  return true;
          }
     }
     return false;
}

bool ATOM::is_hydroxyl(){
     if (symbol_type_ == "O.3" && bond_to_hydrogen() && num_of_bonded_heavy_atom_ == 1)
         for (ATOMVec::iterator aiter = atom_.begin(); aiter != atom_.end(); ++aiter){
             if ((*aiter)->is_hydrogen())
                 continue;
             else if (!(*aiter)->is_carboxyl() && !(*aiter)->is_phosphorus() && !(*aiter)->is_sulfur())
                 return true;
             else
                 return false;
         }
     return false;
}

bool ATOM::is_nitro(){
     if (is_nitrogen() && (num_of_bonded_heavy_atom_ == 3 || num_of_bonded_heavy_atom_ == 2))
     {
         int oxygen_counter = 0;
         for(ATOMVec::iterator aiter = atom_.begin(); aiter != atom_.end(); ++aiter)
         {
              if((*aiter)->is_oxygen())
                  oxygen_counter += 1;
         }
         if (oxygen_counter == 2 || oxygen_counter == 1)
             return true;
         else
             return false;
     }
     return false;
}

bool ATOM::is_nitro_oxygen(){
     if (is_oxygen() && num_of_bonded_heavy_atom_ == 1){
          for (ATOMVec::iterator aiter = atom_.begin(); aiter!=atom_.end(); ++aiter){
              if ((*aiter)->is_nitro())
                  return true;
          }
     }
     return false;
}
 
bool ATOM::is_guanidino()
{
     if(symbol_type_ == "C.cat")
         return true;
     else if (is_carbon() || num_of_bonded_heavy_atom_ == 3)
	 {
         int nitrogen_counter = 0;
         int tri_valence_counter = 0;
         for(ATOMVec::iterator aiter = atom_.begin(); aiter != atom_.end(); ++aiter)
         {
			 if((*aiter)->is_nitrogen() && !(*aiter)->isAromaticRingAtom() && !(*aiter)->is_planar_amide())
			 {
                  nitrogen_counter += 1;
                  if ((*aiter)->get_num_bonded_atom() == 3)
                      tri_valence_counter += 1;
             }
         }
		 // fixed by xfliu, 20090327
         if (nitrogen_counter == 3 && tri_valence_counter >= 2)
            return true;
         else
             return false;
     }
     else
         return false;
}

bool ATOM::is_amidine()
{
     if (is_carbon() || num_of_bonded_heavy_atom_ == 3)
	 {
         int nitrogen_counter = 0, carbon_counter = 0;
         int tri_valence_counter = 0;
         for(ATOMVec::iterator aiter = atom_.begin(); aiter != atom_.end(); ++aiter)
         {
			 if((*aiter)->is_nitrogen() && !(*aiter)->isAromaticRingAtom() && !(*aiter)->is_planar_amide())
			 {
                  nitrogen_counter += 1;
                  if ((*aiter)->get_num_bonded_atom() == 3)
                      tri_valence_counter += 1;
              }
			  else if((*aiter)->is_carbon())
				  carbon_counter += 1;
         }
         if (nitrogen_counter == 2 && tri_valence_counter >= 1 && carbon_counter == 1)
            return true;
         else
             return false;
     }
     else
         return false;
}

bool ATOM::is_amine()
{
	if (is_nitrogen() &&  num_of_bonded_atom_ == 3)
	{
		int counter = 0;
        for(ATOMVec::iterator aiter = atom_.begin(); aiter != atom_.end(); ++aiter)
		{
			if ((*aiter)->get_symbol_type() == "C.3" || (*aiter)->is_hydrogen())
				counter ++ ;
		}
		if (counter == 3)
			return true;
		else
			return false;
	}
	else
		return false;
}

bool ATOM::is_amide()
{
	bool flag = false;
	if(is_nitrogen())
	{
		if(symbol_type_ == "N.am")
			return true;
		for(ATOMVec::iterator it = atom_.begin(); it != atom_.end(); ++it)
		{
			if((*it)->is_carbonyl())
			{
				flag = true;
				break;
			}
		}
	}
	return flag;
}

bool ATOM::is_planar_amide()
{
	bool flag = false;
	if(is_nitrogen())
	{
		if(symbol_type_ == "N.2" || symbol_type_ == "N.pl3" || symbol_type_ == "N.ar")
			flag = true;
	}
	return flag;
}

bool ATOM::is_aromatic_hetero(){
     return (is_hetero() && is_ring && is_aromatic());
}

bool ATOM::is_pos_charged(){
     if(is_nitrogen() || is_phosphorus()){
         if(symbol_type_ == "N.4")
             return true;
         else 
             return false;
         /*else{
              float order = 0.0;
              for(BONDVec::iterator it = get_bond_list().begin(); it != get_bond_list().end();
                  ++it)
                  
                  order += (*it)->get_bond_order();
              if (order == 4.0)
                  return true;
              else return false;
              }*/
         }
     else
         return false;
}

bool ATOM::is_neg_charged(){
     if((is_oxygen() || is_sulfur()) && num_of_bond_ == 1)
         if(bond_[0]->is_single()&&!is_carboxyl_oxygen())
             return true;
         else return false;
     else
         return false;
}

bool ATOM::is_sp3_oxygen_acceptor()
{
	 if(symbol_type_ == "O.3" || symbol_type_ == "S.3")
	 {
		 int pi_counter = 0;
		 for(BONDVec::iterator it = bond_.begin(); it != bond_.end(); ++it)
			 if(!(*it)->is_single())
				 pi_counter += 1;
		 if(pi_counter == 2)
			 return true;
		 else 
			 return false;
	 }
	 else
		 return false;
}

// set rotation mark recursively
void ATOM::set_rotation_mark(ATOM* anchor){
     rotation_mark_ = true;
     for(ATOMVec::iterator aiter = atom_.begin(); aiter != atom_.end(); ++aiter){
         if(*aiter == anchor) continue;
         else if((*aiter)->get_rotation_mark() == true) continue;
         else
             (*aiter)->set_rotation_mark(this);
             
     }
}
unsigned int ATOM::BOSum()
{
	unsigned int bo;
    unsigned int bosum=0;

	for (BONDVec::iterator biter = bond_.begin(); biter != bond_.end(); ++biter)
	{
		
        bo = (*biter)->get_bond_order();
		if(bo==1.5)
		{
			bo=5;
		}
        bosum += (bo < 4) ? 2*bo : 3;
	}

    bosum /= 2;
    return bosum;
  }
 
                            
