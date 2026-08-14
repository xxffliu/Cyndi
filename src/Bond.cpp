#include "../include/Bond.h"
#include <algorithm>

BOND::BOND():
             is_rotor(false),
             id(0),
			 AromaticRingBond_(false),
             type("UND"),
             first_atom(NULL),
             second_atom(NULL),
			 bond_order_(0){}

// we leave the bond atom assignment to the constructor of MOL
BOND::BOND(const BOND& bond0):
			 id(bond0.id),type(bond0.type),is_rotor(bond0.is_rotor), bond_order_(bond0.bond_order_),AromaticRingBond_(bond0.AromaticRingBond_){}

BOND& BOND::operator=(const BOND& bond0){
      id = bond0.id;
      type = bond0.type;
      is_rotor = bond0.is_rotor;
	  bond_order_ = bond0.bond_order_;
	  AromaticRingBond_ = bond0.AromaticRingBond_;
      return *this;
      }

BOND::~BOND(){}

void BOND::clear()
{
     first_atom = NULL;
     second_atom = NULL;
     id = 0;
     type = "UNK";
     is_rotor = false;
	 bond_order_ = 0;
	 AromaticRingBond_ = false;
     }

ATOM* BOND::get_partner(const ATOM* atom)
{
	if (atom == first_atom)
		return second_atom;
	else if (atom == second_atom)
		return first_atom;
	else
	{
		cout<<"Error! BOND::get_partner():this atom doesn't belong to this bond"<<endl;
		exit(1);
	}
}

bool BOND::has_atom(int id)
{
				 if((id == first_atom->get_id()) || id == second_atom->get_id())
					return true;
				 return false;
}

void BOND::swap(BOND& bond1, BOND& bond2)
{
     BOND temp_bond = bond1;
     bond1 = bond2;
     bond2 = temp_bond;
     }

bool BOND::is_in_ring()
{
	if(first_atom->is_ring && second_atom->is_ring)
	{
		for(vector<int>::iterator it = first_atom->ring_id.begin(); it != first_atom->ring_id.end(); ++it)
		{
			if(find(second_atom->ring_id.begin(), second_atom->ring_id.end(), *it) != second_atom->ring_id.end())
				return true;
		}
		return false;
	}
	else
		return false;
}

bool BOND::is_amide()
{
	bool flag = false;
	if(type == "am" || type == "1")
	{
		if(first_atom->is_carbonyl() && second_atom->is_amide())
			flag = true;
		else if(second_atom->is_carbonyl() && first_atom->is_amide())
			flag = true;
	}
	return flag;
}

bool BOND::is_guanidino()
{
	bool flag = false;
	if(type == "1" || type == "2")
	{
		if(first_atom->is_nitrogen() && second_atom->is_guanidino())
		{
			flag = true;
		}
		else if(second_atom->is_nitrogen() && first_atom->is_guanidino())
		{
			flag = true;
		}
	}
	return flag;
}
                                 
