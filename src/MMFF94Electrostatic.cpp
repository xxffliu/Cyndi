#include "../include/MMFF94Electrostatic.h"
#include "../include/utility.h"

using namespace std;

// constructors
MMFF94Ele::MMFF94Ele():FFComponent(),fast_access_chg(),fast_access_pbci(),Ele_data_holder_(),Ele_energy_(0.0)
{
	//set component name
    set_name("MMFF94 Ele");
}
MMFF94Ele::MMFF94Ele(ForceField& ff):FFComponent(ff),fast_access_chg(),fast_access_pbci(),Ele_data_holder_(),Ele_energy_(0.0)
{
	set_name("MMFF94 Ele");
}
//copy constructor
MMFF94Ele::MMFF94Ele(const MMFF94Ele& to_copy):FFComponent(to_copy)
{
	fast_access_chg = to_copy.fast_access_chg;
    fast_access_pbci = to_copy.fast_access_pbci;
    Ele_data_holder_ = to_copy.Ele_data_holder_;
    Ele_energy_ = to_copy.Ele_energy_;
}
// destructor
MMFF94Ele::~MMFF94Ele()
{
	fast_access_chg.clear();
	fast_access_pbci.clear();
    Ele_data_holder_.clear();
}
// extract Ele parameters from the non bonded FFParamter object and build a hashtable for fast access
bool MMFF94Ele::extract_Ele_parameters(FFParameter& ffp)
{
     if(!ffp.is_valid())
         return false;
     // 
     // build a 2 dim array of atom types and loop variables
     //FFParameter::AtomTypes atom_types = ffp.get_atomtypes();
     int num_types = ffp.get_num_types();
     int num_entries = num_types * num_types;
     num_of_atom_types_ = num_types;
     // clear internal data structures
	 fast_access_chg.clear();
	 fast_access_pbci.clear();
     //debug
     //cout<<get_force_field()->get_parameters().params_in_each_section["Ele"].size()<<endl;
     for (vector<vector<string> >::iterator it = force_field_->get_parameters().params_in_each_section["BCI"].begin(); it != force_field_->get_parameters().params_in_each_section["BCI"].end(); ++it)
	 {
		 int type =str2int((*it)[0]);
		 int i=str2int((*it)[1]);
		 int j=str2int((*it)[2]);
		 if (type == -1||i==-1||j==-1)
		 {
			 cout<<"MMFF94Ele::extract_Ele_parameters(): error! no numeric atom type defined for atom "<<type<<endl;
			 return false;
		 }
		 int cxq=mmff94_force_field->GetCXQ(type,i,j);
		 fast_access_chg[cxq].a=i;
		 fast_access_chg[cxq].b=j;
		 fast_access_chg[cxq].bci=str2double((*it)[3]);
		 fast_access_chg[cxq].is_defined=true;
	 }
	 for (vector<vector<string> >::iterator it = force_field_->get_parameters().params_in_each_section["PBCI"].begin(); it != force_field_->get_parameters().params_in_each_section["PBCI"].end(); ++it)
	 {
		 int k=str2int((*it)[1]);
		 if (k==-1)
		 {
			 cout<<"MMFF94Ele::extract_Ele_parameters(): error! no numeric atom type defined for atom "<<k<<endl;
			 return false;
		 }
		 int index=k;
		 fast_access_pbci[index].a=k;
		 fast_access_pbci[index].pbci=str2double((*it)[2]);
		 fast_access_pbci[index].fcadj=str2double((*it)[3]);
		 fast_access_pbci[index].is_defined=true;			 
	 }
     return true;
}

// query a set of Ele parameters has defined for a given combination of atom types
bool MMFF94Ele::has_params_chg(int i, int j,int type)
{
	int cxq=mmff94_force_field->GetCXQ(type,i,j);
    if ((i>=0 && i <= num_of_atom_types_) && (j>=0 && j <= num_of_atom_types_)&&(type==1||type==0||type==4))
	{ 
		return (fast_access_chg[cxq].is_defined);
	}
	return false;
}
bool MMFF94Ele::has_params_pbci(int i)
{
	const int index =i;
	if(i>=0&&i<=num_of_atom_types_)
		return(fast_access_pbci[index].is_defined);
	else 
		return false;
}
// assign the parameters for a given combination of atom types;
bool MMFF94Ele::assign_params( MMFF94Ele::EleForceValues &value,int i, int j)
{
	ATOM * atom1,*atom2;
	atom1=force_field_->get_mol()->get_map(i);
	atom2=force_field_->get_mol()->get_map(j);
	if (has_params_pbci(i)&&has_params_pbci(j))
	{
		value.is_14_interaction=(atom2->is_vicinal_to(atom1));
		value.is_defined=true;
		value.q_i=atom1->get_mmff94_partial_charge();
		value.q_j=atom2->get_mmff94_partial_charge();
		return true;
	}
	else
		return false;
}
// set up method
bool MMFF94Ele::setup()
{
     if (force_field_ == NULL)
	 {
		 cout<<"MMFF94Ele::setup(): force field bound component can not be found"<<endl;
         return false;
	 }
     // clear the Ele parameter container
     Ele_data_holder_.clear();
     // tempararily set this component enabled;
     setenabled(true);
	 bool is_14 = false;
     // extract the L-J Ele parameters
     mmff94_force_field = dynamic_cast<MMFF94*>(force_field_);
	 if (mmff94_force_field == NULL || mmff94_force_field->has_initialized_param())
	 {
		 bool result = extract_Ele_parameters(force_field_->get_parameters());
		 if (!result)
		 {
			 cout << "MMFF94Ele::setup(): can not access L-J Ele section"<<endl;
			 return false;
		 }
	 }
     // now iterate all atom pairs and remove 1-2, 1-3 pairs
	 //MMFF94Ele::EleForceValues value;
	 SetFormalCharges();
	 SetPartialCharges();
	 for (ATOMVec::iterator aiter1 = force_field_->get_atoms().begin();aiter1 != force_field_->get_atoms().end(); ++aiter1)
	 {
		 for (ATOMVec::iterator aiter2 = aiter1 + 1; aiter2 != force_field_->get_atoms().end(); ++aiter2)
		 {
			 bool is_14 = false;
			 // check if is bonded
			 if ((*aiter2)->is_bonded_to(*aiter1))
				 continue;
             // then check if is geminal
             else if ((*aiter2)->is_geminal_to(*aiter1))
				 continue;
			 else
			 {
				 if((*aiter2)->is_vicinal_to(*aiter1))
					 is_14 = true;
				 Ele_data_holder_.push_back(MMFF94Ele::EleData());
				 Ele_data_holder_.back().atom1 = *aiter1;
				 Ele_data_holder_.back().atom2 = *aiter2;                       
				 //Ele_data_holder_.back().value = value; 
				 Ele_data_holder_.back().is_14_interaction = is_14;

			 }
		 }
	 }
	 return true;
}

double MMFF94Ele::update_energy()
{
	energy_ = 0.0;
	Ele_energy_ = 0.0;
	// provisionally set the 14_scale_factor to 1.0 
	float ele_14_scale_factor = 0.75;
	// loop all non bond atom pairs and calculate Ele and ele energies
	for (vector<EleData>::size_type i = 0; i < Ele_data_holder_.size();++i)
	{
		//vector3 rij;
		double r_ij;
		r_ij = Ele_data_holder_[i].atom1->get_position().dist(Ele_data_holder_[i].atom2->get_position());
		//r_ij=rij.length();
		if(r_ij	< 0.0001)
			continue;
		//Ele_data_holder_[i].value.r_ij=r_ij;
		// Unit: Kcal/mol
		double first=332.0716*Ele_data_holder_[i].atom1->get_mmff94_partial_charge()*Ele_data_holder_[i].atom2->get_mmff94_partial_charge();
		double second=1.0/(r_ij+0.05); 
		if(!Ele_data_holder_[i].is_14_interaction)
		{
			energy_ += first*second;
			//cout<<Ele_data_holder_[i].atom1->get_id()<<"-"<<Ele_data_holder_[i].atom2->get_id()<<" "<<first*second<<endl;
		}
		else
		{
			energy_ += ele_14_scale_factor * first * second;
			//cout<<"Short: "<<Ele_data_holder_[i].atom1->get_id()<<"-"<<Ele_data_holder_[i].atom2->get_id()<<" "<<0.75*first*second<<endl;
		}
	}
	return energy_;            
 }
// get Ele energy
double MMFF94Ele::get_Ele_energy()const
{
       return Ele_energy_;
}

// calculate the forces imposed on each atoms by Ele and ele interaction
void MMFF94Ele::update_forces()
{
	if (force_field_ == NULL)
	{
		cout<<"MMFF94Ele::update_force(): error! this component doesn't bond to any force field"<<endl;
		return;
	}
    // iterate all non bond pairs and update forces	
    for (vector<EleData>::size_type i = 0; i< Ele_data_holder_.size();++i)
	{
		vector3 force_1 = Ele_data_holder_[i].atom1->get_force(), force_2 = Ele_data_holder_[i].atom2->get_force();
        vector3 direction(Ele_data_holder_[i].atom1->get_position() - Ele_data_holder_[i].atom2->get_position());
        double distance = direction.length();
        double inverse_distance = 1./distance;
		//direction = direction / distance;
        direction.normalize();
        if (!isNearZero(distance))
		{
			double Ele_factor = 332.0716 * Ele_data_holder_[i].atom1->get_mmff94_partial_charge() * Ele_data_holder_[i].atom2->get_mmff94_partial_charge() / ((distance + 0.05) * (distance + 0.05));
			if(Ele_data_holder_[i].is_14_interaction)
			{
				Ele_factor *= 0.75;
			}
			vector3 force = direction * Ele_factor * FORCE_FACTOR;
			force_1 += force;
			force_2 -= force;
			Ele_data_holder_[i].atom1->set_force(force_1);
            Ele_data_holder_[i].atom2->set_force(force_2);
		}
	}
}

bool MMFF94Ele::SetPartialCharges()
{
	vector<double> charges(( force_field_->get_mol()->get_num_atom()+1), 0);
    double M, Wab, factor, q0a, q0b, Pa, Pb;
	for(ATOMVec::iterator aiter=force_field_->get_atoms().begin();aiter!=force_field_->get_atoms().end();aiter++)
	{
		int type;
		type=(*aiter)->get_mmff94_type();
		switch (type)
		{
			case 32:
			case 35:
			case 72:
				factor = 0.5;
			break;
			case 62:
			case 76:
				factor = 0.25;
			break;
			default:
				factor = 0.0;
				break;
		}
		M = mmff94_force_field->GetCrd(type);
		q0a = (*aiter)->get_mmff94_partial_charge();
		// charge sharing
		if (!factor)
		{
			for(ATOMVec::iterator aiter_1=(*aiter)->get_atom_list().begin();aiter_1!=(*aiter)->get_atom_list().end();aiter_1++)
			{
				if((*aiter_1)->get_mmff94_partial_charge()<0.0)
				{
					q0a+=(*aiter_1)->get_mmff94_partial_charge()/(2.0*(*aiter_1)->get_num_neighbor_bond());
				}
			}
		}
		// needed for SEYWUO, positive charge sharing?
		if (type == 62)
		{
			for(ATOMVec::iterator aiter_1=(*aiter)->get_atom_list().begin();aiter_1!=(*aiter)->get_atom_list().end();aiter_1++)
			{
				if((*aiter_1)->get_mmff94_partial_charge()>0.0)
				{
					q0a-=((*aiter_1)->get_mmff94_partial_charge()/2.0);
				}
			}
		}
		q0b = 0.0;
		Wab = 0.0;
		Pa = Pb = 0.0;
		for(ATOMVec::iterator aiter_1=(*aiter)->get_atom_list().begin();aiter_1!=(*aiter)->get_atom_list().end();aiter_1++)
		{
			int nbr_type=(*aiter_1)->get_mmff94_type();
			q0b+=(*aiter_1)->get_mmff94_partial_charge();
			int bond_type;
			bond_type=mmff94_force_field->GetBondType((*aiter),(*aiter_1));
			int cxq=mmff94_force_field->GetCXQ(bond_type,type,nbr_type);
			int cxq_sym=mmff94_force_field->GetCXQ(bond_type,nbr_type,type);
			if(fast_access_chg[cxq].is_defined)
			{
				Wab += -fast_access_chg[cxq].bci;
			}
			else if(fast_access_chg[cxq_sym].is_defined)
			{
				Wab += fast_access_chg[cxq_sym].bci;
			}
			else
			{
				Pa = fast_access_pbci[type].pbci;
				Pb = fast_access_pbci[nbr_type].pbci;
				Wab += Pa - Pb;
			}
		}
		if (factor)
			charges[(*aiter)->get_id()] = (1.0 - M * factor) * q0a + factor * q0b +Wab;
		else 
			charges[(*aiter)->get_id()] = q0a+Wab;
	}
	for(ATOMVec::iterator aiter_2=force_field_->get_atoms().begin();aiter_2!=force_field_->get_atoms().end();aiter_2++)
	{
		(*aiter_2)->set_mmff94_partial_charge(charges[(*aiter_2)->get_id()]);
	}
	return true;
}
  // we set the the formal charge with SetPartialCharge because formal charges 
  // in MMFF94 are not always and integer
bool MMFF94Ele::SetFormalCharges()
{
	for(ATOMVec::iterator aiter=force_field_->get_atoms().begin();aiter!=force_field_->get_atoms().end();aiter++)
	{
		int type;
		type=mmff94_force_field->GetType((*aiter)).numeric;
		(*aiter)->set_mmff94_partial_charge(0.0);
		bool done = false;
		switch (type)
		{
			case 34:
			case 49:
			case 51:
			case 54:
			case 58:
			case 92:
			case 93:
			case 94:
			case 97:
				(*aiter)->set_mmff94_partial_charge(1.0);
				done = true;
				break;
			case 35:
            case 62:
            case 89:
            case 90:
			case 91:
				(*aiter)->set_mmff94_partial_charge(-1.0);
				done = true;
				break;
            case 55:
				(*aiter)->set_mmff94_partial_charge(0.5);
				done = true;
                break;
			case 56:
				(*aiter)->set_mmff94_partial_charge(1.0/3.0);
				done = true;
                break;
			case 87:
			case 95:
			case 96:
			case 98:
			case 99:
				(*aiter)->set_mmff94_partial_charge(2.0);
				done = true;
				break;
                              //case 98:
							  //(*aiter)->set_mmff94_partial_charge(3.0);
			default:
				break;
		}
		if (done)
			continue;
		if (type == 32)
		{
			int o_count = 0;
			bool sulfonamide = false;
			int s_count = 0; 
			for(ATOMVec::iterator aiter_1=(*aiter)->get_atom_list().begin();aiter_1!=(*aiter)->get_atom_list().end();aiter_1++)
			{
				for(ATOMVec::iterator aiter_2=(*aiter_1)->get_atom_list().begin();aiter_2!=(*aiter_1)->get_atom_list().end();aiter_2++)
				{
					if((*aiter_2)->is_oxygen()&&((*aiter_2)->get_num_neighbor_bond()==1))  
						o_count++;
					if((*aiter_2)->is_sulfur()&&((*aiter_2)->get_num_neighbor_bond()==1))  
						s_count++;
					if((*aiter_2)->is_nitrogen()&&(!(*aiter_2)->is_aromatic()))    
						sulfonamide = true;
				}
				if((*aiter_1)->is_carbon())
					(*aiter)->set_mmff94_partial_charge(-0.5);// O2CM
				if((*aiter_1)->is_nitrogen()&&(o_count == 3))
					(*aiter)->set_mmff94_partial_charge(-1.0 / o_count); // O3N
				if(((*aiter_1)->is_sulfur())&&(!sulfonamide))
				{
					if(((o_count + s_count) == 2)&&(*aiter_1)->get_num_neighbor_bond()==3)
					{
						if (((o_count + s_count) == 2) && ((*aiter_1)->get_num_neighbor_bond() == 3) && ((*aiter_1)->BOSum() == 3))
						{
							(*aiter)->set_mmff94_partial_charge(-0.5);//O2S
						}
						else if ((o_count + s_count) == 3)
							(*aiter)->set_mmff94_partial_charge(-1.0 / 3.0);//O3S
						else if((o_count + s_count) == 4)
							(*aiter)->set_mmff94_partial_charge(-0.5);//O4S
						if((*aiter_1)->is_phosphorus())
						{
							if((o_count + s_count) == 2)
								(*aiter)->set_mmff94_partial_charge(-0.5);//O2P
							else if((o_count + s_count) == 3)
								(*aiter)->set_mmff94_partial_charge(-2.0 / 3.0); // O3P
							else if((o_count + s_count) == 4)
								(*aiter)->set_mmff94_partial_charge(-0.25);//O4P
						}
						if (type == 77)
							(*aiter)->set_mmff94_partial_charge(-0.25);// O4CL
					}
				}
			}
		}
		else if (type == 61)
		{
			bool find_bond=false;
			for(ATOMVec::iterator aiter_1=(*aiter)->get_atom_list().begin();aiter_1!=(*aiter)->get_atom_list().end();aiter_1++)
			{
				for(BONDVec::iterator biter=(*aiter)->get_bond_list().begin();biter!=(*aiter)->get_bond_list().end();biter++)
				{
					if(((*biter)->get_first_atom()==(*aiter)&&(*biter)->get_second_atom()==(*aiter_1))||((*biter)->get_first_atom()==(*aiter_1)&&(*biter)->get_second_atom()==(*aiter)))
					{
						find_bond=true;
						break;
					}
					// fixed by xfliu, 20090310
					if(find_bond == true)
					{
						if((*biter)->is_triple()&&(*aiter_1)->is_nitrogen())
							(*aiter)->set_mmff94_partial_charge(1.0);
					}
				}
			}
		}
		else if (type == 72)
		{
			int s_count = 0;
			for(ATOMVec::iterator aiter_1=(*aiter)->get_atom_list().begin();aiter_1!=(*aiter)->get_atom_list().end();aiter_1++)
			{
				if((*aiter_1)->is_sulfur())
					s_count++;
				if((*aiter_1)->is_phosphorus()||(*aiter_1)->is_sulfur())
				{
					for(ATOMVec::iterator aiter_2=(*aiter_1)->get_atom_list().begin();aiter_2!=(*aiter_1)->get_atom_list().end();aiter_2++)
						if(((*aiter_2)->is_sulfur()||(*aiter_2)->is_oxygen()) && (*aiter_2)->get_num_neighbor_bond()==1 && (*aiter)->get_id()!=(*aiter_2)->get_id())
							((*aiter)->set_mmff94_partial_charge(-0.5));
				}
				else
					(*aiter)->set_mmff94_partial_charge(-1.0);
				if((*aiter_1)->is_carbon())
					for(ATOMVec::iterator aiter_2=(*aiter_1)->get_atom_list().begin();aiter_2!=(*aiter_1)->get_atom_list().end();aiter_2++)
						if((*aiter_2)->is_sulfur()&&((*aiter_2)->get_num_neighbor_bond()==1)&&((*aiter)->get_id()!=(*aiter_2)->get_id()))
							(*aiter)->set_mmff94_partial_charge(-0.5);//SSMO
				if (s_count >= 2)
					(*aiter)->set_mmff94_partial_charge(-0.5); // SSMO
			}
		}
		else if (type == 76)
		{
			vector<RING> vr;
			vr=force_field_->get_mol()->get_ring_vector();
			vector<RING>::iterator ri;
			vector<ATOM*>::iterator rj;
			int n_count;
			for (ri = vr.begin();ri != vr.end();ri++)
			{ // for each ring
				n_count = 0;
				if (((*ri).is_aromatic) && ((*ri).IsMember((*aiter)))&& ((*ri).size== 5))
				{
					for(rj = (*ri).vatom.begin();rj != (*ri).vatom.end();rj++) // for each ring atom
						if ((*rj)->is_nitrogen())
							n_count++;
					if (n_count > 1)
						(*aiter)->set_mmff94_partial_charge(-1.0 / n_count);
				}
			}
		}
		else if (type == 81)
		{
			(*aiter)->set_mmff94_partial_charge(1.0);
			vector<RING>vr;
			vr=force_field_->get_mol()->get_ring_vector();
			vector<RING>::iterator ri;
			vector<ATOM*>::iterator rj;
			for(ri=vr.begin();ri!=vr.end();ri++)// for each ring
				if (((*ri).is_aromatic)&& ((*ri).IsMember((*aiter))) && ((*ri).size == 5))
				{
					int n_count = 0;
					for(rj = (*ri).vatom.begin();rj != (*ri).vatom.end();rj++) // for each ring atom
						if((*rj)->is_nitrogen()&&((*rj)->get_num_neighbor_bond()==3))     
							n_count++;
					(*aiter)->set_mmff94_partial_charge((1.0 / n_count)); // NIM+
					for(ATOMVec::iterator aiter_1=(*aiter)->get_atom_list().begin();aiter_1!=(*aiter)->get_atom_list().end();aiter_1++)
						for(ATOMVec::iterator aiter_2=(*aiter_1)->get_atom_list().begin();aiter_2!=(*aiter_1)->get_atom_list().end();aiter_2++)
							if((*aiter_2)->get_mmff94_type()==56)
								(*aiter)->set_mmff94_partial_charge(1.0 / 3.0);
					for(ATOMVec::iterator aiter_1=(*aiter)->get_atom_list().begin();aiter_1!=(*aiter)->get_atom_list().end();aiter_1++)
						for(ATOMVec::iterator aiter_2=(*aiter_1)->get_atom_list().begin();aiter_2!=(*aiter_1)->get_atom_list().end();aiter_2++)
							if((*aiter_2)->get_mmff94_type()==55)
								(*aiter)->set_mmff94_partial_charge(1.0 / (1.0 + n_count));
				}
		}
	}
	return true;
}
       
            
void MMFF94Ele::update(){;}    