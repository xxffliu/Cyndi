#include "../include/MMFF94OOP.h"
#include "../include/MMFF94Bend.h"
#include "../include/Vector3.h"
#include "../include/utility.h"
using namespace std;
class vector3;
// constructors
/*#ifndef FORCE_FACTOR
#define FORCE_FACTOR 1E10/AVOGADRO
#endif*/
MMFF94OOP::MMFF94OOP():FFComponent(),fast_access(),oop_data_holder_()
{
	set_name("MMFF94 OOP");
}
MMFF94OOP::MMFF94OOP(ForceField& ff):FFComponent(ff),fast_access(),oop_data_holder_()
{
	set_name("MMFF94 OOP");
}

// copy constructor
MMFF94OOP::MMFF94OOP(const MMFF94OOP& to_copy):FFComponent(to_copy)
{
	fast_access = to_copy.fast_access;
    oop_data_holder_ = to_copy.oop_data_holder_;
}
// destructor
MMFF94OOP::~MMFF94OOP()
{
	fast_access.begin();
    oop_data_holder_.clear();
}

// extract angle oop ff parameters from the MMFF94 parameter file
// and build some data structure for fast access these data 
bool MMFF94OOP::extract_OOP_parameters(FFParameter& ffp)
{
     // check the force field is valid
     if (!ffp.is_valid())
	 {
         //debug
         cout<<"MMFF94OOP::extract_OOP_parameters(): ffp is not valid"<<endl;
         return false;
     }
     // build a one dim array of atom types
     //FFParameter::AtomTypes atom_types = ffp.get_atomtypes();
     int num_types = ffp.get_num_types();
     num_of_atom_types_ = num_types;
     // clear internal data structures
  
	 fast_access.clear();
     // start to pack the parameters into the vector of fast access
     int center_type;
	 int i,k,l;
	 int cxo;
	 double k0;

     for (vector<vector<string> >::iterator it = force_field_->get_parameters().params_in_each_section["oopbend"].begin();it != force_field_->get_parameters().params_in_each_section["oopbend"].end(); ++it)
	 {
               i=str2int((*it)[0]);
			   center_type=str2int((*it)[1]);	
			   k=str2int((*it)[2]);
			   l=str2int((*it)[3]);
			   k0=str2double((*it)[4]);
			   if ((center_type == -1)||(i==-1)||(k==-1)||(l==-1))
			   {
				   cout<<"error! these atoms are not existing! "<<endl;
                   return false;
               }
			   cxo= mmff94_force_field->GetCXO(i,center_type,k,l);
			   fast_access[cxo].is_defined=true;
			   fast_access[cxo].type1=i;
			   fast_access[cxo].type_center=center_type;
			   fast_access[cxo].type2=k;
			   fast_access[cxo].type3=l;
			   fast_access[cxo].k0=k0;
			   cxo= mmff94_force_field->GetCXO(i,center_type,l,k);
			   fast_access[cxo].is_defined=true;
			   fast_access[cxo].type1=i;
			   fast_access[cxo].type_center=center_type;
			   fast_access[cxo].type2=l;
			   fast_access[cxo].type3=k;
			   fast_access[cxo].k0=k0;
			   cxo= mmff94_force_field->GetCXO(k,center_type,i,l);
			   fast_access[cxo].is_defined=true;
			   fast_access[cxo].type1=k;
			   fast_access[cxo].type_center=center_type;
			   fast_access[cxo].type2=i;
			   fast_access[cxo].type3=l;
			   fast_access[cxo].k0=k0;
			   cxo= mmff94_force_field->GetCXO(k,center_type,l,i);
			   fast_access[cxo].is_defined=true;
			   fast_access[cxo].type1=k;
			   fast_access[cxo].type_center=center_type;
			   fast_access[cxo].type2=l;
			   fast_access[cxo].type3=i;
			   fast_access[cxo].k0=k0;
			   cxo= mmff94_force_field->GetCXO(l,center_type,i,k);
			   fast_access[cxo].is_defined=true;
			   fast_access[cxo].type1=l;
			   fast_access[cxo].type_center=center_type;
			   fast_access[cxo].type2=i;
			   fast_access[cxo].type3=k;
			   fast_access[cxo].k0=k0;
			   cxo= mmff94_force_field->GetCXO(l,center_type,k,i);
			   fast_access[cxo].is_defined=true;
			   fast_access[cxo].type1=l;
			   fast_access[cxo].type_center=center_type;
			   fast_access[cxo].type2=k;
			   fast_access[cxo].type3=i;
			   fast_access[cxo].k0=k0;
			   //willson angle's value will be given after  the following  caculation
	 }
     return true;
}

// query a center atom parameters has defined for a given center atom type
bool MMFF94OOP::has_params(int i,int j, int k, int l)
{
     if ((i < 0 || i > num_of_atom_types_)||(j<0||j>num_of_atom_types_)||(k<0||k>num_of_atom_types_)||(l<0||l>num_of_atom_types_))
         return false;
	 int cxo = mmff94_force_field->GetCXO(i,j,k,l);
	 return fast_access[cxo].is_defined;
}

// assign the parameters for a given combination of atom types;
bool MMFF94OOP::assign_params(MMFF94OOP::ForceValues& param, int i ,int j, int k ,int l)
{
     if (has_params(i, j, k ,l))
	 {
		 int cxo = mmff94_force_field->GetCXO(i,j,k,l);
		 param.k0 = fast_access[cxo].k0;
         return true;
	 }
     else
         return false;
}
// set up method
bool MMFF94OOP::setup()
{
	if (force_field_ == NULL)
	{
		cout<<"MMFF94OOP::setup(): force field bound component can not be found"<<endl;
		return false;
	}
     // clear the parameter holder
    oop_data_holder_.clear();
    // tempararily set this component enabled;
    setenabled(true);
     
    mmff94_force_field = dynamic_cast<MMFF94*>(force_field_);
    if (mmff94_force_field == NULL || mmff94_force_field->has_initialized_param())
	{
		bool result = extract_OOP_parameters(get_force_field()->get_parameters());
		if (!result)
		{
			cout << "MMFF94OOP::setup(): can not found OOP section"<<endl;
			return false;
		}
	}
                          
   
	//initiating an integer array for storing the center atoms
	int center_atom_array[]={2,3,8,10,17,26,30,37,39,40,41,43,45,49,54,55,56,57,58,63,64,67,69,73,78,80,81,82};
    // retrieve all stretch parameters
    for (ATOMVec::iterator aiter = force_field_->get_atoms().begin(); aiter != force_field_->get_atoms().end(); ++aiter)
	{
		// remove 5 coordinated bonds
		if ((*aiter)->get_num_bonded_atom() > 4)
			continue;
		int center_atom = (*aiter)->get_mmff94_type();
		//converting the array to a vector,and finding whether the atom in the molecular is a center one or not
		vector<int>center_vec (center_atom_array,center_atom_array+28);
		vector<int>::iterator it;
		it=find(center_vec.begin(),center_vec.end(),center_atom);
		if(it==center_vec.end())
			continue;  
		// check if the root atom has 3 or more than 3 neighbor atoms
		if ((*aiter)->get_num_neighbor_bond() < 3)
		{
			continue;
		}
		MMFF94OOP::OOPData data;
		BONDVec OOPBonds = (*aiter)->get_bond_list();
		ATOM* atom1 = OOPBonds[0]->get_partner(*aiter);
		int atom_type_A = atom1->get_mmff94_type();
		ATOM* atom2 = OOPBonds[1]->get_partner(*aiter);
		int atom_type_B = atom2->get_mmff94_type();
		ATOM* atom3 = OOPBonds[2]->get_partner(*aiter);
		int atom_type_C = atom3->get_mmff94_type();
			
		data.atom1 = atom1;
		data.atom2 = atom2;
		data.atom3 = atom3;
		data.center = *aiter;
		if(has_params(atom_type_A, center_atom, atom_type_B ,atom_type_C))
		{
			int cxo = mmff94_force_field->GetCXO(atom_type_A,center_atom,atom_type_B,atom_type_C);
			data.k = fast_access[cxo].k0;
			oop_data_holder_.push_back(data);
			continue;
		}
		else if(has_params(mmff94bend.EqLvl2(atom_type_A),center_atom,mmff94bend.EqLvl2(atom_type_B),mmff94bend.EqLvl2(atom_type_C)))
		{
			int i = mmff94bend.EqLvl2(atom_type_A);
			int k = mmff94bend.EqLvl2(atom_type_B);
			int l = mmff94bend.EqLvl2(atom_type_C);
			int j = center_atom;
			int cxo = mmff94_force_field->GetCXO(i,j,k,l);
			data.k = fast_access[cxo].k0;
			oop_data_holder_.push_back(data);
			continue;
		}
		else if(has_params(mmff94bend.EqLvl3(atom_type_A),center_atom,mmff94bend.EqLvl3(atom_type_B),mmff94bend.EqLvl3(atom_type_C)))
		{
			int i = mmff94bend.EqLvl3(atom_type_A);
			int k = mmff94bend.EqLvl3(atom_type_B);
			int l = mmff94bend.EqLvl3(atom_type_C);
			int j = center_atom;
			int cxo = mmff94_force_field->GetCXO(i,j,k,l);
			data.k = fast_access[cxo].k0;
			oop_data_holder_.push_back(data);
			continue;
		}
		else if(has_params(mmff94bend.EqLvl4(atom_type_A),center_atom,mmff94bend.EqLvl4(atom_type_B),mmff94bend.EqLvl4(atom_type_C)))
		{
			int i = mmff94bend.EqLvl4(atom_type_A);
			int k = mmff94bend.EqLvl4(atom_type_B);
			int l = mmff94bend.EqLvl4(atom_type_C);
			int j = center_atom;
			int cxo = mmff94_force_field->GetCXO(i,j,k,l);
			data.k = fast_access[cxo].k0;
			oop_data_holder_.push_back(data);
			continue;
		}
		else if(has_params(mmff94bend.EqLvl5(atom_type_A),center_atom,mmff94bend.EqLvl5(atom_type_B),mmff94bend.EqLvl5(atom_type_C)))
		{
			int i = mmff94bend.EqLvl5(atom_type_A);
			int k = mmff94bend.EqLvl5(atom_type_B);
			int l = mmff94bend.EqLvl5(atom_type_C);
			int j = center_atom;
			int cxo = mmff94_force_field->GetCXO(i,j,k,l);
			data.k = fast_access[cxo].k0;
			oop_data_holder_.push_back(data);
			continue;
		}
		else
		{
			//this is accroding to the rule of the  page 628 of MMFF.V 's
			int cxo = mmff94_force_field->GetCXO(0,center_atom,0,0);
			data.k = fast_access[cxo].k0;
			oop_data_holder_.push_back(data);
			continue;
		}
#ifdef DEBUG
                     cout<<"MMFF94OOP::setup(): cannot find oop parameters for "<<
                     force_field_->get_parameters().get_type_name(atom_type_root)<<endl;
#endif
	}
	 
// everything goes well
	return true;
}
// update methods
double MMFF94OOP::update_energy()
{
	//energy initializion
    energy_ = 0.0;  
	double part_energy=0.0;
	double r_ij,r_kj,r_lj;
	vector3 ri,rk,rl;
	double wilson[3];
	vector3 ni,nk,nl;
	double product_length_i,product_length_k,product_length_l;   
    // iterate all OOP and summarize the oop energies
	for (vector<OOPData>::size_type i = 0; i<oop_data_holder_.size(); ++i)
	{
		//debug
		//cout<<oop_data_holder_[i].center->get_id()<<"-"<<oop_data_holder_[i].atom1->get_id()<<"-"<<oop_data_holder_[i].atom2->get_id()<<"-"<<oop_data_holder_[i].atom3->get_id()<<"-"<<oop_data_holder_[i].k<<endl;
		bool erro=false;
		ri=oop_data_holder_[i].center->get_position()-oop_data_holder_[i].atom1->get_position();
		rk=oop_data_holder_[i].center->get_position()-oop_data_holder_[i].atom2->get_position();
		rl=oop_data_holder_[i].center->get_position()-oop_data_holder_[i].atom3->get_position();
        ni=rk%rl;
		nk=ri%rl;
		nl=ri%rk;
		product_length_i=ri.length_2()*ni.length_2();
		product_length_k=rk.length_2()*nk.length_2();
		product_length_l=rl.length_2()*nl.length_2();
		r_ij=ri.length();
		r_kj=rk.length();
		r_lj=rl.length();
		//abort the calculation if any bond length is near zero
		if(r_ij<0.0001||r_kj<0.0001||r_lj<0.0001)
		{
			erro=true;
			continue;
		}
		//abort the wilson angle if the angle of the angle of the two bonds in the plane is near zero or near PI
		const double cos_theta1 = (rk*rl)/(r_kj*r_lj);
		const double theta1 = acos(cos_theta1);
		// If theta equals 180 degree or 0 degree
		if (theta1<0.0001||(fabs(theta1- PI))<0.0001)
			wilson[0]=0.0;
		else if(product_length_i>0.0001)
		{
			wilson[0]=asin(abs(ri*ni)/sqrt(product_length_i));
		}
		else wilson[0]=0.0;
		const double cos_theta2 =(ri*rl)/(r_ij*r_lj);
		const double theta2 = acos(cos_theta2);
		// If theta equals 180 degree or 0 degree
		if (theta2<0.0001||(fabs(theta2 - PI))<0.0001)
			wilson[1]=0.0;
		else if(product_length_k>0.0001)
		{
			wilson[1]=asin(abs(rk*nk)/sqrt(product_length_k));
		}
		else wilson[1]=0.0;
		const double cos_theta3 =(rk*ri)/(r_kj*r_ij);
		const double theta3 = acos(cos_theta3);
		// If theta equals 180 degree or 0 degree
		if (theta3<0.0001||(fabs(theta3 - PI))<0.0001)
			wilson[2]=0.0;
		else if(product_length_l>0.0001)
		{
			wilson[2]=asin(abs(rl*nl)/sqrt(product_length_l));
		}
		else wilson[2]=0.0;
		for(int p=0;p<3;p++)
		{
			energy_ += 0.043844/2*oop_data_holder_[i].k*wilson[p]*wilson[p]*RAD_TO_DEG*RAD_TO_DEG;
		}
	}
	//cout<<"oop_energy: "<<energy_<<endl;
    return energy_;
}
// calculate current forces imposed by stretch and add them to the force field;
void MMFF94OOP::update_forces()
{
	if (get_force_field() == 0)
	{
		cout<<"MMFF94OOP::update_force(): error! this component doesn't bond to any force field"<<endl;
		return;
	}
	double K0 = 0.043844 * RAD_TO_DEG * RAD_TO_DEG;
    // iterate all OOP and update forces
    for (vector<OOPData>::size_type i = 0; i < oop_data_holder_.size(); ++i)
	{
		ATOM* i1 = oop_data_holder_[i].atom1, *i2 = oop_data_holder_[i].atom1, *i3 = oop_data_holder_[i].atom2;
		ATOM* k1 = oop_data_holder_[i].atom2, *k2 = oop_data_holder_[i].atom3, *k3 = oop_data_holder_[i].atom3;
		ATOM* l1 = oop_data_holder_[i].atom3, *l2 = oop_data_holder_[i].atom2, *l3 = oop_data_holder_[i].atom1;
		ATOM* center = oop_data_holder_[i].center;
		//vector3 force_cen(oop_data_holder_[i].center->get_force()), force_1(oop_data_holder_[i].atom1->get_force()), force_2(oop_data_holder_[i].atom2->get_force()), force_3(oop_data_holder_[i].atom3->get_force());
		//get the bond vectors from center atom to the other ones
		vector3 vji_1 = i1->get_position() - center->get_position(), vji_2 = i2->get_position() - center->get_position(),vji_3 = i3->get_position() - center->get_position();
		vector3 vjk_1 = k1->get_position() - center->get_position(), vjk_2 = k2->get_position() - center->get_position(),vjk_3 = k3->get_position() - center->get_position();
		vector3 vjl_1 = l1->get_position() - center->get_position(), vjl_2 = l2->get_position() - center->get_position(),vjl_3 = l3->get_position() - center->get_position();
		double r_ji_1 = vji_1.length(), r_ji_2 = vji_2.length(), r_ji_3 = vji_3.length();
		double r_jk_1 = vjk_1.length(), r_jk_2 = vjk_2.length(), r_jk_3 = vjk_3.length();
		double r_jl_1 = vjl_1.length(), r_jl_2 = vjl_2.length(), r_jl_3 = vjl_3.length();
		//vector3 v21 = oop_data_holder_[i].atom2->get_position() - oop_data_holder_[i].center->get_position();
		//vector3 v31 = oop_data_holder_[i].atom3->get_position() - oop_data_holder_[i].center->get_position();
		//double r_ji = v01.length();
		//double r_jk = v21.length();
		//double r_jl = v31.length();
		if(isNearZero(r_ji_1)|| isNearZero(r_jk_1)|| isNearZero(r_jl_1))
		{
			continue;
		}
		else if(isNearZero(r_ji_2)|| isNearZero(r_jk_2)|| isNearZero(r_jl_2))
		{
			continue;
		}
		else if(isNearZero(r_ji_3)|| isNearZero(r_jk_3)|| isNearZero(r_jl_3))
		{
			continue;
		}
		// normalize v01, v21 ,v31
		vji_1.normalize(), vji_2.normalize(), vji_3.normalize();
		vjk_1.normalize(), vjk_2.normalize(), vjk_3.normalize();
		vjl_1.normalize(), vjl_2.normalize(), vjl_3.normalize();
		//vectors for holding the normal planes
		vector3 an_1 = vji_1 % vjk_1, an_2 = vji_2 % vjk_2, an_3 = vji_3 % vjk_3;
		vector3 bn_1 = vjk_1 % vjl_1, bn_2 = vjk_2 % vjl_2, bn_3 = vjk_3 % vjl_3;
		vector3 cn_1 = vjl_1 % vji_1, cn_2 = vjl_2 % vji_2, cn_3 = vjl_3 % vji_3;
		//bond angles between the two bonds of the plane
		const double cos_theta1 = vji_1 * vjk_1;
		const double cos_theta2 = vjk_2 * vjl_2;
		const double cos_theta3 = vjl_3 * vji_3;
	    const double theta1 = acos(cos_theta1);
		const double theta2 = acos(cos_theta2);
		const double theta3 = acos(cos_theta3);
		double wilson1 = asin(an_1 * vjl_1/sin(theta1));
		double wilson2 = asin(an_2 * vjl_2/sin(theta2));
		double wilson3 = asin(an_3 * vjl_3/sin(theta3));
		if(!isNearZero(theta1) && !isNearZero(fabs(theta1-PI)) && !isNearZero(cos(wilson1)))
		{
			double c1 = -wilson1 * K0 * oop_data_holder_[i].k * FORCE_FACTOR;
			double tmp1 = cos(wilson1) / c1;
			const vector3 d_l = (an_1 / sin(theta1) - vjl_1 * sin(wilson1))/(r_jl_1 * tmp1);
			const vector3 d_i = (((bn_1 + (((-vji_1+vjk_1*cos_theta1)*sin(wilson1))/sin(theta1)))/r_ji_1)/tmp1)/sin(theta1);
			const vector3 d_k = (((cn_1+(((-vjk_1+vji_1*cos_theta1)*sin(wilson1))/sin(theta1)))/r_jk_1)/tmp1)/sin(theta1);
			vector3 force_1 = i1->get_force() + d_i, force_2 = k1->get_force() + d_k, force_3 = l1->get_force() + d_l, force_cen = center->get_force() - d_i - d_k - d_l;
			i1->set_force(force_1);
			k1->set_force(force_2);
			l1->set_force(force_3);
			center->set_force(force_cen);

		}
		if(!isNearZero(theta2) && !isNearZero(fabs(theta2-PI)) && !isNearZero(cos(wilson2)))
		{
			double c1 = -wilson2 * K0 * oop_data_holder_[i].k * FORCE_FACTOR;
			double tmp2 = cos(wilson2) / c1;
			const vector3 d_i=((an_2/sin(theta2)-vjl_2*sin(wilson2))/r_jl_2)/tmp2;
			const vector3 d_l=(((bn_2+(((-vji_2+vjk_2*cos_theta2)*sin(wilson2))/sin(theta2)))/r_ji_2)/tmp2)/sin(theta2);
			const vector3 d_k=(((cn_2+(((-vjk_2+vji_2*cos_theta2)*sin(wilson2))/sin(theta2)))/r_jk_2)/tmp2)/sin(theta2);
			vector3 force_1 = i2->get_force() + d_i, force_2 = k2->get_force() + d_k, force_3 = l2->get_force() + d_l, force_cen = center->get_force() - d_i - d_k - d_l;
			i2->set_force(force_1);
			k2->set_force(force_2);
			l2->set_force(force_3);
			center->set_force(force_cen);
		}
		if(!isNearZero(theta3) && !isNearZero(fabs(theta3-PI)) && !isNearZero(wilson3) && isNearZero(fabs(wilson3-PI)) &&!isNearZero(cos(wilson3)))
		{
			double c1 = -wilson3 * K0 * oop_data_holder_[i].k * FORCE_FACTOR;
			double tmp3 = cos(wilson3) / c1;
			const vector3 d_k=((an_3/sin(theta3)-vjl_3*sin(wilson3))/r_jl_3)/tmp3;
			const vector3 d_i=(((bn_3+(((-vji_3+vjk_3*cos_theta3)*sin(wilson3))/sin(theta3)))/r_ji_3)/tmp3)/sin(theta3);
			const vector3 d_l=(((cn_3+(((-vjk_3+vji_3*cos_theta3)*sin(wilson3))/sin(theta3)))/r_jk_3)/tmp3)/sin(theta3);
			vector3 force_1 = i3->get_force() + d_i, force_2 = k3->get_force() + d_k, force_3 = l3->get_force() + d_l, force_cen = center->get_force() - d_i - d_k - d_l;
			i3->set_force(force_1);
			k3->set_force(force_2);
			l3->set_force(force_3);
			center->set_force(force_cen);
		}
		/*oop_data_holder_[i].atom1->set_force(force_1);
		oop_data_holder_[i].atom2->set_force(force_2);
		oop_data_holder_[i].atom3->set_force(force_3);
		oop_data_holder_[i].center->set_force(force_cen);*/
	}
}