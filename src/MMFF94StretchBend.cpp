#include "../include/MMFF94StretchBend.h"
#include "../include/Atom.h"
#include "../include/Bond.h"
#include"../include/Mol.h"
#include"../include/MMFF94Stretch.h"
#include"../include/MMFF94Bend.h"
#include <cmath>
#include "../include/utility.h"


//clear function of Str_BendHashData
MMFF94Str_Bend::MMFF94Str_Bend():FFComponent(),fast_access_strbnd(),Str_Bend_data_holder_()
{
	// set component name
    set_name("MMFF94 Str_Bend");
}
MMFF94Str_Bend::MMFF94Str_Bend(ForceField& ff):FFComponent(ff),fast_access_strbnd(),Str_Bend_data_holder_()
{
	set_name("MMFF94 Str_Bend");
}
// copy constructor
MMFF94Str_Bend::MMFF94Str_Bend(const MMFF94Str_Bend& to_copy):FFComponent(to_copy)
{
	fast_access_strbnd = to_copy.fast_access_strbnd;
    Str_Bend_data_holder_ = to_copy.Str_Bend_data_holder_;
}
// destructor
MMFF94Str_Bend::~MMFF94Str_Bend()
{
	fast_access_strbnd.clear();
    Str_Bend_data_holder_.clear();
}
bool MMFF94Str_Bend::extract_AB_parameters(FFParameter& ffp)
{
     if(!ffp.is_valid())
	 {
         //debug
         cout<<"force field parameters is not valid"<<endl;
         return false;
     }
     // build a 2 dim array of atom types and loop variables
     //FFParameter::AtomTypes atom_types = ffp.get_atomtypes();
     int num_types = ffp.get_num_types();
     int num_entries = num_types * num_types * num_types;
     num_of_atom_types_ = num_types;
     //i is the maximun of the CXS value by gaining from the equation  /CXS = MC * (J * MA^2 + I * MA + K) + STijk
     
     //debug
     //cout<<num_types<<endl;
     // start pack the parameters into the vector fast_access_strbnd

     for (vector<vector<string> >::iterator it = force_field_->get_parameters().params_in_each_section["strbnd"].begin(); it != force_field_->get_parameters().params_in_each_section["strbnd"].end(); ++it)
	 {
			 int Str_Bend_type=str2int((*it)[0]);
			 int type1 = str2int((*it)[1]);
			 int type2 = str2int((*it)[2]);
			 int type3 = str2int((*it)[3]);
			 if ((type1 == -1) || (type2 == -1) || (type3 == -1)||(Str_Bend_type<0||Str_Bend_type>11))
			 {
                   cout<<"error! no numeric atom type defined for atom or no defined Str_Bend class"<<endl;
                   return false;
             }
			 int CXS = mmff94_force_field->GetCXS(Str_Bend_type,type1,type2,type3);
             fast_access_strbnd[CXS].is_defined = true;
			 fast_access_strbnd[CXS].strbnd__type =str2int((*it)[0]);
			 fast_access_strbnd[CXS].type1=str2int((*it)[1]);
			 fast_access_strbnd[CXS].type2=str2int((*it)[2]);
			 fast_access_strbnd[CXS].type3=str2int((*it)[3]);
			 fast_access_strbnd[CXS].k_ijk = str2double((*it)[4]);
			 fast_access_strbnd[CXS].k_kji= str2double((*it)[5]);
			 CXS = mmff94_force_field->GetCXS(Str_Bend_type,type3,type2,type1);
             fast_access_strbnd[CXS].is_defined = true;
			 fast_access_strbnd[CXS].strbnd__type =str2int((*it)[0]);
			 fast_access_strbnd[CXS].type1=str2int((*it)[3]);
			 fast_access_strbnd[CXS].type2=str2int((*it)[2]);
			 fast_access_strbnd[CXS].type3=str2int((*it)[1]);
			 fast_access_strbnd[CXS].k_ijk = str2double((*it)[4]);
			 fast_access_strbnd[CXS].k_kji= str2double((*it)[5]);   
	 }

	 for (vector<vector<string> >::iterator it = force_field_->get_parameters().params_in_each_section["S_B"].begin(); it != force_field_->get_parameters().params_in_each_section["S_B"].end(); ++it)
	 {
			 int type1 = str2int((*it)[0]);
			 int type2 = str2int((*it)[1]);
			 int type3 = str2int((*it)[2]);
			 if ((type1 == -1) || (type2 == -1) || (type3 == -1))
			 {
                   cout<<"error! no numeric atom type defined for atom or no defined S_B class"<<endl;
                   return false;
             }
			 int index = type1 + type2*4 + type3*4*4;
             emp_fast_acess_sb[index].is_defined = true;
			 emp_fast_acess_sb[index].strbnd__type =0;
			 emp_fast_acess_sb[index].type1=str2int((*it)[0]);
			 emp_fast_acess_sb[index].type2=str2int((*it)[1]);
			 emp_fast_acess_sb[index].type3=str2int((*it)[2]);
			 emp_fast_acess_sb[index].k_ijk = str2double((*it)[3]);
			 emp_fast_acess_sb[index].k_kji= str2double((*it)[4]);
			 index = type3 + type2*4 + type1*4*4;
             emp_fast_acess_sb[index].is_defined = true;
			 emp_fast_acess_sb[index].strbnd__type =0;
			 emp_fast_acess_sb[index].type1=str2int((*it)[2]);
			 emp_fast_acess_sb[index].type2=str2int((*it)[1]);
			 emp_fast_acess_sb[index].type3=str2int((*it)[0]);
			 emp_fast_acess_sb[index].k_ijk = str2double((*it)[3]);
			 emp_fast_acess_sb[index].k_kji= str2double((*it)[4]);  
	 }
	 return true;
}

// query a set of Str_Bend Str_Bend parameters has defined for a given combination of atom types
// query a set of parameters has defined for a given combination of atoms;
bool MMFF94Str_Bend::has_params(int type1, int type2,int type3, int Str_Bend_type)
{
	if((type1<0||type1>num_of_atom_types_)||(type2<0||type2>num_of_atom_types_)||(type3<0||type3>num_of_atom_types_))
	{
		cout<<"Warning: MMFF94Str_Bend::has_params: Atoms "<<type1<<" "<<type2<<" "<<type3<<" "<<"haven't been defined in the forcefield"<<endl;
		return  false;
	}
	if(Str_Bend_type<0||Str_Bend_type>11)
		return false;
	int cxs = mmff94_force_field->GetCXS(Str_Bend_type,type1,type2,type3);
	return fast_access_strbnd[cxs].is_defined;
}
     
bool MMFF94Str_Bend::has_params_emperical(int i,int j,int k)
{
	if((i<0||i>4)||(j<0||j>4)||(k<0||k>4))
		return false;
	int index = i + j*4 + k*4*4;
	return emp_fast_acess_sb[index].is_defined;
}

// assign the parameters to a given atom type combination
// if no parameters are define for this 2 combination, return false and nothing is changed
bool  MMFF94Str_Bend::assign_params(MMFF94Str_Bend::ForceValues& param, int i, int j, int k ,int Str_Bend_type) 
{
	
     if (has_params(i,j,k,Str_Bend_type))
	 {
		 param.k_ijk= fast_access_strbnd[mmff94_force_field->GetCXS(Str_Bend_type,i,j,k)].k_ijk;
		 param.k_kji=fast_access_strbnd[mmff94_force_field->GetCXS(Str_Bend_type,i,j,k)].k_kji;
		 //param.is_defined=fast_access_strbnd[mmff94.GetCXS(Str_Bend_type,i,j,k)].is_defined;
		 return true;
	 }
	 return false;
}
bool MMFF94Str_Bend::assign_params_emprical(MMFF94Str_Bend::ForceValues& param, int i, int j, int k)
{
	if(has_params_emperical( i, j, k))
	{
		int index = i+j*4+k*4*4;
		param.k_ijk = emp_fast_acess_sb[index].k_ijk;
		param.k_kji = emp_fast_acess_sb[index].k_kji;
		return true;
	}
	return false;

}

bool  MMFF94Str_Bend::setup()
{
	if (get_force_field() == 0)
	{
		cout<<" MMFF94Str_Bend::setup(): force field bound component can not be found"<<endl;
		return false;
	}
     // clear the parameter container
    Str_Bend_data_holder_.clear();
    // tempararily set this component enabled;
    setenabled(true);
    mmff94_force_field = dynamic_cast< MMFF94*>(force_field_);
    if (mmff94_force_field==0 || mmff94_force_field->has_initialized_param())
	{
		bool result = extract_AB_parameters(get_force_field()->get_parameters());
		if (!result)
		{
			cout << " MMFF94Str_Bend::setup(): can not found Str_Bend Str_Bend section"<<endl;
			return false;
		}
	}
    BONDVec::iterator biter, biter2;
    for (ATOMVec::iterator aiter = force_field_->get_atoms().begin(); aiter != force_field_->get_atoms().end(); ++aiter)
	{
		if((*aiter)->get_symbol_type() == "H")
			continue;
		for (biter = (*aiter)->get_bond_list().begin(); biter != (*aiter)->get_bond_list().end(); ++biter)
		{
			for (biter2 = biter + 1; biter2 != (*aiter)->get_bond_list().end(); ++biter2)
			{
				//cout<<(*biter)->get_first_atom()<<" "<<(*biter)->get_second_atom()<<endl;
				int atom_type_A = (*biter)->get_partner(*aiter)->get_mmff94_type();
				int atom_type_B = (*aiter)->get_mmff94_type();
				int atom_type_C = (*biter2)->get_partner(*aiter)->get_mmff94_type();
				MMFF94Str_Bend::ForceValues value;
				MMFF94Str_Bend::Str_BendData data;
				bool is_defined = false;
				int _Str_Bend_type = mmff94_force_field->GetStrBndType((*biter)->get_partner(*aiter),(*aiter),(*biter2)->get_partner(*aiter));
				if(assign_params(value, atom_type_A, atom_type_B, atom_type_C,_Str_Bend_type))
				{
					//int angle_type;
					data.value = value;
					data.atom1 = (*biter)->get_partner(*aiter);
					data.atom2 = *aiter;
					data.atom3 = (*biter2)->get_partner(*aiter);
					//angle_type=mmff94.GetAngleType(data.atom1,data.atom2,data.atom3);
					data.r_ij=mmff94stretch.GetBondLength(data.atom1,data.atom2);
					data.r_kj=mmff94stretch.GetBondLength(data.atom3,data.atom2);
					data.theta0=mmff94bend.GetBondAngle(data.atom1, data.atom2, data.atom3);
					data.strbnd_type = _Str_Bend_type;
					is_defined = true;
				}
				else
				{
					int rowa, rowb, rowc;
					int bond_type1,bond_type2;
					ATOM *atom_i = (*biter)->get_partner(*aiter);
					ATOM *atom_j = *aiter;
					ATOM *atom_k = (*biter2)->get_partner(*aiter);
	                bond_type1=mmff94_force_field->GetBondType(atom_i,atom_j);
		            bond_type2=mmff94_force_field->GetBondType(atom_k,atom_j);
		            //angle_type=mmff94.GetAngleType(atom_i,atom_j,atom_k);
  
                    // This is not a real empirical rule...
                    rowa = GetElementRow(atom_i);
                    rowb = GetElementRow(atom_j);
                    rowc = GetElementRow(atom_k);
					if (assign_params_emprical(value, rowa, rowb, rowc))
					{
						//int angle_type;
						//angle_type=mmff94.GetAngleType(atom_i,atom_j,atom_k);
						data.r_ij=mmff94stretch.GetBondLength(atom_i,atom_j);
						data.r_kj=mmff94stretch.GetBondLength(atom_k,atom_j);
						data.theta0=mmff94bend.GetBondAngle(atom_i, atom_j, atom_k);
						data.strbnd_type = _Str_Bend_type;
						data.value = value;
						data.atom1 = atom_i;
						data.atom2 = atom_j;
						data.atom3 = atom_k;
						is_defined = true;
					}
				}
				if(is_defined)
				{
					Str_Bend_data_holder_.push_back(data);
				}
				else
				{
					cout<<"Warning: MMFF94Str_Bend::setup: We can't assign proper stretch-bend parameters for angle "<<(*aiter)->get_id()<<"-"<<(*biter)->get_partner(*aiter)->get_id()<<"-"<<(*biter2)->get_partner(*aiter)->get_id()<<endl;
				}
			}
		}
	}
// everything goes well
    return true;
}

// update methods
double MMFF94Str_Bend::update_energy()
{
       //energy initializion
       energy_ = 0.0;
       if (Str_Bend_data_holder_.size() == 0)
           return 0.0;
       // iterate all Str_Bends and summarize the Str_Bend energies
       for (vector<Str_BendData>::size_type i = 0; i<Str_Bend_data_holder_.size(); ++i)
	   {
           //cout<<Str_Bend_data_holder_.size()<<endl;
           vector3 v1, v2;
		   double r_1,r_2;
           v1 = Str_Bend_data_holder_[i].atom1->get_position() - Str_Bend_data_holder_[i].atom2->get_position();
           v2 = Str_Bend_data_holder_[i].atom3->get_position() - Str_Bend_data_holder_[i].atom2->get_position();
           
		   r_1=v1.length();
		   r_2=v2.length();
           
           if((v1.length()<0.0001) && (v2.length()<0.0001))
		   {
			   cout<<v1.length()<<" "<<v2.length()<<endl;
               continue;
           }
		   double costheta = ((v1*v2)/(r_1*r_2));
		   double theta;
		   if (costheta > 1.0)
			   costheta = 1.0;
		   else if (costheta < -1.0) 
			   costheta = -1.0;
		   theta = acos(costheta)*RAD_TO_DEG;
		   int mmff94_type=Str_Bend_data_holder_[i].atom2->get_mmff94_type();
		   double ka_ijk;
		   double ka_kji;
		   ka_ijk=Str_Bend_data_holder_[i].value.k_ijk;
		   ka_kji=Str_Bend_data_holder_[i].value.k_kji;
		   double theta0;
		   theta0=Str_Bend_data_holder_[i].theta0;
		   double r_ij0,r_kj0;
		   
		   r_ij0=Str_Bend_data_holder_[i].r_ij;
		   r_kj0=Str_Bend_data_holder_[i].r_kj;
		   
		   double delta_r1,delta_r2,delta_theta;
		   delta_r1=r_1-r_ij0;
		   delta_r2=r_2-r_kj0;
		   delta_theta=theta-theta0;
		   if (MMFF94_LIN[mmff94_type-1])
		   {
			   energy_+= 0.0;
		   }
		   else 
		   {
			   energy_ += 2.51210 * (ka_ijk * delta_r1 + ka_kji * delta_r2) * delta_theta;
		   }
	   }
	   //cout<<"stretch_bend_energy:"<<energy_<<endl;
	   return energy_;
	   
}

// calculate current forces imposed by Str_Bend and add them to the force field;
void MMFF94Str_Bend::update_forces()
{
	if (get_force_field() == 0)
	{
		cout<<"MMFF94Str_Bend::update_force(): error! this component doesn't bond to any force field"<<endl;
		return;
	}
    // iterate all Str_Bends and update forces
	for (vector<Str_BendData>::size_type i = 0; i < Str_Bend_data_holder_.size(); ++i)
	{
		vector3 force_1 = Str_Bend_data_holder_[i].atom1->get_force(), force_2 = Str_Bend_data_holder_[i].atom2->get_force(), force_3 = Str_Bend_data_holder_[i].atom3->get_force();
        // calculate the vectors between atom1 and atom2, atom3 and atom2 then normalize them
        vector3 v1 = Str_Bend_data_holder_[i].atom1->get_position() - Str_Bend_data_holder_[i].atom2->get_position();
        vector3 v2 = Str_Bend_data_holder_[i].atom3->get_position() - Str_Bend_data_holder_[i].atom2->get_position();
		int mmff94_type=Str_Bend_data_holder_[i].atom2->get_mmff94_type();
        double length = v1.length();
        if (isNearZero(length))
			continue;
        double inverse_v1 = 1/length;
        v1.normalize();
        length = v2.length();
        if (isNearZero(length))
			continue;
        double inverse_v2 = 1/length;
        v2.normalize();
		vector3 cross_product = cross(v1, v2);
		double cross_length = cross_product.length();
		if(isNearZero(cross_length))
			continue;
		cross_product.normalize();
        double costheta = v1 * v2;
		double theta;
		if (costheta > 1.0) 
			costheta = 1.0;
		else if (costheta < -1.0) 
			costheta = -1.0;
        theta = acos(costheta)*RAD_TO_DEG;

		double delta = theta - Str_Bend_data_holder_[i].theta0;
		double d_ij = 1/inverse_v1 - Str_Bend_data_holder_[i].r_ij;
		double d_kj = 1/inverse_v2 - Str_Bend_data_holder_[i].r_kj;

		if (!MMFF94_LIN[mmff94_type-1])
		{
			vector3 r1(VZero), r2(VZero);
			double s1 = Str_Bend_data_holder_[i].value.k_ijk * delta * DEG_TO_RAD;
			double s2 = Str_Bend_data_holder_[i].value.k_kji * delta * DEG_TO_RAD;
			r1 = s1 * v1;
			r2 = s2 * v2;
			double scal = -(Str_Bend_data_holder_[i].value.k_ijk * d_ij + Str_Bend_data_holder_[i].value.k_kji * d_kj);
			vector3 t1 = v1 % cross_product;
			t1.normalize();
			vector3 t2 = v2 % cross_product;
			t2.normalize();
			vector3 n1 = -t1 * inverse_v1;
			vector3 n2 = t2 * inverse_v2;
			r1 += n1 * scal;
			r2 += n2 * scal;
			force_1 -= r1 * 2.51210 * FORCE_FACTOR;
			force_2 += (r1+r2) * 2.51210 * FORCE_FACTOR;
			force_3 -= r2 * 2.51210 * FORCE_FACTOR;
			Str_Bend_data_holder_[i].atom1->set_force(force_1);
            Str_Bend_data_holder_[i].atom2->set_force(force_2);
			Str_Bend_data_holder_[i].atom3->set_force(force_3);
		}
	}
}
  int MMFF94Str_Bend::GetElementRow(ATOM *atom)
  {
    int row;
    
    row = 0;

	if (atom->get_atomic_num() > 2)
      row++;
    if (atom->get_atomic_num() > 10)
      row++;
    if (atom->get_atomic_num() > 18)
      row++;
    if (atom->get_atomic_num() > 36)
      row++;
    if (atom->get_atomic_num() > 54)
      row++;
    if (atom->get_atomic_num() > 86)
      row++;
    
    return row;
  } 
  MMFF94Str_Bend mmff94stretchbend;
