#include "../include/MMFF94Bend.h"
#include "../include/Atom.h"
#include "../include/Bond.h"
#include"../include/Mol.h"
#include"../include/MMFF94Stretch.h"
#include <cmath>
#include "../include/utility.h"

/*#ifndef AVOGADRO
#define AVOGADRO  6.0221367E+23L;
#endif
*/

//using namespace std;
//clear function of BendHashData
MMFF94Bend::MMFF94Bend():FFComponent(),fast_access(),Equal_access(),bend_data_holder_()
{
	// set component name
    set_name("MMFF94 Bend");
}
MMFF94Bend::MMFF94Bend(ForceField& ff):FFComponent(ff),fast_access(),Equal_access(),bend_data_holder_()
{
	set_name("MMFF94 Bend");
}
// copy constructor
MMFF94Bend::MMFF94Bend(const MMFF94Bend& to_copy):FFComponent(to_copy)
{
	fast_access = to_copy.fast_access;
	Equal_access= to_copy.Equal_access;
    bend_data_holder_ = to_copy.bend_data_holder_;
}
// assignment operator
MMFF94Bend& MMFF94Bend::operator =(const MMFF94Bend &bend)
{
	fast_access = bend.fast_access;
	Equal_access = bend.Equal_access;
	bend_data_holder_ = bend.bend_data_holder_;
	num_of_atom_types_ = bend.num_of_atom_types_;
	mmff94_force_field = bend.mmff94_force_field;
	return *this;
}
// destructor
MMFF94Bend::~MMFF94Bend()
{
	fast_access.clear();
	Equal_access.clear();
	bend_data_holder_.clear();
}
bool MMFF94Bend::extract_AB_parameters(FFParameter& ffp)
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
     num_of_atom_types_ = num_types;
     //debug
     //cout<<num_types<<endl;
     // start pack the parameters into the vector fast_access
     int type1, type2, type3, angle_type;
	 unsigned int CXA;
	 int index;
	 
     for (vector<vector<string> >::iterator it = force_field_->get_parameters().params_in_each_section["DEF"].begin(); it != force_field_->get_parameters().params_in_each_section["DEF"].end(); ++it)
	 {
		 index=str2int((*it)[0]);
		 Equal_access[index].equl2=str2int((*it)[1]);
		 Equal_access[index].equl3=str2int((*it)[2]);
		 Equal_access[index].equl4=str2int((*it)[3]);
		 Equal_access[index].equl5=str2int((*it)[4]);	
			 //debug
			 //cout<<Equal_access[index].equl2<<":"<<Equal_access[index].equl3<<";"<<Equal_access[index].equl4<<";"<<Equal_access[index].equl5<<endl;
	 }

     for (vector<vector<string> >::iterator it = force_field_->get_parameters().params_in_each_section["angle"].begin(); it != force_field_->get_parameters().params_in_each_section["angle"].end(); ++it)
	 {
		 angle_type=str2int((*it)[0]);
		 type1 = str2int((*it)[1]);
		 type2 = str2int((*it)[2]);
		 type3 = str2int((*it)[3]);
		 if ((type1 == -1) || (type2 == -1) || (type3 == -1)||(angle_type<0||angle_type>8))
		 {
			 cout<<"error! no numeric atom type defined for atom or no defined angle class"<<endl;
             return false;
		 }
		 CXA = mmff94_force_field->GetCXA(angle_type,type1,type2,type3);
         fast_access[CXA].is_defined = true;
		 fast_access[CXA].angle_type =str2int((*it)[0]);
		 fast_access[CXA].type1=str2int((*it)[1]);
		 fast_access[CXA].type2=str2int((*it)[2]);
		 fast_access[CXA].type3=str2int((*it)[3]);
		 fast_access[CXA].theta0 = str2double((*it)[5]);
         fast_access[CXA].k= str2double((*it)[4]);
		 CXA = mmff94_force_field->GetCXA(angle_type,type3,type2,type1);
         fast_access[CXA].is_defined = true;
		 fast_access[CXA].angle_type =str2int((*it)[0]);
		 fast_access[CXA].type1=str2int((*it)[3]);
		 fast_access[CXA].type2=str2int((*it)[2]);
		 fast_access[CXA].type3=str2int((*it)[1]);
		 fast_access[CXA].theta0 = str2double((*it)[5]);
         fast_access[CXA].k= str2double((*it)[4]);               
         }
	 return true;
}

// query a set of angle bend parameters has defined for a given combination of atom types
// query a set of parameters has defined for a given combination of atoms;
bool MMFF94Bend::has_params(int type1, int type2,int type3, int angle_type)
{

	if((type1<0||type1>num_of_atom_types_)||(type2<0||type2>num_of_atom_types_)||(type3<0||type3>num_of_atom_types_))
	{
		return  false;
		cout<<"Warning: MMFF94Bend::has_params: These atoms haven't defined in the project"<<endl;
	}
	if(angle_type<0||angle_type>8)
		return false;
	int cxa = mmff94_force_field->GetCXA(angle_type,type1,type2,type3);
	return fast_access[cxa].is_defined;
}				 

// assign the parameters to a given atom type combination
// if no parameters are define for this 2 combination, return false and nothing is changed
bool  MMFF94Bend::assign_params(MMFF94Bend::ForceValues& param, int i, int j, int k ,int angle_type)
{
	if((has_params(i,j,k ,angle_type)))
	{
		int cxa = mmff94_force_field->GetCXA(angle_type,i,j,k);
		param.k = fast_access[cxa].k;
		param.theta0 = fast_access[cxa].theta0;
		return true;
	}
	else 
		return false;
		

}

// get refernece bond angle
double MMFF94Bend::GetBondAngle(ATOM* a, ATOM* b, ATOM* c)
{
	int type1 = a->get_mmff94_type(), type2 = b->get_mmff94_type(), type3 = c->get_mmff94_type();
	int angle_type = mmff94_force_field->GetAngleType(a, b, c);
	int cxa = mmff94_force_field->GetCXA(angle_type, type1, type2, type3);
	if(has_params(type1, type2, type3, angle_type))
	{
		return fast_access[cxa].theta0;
	}
	else
		return 120.0;
}
// set up method
bool  MMFF94Bend::setup()
{
	if (get_force_field() == 0)
	{
		cout<<" MMFF94Bend::setup(): force field bound component can not be found"<<endl;
		return false;
	}
	// clear the parameter container
	bend_data_holder_.clear();
	// tempararily set this component enabled;
	setenabled(true);

	mmff94_force_field = dynamic_cast< MMFF94*>(force_field_);
	if (mmff94_force_field== 0 || mmff94_force_field->has_initialized_param())
	{
		bool result = extract_AB_parameters(get_force_field()->get_parameters());
		if (!result)
		{
			cout << " MMFF94Bend::setup(): can not found angle bend section"<<endl;
			return false;
		}
	}
	//mmff94stretch.fast_access = stretch_fast_access;
	//mmff94stretch.fast_access_emperical = stretch_fast_access_emperical;
	//mmff94stretch.set_num_atom_types(num_of_atom_types_);
	// retrieve all bend parameters
	BONDVec::iterator biter, biter2;
	for (ATOMVec::iterator aiter = force_field_->get_atoms().begin();aiter != force_field_->get_atoms().end(); ++aiter)
	{
		if((*aiter)->get_symbol_type() == "H")
			continue;
		for (biter = (*aiter)->get_bond_list().begin(); biter != (*aiter)->get_bond_list().end(); ++biter)
		{
			for (biter2 = biter, ++biter2; biter2 != (*aiter)->get_bond_list().end(); ++biter2)
			{
				MMFF94Bend::ForceValues value;
				ATOM* atom1 = (*biter)->get_partner(*aiter);
				ATOM* atom2 = *aiter;
				ATOM* atom3 = (*biter2)->get_partner(*aiter);
				int atom_type_A = atom1->get_mmff94_type();
				int atom_type_B = atom2->get_mmff94_type();
				int atom_type_C = atom3->get_mmff94_type();
				int _angle_type = mmff94_force_field->GetAngleType((*biter)->get_partner(*aiter),(*aiter),(*biter2)->get_partner(*aiter));
				//bend_data_holder_.push_back(MMFF94Bend::BendData());
				//angle_type=_angle_type;
				if (has_params(atom_type_A, atom_type_B,atom_type_C,_angle_type))
				{
					value.k = fast_access[mmff94_force_field->GetCXA(_angle_type,atom_type_A, atom_type_B,atom_type_C)].k;
					value.theta0= fast_access[mmff94_force_field->GetCXA(_angle_type,atom_type_A, atom_type_B,atom_type_C)].theta0;
				}
				else if(has_params(EqLvl2(atom_type_A),atom_type_B,EqLvl2(atom_type_C),_angle_type))
				{
                    value.k = fast_access[mmff94_force_field->GetCXA(_angle_type,EqLvl2(atom_type_A),atom_type_B,EqLvl2(atom_type_C))].k;		
					value.theta0 = fast_access[mmff94_force_field->GetCXA(_angle_type,EqLvl2(atom_type_A),atom_type_B,EqLvl2(atom_type_C))].theta0;		
				}
				else if (has_params(EqLvl3(atom_type_A),atom_type_B,EqLvl3(atom_type_C),_angle_type))
				{
					value.k = fast_access[mmff94_force_field->GetCXA(_angle_type,EqLvl3(atom_type_A),atom_type_B,EqLvl3(atom_type_C))].k;		
					value.theta0 = fast_access[mmff94_force_field->GetCXA(_angle_type,EqLvl3(atom_type_A),atom_type_B,EqLvl3(atom_type_C))].theta0;		
				}
				else if (has_params(EqLvl4(atom_type_A),atom_type_B,EqLvl4(atom_type_C),_angle_type))
				{
					value.k  = fast_access[mmff94_force_field->GetCXA(_angle_type,EqLvl4(atom_type_A),atom_type_B,EqLvl4(atom_type_C))].k;
					value.theta0 = fast_access[mmff94_force_field->GetCXA(_angle_type,EqLvl4(atom_type_A),atom_type_B,EqLvl4(atom_type_C))].theta0;
				}
				else if (has_params(EqLvl5(atom_type_A),atom_type_B,EqLvl5(atom_type_C),_angle_type))
				{
					value.k = fast_access[mmff94_force_field->GetCXA(_angle_type,EqLvl5(atom_type_A),atom_type_B,EqLvl5(atom_type_C))].k;
					value.theta0= fast_access[mmff94_force_field->GetCXA(_angle_type,EqLvl5(atom_type_A),atom_type_B,EqLvl5(atom_type_C))].theta0;
				}
				else
				{
					value.k = 0.0;
					value.theta0 = 120.0;
					if (mmff94_force_field->GetCrd(atom_type_B) == 4)
						value.theta0 = 109.45;
					if ((mmff94_force_field->GetCrd(atom_type_B) == 2) &&((*aiter)->is_oxygen()))
						value.theta0= 105.0;
					if ((*aiter)->get_atomic_num() > 10)
						value.theta0 = 95.0;
					if (mmff94_force_field->HasLinSet(atom_type_B))
						value.theta0 = 180.0;
					if ((mmff94_force_field->GetCrd(atom_type_B) == 3) && (mmff94_force_field->GetVal(atom_type_B) == 3) && !(mmff94_force_field->GetMltb(atom_type_B)))
					{
						if ((*aiter)->is_nitrogen())
						{
							value.theta0= 107.0;
						}
						else
						{
							value.theta0 = 92.0;
						}
					}
					if (force_field_->get_mol()->IsInRingSize(atom1 ,3) && force_field_->get_mol()->IsInRingSize(atom2,3) && force_field_->get_mol()->IsInRingSize(atom3,3) && force_field_->In_the_sameRing(atom1 ,atom2))
						value.theta0 = 60.0;
					if (force_field_->get_mol()->IsInRingSize(atom1 ,4) && force_field_->get_mol()->IsInRingSize(atom2,4) && force_field_->get_mol()->IsInRingSize(atom3,4) && force_field_->In_the_sameRing(atom1 ,atom2))
						value.theta0= 90.0;    
				}
				//emperical rules
				if (value.k==0.0 )
				{
					double beta, Za, Zc, Cb, r0ab, r0bc, theta, theta2, D, rr, rr2;
					Za =mmff94_force_field->GetZParam(atom1);
					Cb =mmff94_force_field->GetZParam(atom2);
					Zc =mmff94_force_field->GetZParam(atom3);
					int bond_type_ab=mmff94_force_field->GetBondType(atom1 ,atom2);
					int bond_type_bc=mmff94_force_field->GetBondType(atom2, atom3);
					unsigned int cxb_ab=mmff94_force_field->GetCXB(bond_type_ab,atom_type_A,atom_type_B);
					r0ab=mmff94stretch.GetBondLength(atom1,atom2);
					r0bc=mmff94stretch.GetBondLength(atom2,atom3);
					rr = r0ab + r0bc;
					rr2 = rr * rr;
					D = (r0ab - r0bc) / rr2;
					theta =  value.theta0;
					theta2 = theta * theta;
					beta = 1.75;
					if (force_field_->get_mol()->IsInRingSize(atom1,4)&& force_field_->get_mol()->IsInRingSize(atom2,4) && force_field_->get_mol()->IsInRingSize(atom3,4) && mmff94_force_field->In_the_sameRing(atom1,atom3))
						beta = 0.85 * beta;
					if (force_field_->get_mol()->IsInRingSize(atom1,3) && force_field_->get_mol()->IsInRingSize(atom2,3)&& force_field_->get_mol()->IsInRingSize(atom3,3)&& mmff94_force_field->In_the_sameRing(atom1,atom3))
						beta = 0.05 * beta;
					value.k = (beta * Za * Cb * Zc * exp(-2 * D)) / (rr * theta2);
				}
				MMFF94Bend::BendData data;
				data.atom1 = atom1;
				data.atom2 = atom2;
				data.atom3 = atom3;
				data.value = value;
				data.angle_type = _angle_type;
				bend_data_holder_.push_back(data);
			}
		}
	}
	// everything goes well
	mmff94bend = *this;
	return true;
}

// update methods
double MMFF94Bend::update_energy()
{
       //energy initializion
       energy_ = 0.0;
       if (bend_data_holder_.size() == 0)
           return 0.0;
       // iterate all angles and summarize the Bend energies
       for (vector<BendData>::size_type i = 0; i<bend_data_holder_.size(); ++i)
	   {
           //cout<<bend_data_holder_.size()<<endl;
           vector3 v1, v2;
		   double v1_length,v2_length;
           v1 = bend_data_holder_[i].atom1->get_position() - bend_data_holder_[i].atom2->get_position();
           v2 = bend_data_holder_[i].atom3->get_position() - bend_data_holder_[i].atom2->get_position();
           
		   v1_length=v1.length();
		   v2_length=v2.length();
           if((v1.length()<1e-3) && (v2.length()<1e-3))
		   {
                               cout<<v1.length()<<" "<<v2.length()<<endl;
                               exit(1);
		   }
		   double costheta = (v1*v2)/(v1_length*v2_length);
		   double theta;
		   if (costheta > 1.0)
			   costheta = 1.0;
		   else if (costheta < -1.0) 
			   costheta = -1.0;
		   theta = acos(costheta)*180.0/PI;
		   int mmff94_type=bend_data_holder_[i].atom2->get_mmff94_type();
		   double ka;
		   ka=bend_data_holder_[i].value.k;
		   double theta0;
		   theta0=bend_data_holder_[i].value.theta0;
		   if (MMFF94_LIN[mmff94_type-1])
		   {
			   //md*a->KJ/mol==md*a*143.9325
			   energy_+= 143.9325 * ka * (1.0 + cos(theta * DEG_TO_RAD));
		   }
		   else 
		   {
			   double delta,delta2;
			   delta=theta-theta0;
			   delta2=delta*delta;
			   energy_+= 0.043844 * 0.5 * ka * delta2 *(1.0 - 0.007 * delta);
		   }

	   }
	   //debug
	   //cout<<"angle_bend_energy:"<<energy_<<endl;
	   return energy_;
}

// calculate current forces imposed by bend and add them to the force field;
void MMFF94Bend::update_forces()
{
	if (get_force_field() == 0)
	{
		cout<<"MMFF94Bend::update_force(): error! this component doesn't bond to any force field"<<endl;
		return;
	}
    // iterate all angles and update forces
	for (vector<BendData>::size_type i = 0; i < bend_data_holder_.size(); ++i)
	{
		vector3 force_1 = bend_data_holder_[i].atom1->get_force(), force_2 = bend_data_holder_[i].atom2->get_force(), force_3 = bend_data_holder_[i].atom3->get_force();
        // calculate the vectors between atom1 and atom2, atom3 and atom2 then normalize them
        vector3 v1 = bend_data_holder_[i].atom1->get_position() - bend_data_holder_[i].atom2->get_position();
        vector3 v2 = bend_data_holder_[i].atom3->get_position() - bend_data_holder_[i].atom2->get_position();
		
		int mmff94_type=bend_data_holder_[i].atom2->get_mmff94_type();
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
	
		vector3 cross = v1 % v2;
		double cross_length = cross.length();
		if(isNearZero(cross_length))
			continue;
        double angle = vectorAngle(v1, v2);
        double costheta = v1 * v2;
		if (costheta > 1.0) 
			costheta = 1.0;
		else if (costheta < -1.0) 
			costheta = -1.0;
        double theta = acos(costheta) * RAD_TO_DEG;
		vector3 n1(VZero),n2(VZero);

		double ka = bend_data_holder_[i].value.k * RAD_TO_DEG;
	 	if (MMFF94_LIN[mmff94_type-1])
		{
			n1 = (v2 - v1 * costheta) * inverse_v1;
			n2 = (v1 + v2 * costheta) * inverse_v2;
		}
		double delta = theta - bend_data_holder_[i].value.theta0;
		if(!MMFF94_LIN[mmff94_type-1])
		{
			vector3 t1 = v1 % cross;
			t1.normalize();
			vector3 t2 = v2 % cross ;
			t2.normalize();
			n1 =  -t1 * inverse_v1;
			n2 =  t2 * inverse_v2;
		}
		double factor(0.);
		// convert the units from kcal/mol A to N
        // kcal -> J: 1e3 * 4.2
        // A -> m: 1e-10
        // J/mol -> J: Avogadro 
		if(MMFF94_LIN[mmff94_type-1])
		{
			factor= -143.9325 * FORCE_FACTOR * ka;
		} 
		else
		{
			factor = 0.043844/2 * FORCE_FACTOR * ka * delta * (2.0 + 3.0 * 0.007 * delta);
		}
		n1 *= factor;
		n2 *= factor;
        force_1 += n1;
        force_3 += n2;
        force_2 -= (n1 + n2);
        bend_data_holder_[i].atom1->set_force(force_1);
        bend_data_holder_[i].atom2->set_force(force_2);
        bend_data_holder_[i].atom3->set_force(force_3);
	}
	return;
}
                      
 // MMFF part I - TABLE IV
  int MMFF94Bend ::EqLvl2(int type)
  {
	  if(type>0&&type<100)
		  return (Equal_access[type].equl2);
	  else {
	//	  cout<<"the atom type is not exist!"<<endl;
		  return type;
	  }
	    
  }
  
  // MMFF part I - TABLE IV
  int MMFF94Bend ::EqLvl3(int type)
  {
	  if(type>0&&type<100){
	  //debug
		 // cout<<Equal_access[type].equl3<<endl;
		  return  (Equal_access[type].equl3);}
	  else{
		  cout<<"the atom type is not exist!"<<endl;
		  return type;
	  }
  }
  
  // MMFF part I - TABLE IV
  int MMFF94Bend ::EqLvl4(int type)
  {
	  if(type>0&&type<100){
		  //debug
		//  cout<<Equal_access[type].equl4<<endl;
		  return (Equal_access[type].equl4);}
	  else{
		  cout<<"the atom type is not exist!"<<endl;
		  return type;
	  }
   
  }

  // MMFF part I - TABLE IV
  int MMFF94Bend ::EqLvl5(int type)
  {
	  if(type>0&&type<100){
		  //debug
		//  cout<<Equal_access[type].equl5<<endl;
		  return (Equal_access[type].equl5);}
	  else{
		  cout<<"the atom type is not exist!"<<endl;
		  return type;
	  }
  }
  MMFF94Bend mmff94bend;
  //map<int,MMFF94Bend::BendHashData> bend_fast_access;