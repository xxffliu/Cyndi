#include "../include/MMFF94Stretch.h"
#include "../include/MMFF94.h"
#include "../include/utility.h"
#include <cmath>
using namespace std;

// constructors
MMFF94Stretch::MMFF94Stretch():FFComponent(),fast_access(),fast_access_emperical(),stretch_data_holder_()
{
	// set component name
	set_name("MMFF94 Stretch");
}
MMFF94Stretch::MMFF94Stretch(ForceField& ff):FFComponent(ff),fast_access(),fast_access_emperical(),stretch_data_holder_()
{
	set_name("MMFF94 Stretch");
}
// copy constructor
MMFF94Stretch::MMFF94Stretch(const MMFF94Stretch& to_copy):FFComponent(to_copy)
{
	fast_access = to_copy.fast_access;
    stretch_data_holder_ = to_copy.stretch_data_holder_;
	fast_access_emperical=to_copy.fast_access_emperical;
}
MMFF94Stretch& MMFF94Stretch::operator =(const MMFF94Stretch &str)
{
	fast_access = str.fast_access;
	fast_access_emperical = str.fast_access_emperical;
	stretch_data_holder_ = str.stretch_data_holder_;
	mmff94_force_field = str.mmff94_force_field;
	num_of_atom_types_ = str.num_of_atom_types_;
	return *this;
}
// destructor
MMFF94Stretch::~MMFF94Stretch()
{
	fast_access.clear();
    stretch_data_holder_.clear();
	fast_access_emperical.clear();
}

// extract bond stretch parameters from FFParameter object bonded to the force field
// and establish a hash table for fast access
bool MMFF94Stretch::extract_BS_parameters(FFParameter& ffp)
{
     if (!ffp.is_valid())
         return false;
     // build a two dim array of atom types and loop variable
     //FFParameter::AtomTypes& atom_types = ffp.get_atomtypes();
     int num_types = ffp.get_num_types();
     num_of_atom_types_ = num_types;
	 //i is the maximun of the CXB value by gaining from the equation CXB = MC * (I * MA + J) + BTij 
     int bond_type,type1, type2;
	 unsigned int CXB;
	 for (vector<vector<string> >::iterator it = force_field_->get_parameters().params_in_each_section["BEND"].begin();it != force_field_->get_parameters().params_in_each_section["BEND"].end(); ++it)
	 {
		 type1=str2int((*it)[0]);
		 type2=str2int((*it)[1]);
		 //53 is the largest atomic number of the atom in this emperical rule table;
		 int index =(type1-1) + (type2-1) * 53;
		 // cout<<"emp_index__"<<index<<endl;
		 fast_access_emperical[index].is_defined=true;
		 fast_access_emperical[index].type1=type1;
		 fast_access_emperical[index].type2=type2;
		 fast_access_emperical[index].r0=str2double((*it)[2]);
		 fast_access_emperical[index].k=str2double((*it)[3]);
		 index =(type2-1) + (type1-1) * 53;
		 // cout<<"emp_index__"<<index<<endl;
		 fast_access_emperical[index].is_defined=true;
		 fast_access_emperical[index].type1=type1;
		 fast_access_emperical[index].type2=type2;
		 fast_access_emperical[index].r0=str2double((*it)[2]);
		 fast_access_emperical[index].k=str2double((*it)[3]);
	 }
	 type1=0;
	 type2=0;
	 //string name_type1, name_type2;  (Cancel!!)
     for (vector<vector<string> >::iterator it = force_field_->get_parameters().params_in_each_section["bond"].begin();it != force_field_->get_parameters().params_in_each_section["bond"].end(); ++it)
	 {
		 bond_type =str2int((*it)[0]);
		 type1 = str2int((*it)[1]);
		 type2 = str2int((*it)[2]);
		 CXB = mmff94_force_field-> GetCXB(bond_type,type1,type2);
		 //CXB equal the returned value by function of GetCXB which is used for caculating the canonical bond index
		 fast_access[CXB].is_defined = true;
		 fast_access[CXB].bond_type =str2int((*it)[0]);
		 fast_access[CXB].type1=str2int((*it)[1]);
		 fast_access[CXB].type2=str2int((*it)[2]);
		 fast_access[CXB].k = str2double((*it)[3]);
		 fast_access[CXB].r0 = str2double((*it)[4]);
		 CXB = mmff94_force_field-> GetCXB(bond_type,type2,type1);
		 fast_access[CXB].is_defined = true;
		 fast_access[CXB].bond_type =str2int((*it)[0]);
		 fast_access[CXB].type1=str2int((*it)[2]);
		 fast_access[CXB].type2=str2int((*it)[1]);
		 fast_access[CXB].k = str2double((*it)[3]);
		 fast_access[CXB].r0 = str2double((*it)[4]);     
     }
	 //stretch_fast_access = fast_access;
	 //stretch_fast_access_emperical = fast_access_emperical;
     return true;
}
  
// query a set of parameters has defined for a given combination of atoms;
bool MMFF94Stretch::has_params(int type1, int type2, int bond_type)
{
	int cxb = mmff94_force_field->GetCXB(bond_type,type1,type2);
	if((type1<0||type1>num_of_atom_types_)||(type2<0||type2>num_of_atom_types_))
	{
		cout<<"Warning: MMFF94Stretch::has_params: Atoms "<<type1<<" "<<type2<<" "<<"haven't been defined in the forcefield"<<endl;
		return  false;
	}
	if(!(bond_type==0||bond_type==1))
		return false;
	return fast_access[cxb].is_defined;
}
     
bool MMFF94Stretch::has_params_emprical(int i,int j)
{
	if(i<0||j>53)
		return false;
	return fast_access_emperical[(i-1)+(j-1)*53].is_defined;
	
}

// assign the parameters to a given atom type combination
// if no parameters are define for this 2 combination, return false and nothing is changed
bool MMFF94Stretch::assign_params(MMFF94Stretch::ForceValues& param, int i, int j, int bond_type)
{
     if (has_params(i, j,bond_type))
	 {
		 param.k = fast_access[mmff94_force_field->GetCXB(bond_type,i,j)].k;
		 param.r0 = fast_access[mmff94_force_field->GetCXB(bond_type,i,j)].r0;
		 return true;
     }
     else
         return false;
}


bool MMFF94Stretch::assign_params_emperical(MMFF94Stretch::ForceValues& param, int i, int j)
{
	if(has_params_emprical(i,j))
	{
		param.k=fast_access_emperical[(i-1)+(j-1)*53].k;
		param.r0=fast_access_emperical[(i-1)+(j-1)*53].r0;
		return true;
	}
	else
		return false;
}

double MMFF94Stretch::GetBondLength(ATOM* a,ATOM* b)
{
	
	int bond_type;
	int cxb,type1,type2;
	type1=a->get_mmff94_type();
	type2=b->get_mmff94_type();
	bond_type=mmff94_force_field->GetBondType(a,b);
	cxb=mmff94_force_field->GetCXB(bond_type,type1,type2);
	if(has_params(type1,type2, bond_type))
	{
		return fast_access[cxb].r0;
	}
	else if(has_params_emprical(a->get_atomic_num(),b->get_atomic_num()))
	{
		return fast_access_emperical[(a->get_atomic_num()-1)+(b->get_atomic_num()-1)*53].r0;
	}
	else 
		return (mmff94_force_field->GetRuleBondLength(a,b));

  
}
// set up method
bool MMFF94Stretch::setup()
{
     if (force_field_ == NULL)
	 {
		 cout<<"MMFF94Stretch::setup(): force field bound component can not be found"<<endl;
		 return false;
	 }
     // clear the parameter holder
     stretch_data_holder_.clear();
     // tempararily set this component enabled;
     setenabled(true);
     
     mmff94_force_field = dynamic_cast<MMFF94*>(force_field_);
     if ((mmff94_force_field == NULL) || mmff94_force_field->has_initialized_param())
	 {
		 bool result = extract_BS_parameters(force_field_->get_parameters());
		 if (!result)
		 {
			 cout << "MMFF94Stretch::setup(): can not found bond stretch section"<<endl;
			 return false;
		 }
	 }                    
                          
     // retrieve all stretch parameters
	 //int swap;
	 //ATOM* Swap;
     for (BONDVec::iterator biter = force_field_->get_bonds().begin();biter != force_field_->get_bonds().end(); ++biter)
	 {
		 //swap=0;
		 MMFF94Stretch::ForceValues value;
		 ATOM* atom_type_A = (*biter)->get_first_atom();
		 ATOM* atom_type_B = (*biter)->get_second_atom();
		 int bond_type =mmff94_force_field->GetBondType(atom_type_A,atom_type_B) ;
		 //cout<<atom_type_A->get_id()<<"-"<<atom_type_B->get_id()<<"-"<<bond_type<<endl;
		 int atom_1=atom_type_A->get_mmff94_type();
		 int atom_2=atom_type_B->get_mmff94_type();              
		 if(assign_params(value, atom_1,atom_2, bond_type))
		 {
			 //cout<<atom_1<<"-"<<atom_2<<": "<<value.k<<" "<<value.r0<<endl;
			 MMFF94Stretch::StretchData data;
			 data.value = value;
			 data.atom1 = atom_type_A;
			 data.atom2 = atom_type_B;
			 data.bond_type = bond_type;
			 stretch_data_holder_.push_back(data);
		 }
		 else if(assign_params_emperical(value,atom_type_A->get_atomic_num(),atom_type_B->get_atomic_num()))
		 {	
			 // if we cannot assign proper type for the bond, using the emperical rules according the pa
			 double rr, rr2, rr4, rr6;
			 MMFF94Stretch::StretchData data;
			 //stretch_data_holder_.push_back(MMFF94Stretch::StretchData());
			 data.atom1 =atom_type_A;
			 data.atom2 = atom_type_B;
			 data.bond_type = bond_type;    
			 value.r0= mmff94_force_field->GetRuleBondLength(atom_type_A, atom_type_B); 
			 rr = fast_access_emperical[(atom_1-1)+(atom_2-1)*53].r0/value.r0;
			 rr2 = rr * rr;
			 rr4 = rr2 * rr2;
			 rr6 = rr4 * rr2;
			 value.k=fast_access_emperical[(atom_type_A->get_atomic_num()-1)+(atom_type_B->get_atomic_num()-1)*53].k*rr6;
			 data.value = value;
			 stretch_data_holder_.push_back(data);
		 }
		 else
		 {
			 cout<<"Error: MMFF94Stretch::setup: We can't assign proper parameters for bond "<<atom_type_A->get_id()<<"-"<<atom_type_B->get_id()<<endl;
			 exit(1);
		 }
	 }
// everything goes well
	 mmff94stretch = *this;
     return true;
}
// update methods
double MMFF94Stretch::update_energy()
{
       //energy initializion
       energy_ = 0.0;
       //debug
       //cout<<stretch_data_holder_.size()<<endl;
       // iterate all bonds and summarize the stretch energies
       for (vector<StretchData>::size_type i = 0; i<stretch_data_holder_.size(); ++i)
	   {
           vector3  r_12 = stretch_data_holder_[i].atom1->get_position()-stretch_data_holder_[i].atom2->get_position();
		   double r12=r_12.length();
           //debug
		   //cout<<stretch_data_holder_[i].atom1->get_mmff94_symbol_type()<<"-"<<stretch_data_holder_[i].atom2->get_mmff94_symbol_type()<<" "<<stretch_data_holder_[i].value.k<<" "<<stretch_data_holder_[i].value.r0<<" ";
		   //delta equal the difference between the actual and reference bond lenghth
		   double delta= r12 - stretch_data_holder_[i].value.r0;
		   double delta2=delta *delta;
		   double kb_ij=stretch_data_holder_[i].value.k;
		   //"cubic-stretch" -cs is a constant which was given -2 A¨B
		   //                   kb_ij                              7
		   // EB_ij = 143.9325 ------- /\r_ij^2 (1 + cs /\_rij + ---- cs^2 r_ij^2)
           //                     2                               12
		   //md*a->Kcal/mol==md*a*143.9325
		   //cout<<143.9325*kb_ij/2*delta2*(1-2*delta+7/3*delta2)<<endl;
           energy_ += 143.9325*kb_ij/2*delta2*(1-2*delta+7/3*delta2);
       }
	   //debug
	   //cout<<"Stretch_energy: "<<energy_<<endl;
       return energy_;
}
// calculate current forces imposed by stretch and add them to the force field;
void MMFF94Stretch::update_forces()
{
	if (get_force_field() == 0)
	{
		cout<<"MMFF94Stretch::update_force(): error! this component doesn't bond to any force field"<<endl;
		return;
	}    
    // iterate all bonds and update forces
    for (vector<StretchData>::size_type i = 0; i<stretch_data_holder_.size(); ++i)
	{
		vector3 force_1(stretch_data_holder_[i].atom1->get_force()), force_2(stretch_data_holder_[i].atom2->get_force());
		vector3 direction(stretch_data_holder_[i].atom1->get_position() - stretch_data_holder_[i].atom2->get_position());
		double distance = direction.length();
		direction.normalize();
		double delta = distance - stretch_data_holder_[i].value.r0;
		double delta2=delta * delta;
		// convert the units from kcal/mol A to N            
        // A -> m: 1e-10
        // J/mol -> J: Avogadro
		double  kb = 143.9325 * FORCE_FACTOR / 2 * stretch_data_holder_[i].value.k * delta;
		direction =direction * kb * (2 - 3 * 2.0 * delta + 28.0/3.0 * delta2);
        force_1 -= direction;
        force_2 += direction;
                 
        stretch_data_holder_[i].atom1->set_force(force_1);
        stretch_data_holder_[i].atom2->set_force(force_2);
	}
}

MMFF94Stretch mmff94stretch;
//map<int,MMFF94Stretch::StretchHashData> stretch_fast_access;
//map<int,MMFF94Stretch::StretchHashData_emprical> stretch_fast_access_emperical;