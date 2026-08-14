#include "../include/TAFF.h"
#include "../include/TAFFStretch.h"
#include "../include/TAFFBend.h"
#include "../include/TAFFOOP.h"
#include "../include/TAFFTorsion.h"
#include "../include/TAFFVDW.h"
#include "../include/TAFFElectrostatic.h"
#include <sstream>
using namespace std;

// default constructor
TAFF::TAFF():
param_is_initialized_(false),
file_name_(DEFAULT_TAFF_PARAM_FILE)
{
	// set force name
	set_ff_name("TAFF ["+file_name_+"]");
	//creat the component list
	insert_component(new TAFFStretch(*this));
	insert_component(new TAFFBend(*this));
	insert_component(new TAFFOOP(*this));
	insert_component(new TAFFTorsion(*this));
	insert_component(new TAFFVDW(*this));
	insert_component(new TAFFEle(*this));
}
// constructor with a MOL
TAFF::TAFF(MOL& mol):ForceField(),
param_is_initialized_(false),
file_name_(DEFAULT_TAFF_PARAM_FILE)
{
	// set force name
	set_ff_name("TAFF ["+file_name_+"]");
	if(!mol.is_initialized())
		return;
	//creat the component list
	insert_component(new TAFFStretch(*this));
	insert_component(new TAFFBend(*this));
	insert_component(new TAFFOOP(*this));
	insert_component(new TAFFTorsion(*this));
	insert_component(new TAFFVDW(*this));
	insert_component(new TAFFEle(*this));
	//
	bool result = setup(mol);
	if(!result){
		cout<<"force field setup failed!"<<endl;
		valid_ = false;
	}
}
// special constructor with a MOL and only specified energy terms are included
TAFF::TAFF(MOL& mol, vector<TAFFCOMPONENT> energy_terms_list):ForceField(),
param_is_initialized_(false),
file_name_(DEFAULT_TAFF_PARAM_FILE)
{
	// set force name
	set_ff_name("TAFF ["+file_name_+"]");
	// loop the energy terms list and creat corresponding terms
	for(vector<TAFFCOMPONENT>::iterator it = energy_terms_list.begin();
		it != energy_terms_list.end(); ++it)
	switch(*it)
	{
		case TAFF_BOND_STRETCH:
			insert_component(new TAFFStretch(*this));
			break;
		case TAFF_ANGLE_BEND:
			insert_component(new TAFFBend(*this));
			break;
		case TAFF_DIHEDRAL_TORSION:
			insert_component(new TAFFTorsion(*this));
			break;
		case TAFF_OOP_BEND:
			insert_component(new TAFFOOP(*this));
			break;
		case TAFF_VDW:
			insert_component(new TAFFVDW(*this));
			break;
		case TAFF_ELE:
			insert_component(new TAFFEle(*this));
			break;
	}
	//
	bool result = setup(mol);
	if(!result){
		cout<<"force field setup failed!"<<endl;
		valid_ = false;
	}
}
// copy constructor
TAFF::TAFF(const TAFF& taff):
ForceField(taff),
param_is_initialized_(taff.param_is_initialized_),
file_name_(taff.file_name_){}
// destructor
TAFF::~TAFF(){
}

void TAFF::clear(){
	ForceField::clear();
	file_name_ = DEFAULT_TAFF_PARAM_FILE;
	param_is_initialized_ = false;
}

// assignment operator
const TAFF& TAFF::operator=(const TAFF& taff){
	if (this != &taff){
		ForceField::operator=(taff);
		file_name_ = taff.file_name_;
		param_is_initialized_ = taff.param_is_initialized_;
	}
	return *this;
}

bool TAFF::specific_setup(){
	// check whether the molecule is aasigned
	if (get_mol() == 0)
		return false;
	//debug
#ifdef DEBUG
	cout<<"Atom Typering..."<<endl;
#endif
	// bond all ff params with current FFParameter object 
	// assign the atom numeric type according to symbolic type implemented in mol2 file from the FFParameter object bonded with the FF
	parameter_.read_parameter(DEFAULT_TAFF_PARAM_FILE);
	if ((!parameter_.is_initialized()) && (!parameter_.is_valid())){
		parameter_.clear();
		param_is_initialized_ = false;
		return false;
	}
	for (ATOMVec::iterator aiter = get_atoms().begin(); aiter!=get_atoms().end(); ++aiter){
		bool has_found = false;
		for (vector<vector<string> >::iterator it = get_parameters().params_in_each_section["atom"].begin();
			it != get_parameters().params_in_each_section["atom"].end(); ++it)

			if ((*aiter)->get_symbol_type() == (*it)[1]){
				has_found = true;
				(*aiter)->set_type(str2int((*it)[0]));
				break;
			}

			if(has_found == false)
				cout<<"TAFF::specific_setup(): No atom type assigned for :"<<(*aiter)->get_symbol_type()<<endl;
	}
	param_is_initialized_ = true;
	return true;
}

int TAFF::get_update_frequency() const{
	return DEFAULT_UPDATEFREQUENCY;
}

double TAFF::get_stretch_energy() const{
	FFComponent* component = get_component("TAFF Stretch");
	if (component != 0){
		//debug
		//cout<<"stretch flag"<<endl;
		return component->get_energy();
	}
	else
		return 0.0;
}

double TAFF::get_bend_energy() const{
	FFComponent* component = get_component("TAFF Bend");
	if (component != 0){
		return component->get_energy();
	}
	else
		return 0.0;
}

double TAFF::get_torsion_energy() const{
	FFComponent* component = get_component("TAFF Torsion");
	if (component != 0){
		//debug
		//cout<<"torsion flag"<<endl;
		return component->get_energy();
	}
	else
		return 0.0;
}

double TAFF::get_oop_energy() const{
	FFComponent* component = get_component("TAFF OOP");
	if (component != 0){
		//debug
		//cout<<"oop flag"<<endl;
		return component->get_energy();
	}
	else
		return 0.0;
}

double TAFF::get_vdw_energy() const{
	FFComponent* component = get_component("TAFF VDW");
	if (component != 0){
		//debug
		//cout<<"vdw flag"<<endl;
		return component->get_energy();
	}
	return 0.0;
}

double TAFF::get_ele_energy() const{
	FFComponent* component = get_component("TAFF ELE");
	if (component != 0){
		//debug
		//cout<<"ele flag"<<endl;
		return component->get_energy();
	}
	return 0.0;
}

bool TAFF::has_initialized_param() const{
	return param_is_initialized_;
}


string TAFF::get_results() const{
	ostringstream os;

	os<<"\n"
		<<"TAFF Energy:\n"
		<<" - electrostatic     : " <<get_ele_energy()<<  " kcal/mol\n" 
		<<" - van der Waals     : " <<get_vdw_energy()<<  " kcal/mol\n"
		<< " - bond stretch      : " <<get_stretch_energy()<<  " kcal/mol\n"
		<< " - angle bend        : " <<get_bend_energy()<<  " kcal/mol\n" 
		<< " - torsion           : " <<get_torsion_energy()<<  " kcal/mol\n"
		<< " - oop bend          : " <<get_oop_energy()<< " kcal/mol\n" 
		<< "---------------------------------------\n" 
		<< "  total energy       : " <<get_energy()<< " kcal/mol\n";

	return os.str();
}

TAFF taff;                 

