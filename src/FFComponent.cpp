#include "../include/FFComponent.h"

// default constructor 
FFComponent::FFComponent()
		: force_field_(NULL),
		  energy_(0),
		  num_of_atom_types_(0),
		  name_("GenericForceFieldComponent"),
		  enabled_(true){}

// constructor 
FFComponent::FFComponent(ForceField& force_field)
		: force_field_(&force_field),
		  energy_(0),
		  num_of_atom_types_(0),
		  name_("GenericForceFieldComponent"),
		  enabled_(true){}

// copy constructor 
FFComponent::FFComponent(const FFComponent& force_field_component)
		: force_field_(force_field_component.force_field_),
		  energy_(force_field_component.energy_),
		  num_of_atom_types_(force_field_component.num_of_atom_types_),
		  name_(force_field_component.name_),
		  enabled_(force_field_component.enabled_){}

// destructor
FFComponent::~FFComponent(){}

// setup
bool FFComponent::setup(){
     return true;
     }
// update pair lists - empty
//void FFComponent::update(){}

// set name of the component
void FFComponent::set_name(const string& name){
     name_ = name;
     }
     
// get name of the component
string FFComponent::get_name() const{
       return name_;
       }

//return a pointer to the force field
ForceField* FFComponent::get_force_field() const{
            return force_field_;
}

//set the forcefield
void FFComponent::set_force_field(ForceField& ff){
     force_field_ = &ff;
     }

double FFComponent::get_energy() const{
       return energy_;
       }

double FFComponent::update_energy(){
       return 0.0;
       }
void FFComponent::update_forces(){}

void FFComponent::update(){}
            
