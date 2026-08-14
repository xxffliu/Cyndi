#include "../include/ForceField.h"
#include "../include/FFComponent.h"
//#include "utility.h"
using namespace std;

//default constructor
ForceField::ForceField():
mol_(0),
atom_vec_(),
bond_vec_(),
parameter_(),
valid_(false),
ff_name_("Unassigned"),
energy_(0.0),
unassignedAtoms_(),
component_(){}
// copy constructor
ForceField::ForceField(const ForceField& ff):
mol_(ff.mol_),
atom_vec_(ff.atom_vec_),
bond_vec_(ff.bond_vec_),
parameter_(ff.parameter_),
valid_(ff.valid_),
ff_name_(ff.ff_name_),
energy_(ff.energy_)
{
	for (vector<FFComponent*>::size_type i = 0; i < ff.component_.size(); i++)
	{
		//FFComponent* temp(*ffciter);
		component_.push_back((FFComponent*)ff.component_[i]);
		//delete temp;
	}
}
// assignment operator
ForceField& ForceField::operator=(const ForceField& ff)
{
	if (&ff != this)
	{
		mol_ = ff.mol_;
		atom_vec_.clear();
		bond_vec_.clear();
		atom_vec_ = ff.atom_vec_;
		bond_vec_ = ff.bond_vec_;
		parameter_ = ff.parameter_;
		valid_ = ff.valid_;
		ff_name_ = ff.ff_name_;
		energy_ = ff.energy_;
		for(vector<FFComponent*>::size_type i = 0; i< component_.size(); i++)
			delete component_[i];
		component_.clear();
		for (vector<FFComponent*>::size_type i = 0; i < ff.component_.size(); i++)
		{
			//FFComponent* temp(*ffciter);
			component_.push_back((FFComponent*)ff.component_[i]);
			//delete temp;
		}
		return (*this);
	}
	return *this;
}

// clear method called by destructor

void ForceField::clear()
{
	mol_ = 0;
	/*for(ATOMVec::iterator it = atom_vec_.begin(); it != atom_vec_.end(); ++it)
	delete *it;
	for(BONDVec::iterator it = bond_vec_.begin(); it != bond_vec_.end(); ++it)
	delete *it;*/
	atom_vec_.clear();
	bond_vec_.clear();
	parameter_.clear();
	valid_ = false;
	ff_name_ = "Unassigned";
	energy_ = 0.0;
	for(vector<FFComponent*>::size_type i = 0; i< component_.size(); i++)
		delete component_[i];
	component_.clear();
	/*for(ATOMVec::iterator it = unassignedAtoms_.begin(); it != unassignedAtoms_.end(); ++it)
		delete *it;*/
	unassignedAtoms_.clear();
}

// constructor initiated with a mol;
ForceField::ForceField(MOL& mol):
parameter_(),
valid_(false),
ff_name_("Unassigned"),
energy_(0.0),
unassignedAtoms_()
{
	bool result = setup(mol);
	if (!result)
	{
		cout<<"Force field set up failed!"<<endl;
		valid_ = false;
	}
}

// destructor
ForceField::~ForceField()
{
	clear();
	valid_ = false;
}

bool ForceField::isValid() const
{
	return valid_;
}

// setup methods
bool ForceField::setup(MOL& mol)
{
	unassignedAtoms_.clear();
	mol_ = &mol;
	mol_->init_neighbor_list();
	atom_vec_ = mol_->get_atom_vector();
	bond_vec_ = mol_->get_bond_vector();
	//cout<<&atom_vec_[0]<<" "<<&mol_->get_atom_vector()[0]<<endl; //debug
	// construct neighbor list;
	for(ATOMVec::iterator aiter = atom_vec_.begin(); aiter!=atom_vec_.end(); ++aiter)
	{
		(*aiter)->clear_neighbor_bond_list();
		(*aiter)->clear_neighbor_atom_list();
		int i = 0;
		//cout<<(*aiter)->get_id()<<" "<<(*aiter)<<endl;
		for(BONDVec::iterator biter = bond_vec_.begin(); biter!=bond_vec_.end(); ++biter)
		{
			//debug
			//cout<<(*biter)->get_first_atom()->get_id()<<" "<<(*aiter)->get_id()<<" "<<(*biter)->get_second_atom()->get_id()<<endl;
			if(((*biter)->get_first_atom()->get_id() == (*aiter)->get_id()) || ((*biter)->get_second_atom()->get_id() == (*aiter)->get_id()))
			{
				(*aiter)->add_neighbor_bond_list(*biter);
				(*aiter)->add_neighbor_atom_list((*biter)->get_partner(*aiter));
				i++;
			}
			(*aiter)->set_num_neighbor_bond(i);
			(*aiter)->set_num_neighbor_atom(i);
		}
	}

	bool success = false;
	// specific FF setup

	success = specific_setup();
	//unassignedAtoms_.clear();

	if (!success)
	{
		cout<<" Force Field specific setup faild "<<endl;
		return false;
	}
	//debug
#ifdef DEBUG
	cout<<"Assigning force field parameters..."<<endl;
#endif
	// Call the setup method for each force field component.
	vector<FFComponent*>::iterator  it;
	for (it = component_.begin(); (it != component_.end()) && success; ++it)
	{
		success = false;
		success = (*it)->setup();
		if (!success)
			cout << "Force Field Component setup of " << (*it)->name_ <<  " failed!" << endl;
	}
	valid_ = success;
	return success;
}

bool ForceField::specific_setup()
{
	return true;
}


FFParameter& ForceField::get_parameters()
{
	return parameter_;
}
string ForceField::get_ff_name()
{
	return ff_name_;
}
void ForceField::set_ff_name(const string& name)
{
	ff_name_ = name;
}
int ForceField::get_num_atoms()
{
	if (!isValid())
		return 0;
	return int(atom_vec_.size());
}

void ForceField::update_forces()
{
	// check for validity of the forcefield
	if (!isValid())
		return;
	// set forces to zero
	for (ATOMVec::iterator aiter = atom_vec_.begin(); aiter != atom_vec_.end(); ++aiter){
		(*aiter)->set_force(0.0,0.0,0.0);
	}
	// call each FF component to add their individual forces to the atoms
	for (vector<FFComponent*>::iterator ffciter = component_.begin();
		ffciter != component_.end(); ++ ffciter){
			if (!((*ffciter)->isenabled()))
				continue;
			(*ffciter)->update_forces();
	}
}

// calculate the RMS of the gradient
double ForceField::get_rms_gradient() const
{
	double sum = 0.0;
	for (ATOMVec::const_iterator const_aiter = atom_vec_.begin(); const_aiter != atom_vec_.end(); ++const_aiter)
		sum += (*const_aiter)->get_force().length_2();
	sum = sqrt(sum/(3 * (double)atom_vec_.size()));
	sum *= AVOGADRO / 1e13;
	return sum;
}

double ForceField::get_energy() const
{
	return energy_;
}
double ForceField::update_energy()
{
	// check for validity of the FF
	if (!isValid())
		return 0;
	//clear the total energy;
	energy_ = 0.0;
	//call each FF component and update their energies
	for (vector<FFComponent*>::iterator ffciter = component_.begin(); ffciter != component_.end(); ++ ffciter)
	{
		if(!(*ffciter)->isenabled())
			continue;
		else
		{
			//cout<<(*ffciter)->get_name()<<" "<<(*ffciter)->update_energy()<<endl;
			energy_ += (*ffciter)->update_energy();
		}
	}
	mol_->set_energy(energy_);
	return energy_;
}

int ForceField::get_update_frequency() const
{
	return 1;
}

void ForceField::update()
{
	// check for validity of the FF
	if (!isValid())
		return ;
	// iterates all of the components and call their update methods
	for (vector<FFComponent*>::iterator ffciter = component_.begin();
		ffciter != component_.end(); ++ ffciter){
			if (!((*ffciter)->isenabled()))
				continue;
			(*ffciter)->update();
	}
}
//insert a new component
void ForceField::insert_component(FFComponent* ffp)
{
	component_.push_back(ffp);
	ffp->set_force_field(*this);
}
// remove a component
void ForceField::remove_component(const FFComponent* ffp)
{
	// iterate all component, test if equal, then remove
	for (vector<FFComponent*>::iterator ffciter = component_.begin();
		ffciter != component_.end(); ++ ffciter){
			if (*ffciter == ffp){
				delete *ffciter;
				component_.erase(ffciter);
				break;
			}
	}
}
void ForceField::remove_component(const string& name)
{
	for (vector<FFComponent*>::iterator ffciter = component_.begin();
		ffciter != component_.end(); ++ ffciter)
	{
			if ((*ffciter)->get_name() == name)
			{
				delete *ffciter;
				component_.erase(ffciter);
				break;
			}
	}
}

// return FF component by its index
FFComponent* ForceField::get_component(const int id) const
{
	if (id > component_.size())
		return 0;
	return component_[id];
}

// return FF component by its name
FFComponent* ForceField::get_component(const string& component_name) const
{
	for (vector<FFComponent*>::size_type i = 0; i < component_.size(); i++)
	{
		if (component_[i]->get_name() == component_name)
			return (component_[i]);
	}
	return 0;
}             

void ForceField::add_unassigned_atom(ATOM* unassigned_atom_ptr)
{
	if(find(unassignedAtoms_.begin(), unassignedAtoms_.end(), unassigned_atom_ptr) != unassignedAtoms_.end())
		unassignedAtoms_.push_back(unassigned_atom_ptr);
	return;
}
ATOMVec& ForceField::get_unassigned_atoms()
{
	return unassignedAtoms_;
}

bool  ForceField::In_the_sameRing(ATOM* a,ATOM* b)
{
	bool InSame=false;
	vector<int>::iterator a_ring_id;
	vector<int>::iterator b_ring_id;
	for(a_ring_id=a->ring_id.begin();a_ring_id!=a->ring_id.end();a_ring_id++)
	{ 
		b_ring_id=find(b->ring_id.begin(),b->ring_id.end(),(*a_ring_id));
		if(b_ring_id!=b->ring_id.end())
		{
			InSame=true;
			break;
		}
	}
	return InSame;
}  




