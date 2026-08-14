#include "../include/Mol.h"
#include "../include/utility.h"
#include <algorithm>
using namespace std;

// fixed by xfliu, 20060226
bool operator<(vector<int> lhs, vector<int> rhs)
{
	return lhs.size() < rhs.size();
}
// member fuctions of class RING
void RING::clear(){
	vatom.clear();
	//vbond.clear();
	size = 0;
	is_aromatic = false;
}
RING::RING():vatom(),size(0),is_aromatic(false){}
RING::RING(const RING& to_copy):vatom(to_copy.vatom), size(to_copy.size), is_aromatic(to_copy.is_aromatic){}
RING& RING::operator = (const RING& rhs)
{
	vatom = rhs.vatom;
	size = rhs.size;
	is_aromatic = rhs.is_aromatic;
	return *this;
}
RING::~RING()
{
	clear();
}
vector3 RING::get_center(){
	vector3 center(VZero);
	if(size != 0){
		for(ATOMVec::iterator it = vatom.begin(); it != vatom.end(); ++it)
			center += (*it)->get_position();
		center /= size;
	}
	return center;
}
bool RING::IsMember(ATOM *a){
	bool _ismember=false;
	vector<ATOM*>::iterator atom;
	atom=find(vatom.begin(),vatom.end(),a);
	if(atom!=vatom.end())
		_ismember=true;
	return _ismember;
}
int RING::NumPIElectrons()
{
	int numPI(0);
	for(ATOMVec::iterator it = vatom.begin(); it != vatom.end(); ++it)
		numPI += (*it)->get_num_PI_electrons();
	return numPI;
}
bool operator<(const MOL& lhs, const MOL& rhs){
	return (lhs.fitness_ < rhs.fitness_);
}
bool GreaterMol(const MOL& lhs, const MOL& rhs)
{
	return (lhs.fitness_ >= rhs.fitness_);
}
MOL::MOL():comment_(),pharmacophore_name_("UNK"),name("MOLECULE"),num_atom(0),num_bond(0),num_rotatable_bonds_(0),
			num_fragment(0),is_fragmentized_(false),fitness_(0.0),energy_(0.),rmsd_(0.),is_initialized_(false),
			_vatom(),_vbond(),barycenter_(),_vfragment(), MW_(0.0)
			{
				_vring.clear();
			}

MOL::MOL(const MOL& mol0):name(mol0.name),num_atom(mol0.num_atom),num_bond(mol0.num_bond),
						  num_fragment(mol0.num_fragment),num_rotatable_bonds_(mol0.num_rotatable_bonds_),
						  pharmacophore_name_(mol0.pharmacophore_name_),comment_(mol0.comment_),
						  _vfragment(mol0._vfragment),// Warning, fragment are also stored in pointers
						  barycenter_(mol0.barycenter_),
						  energy_(mol0.energy_),fitness_(mol0.fitness_), MW_(mol0.MW_)
	{
	//rotation_map(mol0.rotation_map),
	//rotor_id_map(mol0.rotor_id_map){
	// notice that in our MOL data structure, atoms and bonds are stored 
	// in terms of pointer. so we have to create new pointers and assign
	// them with the values from the to-copy mol.
	for(int i = 0; i< mol0.num_atom; i++){
		ATOM* atom = new ATOM;
		*atom = *(mol0._vatom[i]);
		_vatom.push_back(atom);

	}
	for(int i = 0; i < mol0.num_bond; i++){
		BOND* bond = new BOND;
		int atom1 = mol0._vbond[i]->get_first_atom()->get_id() -1;
		int atom2 = mol0._vbond[i]->get_second_atom()->get_id() -1;
		bond->set_first_atom(_vatom[atom1]);
		bond->set_second_atom(_vatom[atom2]);
		*bond = *(mol0._vbond[i]);
		_vbond.push_back(bond);
	}
	//initialize();
}

MOL& MOL::operator=(const MOL& mol0){
	clear();
	name = mol0.name;
	num_atom = mol0.num_atom;
	num_bond = mol0.num_bond;
	num_fragment = mol0.num_fragment;
	num_rotatable_bonds_ = mol0.num_rotatable_bonds_;
	pharmacophore_name_ = mol0.pharmacophore_name_;
	comment_ = mol0.comment_;
	energy_ = mol0.energy_;
	fitness_ = mol0.fitness_;
	//is_fragmentized_ = mol0.is_fragmentized_;
	_vfragment = mol0._vfragment; // Warning, fragment are also stored in pointers
	//_vring = mol0._vring;
	barycenter_ = mol0.barycenter_;
	MW_ = mol0.MW_;
	//rotation_map = mol0.rotation_map;
	//rotor_id_map = mol0.rotor_id_map;

	for(int i = 0; i< mol0.num_atom; i++){
		ATOM* atom = new ATOM;
		*atom = *(mol0._vatom[i]);
		_vatom.push_back(atom);
	}
	for(int i = 0; i < mol0.num_bond; i++){
		BOND* bond = new BOND;
		int atom1 = mol0._vbond[i]->get_first_atom()->get_id() -1;
		int atom2 = mol0._vbond[i]->get_second_atom()->get_id() -1;
		bond->set_first_atom(_vatom[atom1]);
		bond->set_second_atom(_vatom[atom2]);
		*bond = *(mol0._vbond[i]);
		_vbond.push_back(bond);
	}
	//if(!mol0.is_initialized())
	//initialize();
	return *this;
}

MOL::~MOL(){
	clear();
}

void MOL::clear(){
	comment_ = "";
	pharmacophore_name_ = "UNK";
	name = "MOLECULE";
	num_atom = 0;
	num_bond = 0;
	num_fragment = 0;
	num_rotatable_bonds_ = 0;
	fitness_ = 0.;
	energy_ = 0.;
	rmsd_ = 0.;
	MW_ = 0.0;
	is_fragmentized_ = false;
	is_initialized_ = false;
	//bk_pos_.clear();
	//_vatom.clear();
	//_vbond.clear();
	//_vfragment.clear();

	rings_.clear();
	vertex_.clear();
	edge_.clear();
	barycenter_ = VZero;
	clear_rotation_mark();
	// noe delete every pointer in _vatom, _vbond, _vfragment
	for(ATOMVec::iterator it = _vatom.begin(); it != _vatom.end(); ++it)
		delete *it;
	_vatom.clear();
	for(BONDVec::iterator it = _vbond.begin(); it != _vbond.end(); ++it)
		delete *it;
	_vbond.clear();
	for(FRAGVec::iterator it = _vfragment.begin(); it != _vfragment.end(); ++it)
		delete *it;
	_vfragment.clear();
	_vring.clear();
	return;
}


string MOL::get_comment()
{
	return comment_;
}

string MOL::get_pharmacophore_name()
{
	return pharmacophore_name_;
}
string MOL::get_name(){
	return name;
}
string MOL::get_charge_type(){
	return charge_type_;
}
int MOL::get_num_atom(){
	return num_atom;
}
int MOL::get_num_bond(){
	return num_bond;
}
int MOL::get_num_fragment(){
	return num_fragment;
}
int MOL::get_num_ring(){
	return _vring.size();
}
int MOL::get_num_of_rot_bonds(){
	return num_rotatable_bonds_;
}
ATOMVec& MOL::get_atom_vector(){
	return _vatom;
}
vector<vector3> MOL::get_coordinates(){
	vector<vector3> tmp;
	for(ATOMVec::iterator it = _vatom.begin(); it != _vatom.end(); ++it)
		tmp.push_back((*it)->get_position());
	return tmp;
}
vector<vector3> MOL::get_heavy_coordinates(){
	vector<vector3> tmp;
	for(ATOMVec::iterator it = _vatom.begin(); it != _vatom.end(); ++it)
	{
		if(!(*it)->is_hydrogen())
			tmp.push_back((*it)->get_position());
	}
	return tmp;
}
ATOM* MOL::get_atom(int id){
	if (id>0 && id<=num_atom) 
		return _vatom[id-1];
	else
		return (ATOM*)0;
}
ATOM* MOL::get_first_atom(){
	return (_vatom.empty()? NULL : _vatom[0]);
}
BONDVec& MOL::get_bond_vector(){
	return _vbond;
}
BOND* MOL::get_bond(int id){
	if (id>0 && id<=num_bond)
		return _vbond[id-1];
	else
		return NULL;
}
BOND* MOL::get_first_bond(){
	return (_vbond.empty() ? NULL : _vbond[0]);
}

FRAGMENT* MOL::get_fragment(int id){
	if (id>0 && id<=num_fragment)
		return (_vfragment[id-1]);
	else
		return NULL;
}

vector<RING>& MOL::get_ring_vector(){
	return _vring;
} 

FRAGVec& MOL::get_fragment_vector(){
	return _vfragment;
}

double MOL::updateRMSD(){
	double rmsd = 0.0;
	for (ATOMVec::iterator aiter = _vatom.begin(); aiter != _vatom.end(); ++aiter){
		if((*aiter)->is_hydrogen())
			continue;
		vector3 v = (*aiter)->get_position();
		vector3 u = (*aiter)->get_orig_position();
		rmsd += (v-u).length_2();
	}
	rmsd_ = rmsd/num_atom;
	return rmsd_;
}

vector3 MOL::get_barycenter(){
	return barycenter_;
}

float MOL::get_fitness(){
	return fitness_;
}

float MOL::get_energy(){
	return energy_;
}

double MOL::get_rmsd(){
	return rmsd_;
}

double MOL::get_RMSD(MOL& mol){
	//first check the 2 molecules have the same number of atoms
	if(num_atom != mol.get_num_atom()){
		cout<<"MOL::get_RMSD(): Error: the two molecules are diffenrent in size"<<endl;
		exit(1);
	}
	else{  
		double rmsd = 0.0;
		for (int i = 0; i < num_atom; ++i){
			vector3 v = _vatom[i]->get_position();
			vector3 u = mol.get_atom_vector()[i]->get_position();
			rmsd += (v - u).length_2();
		}
		return sqrt(rmsd/num_atom);
	}
}

void MOL::init_neighbor_list(){
	//cout<<"Initiating neighbor list..."<<endl;

	for (ATOMVec::iterator aiter = _vatom.begin(); aiter != _vatom.end(); ++aiter){
		(*aiter)->clear_neighbor_atom_list();
		(*aiter)->clear_neighbor_bond_list();
		int i = 0, k = 0;
		for (BONDVec::size_type j = 0; j < num_bond; j++)
			if ((_vbond[j]->get_first_atom() == *aiter) || (_vbond[j]->get_second_atom() == *aiter)){
				(*aiter)->add_neighbor_bond_list(_vbond[j]);
				(*aiter)->add_neighbor_atom_list(_vbond[j]->get_partner(*aiter));
				i++;
				if(!(_vbond[j]->get_partner(*aiter)->is_hydrogen()))
					k++;
			}
			(*aiter)->set_num_neighbor_bond(i);
			(*aiter)->set_num_neighbor_atom(i);
			if(k == 0)
				(*aiter)->set_num_neighbor_heavy_atom(i);
			else
				(*aiter)->set_num_neighbor_heavy_atom(k);
	}
}

bool MOL::is_initialized(){
	return is_initialized_;
}
//

void MOL::set_comment(string comment)
{
	comment_ = comment;
	return;
}
void MOL::set_pharmacophore_name(const string& name)
{
	pharmacophore_name_ = name;
	return;
}
void MOL::set_name(const string& name0){
	name = name0;
}
void MOL::set_charge_type(string type){
	charge_type_ = type;
}
void MOL::set_num_atom(int num_atom0){
	num_atom = num_atom0;
}
void MOL::set_num_bond(int num_bond0){
	num_bond = num_bond0;
}
void MOL::set_num_fragment(int num_frag){
	num_fragment = num_frag;
}
void MOL::set_atom_vector(ATOMVec& atom_ptr){
	_vatom.clear();
	_vatom = atom_ptr;
}
void MOL::set_bond_vector(BONDVec& _vbond0){
	_vbond.clear();
	_vbond = _vbond0;
}
void MOL::set_fragment_vector(FRAGVec& _vfragment0){
	_vfragment.clear();
	_vfragment = _vfragment0;
}
void MOL::set_fitness(float fit){
	fitness_ = fit;
}

void MOL::set_energy(float energy){
	energy_ = energy;
}
void MOL::set_rmsd(float rmsd){
	rmsd_ = rmsd;
}

void MOL::set_coordinates(vector<vector3> coord){
	if(coord.size() != num_atom){
		cout<<"MOL::set_coordinates: Error: The size of atoms and target coordinates con't match"<<endl;
		exit(-1);
	}
	for(int i = 0; i < num_atom; i++)
		_vatom[i]->set_position(coord[i]);
}     

//
void MOL::swap(MOL& mol1, MOL& mol2){
	MOL temp_mol = mol1;
	mol1 = mol2;
	mol2 = temp_mol;
}

// identity ring atoms and bonds
void MOL::find_rings(){
	//cout<<"Finding ring sets..."<<endl;
	RING ring;
	InitiateGraph_();

	PruneEndBonds_();

	while(!vertex_.empty()){
		sort(vertex_.begin(), vertex_.end());
		VERTEX curr_vertex(vertex_.front());
		CollapsRing_(curr_vertex);
		vertex_.pop_front();
		// update the vertex neighboring list
		for(deque<VERTEX>::iterator viter = vertex_.begin(); viter != vertex_.end(); ++viter)
		{
			viter->num_neighbor = 0;
			for(deque<EDGE>::iterator eiter = edge_.begin(); eiter != edge_.end(); ++eiter)
			{
				if(eiter->has_vertex(*viter))
					viter->num_neighbor += 1;
			}
		}
	}
	// Here we merge the rings bridged by more than 2 atoms into one ring
	/*vector<vector<int> > tmp_ring, common;
	vector<int> boolian(rings_.size(), 1);
	for (vector<vector<int> >::iterator riter = rings_.begin(); riter != rings_.end(); ++riter)
	{
	for (vector<vector<int> >::iterator riter1 = riter+1; riter1 != rings_.end(); ++riter1)
	{
	for(vector<int>::iterator it = riter1->begin(); it != riter1->end(); ++it)
	{
	if(find(riter->begin(), riter->end(), *it) != riter->end())
	{
	tmp.push_back(*it);
	}
	}
	if(tmp.size() > 2 && find(common.begin(), common.end(), tmp) == common.end())
	common.push(tmp);*/

	//the redundency rmoval methods used here sucks...anyway, it works
	/*for (vector<vector<int> >::iterator riter = rings_.begin(); riter != rings_.end(); ++riter){
	sort(riter->begin(),riter->end());
	vector<int>::iterator iter = unique(riter->begin(),riter->end());
	riter->erase(iter,riter->end());
	}
	vector<vector<int> >::iterator iter = unique(rings_.begin(),rings_.end());
	rings_.erase(iter,rings_.end());*/

	// to eliminate the bridged rings, we discard any rings contaning more than 8 atoms
	rings_.erase(remove_if(rings_.begin(), rings_.end(),redundantRing),rings_.end());
	// now we fuse the fused rings into 1 ring, if two of the atoms of one ring are shared by other 2 or 3 rings, 
	// we fuse all the atoms in these rings into one
	//CombineFusedRings_();
	_vring.resize(rings_.size());
	// fixed by xfliu, 20060226
	// first we want to sort the ring vector assendently according to ring size
	sort(rings_.begin(), rings_.end());
	//debug
#ifdef DEBUG
	cout<<" rings size: "<<rings_.size()<<endl;
	for (vector<vector<int> >::iterator riter = rings_.begin(); riter != rings_.end(); ++riter){
		for(vector<int>::iterator iter = riter->begin(); iter!=riter->end(); ++iter)
			cout<<*iter<<" ";
		cout<<endl;
	}
#endif
	int ring_counter = 0;
	for (vector<vector<int> >::iterator riter = rings_.begin(); riter != rings_.end(); ++riter){
		int aromatic_counter = 0;
		for(vector<int>::iterator iter = riter->begin(); iter!=riter->end(); ++iter)
			for(ATOMVec::iterator aiter = _vatom.begin(); aiter != _vatom.end(); ++aiter)
			{
				if (*iter == (*aiter)->get_id()){
					(*aiter)->is_ring = true;
					(*aiter)->ring_id.push_back(ring_counter);
					_vring[ring_counter].vatom.push_back(*aiter);
					//if((*aiter)->is_aromatic())
					//aromatic_counter += (*aiter)->get_aromacity();
				}
			}
			_vring[ring_counter].size = _vring[ring_counter].vatom.size();
			// here we use simple 4n+2 rule to check if the ring is aromatic, here n >= 1, ring size cannot be larger than 6 member
			// additional criteria: 
			/*if(aromatic_counter > 2 && aromatic_counter%4 == 2 && _vring[ring_counter].size <= 6)
			{
				_vring[ring_counter].is_aromatic = true;
				//debug
#ifdef DEBUG
				cout<<"aromatic ";
				for(ATOMVec::iterator it = _vring[ring_counter].vatom.begin(); it != _vring[ring_counter].vatom.end(); ++it)
					cout<<(*it)->get_id()<<" ";
				cout<<endl;
#endif
			}*/
			ring_counter += 1;
	}
	// dump the bonds to  rings respectively
	for(vector<RING>::iterator it = _vring.begin(); it != _vring.end(); ++it)
	{
		for(ATOMVec::iterator ait = it->vatom.begin(); ait != it->vatom.end(); ++ait)
		{
			for(ATOMVec::iterator ait1 = ait + 1; ait1 != it->vatom.end(); ++ait1)
			{
				if((*ait)->is_bonded_to(*ait1))
					it->vbond.push_back(GetBond(*ait, *ait1));
			}
		}
	}
}
// perceive aromatic rings
void MOL::PerceiveAromaticRings()
{
	for(vector<RING>::iterator rit = _vring.begin(); rit != _vring.end(); ++rit)
	{
		bool isAromatic(false);
		// here we use simple 4n+2 rule to check if the ring is aromatic, here n >= 1, ring size cannot be larger than 6 member
		int numDelocateBonds = 0;
		int numDoubleBonds = 0;
		int PI_Electron = 0;
		for(BONDVec::iterator biter = rit->vbond.begin(); biter != rit->vbond.end(); ++biter)
		{
			// stronger criteria: all the double bonds must be in the same ring alternatively
			if((*biter)->is_aromatic())
				numDelocateBonds += 1;
			else if((*biter)->is_double())
				numDoubleBonds += 1;
		}
		// weaker criteria: the ring atoms may be connected with double bonds in the other rings
		int numDummyDoubleBonds(0);
		for(ATOMVec::iterator ait = rit->vatom.begin(); ait != rit->vatom.end(); ++ait)
		{
			for(BONDVec::iterator it1 = (*ait)->get_bond_list().begin(); it1 != (*ait)->get_bond_list().end(); ++it1)
				if((*it1)->is_double())
				{
					ATOM* a = (*it1)->get_partner(*ait);
					// fixed by xfliu, 20060226
					if(a->isAromaticRingAtom() && !rit->IsMember(a))
						numDummyDoubleBonds += 1;
				}
		}
		// empericailly, this ring only possess half of such double bonds
		numDoubleBonds += numDummyDoubleBonds / 2;
		PI_Electron = rit->NumPIElectrons();
		if(PI_Electron % 4 == 2)
		{
			if(rit->size == 5)
			{
				if(numDelocateBonds + numDoubleBonds >= 2 && numDoubleBonds < 3 && numDelocateBonds <= 5)
					isAromatic = true;
			}
			else if(rit->size == 6)
			{
				if(numDelocateBonds + numDoubleBonds >= 3 && numDoubleBonds < 4 && numDelocateBonds <= 6)
					isAromatic = true;
			}
		}
		if(isAromatic)
		{
			rit->is_aromatic = true;
			for(ATOMVec::iterator atom_in_ring= rit->vatom.begin();atom_in_ring!=rit->vatom.end();atom_in_ring++)
			{
				if(!((*atom_in_ring)->isAromaticRingAtom()))
				{
					(*atom_in_ring)->setAromaticRingAtom();
					//cout<<(*atom_in_ring)->get_id()<<" ";
				}
			}
			//cout<<endl;
			for(BONDVec::iterator bond_in_ring=rit->vbond.begin();bond_in_ring!=rit->vbond.end();bond_in_ring++)
			{
				if(!((*bond_in_ring)->isAromaticRingBond()))
					(*bond_in_ring)->setAromaticRingBond();
			
			}
		}
	}
	return;
}

// initiate path graph
void MOL::InitiateGraph_(){
	for(ATOMVec::iterator aiter = _vatom.begin(); aiter!=_vatom.end(); ++aiter){
		VERTEX temp;
		if(!(*aiter)->is_hydrogen()){
			temp.id = (*aiter)->get_id();
			temp.num_neighbor = (*aiter)->get_num_bonded_heavy_atom();
			//cout<<temp.id<<" "<<temp.num_neighbor<<endl;
			temp.set_is_deleted(false);
			vertex_.push_back(temp);
		}
	}
	for(BONDVec::iterator biter = _vbond.begin(); biter!=_vbond.end(); ++biter){
		EDGE temp;
		if((*biter)->get_first_atom()->is_hydrogen() || (*biter)->get_second_atom()->is_hydrogen())
			continue;
		temp.first = (*biter)->get_first_atom()->get_id();
		temp.last = (*biter)->get_second_atom()->get_id();
		temp.path.clear();
		temp.is_deleted = false;
		//cout<<temp.first<<" "<<temp.last<<endl;
		edge_.push_back(temp);
	}
	//debug
	//cout<<vertex_.size()<<" "<<edge_.size()<<endl;
}

// prune the end bonds
void MOL::PruneEndBonds_(){
	sort(vertex_.begin(), vertex_.end());
	for(deque<VERTEX>::iterator viter = vertex_.begin(); viter != vertex_.end(); ++viter){
		//cout<<viter->id<<" "<<viter->num_neighbor<<endl;
		if(viter->num_neighbor == 1){
			for (deque<EDGE>::iterator eiter = edge_.begin(); eiter != edge_.end(); ++eiter)
			{
				if (eiter->first == viter->id)
				{

					for (deque<VERTEX>::iterator viter1 = vertex_.begin(); viter1 != vertex_.end(); ++viter1)
					{
						if (viter1->id == eiter->last)
						{
							viter1->num_neighbor -= 1;
						}
					}
					eiter->is_deleted = true;
				}
				else if (eiter->last == viter->id)
				{

					for (deque<VERTEX>::iterator viter1 = vertex_.begin(); viter1 != vertex_.end(); ++viter1)
					{
						if (viter1->id == eiter->first)
						{
							viter1->num_neighbor -= 1;
						}
					}
					eiter->is_deleted = true;
				}
			}
			viter->set_is_deleted(true);
			sort(vertex_.begin(), vertex_.end());
			//cout<<viter->id<<" "<<viter->num_neighbor<<endl;
		}
	}
	// remove the prunned vertex and edge
	for(deque<VERTEX>::iterator viter = vertex_.begin(); viter != vertex_.end();)
		if(viter->is_deleted()){
			//cout<<"prunning: "<<viter->id<<" "<<viter->num_neighbor<<endl;
			viter = vertex_.erase(viter);   
		}
		else
			++viter;
	for(deque<EDGE>::iterator eiter = edge_.begin(); eiter != edge_.end();)
		if(eiter->is_deleted){
			//cout<<"prunning: "<<viter->id<<" "<<viter->num_neighbor<<endl;
			eiter = edge_.erase(eiter);   
		}
		else
			++eiter;       

}

void MOL::CollapsRing_(VERTEX& curr_vertex){
	EDGE new_edge;
	deque<EDGE> new_edge_set;
	for (deque<EDGE>::iterator eiter = edge_.begin(); eiter != edge_.end(); ++eiter)
	{
		if (eiter->has_vertex(curr_vertex))
		{
			if(eiter->first == eiter->last)
				continue;
			for (deque<EDGE>::iterator eiter1 = eiter + 1; eiter1 != edge_.end(); ++eiter1)
			{
				if(eiter1->has_vertex(curr_vertex))
				{
					if(eiter1->first == eiter1->last)
						continue;
					//debug
					//cout<<"before collaps "<<eiter->first<<" "<<eiter->last<<" "<<eiter1->first<<" "<<eiter1->last<<endl;
					new_edge = contatenate(*eiter, *eiter1, curr_vertex);
					//cout<<"after collaps "<<new_edge.first<<" "<<new_edge.last<<endl;
					new_edge_set.push_back(new_edge);

				}
			}
		}
	}

	for(deque<EDGE>::iterator it = new_edge_set.begin(); it != new_edge_set.end(); ++it)
		edge_.push_back(*it);
	for (deque<EDGE>::iterator it = edge_.begin(); it != edge_.end();)
	{
		vector<int> tmp = it->path;
		sort(tmp.begin(), tmp.end());
		if(unique(tmp.begin(), tmp.end()) != tmp.end())
		{
			it = edge_.erase(it);
			continue;
		}
		if(it->has_vertex(curr_vertex))
		{
			if (it->first == it->last)
			{   
				if(find(it->path.begin(),it->path.end(),curr_vertex.id) == it->path.end())
					it->path.push_back(it->first);            
				rings_.push_back(it->path);
				//it = edge_.erase(it);
			}
			it = edge_.erase(it);
		}
		else
			++it;
	}
}

void MOL::CombineFusedRings_()
{
	if(rings_.size() < 2)
		return;
	else
	{
		vector<vector<int> > tobeFused, candidate, atoms(num_atom);
		int ring_id = 0;
		for(vector<vector<int> >::iterator rit = rings_.begin(); rit != rings_.end(); ++rit)
		{
			for(vector<int>::iterator ait = rit->begin(); ait != rit->end(); ++ait)
			{
				atoms[*ait-1].push_back(ring_id);
			}
			ring_id += 1;
		}
		for(vector<vector<int> >::iterator it = atoms.begin(); it != atoms.end(); ++it)
			if(it->size() == 2)
				candidate.push_back(*it);
		if(candidate.size() <= 1)
			return;
		else
			for(vector<vector<int> >::iterator it = candidate.begin(); it != candidate.end(); ++it)
				for(vector<vector<int> >::iterator it1 = it + 1; it1 != candidate.end(); ++it1)
				{
					sort(it->begin(), it->end());
					sort(it1->begin(), it1->end());
					if(equal(it->begin(),it->end(), it1->begin()))
						tobeFused.push_back(*it);
				}
				//now fuse the rings and remove every single ring
				vector<int> fused;
				ring_id = 0;
				for(vector<vector<int> >::iterator fuse_it = tobeFused.begin(); fuse_it != tobeFused.end(); ++fuse_it)
					for(vector<vector<int> >::iterator rit = rings_.begin(); rit != rings_.end();)
					{
						if(find(fuse_it->begin(), fuse_it->end(), ring_id) != fuse_it->end())
						{
							fused.insert(fused.end(), rit->begin(), rit->end());
							rit = rings_.erase(rit);
						}
						else
							++rit;
						ring_id += 1;
					}
					// remove the redundent atoms and add the fused ring to rings_;
					sort(fused.begin(), fused.end());
					fused.erase(unique(fused.begin(), fused.end()), fused.end());
					rings_.push_back(fused);
					return;
	}
}

// End of ring searching

void MOL::find_rotatable_bonds()
{
	//cout<<"Finding rotatable bonds..."<<endl;
	num_rotatable_bonds_ = 0;
	vector<int> numRotors(num_atom, 0);
	for (BONDVec::iterator biter = _vbond.begin(); biter != _vbond.end(); ++biter)
	{
		// first excluding ring atoms and bonds
		if ((*biter)->is_in_ring())
		{
			//debug
			//cout<<"ring: "<<(*biter)->get_first_atom()->get_id()<<" "<<(*biter)->get_second_atom()->get_id()<<endl;
			continue;
		}
		// Then excluding bonds bearing mutiple valences
		if ((*biter)->is_double() || (*biter)->is_triple() || (*biter)->is_aromatic() || (*biter)->is_amide() || (*biter)->is_guanidino())
		{
			//debug
			//cout<<"multiple valence: "<<(*biter)->get_first_atom()->get_id()<<" "<<(*biter)->get_second_atom()->get_id()<<endl;
			continue;
		}
		// excluding hydrogen bonds
		if((*biter)->get_first_atom()->is_hydrogen() || (*biter)->get_second_atom()->is_hydrogen())
			continue;
		// excluding terminal planar -NH2, charged hydroxyl
		// fixed by xfliu, 20090410
		if ((*biter)->get_first_atom()->get_num_bonded_heavy_atom() == 1 || (*biter)->get_second_atom()->get_num_bonded_heavy_atom() == 1)
		{
			if((*biter)->get_first_atom()->is_nitrogen() || (*biter)->get_second_atom()->is_nitrogen())
			//debug
			//cout<<"hydrogen bond: "<<(*biter)->get_first_atom()->get_id()<<" "<<(*biter)->get_second_atom()->get_id()<<endl;
				continue;
			else if((*biter)->get_first_atom()->get_num_bonded_atom() == 1 || (*biter)->get_second_atom()->get_num_bonded_atom() == 1)
				continue;
		}
		// then excluding -CH3 and -CX3
		if ((*biter)->get_first_atom()->get_symbol_type() == "C.3")
		{
			vector<int> neighbors;
			for(ATOMVec::iterator it = (*biter)->get_first_atom()->get_atom_list().begin(); it != (*biter)->get_first_atom()->get_atom_list().end(); ++it)
				if((*it)->is_hydrogen() || (*it)->is_halogen())
					neighbors.push_back((*it)->get_type());
			if(neighbors.size() == 3 && (neighbors[0] == neighbors[1] && neighbors[1] == neighbors[2] && neighbors[0] == neighbors[2]))
				continue;
		}
		if ((*biter)->get_second_atom()->get_symbol_type() == "C.3")
		{
			vector<int> neighbors;
			for(ATOMVec::iterator it = (*biter)->get_second_atom()->get_atom_list().begin(); it != (*biter)->get_second_atom()->get_atom_list().end(); ++it)
				if((*it)->is_hydrogen() || (*it)->is_halogen())
					neighbors.push_back((*it)->get_type());
			if(neighbors.size() == 3 && (neighbors[0] == neighbors[1] && neighbors[1] == neighbors[2] && neighbors[0] == neighbors[2]))
				continue;
		}
		// then exclude any single bonds neighbor to triple bonds, C=-CH, C=-CX or C=-N
		if ((*biter)->get_first_atom()->get_symbol_type() == "C.1")
		{
			bool flag = false;
			for(ATOMVec::iterator it = (*biter)->get_first_atom()->get_atom_list().begin(); it != (*biter)->get_first_atom()->get_atom_list().end(); ++it)
			{
				if((*it)->get_symbol_type() == "N.1")
				{
					flag = true;
					break;
				}
				else if((*it)->get_symbol_type() == "C.1")
				{
					for(ATOMVec::iterator it1 = (*it)->get_atom_list().begin(); it1 != (*it)->get_atom_list().end(); ++it1)
					{
						if((*it1)->is_hydrogen() || (*it1)->is_halogen())
						{
							flag = true;
							break;
						}
					}
				}
			}
			if(flag)
				continue;
		}
		if ((*biter)->get_second_atom()->get_symbol_type() == "C.1")
		{
			bool flag = false;
			for(ATOMVec::iterator it = (*biter)->get_second_atom()->get_atom_list().begin(); it != (*biter)->get_second_atom()->get_atom_list().end(); ++it)
			{
				if((*it)->get_symbol_type() == "N.1")
				{
					flag = true;
					break;
				}
				else if((*it)->get_symbol_type() == "C.1")
				{
					for(ATOMVec::iterator it1 = (*it)->get_atom_list().begin(); it1 != (*it)->get_atom_list().end(); ++it1)
					{
						if((*it1)->is_hydrogen() || (*it1)->is_halogen())
						{
							flag = true;
							break;
						}
					}
				}
			}
			if(flag)
				continue;
		}
		// debug
#ifdef DEBUG
		cout<<(*biter)->get_first_atom()->get_id()<<"----"<<(*biter)->get_second_atom()->get_id()<<endl;;
#endif
		(*biter)->is_rotor = true;
		(*biter)->get_first_atom()->is_root_node = true;
		numRotors[(*biter)->get_first_atom()->get_id()-1] += 1;
		(*biter)->get_second_atom()->is_root_node = true;
		numRotors[(*biter)->get_second_atom()->get_id()-1] += 1;
		num_rotatable_bonds_ += 1;     
	}
	for(ATOMVec::iterator ait = _vatom.begin(); ait != _vatom.end(); ++ait)
	{
		if((*ait)->is_hydrogen() || !(*ait)->is_root_node)
			continue;
		(*ait)->set_num_rotor_bond(numRotors[(*ait)->get_id() - 1]);
	}
}


void MOL::set_rotation_mark()
{
	int order = 0;
	for (BONDVec::iterator biter = _vbond.begin(); biter != _vbond.end(); ++biter){
		if(!(*biter)->is_rotor)
			continue;
		int rotor_id = (*biter)->get_id();
		clear_rotation_mark();
		(*biter)->get_second_atom()->set_rotation_mark((*biter)->get_first_atom());
		ATOMVec tmp1, tmp2;
		int count1 = 0, count2 = 0;
		for(ATOMVec::iterator aiter = _vatom.begin(); aiter != _vatom.end(); ++aiter){
			if((*aiter)->get_rotation_mark())
			{
				//cout<<(*aiter)->get_id()<<endl;
				//rotation_map[(*biter)->get_id()].push_back(*aiter);
				tmp1.push_back(*aiter);
				count1 += (*aiter)->get_num_rotor_bond();
			}
			else
			{
				tmp2.push_back(*aiter);
				count2 += (*aiter)->get_num_rotor_bond();
			}
		}
		// fixed by xfliu, 20090831
		// the more flexible fragment are set as the rotatable fragments
		if(count1 <= count2)
		{
			//debug
			/*cout<<rotor_id<<"---->";
			for(ATOMVec::iterator it = tmp1.begin(); it != tmp1.end(); ++it)
			cout<<(*it)->get_id()<<" ";
			cout<<endl;*/
			rotation_map[rotor_id] = tmp1;
		}
		else
		{
			/*cout<<rotor_id<<"---->";
			for(ATOMVec::iterator it = tmp2.begin(); it != tmp2.end(); ++it)
			cout<<(*it)->get_id()<<" ";
			cout<<endl;*/
			rotation_map[rotor_id] = tmp2;
		}
		rotor_id_map[order] = rotor_id;
		order += 1;
	}
}

void MOL::clear_rotation_mark()
{
	for(ATOMVec::iterator aiter = _vatom.begin(); aiter != _vatom.end(); ++aiter)
		(*aiter)->clear_rotation_mark();
}

//backup the position
void MOL::bk_position(){
	bk_pos_.clear();
	for(ATOMVec::iterator aiter = _vatom.begin(); aiter != _vatom.end(); ++aiter){
		vector3 tmp;
		tmp.Set((*aiter)->get_position().x(),(*aiter)->get_position().y(),(*aiter)->get_position().z());
		bk_pos_.push_back(tmp);
	}
}
//reset the position
void MOL::reset(){
	if(bk_pos_.size() == num_atom)
		for(int i = 0; i < num_atom; ++i)
			_vatom[i]->set_position(bk_pos_[i]);
}

// transform the molecule coordinates according to a transform matrix
void MOL::transform(matrix3x3& mat){
	for(ATOMVec::iterator aiter = _vatom.begin(); aiter != _vatom.end(); ++aiter){
		vector3 coord((*aiter)->get_position());
		(*aiter)->set_position(mat*coord);
	}
}
void MOL::transform(double (*mat)[4]){
	for(ATOMVec::iterator it = _vatom.begin(); it != _vatom.end(); ++it){
		double tmpx = (*it)->get_position().x(),tmpy = (*it)->get_position().y(),
			tmpz = (*it)->get_position().z(),tmpw = 1.0;
		double x = tmpx*mat[0][0]+tmpy*mat[1][0]+tmpz*mat[2][0]+tmpw*mat[3][0];
		double y = tmpx*mat[0][1]+tmpy*mat[1][1]+tmpz*mat[2][1]+tmpw*mat[3][1];
		double z = tmpx*mat[0][2]+tmpy*mat[1][2]+tmpz*mat[2][2]+tmpw*mat[3][2];
		vector3 v(x, y, z);
		(*it)->set_position(v);
	}
}

// rotate part of the molecule around the specified rotor by an angle value
void MOL::rotate(int rotor_id, double angle){
	BOND* b = get_bond(rotor_id);
	if(b == NULL){
		cout<<"MOL::rotate: Can't retrieve bond with id "<<rotor_id<<endl;
		exit(1);
	}
	vector3 origin = b->get_second_atom()->get_position();
	vector3 v = b->get_first_atom()->get_position() - origin;
	matrix3x3 mat;
	mat.RotAboutAxisByAngle(v, angle);   
	for(ATOMVec::iterator aiter = rotation_map[rotor_id].begin(); 
		aiter != rotation_map[rotor_id].end(); ++aiter){
			vector3 tmp = (*aiter)->get_position()-origin;
			(*aiter)->set_position(mat*tmp + origin);
	}
}

void MOL::apply_rotor(int order, double angle){
	int id = rotor_id_map[order];
	if (id == 0){
		cout<<"MOL::apply_rotor: Error! No rotor id was found with order "<<order<<endl;
		exit(1);
	}
	rotate(id, angle);
}
// compute the barycenter
void MOL::compute_center(){
	barycenter_ = VZero;
	for(ATOMVec::iterator aiter = _vatom.begin(); aiter != _vatom.end(); ++aiter)
		barycenter_ += (*aiter)->get_position();
	barycenter_ /= num_atom;
}

vector3 MOL::update_center(){
	barycenter_ = VZero;
	for(ATOMVec::iterator aiter = _vatom.begin(); aiter != _vatom.end(); ++aiter)
		barycenter_ += (*aiter)->get_position();
	barycenter_ /= num_atom;
	return barycenter_;
}

// transform the mol to the target origin and apply 3 base-transformation and 3 base-rotation
void MOL::apply_trans_rot(const vector3& center, const vector3& target_center, 
						  const vector3& rot, const vector3& trans){
							  vector3 vx(1,0,0), vy(0,1,0), vz(0,0,1);
							  matrix3x3 mat, mat_x, mat_y, mat_z;
							  mat_x.RotAboutAxisByAngle(vx, rot.x()), mat_y.RotAboutAxisByAngle(vy, rot.y()),mat_z.RotAboutAxisByAngle(vz, rot.z());
							  mat = mat_x * mat_y * mat_z;
							  for(ATOMVec::iterator aiter = _vatom.begin(); aiter != _vatom.end(); ++aiter)
								  (*aiter)->set_position(mat * ((*aiter)->get_position()- center) + target_center + trans);
							  return;
}

// center the molecule
/*void MOL::center_pos(){
compute_center();
for(ATOMVec::iterator aiter = _vatom.begin(); aiter != _vatom.end(); ++aiter){
//first record the original position
(*aiter)->set_orig_position((*aiter)->get_position());
// then transform the mol to center
(*aiter)->set_position((*aiter)->get_position() - barycenter_);
}
}
*/               
// minimize rmsd between current coordinates and original ones
double MOL::minimizeRMSD(){
	// compute the barycenters of current position and original position
	vector3 curr_cen, orig_cen;
	int num_heavy_atom = 0;
	for(int i = 0; i< num_atom; ++i){
		//debug
		//cout<<_vatom[i]->get_position()<<" "<<bk_pos_[i]<<endl;
		if(_vatom[i]->is_hydrogen())
			continue;
		curr_cen += _vatom[i]->get_position();
		orig_cen += bk_pos_[i];
		num_heavy_atom += 1;
	}
	curr_cen /= num_heavy_atom;
	orig_cen /= num_heavy_atom;
	//debug
	//cout<<curr_cen<<" "<<orig_cen<<endl;
	//now set both the coordinates systenms to their barycenters;
	vector<vector3> v, u;
	//vector<float> weight;
	for(int i = 0; i< num_atom; ++i){
		if(_vatom[i]->is_hydrogen())
		{
			//weight.push_back(0.0);
			continue;
		}
		v.push_back(_vatom[i]->get_position()-curr_cen);
		u.push_back(bk_pos_[i]-orig_cen);
		//weight.push_back(1.0);
	}
	matrix3x3 R;
	// multiply uk and vk by the weights k
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			R.Set(i,j,0.0);
			double tmp = 0.0;
			for(int k = 0; k < u.size(); k++)
				tmp += v[k][i]*u[k][j];
			R.Set(i,j,tmp);
		}
	}
	// Compute the residual matrix F (see original paper for details, names should
	// be identical). 
	// F is a symmetric matrix
	double F[4][4] = {{R.Get(0,0)+R.Get(1,1)+R.Get(2,2),R.Get(1,2)-R.Get(2,1),R.Get(2,0)-R.Get(0,2),R.Get(0,1)-R.Get(1,0)},
	{R.Get(1,2)-R.Get(2,1),R.Get(0,0)-R.Get(1,1)-R.Get(2,2), R.Get(0,1)+R.Get(1,0),R.Get(0,2)+R.Get(2,0)},
	{R.Get(2,0)-R.Get(0,2),R.Get(0,1)+R.Get(1,0),-R.Get(0,0)+R.Get(1,1)-R.Get(2,2),R.Get(1,2)+R.Get(2,1)},
	{R.Get(0,1)-R.Get(1,0),R.Get(0,2)+R.Get(2,0),R.Get(1,2)+R.Get(2,1),-R.Get(0,0)-R.Get(1,1)+R.Get(2,2)}};
	double eigenmatrix[4][4], eigenvals[4];

	// now calculate the eigenvalues and corresponding eigen vectors of F
	//vector3 eigenvals(VZero);
	//matrix3x3 eigenmatrix = F.findEigenvectorsIfSymmetric(eigenvals);
	matrix3x3::jacobi(4, F[0], eigenvals, eigenmatrix[0]);
	// get the max eigen value and corresponding eigen vector
	double eval_max = eigenvals[3];
	if(eval_max == 0.0){
		cout<<"Error! eval_max has zero value!"<<endl;
		exit(1);
	}
	//vector3 evec_max = eigenmatrix.GetColumn(2);
	// Compute final RMSD
	double sum_of_squares = 0.0;
	//double sum_weight = 0.0;
	for (int i = 0; i < u.size(); ++i)
	{
		sum_of_squares += (v[i].length_2() + u[i].length_2());
		//sum_weight += weight[i];
	}
	double rmsd = sqrt(fabs((sum_of_squares - 2.0 * eval_max)) / double(u.size()));
	rmsd_ = rmsd;
	return rmsd;
}

// compute the minimized rmsd value between two conformers
// the target conformer is fixed during this process
double MOL::minimizeRMSD(MOL& target){
	//first check two conformers are the same size
	if(num_atom != target.get_num_atom()){
		cout<<"MOL::minimizeRMSD: Error! 2 Coordinates sets are incompatible in size"<<endl;
		exit(1);
	}

	compute_center();
	target.compute_center();
	vector<vector3> v, u;//coordinates container for this mol and target mol
	//vector<float> weight; // weight vector, heavy atom set to 1, hydrogen set to 0
	for(ATOMVec::iterator it1 = _vatom.begin(), it2 = target.get_atom_vector().begin();
		it1 != _vatom.end() && it2 != target.get_atom_vector().end(); ++it1, ++it2){
			// skip if EITHER atom is hydrogen so v and u stay aligned heavy-atom sets
			if((*it1)->is_hydrogen() || (*it2)->is_hydrogen())
				continue;
			v.push_back((*it1)->get_position()-barycenter_);
			u.push_back((*it2)->get_position()-target.get_barycenter());
	}
	//Match vector sets v to vector sets u, two sets should 
	//have same dimension,3D, same size of vector

	matrix3x3 R;
	// multiply uk and vk by the weights k
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			R.Set(i,j,0.0);
			double tmp = 0.0;
			for(int k = 0; k < u.size(); k++)
				tmp += v[k][i]*u[k][j];
			R.Set(i,j,tmp);
		}
	}
	// Compute the residual matrix F (see original paper for details, names should
	// be identical). 
	// F is a symmetric matrix
	double F[4][4] = {{R.Get(0,0)+R.Get(1,1)+R.Get(2,2),R.Get(1,2)-R.Get(2,1),R.Get(2,0)-R.Get(0,2),R.Get(0,1)-R.Get(1,0)},
	{R.Get(1,2)-R.Get(2,1),R.Get(0,0)-R.Get(1,1)-R.Get(2,2), R.Get(0,1)+R.Get(1,0),R.Get(0,2)+R.Get(2,0)},
	{R.Get(2,0)-R.Get(0,2),R.Get(0,1)+R.Get(1,0),-R.Get(0,0)+R.Get(1,1)-R.Get(2,2),R.Get(1,2)+R.Get(2,1)},
	{R.Get(0,1)-R.Get(1,0),R.Get(0,2)+R.Get(2,0),R.Get(1,2)+R.Get(2,1),-R.Get(0,0)-R.Get(1,1)+R.Get(2,2)}};
	double eigenmatrix[4][4], eigenvals[4];

	// now calculate the eigenvalues and corresponding eigen vectors of F
	//vector3 eigenvals(VZero);
	//matrix3x3 eigenmatrix = F.findEigenvectorsIfSymmetric(eigenvals);
	matrix3x3::jacobi(4, F[0], eigenvals, eigenmatrix[0]);
	// get the max eigen value and corresponding eigen vector
	double eval_max = eigenvals[3];
	if(eval_max == 0.0){
		cout<<"Error! eval_max has zero value!"<<endl;
		exit(1);
	}
	//vector3 evec_max = eigenmatrix.GetColumn(2);
	// Compute final RMSD
	double sum_of_squares = 0.0;
	//double sum_weight = 0.0;
	for (int i = 0; i < u.size(); ++i)
	{
		sum_of_squares += (v[i].length_2() + u[i].length_2());
		//sum_weight += weight[i];
	}
	double rmsd = sqrt(fabs((sum_of_squares - 2.0 * eval_max)) / double(u.size()));

	/*double q[4] = {eigenmatrix[3][0], eigenmatrix[3][1], eigenmatrix[3][2], eigenmatrix[3][3]};
	vector3 c1(q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3],2.0*(q[1]*q[2]-q[0]*q[3]),2.0*(q[1]*q[3]+q[0]*q[2]));
	vector3 c2(2.0*(q[2]*q[1]+q[0]*q[3]),q[0]*q[0]-q[1]*q[1]+q[2]*q[2]-q[3]*q[3],2.0*(q[2]*q[3]-q[0]*q[1]));
	vector3 c3(2.0*(q[3]*q[1]-q[0]*q[2]),2.0*(q[3]*q[2]+q[0]*q[1]),q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3]);
	//vector3 c1(1.0-2.0*(q[2]*q[2]+q[3]*q[3]),2.0*(q[1]*q[2]-q[0]*q[3]),2.0*(q[1]*q[3]+q[0]*q[2]));
	//vector3 c2(2.0*(q[2]*q[1]+q[0]*q[3]),1.0-2.0*(q[1]*q[1]+q[3]*q[3]),2.0*(q[2]*q[3]-q[0]*q[1]));
	//vector3 c3(2.0*(q[3]*q[1]-q[0]*q[2]),2.0*(q[3]*q[2]+q[0]*q[1]),1.0-2.0*(q[2]*q[2]+q[1]*q[1]));
	matrix3x3 rot;
	rot.SetRow(1,c1), rot.SetRow(2,c2), rot.SetRow(3,c3);
	transform(rot);*/
	rmsd_ = rmsd;
	return rmsd;

	/*double  v1[num_atom * 3], u1[num_atom * 3];
	vector<vector3> v, u;//coordinates container for this mol and target mol
	vector<float> weight; // weight vector, heavy atom set to 1, hydrogen set to 0
	for(ATOMVec::iterator it1 = _vatom.begin(), it2 = target.get_atom_vector().begin();
	it1 != _vatom.end(), it2 != target.get_atom_vector().end(); ++it1, ++it2){
	v.push_back((*it1)->get_position());
	u.push_back((*it2)->get_position());
	if((*it1)->is_hydrogen())
	weight.push_back(0.0);
	else
	weight.push_back(1.0);
	}
	double sum_of_squares = 0.0;
	double sum_weight = 0.0;
	for(int i = 0; i < num_atom; i++){
	v1[i*3] = v[i].x(), v1[i*3+1] = v[i].y(), v1[i*3+2] = v[i].z();
	u1[i*3] = u[i].x(), u1[i*3+1] = u[i].y(), u1[i*3+2] = u[i].z();
	sum_of_squares += v[i].length_2() + u[i].length_2();
	sum_weight += weight[i];
	}
	double r[3][3]; 
	double eval_max = qtrfit(v1, u1, num_atom, r);
	double rmsd = sqrt(fabs((sum_of_squares - 2.0 * eval_max)) / sum_weight);
	matrix3x3 rot(r);
	transform(rot);
	//double rmsd = get_RMSD(target);
	return rmsd;*/
}

// compute the atomic weighted gyration radius
float MOL::ComputGyrationRadius()
{
	compute_center();
	vector3 v;
	float r2(0.0);
	for(ATOMVec::iterator it = _vatom.begin(); it != _vatom.end(); ++it)
	{
		//if((*it)->is_hydrogen())
		//	continue;
		v = (*it)->get_position() - barycenter_;
		r2 += (*it)->get_AW() * v.length_2();
	}
	return sqrt(r2/MW_);
}

// fragmentize the molecule if no set information is read from mol2 or some other file format (is_fragmentize_ == false)
// typically the ring system, unsaturated alkane system and amide group are extracted

bool MOL::fragmentize(){return false;}

// initialize the mol
void MOL::initialize(){
	init_neighbor_list();
	find_rings();
	// patch the carboxyl and guanidino group if exists
	/*for(ATOMVec::iterator aiter = _vatom.begin(); aiter != _vatom.end(); ++aiter)
	{
		if((*aiter)->is_carboxyl_oxygen())
			if ((*aiter)->get_num_neighbor_bond() == 1)
			{
				if ((*aiter)->get_symbol_type() != "O.co2")
				{
					(*aiter)->set_symbol_type("O.co2");
				}
				//debug
				//cout<<(*aiter)->get_id()<<" "<<(*aiter)->get_bond_list().size()<<endl;
				BONDVec::iterator biter = (*aiter)->get_bond_list().begin();
				(*biter)->set_type("ar");
			}
			if((*aiter)->is_guanidino()){
				(*aiter)->set_symbol_type("C.cat");
				for(ATOMVec::iterator it = (*aiter)->get_atom_list().begin();
					it != (*aiter)->get_atom_list().end(); ++it){
						(*it)->set_symbol_type("N.pl3");
				}
				for (BONDVec::iterator it = (*aiter)->get_bond_list().begin();
					it != (*aiter)->get_bond_list().end(); ++it){
						(*it)->set_type("ar");
				}
			}
	}*/
	compute_center();
	find_rotatable_bonds();
	PerceiveAromaticRings();
	set_rotation_mark();
	is_initialized_ = true;
}
//function of  GasteigerCharge caculation
bool MOL::calculateGasteigerCharge()
{
	_gsv.resize(num_atom + 1);
	double a,b,c;
	ATOMVec::iterator aiter ;
	InitialPartialCharges();
	for(aiter  = _vatom.begin(); aiter != _vatom.end(); ++aiter)
	{
		if(!GasteigerSigmaChi(*aiter,a,b,c))
			return(false);
		_gsv[(*aiter)->get_id()].SetValues(a,b,c,(*aiter)->get_charge());
	}

	double alpha,charge,denom;
	int iter;
	alpha =1.0;
	ATOM *src,*dst;
	int j;
	for(iter = 0;iter < GASTEIGER_ITERS;++iter)
	{
		alpha *=GASTEIGER_DAMP;
		for(j=1;j < _gsv.size();++j)
		{
			charge=_gsv[j].q;
		//cout<<src->get_id()<<" "<<dst->get_id()<<" "<<charge<<" "<<denom<<" "<<endl;
		_gsv[j].chi=(_gsv[j].c*charge+_gsv[j].b)*charge+_gsv[j].a;
		}

	for(BONDVec::iterator biter = _vbond.begin() ;biter!=_vbond.end();++biter)
	{ 
		src = (*biter)->get_first_atom();
		dst = (*biter)->get_second_atom();

		if (_gsv[src->get_id()].chi >=_gsv[dst->get_id()].chi)
		{
			if (dst->is_hydrogen())
				denom = double(GASTEIGER_DENOM);
			else
				denom =_gsv[dst->get_id()].denom;

		}
		else 
		{
			if (src->is_hydrogen())
				denom=double(GASTEIGER_DENOM);
			else
				denom=_gsv[src->get_id()].denom;
		} 
		charge = (_gsv[src->get_id()].chi - _gsv[dst->get_id()].chi)/denom;
		_gsv[src->get_id()].q -=alpha*charge;
		_gsv[dst->get_id()].q +=alpha*charge;
	}
	}

		for (ATOMVec::iterator aiter = _vatom.begin(); aiter!=_vatom.end() ; ++aiter)
			(*aiter)->set_charge(_gsv[(*aiter)->get_id()].q);

		return(true);
}
   //set initial charges in @p MOL
   // Carbonyl O => -0.5
   // Phosphate O => -0.666
   // Sulfate O => -0.5
   // All other atoms are set to have their initial charge from their formal charge
  void MOL:: InitialPartialCharges()
 {
	 for(ATOMVec::iterator aiter = _vatom.begin(); aiter!=_vatom.end(); ++aiter)
    {
		if ((*aiter)->is_carboxyl_oxygen())
			(*aiter)->set_charge(-0.500);
		else if ((*aiter)->is_phosphate_oxygen()&&
			(*aiter)->get_num_bonded_heavy_atom()==1)
			(*aiter)->set_charge(-0.666);
		else if ((*aiter)->is_sulphate_oxygen())
			(*aiter)->set_charge(-0.500);
		else
			(*aiter)->set_charge((double)(*aiter)->get_charge());
       }
 }

 
  //get a,b,c values
  bool MOL::GasteigerSigmaChi(ATOM *atom,double &a,double &b,double &c)
  {

	  double val[3] = {0.0,0.0,0.0};
	  if(atom->get_element()=="H")
	  {
		  val[0] = 0.37;
		  val[1] = 7.17;
		  val[2] = 12.85;}
	  else if(atom->get_element()=="C")
	  {
		  if(atom->get_symbol_type()=="C.3")
	  { 
		  val[0] = 0.68;
		  val[1] = 7.98;
		  val[2] = 19.04;}
	  else if(atom->get_symbol_type()=="C.2" || atom->get_symbol_type() == "C.ar")  
	  {
		  val[0] = 0.98;
		  val[1] = 8.79;
		  val[2] = 19.62;
	  }
	  else if (atom->get_symbol_type()=="C.1") 
	  { 
		  val[0] = 1.67;
		  val[1] = 10.39;
		  val[2] = 20.57;

	  } 

	  }

	  if(atom->get_element()=="N") 
	  {
		  if(atom->get_symbol_type()=="N.3")
		  {
			  if(atom->get_num_neighbor_bond()==4|| atom->get_charge())
			  {
				  val[0] = 0.0;
				  val[1] = 0.0;
				  val[2] = 23.72;
			  }
			  else
			  {
				  val[0] = 2.08;
				  val[1] = 11.54;
				  val[2] = 23.72;
			  } 

		  }
		  if (atom->get_symbol_type()=="N.2" || atom->get_symbol_type() == "N.ar")
		  {

			  val[0] = 2.57;
			  val[1] = 12.87;
			  val[2] = 24.87;
		  } 
		  if(atom->get_symbol_type()=="N.am"||atom->get_symbol_type()=="N.pl3")
		  {
			  val[0] = 2.46;
			  val[1] = 12.32;
			  val[2] = 24.86;
		  }


		  if (atom->get_symbol_type()=="N.1" )
		  {
			  val[0] = 3.71;
			  val[1] = 15.68;
			  val[2] = 27.11;
		  } 
	  }
	  else if (atom->get_element()=="O")
	  {
		  if (atom->get_symbol_type()=="O.3") 
		  {
			  val[0] = 2.65;
			  val[1] = 14.18;
			  val[2] = 28.49;
		  } 
		  else if(atom->get_symbol_type()=="O.2")
		  { 
			  val[0] = 3.75;
			  val[1] = 17.07;
			  val[2] = 31.33;
		  }
	  } 
	  else if(atom->get_element()=="F")
	  {
		  val[0] = 3.12;
		  val[1] = 14.66;
		  val[2] = 30.82;
	  }
	  else if (atom->get_element()=="P")
	  {
		  val[0] = 1.62;
		  val[1] = 8.90;
		  val[2] = 18.10;
	  }
	  else if(atom->get_element()=="S")
	  {     
		  if(atom->get_symbol_type()=="S.2"||atom->get_symbol_type()=="S.3"||atom->get_symbol_type()=="S.O")
			{
		  val[0] = 2.39;
		  val[1] = 10.14;
		  val[2] = 20.65;
		  }
		  else
		{
		  val[0] = 2.39;
		  val[1] = 12.00;
		  val[2] = 24.00;
		}   
	  }
	  else if(atom->get_element()=="Cl") 
	  {
		  val[0] = 2.66;
		  val[1] = 11.00;
		  val[2] = 22.04;
	  }
	  else if(atom->get_element()=="Br")
	  { val[0] = 2.77;
	  val[1] = 10.08;
	  val[2] = 19.71;}
	  else if(atom->get_element()=="I")
	  {
		  val[0] = 2.90;
		  val[1] = 9.90;
		  val[2] = 18.82;
	  }
	  else if(atom->get_element()=="Al")
	  { 
		  val[0] = 1.06;
		  val[1] = 5.47;
		  val[2] = 11.65;   
	  }

	  a=val[1];
	  b = (val[2]-val[0])/2;
	  c = (val[2]+val[0])/2 - val[1];
	  return (true);

  }
bool MOL::IsInAromaticRingSize(ATOM* atom, int size)
{
	bool isIn = false;
	for(vector<RING>::iterator rit = _vring.begin(); rit != _vring.end(); ++rit)
	{
		vector<ATOM*> ::iterator a = find(rit->vatom.begin(), rit->vatom.end(), atom);
		if(a != rit->vatom.end() && rit->size == size && rit->is_aromatic)
			isIn = true;
	}
	return isIn;
}
   bool MOL::IsInRingSize(ATOM *atom,int size) 
   {  
	   bool IsIn=false;
	   vector<RING>::iterator r;
	   vector<ATOM*>::iterator a;
	   for(r=_vring.begin();r!=_vring.end();r++){
		   a=find((*r).vatom.begin(),(*r).vatom.end(),atom);
		   if(a!=(*r).vatom.end()&&(*r).size==size){
			   IsIn=true;
		   return IsIn;}
	  }
	   return IsIn;

		} 

	 BOND *MOL::GetBond(ATOM *bgn,ATOM *end) const
  {
	  if(bgn->is_bonded_to(end)) { 
     vector<BOND*>::iterator i;

	 for (i=bgn->get_bond_list().begin();i!=bgn->get_bond_list().end();i++ ){ 
		 if((*i)->get_partner(bgn)==end)
			 return (*i);} 

    return(NULL); //just to keep the SGI compiler happy
     }
	  return (NULL);
   }
 ATOM *MOL::get_map(int type) {
	int index;
	vector<ATOM*>::iterator atom;
	for(atom= _vatom.begin();atom!=_vatom.end();atom++){
		index=(*atom)->get_mmff94_type();
		map_atom[index]=(*atom);
	}
	return  (map_atom[type]);

}
 	
	