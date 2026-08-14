#include "../include/TAFFTorsion.h"
#include "../include/TAFF.h"
#include "../include/utility.h"
#include "../include/Vector3.h"
#include <cmath>
using namespace std;

// constructors

TAFFTorsion::TAFFTorsion():
FFComponent(),
fast_access(),
num_of_atom_types_(0),
torsion_data_holder_(){
	// set component name
	set_name("TAFF Torsion");
}
TAFFTorsion::TAFFTorsion(ForceField& ff):
FFComponent(ff),
fast_access(),
num_of_atom_types_(0),
torsion_data_holder_(){
	set_name("TAFF Torsion");
}
// copy constructor
TAFFTorsion::TAFFTorsion(const TAFFTorsion& to_copy):
FFComponent(to_copy){
	fast_access = to_copy.fast_access;
	num_of_atom_types_ = to_copy.num_of_atom_types_;
	torsion_data_holder_ = to_copy.torsion_data_holder_;
}
//destructor
TAFFTorsion::~TAFFTorsion(){
	fast_access.clear();
	torsion_data_holder_.clear();
}

// extract torsion parameters from FFParameter object bonded to the force field
// and establish a hash table for fast access
bool TAFFTorsion::extract_TOR_parameters(FFParameter& ffp){
	if (!ffp.is_valid())
		return false;
	// build a two dim array of atom types and loop variable
	//FFParameter::AtomTypes& atom_types = ffp.get_atomtypes();
	int num_types = ffp.get_num_types();
	int num_entries = num_types * num_types * num_types * num_types;
	num_of_atom_types_ = num_types;
	fast_access.clear();
	for (int i = 1; i <=5; i++){
		//fast_access[i].num_of_atom_types = num_types;
		// allocate two onedimensional fields for the two parameters
		// and a two dimensional field of boolean variables.
		// we might have to delete old stuff lying around
		/*if (fast_access[i].V)
		delete [] fast_acess[i].k;
		if (fast_access[i].n)
		delete [] fast_access[i].r0;
		if (fast_access[i].is_defined)
		delete [] fast_access[i].is_defined;
		fast_access[i].V  = new float[num_entries];
		fast_access[i].n = new float[num_entries];
		fast_access[i].is_defined = new bool[num_entries];*/
		//fast_access[i].clear();
		fast_access[i].V.resize(num_entries);
		fast_access[i].n.resize(num_entries);
		fast_access[i].is_defined.resize(num_entries);
		for (int j = 0; j < num_entries; j++) 
			fast_access[i].is_defined[j] = false;
	}
	int type1, type2, type3, type4, index, bond_order;
	string name_type1, name_type2, name_type3, name_type4;
	for (vector<vector<string> >::iterator it = force_field_->get_parameters().params_in_each_section["torsion"].begin();
		it != force_field_->get_parameters().params_in_each_section["torsion"].end(); ++it){
			name_type1 = (*it)[0];
			name_type2 = (*it)[1];
			name_type3 = (*it)[2];
			name_type4 = (*it)[3];
			bond_order = get_bond_order((*it)[4]);
			type1 = ffp.get_type(name_type1);
			type2 = ffp.get_type(name_type2);
			type3 = ffp.get_type(name_type3);
			type4 = ffp.get_type(name_type4);
			//debug
			//cout<<type1<<" "<<type2<<" "<<type3<<" "<<type4<<endl;
			if ((type1 == -1) || (type2 == -1) || (type3 == -1) || (type4 == -1)){
				cout<<"error! no numeric atom type defined for atom "<<type1<<" or "<<type2<<" or "<<type3<<" or "<<type4<<endl;
				return false;
			}
			index = (type1-1) + (type2-1)* num_types + (type3-1)*num_types*num_types + (type4-1)*num_types*num_types*num_types;
			fast_access[bond_order].is_defined[index] = true;
			fast_access[bond_order].V[index] = str2double((*it)[5]);
			fast_access[bond_order].n[index] = str2double((*it)[6]);
			//build the symmetric index
			index = (type4-1) + (type3-1)* num_types + (type2-1)*num_types*num_types + (type1-1)*num_types*num_types*num_types;
			fast_access[bond_order].is_defined[index] = true;
			fast_access[bond_order].V[index] = str2double((*it)[5]);
			fast_access[bond_order].n[index] = str2double((*it)[6]);
	}
	return true;
}
// query a set of torsion parameters has defined for a given combination of atom types and bond order of the bond between atom j and k
bool TAFFTorsion::has_params(int i, int j, int k, int l, int bond_order){
	//debug
	//cout<<bond_order<<endl;
	if ((i > 0 && i <= num_of_atom_types_) && (j>0 && j <= num_of_atom_types_) && (k>0 && k <= num_of_atom_types_) && (l>0 && l <= num_of_atom_types_) 
		&& (bond_order > 0 && bond_order <= 5)){
			//debug
			//cout<<i<<" "<<j<<" "<<k<<" "<<l<<" is defined"<<endl;
			return (fast_access[bond_order].is_defined[(i-1) + (j-1) * num_of_atom_types_ + (k-1) * num_of_atom_types_*num_of_atom_types_ + (l-1)*num_of_atom_types_*num_of_atom_types_*num_of_atom_types_]
			|| fast_access[bond_order].is_defined[(l-1) + (k-1) * num_of_atom_types_ + (j-1) * num_of_atom_types_*num_of_atom_types_ + (i-1)*num_of_atom_types_*num_of_atom_types_*num_of_atom_types_]);
	}
	else
		return false;
}
// return the parameters for a given combination of atom types and bond order of the bond between atom j and k
TAFFTorsion::ForceValues TAFFTorsion::get_params(int i, int j, int k, int l ,int bond_order){
	TAFFTorsion::ForceValues value;
	if (assign_params(value, i, j, k, l ,bond_order))
		return value;
	else{
		cout<<"TAFFBend::get_params():Warning: no angle bend parameters are assigned for atom "<<
			i<<" and "<<j<<" and "<<k<<" and "<<l<<" with bond order "<<bond_order<<endl;
		value.V = 0.0;
		value.n = 0.0;
		return value;
	}
}
// assign the parameters for a given combination of atom types and bond order of the bond between atom j and k;
bool TAFFTorsion::assign_params(TAFFTorsion::ForceValues& param, int i, int j, int k, int l, int bond_order){
	if (has_params(i, j, k, l, bond_order)){
		int t2 = num_of_atom_types_ * num_of_atom_types_;
		int t3 = num_of_atom_types_ * num_of_atom_types_ * num_of_atom_types_;
		param.V = fast_access[bond_order].V[(i-1) + (j-1)*num_of_atom_types_ + (k-1)*t2 + (l-1)*t3];
		param.n = fast_access[bond_order].n[(i-1) + (j-1)*num_of_atom_types_ + (k-1)*t2 + (l-1)*t3];
		return true;
	}
	else
		return false;
}

// set up method
bool TAFFTorsion::setup(){
	if (force_field_ == NULL){
		cout<<"TAFFTorsion::setup(): force field bound component can not be found"<<endl;
		return false;
	}
	// clear the parameter container
	torsion_data_holder_.clear();
	// tempararily set this component enabled;
	setenabled(true);

	TAFF* taff_force_field = dynamic_cast<TAFF*>(force_field_);
	if (taff_force_field == NULL || taff_force_field->has_initialized_param()){
		bool result = extract_TOR_parameters(force_field_->get_parameters());
		if (!result){
			cout << "TAFFTorsion::setup(): can not found angle Torsion section"<<endl;
			return false;
		}
	}
	TAFFTorsion::ForceValues value;
	// retrieve all torsion parameters
	BONDVec::iterator biter1, biter2, biter3;
	ATOM* a1;
	ATOM* a2;
	ATOM* a3;
	ATOM* a4;
	// loop all atoms in the current molecule and proper torsions will be packed into the vector
	for (ATOMVec::iterator aiter = force_field_->get_atoms().begin();
		aiter != force_field_->get_atoms().end(); ++aiter){
			for (biter1 = (*aiter)->get_bond_list().begin(); biter1 != (*aiter)->get_bond_list().end(); ++biter1){
				// central atoms
				a2 = *aiter;
				if ((*biter1)->get_second_atom() == a2)
					a3 = const_cast<ATOM*>((*biter1)->get_first_atom());
				else if((*biter1)->get_first_atom() == a2)
					a3 = const_cast<ATOM*>((*biter1)->get_second_atom());
				for (biter2 = (*aiter)->get_bond_list().begin(); biter2 != (*aiter)->get_bond_list().end(); ++biter2){
					if ((*biter2)->get_second_atom() != (*biter1)->get_second_atom()){
						// determine the first atom
						if ((*biter2)->get_first_atom() == *aiter)
							a1 = const_cast<ATOM*>((*biter2)->get_second_atom());
						else
							a1 = const_cast<ATOM*>((*biter2)->get_first_atom());
						for (biter3 = const_cast<ATOM*>((*biter1)->get_second_atom())->get_bond_list().begin();
							biter3 != const_cast<ATOM*>((*biter1)->get_second_atom())->get_bond_list().end() ; ++biter3){
								if (!(*biter3)->has_atom(a2->get_id())){
									// determine the fourth atom a4
									if ((*biter3)->get_first_atom() == a3)
										a4 = const_cast<ATOM*>((*biter3)->get_second_atom());
									else
										a4 = const_cast<ATOM*>((*biter3)->get_first_atom());
									//search torsion parameters for a1, a2, a3, a4
									int atom_type_1 = a1->get_type();
									int atom_type_2 = a2->get_type();
									int atom_type_3 = a3->get_type();
									int atom_type_4 = a4->get_type();
									//debug
									//cout<<a1->get_id()<<" "<<a2->get_id()<<" "<<a3->get_id()<<" "<<a4->get_id()<<endl;
									int bond_order = get_bond_order((*biter1)->get_type());
									torsion_data_holder_.push_back(TAFFTorsion::TorsionData());
									torsion_data_holder_.back().atom1 = a1;
									torsion_data_holder_.back().atom2 = a2;
									torsion_data_holder_.back().atom3 = a3;
									torsion_data_holder_.back().atom4 = a4;
									//debug
									//cout<<atom_type_1<<" "<<atom_type_2<<" "<<atom_type_3<<" "<<atom_type_4<<endl;
									if(!assign_params(value, atom_type_1, atom_type_2, atom_type_3, atom_type_4, bond_order))
										if(!assign_params(value, atom_type_4, atom_type_3, atom_type_2, atom_type_1, bond_order))
											if(!assign_params(value, ATOM::ANY_TYPE, atom_type_2, atom_type_3, atom_type_4, bond_order))
												if(!assign_params(value, atom_type_4, atom_type_3, atom_type_2, ATOM::ANY_TYPE, bond_order))
													if(!assign_params(value, atom_type_1, atom_type_2, atom_type_3, ATOM::ANY_TYPE, bond_order))
														if(!assign_params(value, ATOM::ANY_TYPE, atom_type_3, atom_type_2, atom_type_1, bond_order))
															if(!assign_params(value, ATOM::ANY_TYPE, atom_type_2, atom_type_3, ATOM::ANY_TYPE, bond_order))
																if(!assign_params(value, ATOM::ANY_TYPE, atom_type_3, atom_type_2, ATOM::ANY_TYPE, bond_order))
																{
#ifdef DEBUG
																	cout<<"TAFFTorsion::setup(): cannot find torsion parameters for "<<
																		force_field_->get_parameters().get_type_name(atom_type_1)<<"-"<<
																		force_field_->get_parameters().get_type_name(atom_type_2)<<"-"<<
																		force_field_->get_parameters().get_type_name(atom_type_3)<<"-"<<
																		force_field_->get_parameters().get_type_name(atom_type_4)<<
																		" with bond order "<<bond_order<<" and default parameters are assigned"<<endl;
#endif
																	value.V = 1.0;
																	value.n = 0.5;
																	get_force_field()->add_unassigned_atom(a1);
																	get_force_field()->add_unassigned_atom(a2);
																	get_force_field()->add_unassigned_atom(a3);
																	get_force_field()->add_unassigned_atom(a4);
																	//continue;
																}
																torsion_data_holder_.back().value = value;
								}
						}
					}
				}
			}
	}
	// everything goes well
	return true;
}

// update methods
double TAFFTorsion::update_energy(){
	//energy initializion
	energy_ = 0.0;
	if (torsion_data_holder_.size() == 0)
		return 0.0;
	// iterate all torsions and summarize the torsion energies
	for (vector<TorsionData>::size_type i = 0; i<torsion_data_holder_.size(); ++i){
		vector3 v1, v2, v3, v4;
		v1 = torsion_data_holder_[i].atom1->get_position();
		v2 = torsion_data_holder_[i].atom2->get_position();
		v3 = torsion_data_holder_[i].atom3->get_position();
		v4 = torsion_data_holder_[i].atom4->get_position();
		double torsion = CalcTorsionAngle(v1, v2, v3, v4);
		//convert angle unit from rad to deg
		torsion *= DEG_TO_RAD;


		if (torsion_data_holder_[i].value.n < 0)   
			energy_ += torsion_data_holder_[i].value.V/2*(1 - cos(fabs(torsion_data_holder_[i].value.n)*torsion));

		else
			energy_ += torsion_data_holder_[i].value.V/2*(1 + cos(torsion_data_holder_[i].value.n*torsion));

	}
	return energy_;
}

// calculate current forces imposed by torsion and add them to the force field;
void TAFFTorsion::update_forces(){
	if (get_force_field() == 0){
		cout<<"TAFFTorsion::update_force(): error! this component doesn't bond to any force field"<<endl;
		return;
	}
	// iterate all torsions and update forces
	for (vector<TorsionData>::size_type i = 0; i<torsion_data_holder_.size(); ++i){
		vector3 force_1 = torsion_data_holder_[i].atom1->get_force(), force_2 = torsion_data_holder_[i].atom2->get_force(), force_3 = torsion_data_holder_[i].atom3->get_force(),
			force_4 = torsion_data_holder_[i].atom4->get_force();
		// calculate the vectors between atom1 and atom2, atom2 and atom3, atom3 and atom4
		vector3 ab = torsion_data_holder_[i].atom1->get_position() - torsion_data_holder_[i].atom2->get_position();
		vector3 ba = torsion_data_holder_[i].atom2->get_position() - torsion_data_holder_[i].atom1->get_position();
		vector3 cb = torsion_data_holder_[i].atom3->get_position() - torsion_data_holder_[i].atom2->get_position();
		vector3 dc = torsion_data_holder_[i].atom4->get_position() - torsion_data_holder_[i].atom3->get_position();
		double length_ab = ab.length();
		if (length_ab == 0.0) continue;
		double length_cb = cb.length();
		if (length_cb == 0.0) continue;
		double length_dc = dc.length();
		if (length_dc == 0.0) continue;
		// calculate the cross product of v1 and v2, v3 and v2
		vector3 t = ba%cb;
		if (t.length() == 0) continue;
		vector3 u = cb%dc;
		if (u.length() == 0) continue;


		/*double torsion = CalcTorsionAngle(torsion_data_holder_[i].atom1->get_position(),torsion_data_holder_[i].atom2->get_position(),
		torsion_data_holder_[i].atom3->get_position(),torsion_data_holder_[i].atom4->get_position());
		//convert angle unit from rad to deg
		torsion *= DEG_TO_RAD;*/
		double cosphi = (t*u)/(t.length()*u.length());
		if(cosphi>1.0)
			cosphi = 1.0;
		if(cosphi<-1.0)
			cosphi = -1.0;
		double torsion = acos(cosphi);
		// convert the units from kcal/mol A to N
		// kcal -> J: 1e3 * 4.2
		// A -> m: 1e-10
		// J/mol -> J: Avogadro
		double dEdphi;
		if (torsion_data_holder_[i].value.n < 0)
			dEdphi = 4.184*1e13/AVOGADRO*torsion_data_holder_[i].value.V/2*fabs(torsion_data_holder_[i].value.n)*sin(fabs(torsion_data_holder_[i].value.n)*torsion);
		else
			dEdphi = -4.184*1e13/AVOGADRO*torsion_data_holder_[i].value.V/2*torsion_data_holder_[i].value.n*sin(torsion_data_holder_[i].value.n*torsion);
		if ((t%u) * cb > 0.0)
			dEdphi = -dEdphi;
		// calculate vectors between v3 and v1, v4 and v2
		vector3 ca = torsion_data_holder_[i].atom3->get_position() - torsion_data_holder_[i].atom1->get_position();
		vector3 db = torsion_data_holder_[i].atom4->get_position() - torsion_data_holder_[i].atom2->get_position();
		vector3 dEdt = (double)(dEdphi/(t.length_2()*length_cb)) * (t%cb);
		vector3 dEdu = -(double)(dEdphi/(u.length_2()*length_cb)) * (u%cb);
		// calculate the forces on atom 1-4
		force_1 += dEdt % cb;
		force_2 += ca % dEdt + dEdu % dc;
		force_3 += dEdt % ba + db %dEdu;
		force_4 += dEdu % cb;
		/*force_1 += cb%dEdt;
		force_2 += ca%dEdt + dc%dEdu;
		force_3 += ba%dEdt + db%dEdu;
		force_4 += cb%dEdu;*/
		// update the forces
		torsion_data_holder_[i].atom1->set_force(force_1);
		torsion_data_holder_[i].atom2->set_force(force_2);
		torsion_data_holder_[i].atom3->set_force(force_3);
		torsion_data_holder_[i].atom4->set_force(force_4);
	}
}
