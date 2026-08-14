#include "../include/MMFF94Torsion.h"
#include "../include/MMFF94Bend.h"
#include "../include/utility.h"
#include "../include/Vector3.h"
#include "../include/MMFF94StretchBend.h"
#include <cmath>
#include<fstream>
#include<iostream>
using namespace std;

// constructors

MMFF94Torsion::MMFF94Torsion():FFComponent(),fast_access_tor(),torsion_data_holder_()
{
	// set component name
	set_name("MMFF94 Torsion");
}
MMFF94Torsion::MMFF94Torsion(ForceField& ff):FFComponent(ff),fast_access_tor(),torsion_data_holder_()
{
	set_name("MMFF94 Torsion");
}
// copy constructor
MMFF94Torsion::MMFF94Torsion(const MMFF94Torsion& to_copy):FFComponent(to_copy)
{
	fast_access_tor = to_copy.fast_access_tor;
	//num_of_atom_types_ = to_copy.num_of_atom_types_;
	torsion_data_holder_ = to_copy.torsion_data_holder_;
}
//destructor
MMFF94Torsion::~MMFF94Torsion()
{
	fast_access_tor.clear();
	torsion_data_holder_.clear();
}

// extract torsion parameters from FFParameter object bonded to the force field
// and establish a hash table for fast access
bool MMFF94Torsion::extract_TOR_parameters(FFParameter& ffp)
{
	if (!ffp.is_valid())
		return false;
	// build a two dim array of atom types and loop variable
	//FFParameter::AtomTypes& atom_types = ffp.get_atomtypes();
	int num_types = ffp.get_num_types();
	num_of_atom_types_ = num_types;
	fast_access_tor.clear();
	int type1, type2, type3, type4,  torsion_type;
	unsigned cxt;
	string name_type1, name_type2, name_type3, name_type4;
	int type_array[]={0,1,2,4,5};
	vector<int>type_vector (type_array,type_array+5);
	vector<int>::iterator t;
	for (vector<vector<string> >::iterator it = force_field_->get_parameters().params_in_each_section["torsion"].begin();it != force_field_->get_parameters().params_in_each_section["torsion"].end(); ++it)
	{
		torsion_type=str2int((*it)[0]);
		type1=str2int((*it)[1]);
		type2=str2int((*it)[2]);
		type3=str2int((*it)[3]);
		type4=str2int((*it)[4]);
		//debug
		//cout<<type1<<" "<<type2<<" "<<type3<<" "<<type4<<endl;
		bool type_found=false;
		t=find(type_vector.begin(),type_vector.end(),torsion_type);
		if(t!=type_vector.end())
		{
			type_found=true;
		}
		if ((type1 == -1) || (type2 == -1) || (type3 == -1) || (type4 == -1)||(!type_found))
		{
			cout<<"error! no numeric atom type defined for atom "<<type1<<" or "<<type2<<" or "<<type3<<" or "<<type4<<"or"<<torsion_type<<endl;
			return false;
		}
		cxt=mmff94_force_field->GetCXT(torsion_type,type1,type2,type3,type4);
		fast_access_tor[cxt].is_defined = true;
		fast_access_tor[cxt].V1 = str2double((*it)[5]);
		fast_access_tor[cxt].V2= str2double((*it)[6]);
		fast_access_tor[cxt].V3= str2double((*it)[7]);
		fast_access_tor[cxt].torsion_type= str2int((*it)[0]);
		cxt=mmff94_force_field->GetCXT(torsion_type,type4,type3,type2,type1);
		fast_access_tor[cxt].is_defined = true;
		fast_access_tor[cxt].V1 = str2double((*it)[5]);
		fast_access_tor[cxt].V2= str2double((*it)[6]);
		fast_access_tor[cxt].V3= str2double((*it)[7]);
		fast_access_tor[cxt].torsion_type= str2int((*it)[0]);
	}
	return true;
}
// query a set of torsion parameters has defined for a given combination of atom types and bond order of the bond between atom j and k
bool MMFF94Torsion::has_params(int i, int j, int k, int l, int torsion_type)
{
	unsigned cxt = mmff94_force_field->GetCXT(torsion_type,i,j,k,l);
	if ((i >= 0 && i <= num_of_atom_types_) && (j>=0 && j <= num_of_atom_types_) && (k>=0 && k <= num_of_atom_types_) && (l>=0 && l <= num_of_atom_types_) && (torsion_type >= 0 && torsion_type <= 5 && torsion_type!=3))
	{
		return fast_access_tor[cxt].is_defined;
	}
	else
		return false;
}
// assign the parameters for a given combination of atom types and torsion type 
bool MMFF94Torsion::assign_params(MMFF94Torsion::ForceValues& param, int i, int j, int k, int l, int torsion_type){
	unsigned cxt=mmff94_force_field->GetCXT(torsion_type,i,j,k,l);
	if (has_params(i, j, k, l, torsion_type))
	{
		param.V1=fast_access_tor[cxt].V1;
		param.V2=fast_access_tor[cxt].V2;
		param.V3=fast_access_tor[cxt].V3;
		return true;
	}
	else
		return false;
}

// set up method
bool MMFF94Torsion::setup()
{
	if (force_field_ == NULL)
	{
		cout<<"MMFF94Torsion::setup(): force field bound component can not be found"<<endl;
		return false;
	}
	// clear the parameter container
	torsion_data_holder_.clear();
	// tempararily set this component enabled;
	setenabled(true);

	mmff94_force_field = dynamic_cast<MMFF94*>(force_field_);
	if (mmff94_force_field == NULL || mmff94_force_field->has_initialized_param())
	{
		bool result = extract_TOR_parameters(force_field_->get_parameters());
		if (!result)
		{
			cout << "MMFF94Torsion::setup(): can not found angle Torsion section"<<endl;
			return false;
		}
	}
	//  MMFF94Torsion::ForceValues value;
	// retrieve all torsion parameters
	BONDVec::iterator biter1, biter2, biter3;
	ATOM* a1;
	ATOM* a2;
	ATOM* a3;
	ATOM* a4;
	// loop all atoms in the current molecule and proper torsions will be packed into the vector
	for (ATOMVec::iterator aiter = force_field_->get_atoms().begin();aiter != force_field_->get_atoms().end(); ++aiter)
	{
		for (biter1 = (*aiter)->get_bond_list().begin(); biter1 != (*aiter)->get_bond_list().end(); ++biter1)
		{
			// central atoms
			a2 = *aiter;
			if ((*biter1)->get_second_atom() == a2)
				a3 = const_cast<ATOM*>((*biter1)->get_first_atom());
			else if((*biter1)->get_first_atom() == a2)
				a3 = const_cast<ATOM*>((*biter1)->get_second_atom());
			for (biter2 = (*aiter)->get_bond_list().begin(); biter2 != (*aiter)->get_bond_list().end(); ++biter2)
			{
				if ((*biter2)->get_second_atom() != (*biter1)->get_second_atom())
				{
					// determine the first atom
					if ((*biter2)->get_first_atom() == *aiter)
						a1 = const_cast<ATOM*>((*biter2)->get_second_atom());
					else
						a1 = const_cast<ATOM*>((*biter2)->get_first_atom());
					for (biter3 = const_cast<ATOM*>((*biter1)->get_second_atom())->get_bond_list().begin();biter3 != const_cast<ATOM*>((*biter1)->get_second_atom())->get_bond_list().end() ; ++biter3)
					{
						if (!(*biter3)->has_atom(a2->get_id()))
						{
							// determine the fourth atom a4
							if ((*biter3)->get_first_atom() == a3)
								a4 = const_cast<ATOM*>((*biter3)->get_second_atom());
							else
								a4 = const_cast<ATOM*>((*biter3)->get_first_atom());
							//search torsion parameters for a1, a2, a3, a4
							int atom_type_1 = a1->get_mmff94_type();
							int atom_type_2 = a2->get_mmff94_type();
							int atom_type_3 = a3->get_mmff94_type();
							int atom_type_4 = a4->get_mmff94_type();
							// a provisional ForceValues container
							MMFF94Torsion::ForceValues value;
							//debug
							//cout<<a1->get_id()<<" "<<a2->get_id()<<" "<<a3->get_id()<<" "<<a4->get_id()<<endl;
						    int torsion_type = mmff94_force_field->GetTorsionType(a1,a2,a3,a4);
							//torsion_data_holder_.push_back(MMFF94Torsion::TorsionData());
							int cxt = mmff94_force_field->GetCXT(torsion_type,atom_type_1,atom_type_2,atom_type_3,atom_type_4);
							bool match=false; 
							//debug
							//cout<<atom_type_1<<" "<<atom_type_2<<" "<<atom_type_3<<" "<<atom_type_4<<endl;

							// try exact match
							if(has_params(atom_type_1,atom_type_2,atom_type_3,atom_type_4,torsion_type))
							{
								unsigned cxt = mmff94_force_field->GetCXT(torsion_type,atom_type_1,atom_type_2,atom_type_3,atom_type_4);
								value.V1=fast_access_tor[cxt].V1;
								value.V2=fast_access_tor[cxt].V2;
								value.V3=fast_access_tor[cxt].V3;
								value.hasDefined = true;
								match=true;
							}
							// try 3-2-2-5
							else if(has_params(mmff94bend.EqLvl3(atom_type_1),atom_type_2,atom_type_3,mmff94bend.EqLvl5(atom_type_4),torsion_type))
							{
								unsigned cxt = mmff94_force_field->GetCXT(torsion_type,mmff94bend.EqLvl3(atom_type_1),atom_type_2,atom_type_3,mmff94bend.EqLvl5(atom_type_4));
								value.V1=fast_access_tor[cxt].V1;
								value.V2=fast_access_tor[cxt].V2;
								value.V3=fast_access_tor[cxt].V3;
								value.hasDefined = true;
								match=true;
							}
							// try 5-2-2-3
							else if(has_params(mmff94bend.EqLvl5(atom_type_1),atom_type_2,atom_type_3,mmff94bend.EqLvl3(atom_type_4),torsion_type))
							{
								unsigned cxt = mmff94_force_field->GetCXT(torsion_type,(mmff94bend.EqLvl5(atom_type_1)),atom_type_2,atom_type_3,(mmff94bend.EqLvl3(atom_type_4)));
								value.V1=fast_access_tor[cxt].V1;
								value.V2=fast_access_tor[cxt].V2;
								value.V3=fast_access_tor[cxt].V3;
								value.hasDefined = true;
								match=true;
							}
							// try 5-2-2-5
							else if(has_params(mmff94bend.EqLvl5(atom_type_1),atom_type_2,atom_type_3,mmff94bend.EqLvl5(atom_type_4),torsion_type))
							{
								unsigned cxt = mmff94_force_field->GetCXT(torsion_type,(mmff94bend.EqLvl5(atom_type_1)),atom_type_2,atom_type_3,(mmff94bend.EqLvl5(atom_type_4)));
								value.V1=fast_access_tor[cxt].V1;
								value.V2=fast_access_tor[cxt].V2;
								value.V3=fast_access_tor[cxt].V3;
								value.hasDefined = true;
								match=true;
							}
									
							if(!match)
							{
								bool found_rule=false;
								// rule (a) page 631
								if (mmff94_force_field->HasLinSet(atom_type_2) ||mmff94_force_field-> HasLinSet(atom_type_3))
								// rule (b) page 631
								if(mmff94_force_field->get_mol()->GetBond(a2,a3)->is_aromatic())
								{
									double Ub, Uc, pi_bc, beta;
									Ub = mmff94_force_field->GetUParam(a2);
									Uc =mmff94_force_field->GetUParam(a3);
									if ((!mmff94_force_field->HasPilpSet(atom_type_2)) && ((!mmff94_force_field->HasPilpSet(atom_type_3))))
										pi_bc = 0.5;
									else
										pi_bc = 0.3;
									if (((mmff94_force_field->GetVal(atom_type_2) == 3) && (mmff94_force_field->GetVal(atom_type_3) == 4)) || ((mmff94_force_field->GetVal(atom_type_2) == 4) && (mmff94_force_field->GetVal(atom_type_3) == 3)))
										beta = 3.0;
									else
										beta = 6.0;
									value.V1 = 0.0;
									value.V2= beta * pi_bc * sqrt(Ub * Uc);
									value.V3=0.0;
									value.hasDefined = true;
									found_rule = true;
								}
								else
								{
									// rule (c) page 631	
									double Ub, Uc, pi_bc, beta;
									Ub = mmff94_force_field->GetUParam(a2);
									Uc = mmff94_force_field->GetUParam(a3);
									if (((mmff94_force_field->GetMltb(atom_type_2) == 2) && (mmff94_force_field->GetMltb(atom_type_3) == 2)) && (mmff94_force_field->get_mol()->GetBond(a1,a2)->is_double()))
										pi_bc = 1.0;
									else
										pi_bc = 0.4;
									beta = 6.0;
									value.V1=0.0;
									value.V2=beta * pi_bc * sqrt(Ub * Uc);
									value.V3=0.0;
									value.hasDefined = true;
									found_rule=true;
								}
								// rule (d) page 632
								if (!found_rule)
									if (((mmff94_force_field->GetCrd(atom_type_2) == 4) && (mmff94_force_field->GetCrd(atom_type_3) == 4)))
									{
										double Vb, Vc;
										Vb =mmff94_force_field->GetVParam(a2);
										Vc = mmff94_force_field->GetVParam(a3);
										value.V1=0.0;
										value.V2=0.0;
										value.V3=sqrt(Vb * Vc) / 9.0;
										value.hasDefined = true;
										found_rule=true;
									}
									// rule (e) page 632
									if (!found_rule)
										if (((mmff94_force_field->GetCrd(atom_type_2) == 4) && (mmff94_force_field->GetCrd(atom_type_3) != 4)))
										{
											if (mmff94_force_field->GetCrd(atom_type_3) == 3) // case (1)
												if ((mmff94_force_field->GetVal(atom_type_3) == 4) || (mmff94_force_field->GetVal(atom_type_3) == 34) || (mmff94_force_field->GetMltb(atom_type_3) != 0))
													continue;
											if (mmff94_force_field->GetCrd(atom_type_3) == 2) // case (2)
												if ((mmff94_force_field->GetVal(atom_type_3) == 3) || (mmff94_force_field->GetMltb(atom_type_3) != 0))
													continue;
										// case (3) saturated bonds -- see rule (h)
										}
									// rule (f) page 632
									if (!found_rule)
										if (((mmff94_force_field->GetCrd(atom_type_2) != 4) && (mmff94_force_field->GetCrd(atom_type_3) == 4)))
										{
											if (mmff94_force_field->GetCrd(atom_type_2) == 3) // case (1)
												if ((mmff94_force_field->GetVal(atom_type_2) == 4) || (mmff94_force_field->GetVal(atom_type_2) == 34) || (mmff94_force_field->GetMltb(atom_type_2) != 0))
													continue;
											if (mmff94_force_field->GetCrd(atom_type_2) == 2) // case (2)
												if ((mmff94_force_field->GetVal(atom_type_2) == 3) || (mmff94_force_field->GetMltb(atom_type_2) != 0))
													continue;
											// case (3) saturated bonds
										}
									// rule (g) page 632
									if (!found_rule)
										if (mmff94_force_field->get_mol()->GetBond(a2,a3)->is_single()&&((mmff94_force_field->GetMltb(atom_type_2) && mmff94_force_field->GetMltb(atom_type_3)) || (mmff94_force_field->GetMltb(atom_type_2) && mmff94_force_field->HasPilpSet(atom_type_3)) ||(mmff94_force_field->GetMltb(atom_type_3) && mmff94_force_field->HasPilpSet(atom_type_2))))
										{
											if (mmff94_force_field->HasPilpSet(atom_type_2) && mmff94_force_field->HasPilpSet(atom_type_3)) // case (1)
												continue;
											double Ub, Uc, pi_bc, beta;
											Ub = mmff94_force_field->GetUParam(a2);
											Uc = mmff94_force_field->GetUParam(a3);
											beta = 6.0;
											if (mmff94_force_field->HasPilpSet(atom_type_2) && mmff94_force_field->GetMltb(atom_type_3))
											{ // case (2)
												if (mmff94_force_field->GetMltb(atom_type_3) == 1)
													pi_bc = 0.5;
												else if ((mmff94stretchbend.GetElementRow(a2) == 1) && (mmff94stretchbend.GetElementRow(a3) == 1))
													pi_bc = 0.3;
												else
													pi_bc = 0.15;
												found_rule = true;
											}
											if (mmff94_force_field->HasPilpSet(atom_type_3) && mmff94_force_field->GetMltb(atom_type_2))
											{// case (3)
												if (mmff94_force_field->GetMltb(atom_type_2) == 1)
													pi_bc = 0.5;
												else if ((mmff94stretchbend.GetElementRow(a2) == 1) && (mmff94stretchbend.GetElementRow(a3) == 1))
												pi_bc = 0.3;
												else
													pi_bc = 0.15;
												found_rule = true;
											}
											if (!found_rule)
												if (((mmff94_force_field->GetMltb(atom_type_2) == 1) || (mmff94_force_field->GetMltb(atom_type_3) == 1)) && (!a2->is_carbon() || !a3->is_carbon()))
											{
												pi_bc = 0.4;
												found_rule = true;
											}	    
											if (!found_rule)
												pi_bc = 0.15;
											value.V1=0.0;
											value.V2= beta * pi_bc * sqrt(Ub * Uc);
											value.V3=0.0;
											value.hasDefined = true;
											found_rule=true;
										}	
										// rule (h) page 632
										if (!found_rule)
											if ((a2->is_oxygen() || a2->is_sulfur()) && (a3->is_oxygen() || a3->is_sulfur()))
											{
												double Wb, Wc;
												if (a2->is_oxygen())
												{
													Wb = 2.0;
												}
												else
												{
													Wb = 8.0;
												}
												if (a3->is_oxygen())
												{
													Wc = 2.0;
												}
												else
												{
													Wc = 8.0;
												}
												value.V1=0.0;
												value.V2=-sqrt(Wb * Wc);
												value.V3=0.0;
												value.hasDefined = true;
											}
											else
											{
												double Vb, Vc, Nbc;
												Vb = mmff94_force_field->GetVParam(a2);
												Vc = mmff94_force_field->GetVParam(a3);
												Nbc =mmff94_force_field-> GetCrd(atom_type_2) *mmff94_force_field-> GetCrd(atom_type_3);
												value.V1=0.0;
												value.V3=sqrt(Vb * Vc) / Nbc;
												value.V2=0.0;
												value.hasDefined = true;
											}
							}
							if(value.hasDefined)
							{
								MMFF94Torsion::TorsionData data;
								data.atom1 = a1;
								data.atom2 = a2;
								data.atom3 = a3;
								data.atom4 = a4;
								data.torsion_type = torsion_type;
								data.value = value;
								torsion_data_holder_.push_back(data);
							}
							else
							{
								cout<<"Warning: MMFF94Torsion::setup: We can't assign proper parameters for dihedral angle "<<a1->get_id()<<"-"<<a2->get_id()<<"-"<<a3->get_id()<<"-"<<a4->get_id()<<endl;
							}
						}
					}
				}  
			}
		}
	}	// everything goes well
	return true;
}

// update methods
double MMFF94Torsion::update_energy()
{
	//energy initializion
	energy_ = 0.0;
	if (torsion_data_holder_.size() == 0)
		return 0.0;
	// iterate all torsions and summarize the torsion energies
	//ofstream out ;
	//out.open("Torsion_test_log.log");
	//out<<"Atom types    "<<"FF Class  "<<"Torsion Angle   "<<"Steric Energy   "<<"Force Constant  "<<"V1  "<<"  V2 "" V3 "<<endl;
	for (vector<TorsionData>::size_type i = 0; i<torsion_data_holder_.size(); ++i)
	{

		vector3 v1, v2, v3, v4;
		vector3 a21,a23,a34;
		vector3 cross2321,cross2334;
		//debug
		//cout<<torsion_data_holder_[i].atom1->get_id()<<"-"<<torsion_data_holder_[i].atom2->get_id()<<"-"<<torsion_data_holder_[i].atom3->get_id()<<"-"<<torsion_data_holder_[i].atom4->get_id()<<"-"<<torsion_data_holder_[i].value.V1<<"-"<<torsion_data_holder_[i].value.V2<<"-"<<torsion_data_holder_[i].value.V3<<endl;
		v1 = torsion_data_holder_[i].atom1->get_position();
		v2 = torsion_data_holder_[i].atom2->get_position();
		v3 = torsion_data_holder_[i].atom3->get_position();
		v4 = torsion_data_holder_[i].atom4->get_position();
		double torsion_angle = CalcTorsionAngle(v1, v2, v3, v4);
        a21=v1-v2;
		a23=v3-v2;
		a34=v4-v3;
		cross2321 = a23 % a21;
		cross2334 = a23 % a34;
		double cosphi;
		double length_cross2321=cross2321.length();
		double length_cross2334=cross2334.length();
		if(length_cross2321 != 0 && length_cross2334 != 0)
		{
			cross2321 /= length_cross2321;
     		cross2334 /= length_cross2334;
			cosphi = cross2321 * cross2334;
			if (cosphi > 1.0)
			{
				cosphi = 1.0;
			}
			if (cosphi < -1.0)
			{
				cosphi = -1.0;
			}

			const double phi = torsion_angle*DEG_TO_RAD;
			double es = 0.0;
			es = 0.5 * (torsion_data_holder_[i].value.V1 * (1.0 + cos(phi)) +torsion_data_holder_[i].value.V2 * (1.0 - cos(phi * 2.0)) +	torsion_data_holder_[i].value.V3 * (1.0 + cos(phi * 3.0)));
			energy_ += es;
					
		}
	}
	return energy_;
}

// calculate current forces imposed by torsion and add them to the force field;
void MMFF94Torsion::update_forces()
{
	if (get_force_field() == 0)
	{
		cout<<"MMFF94Torsion::update_force(): error! this component doesn't bond to any force field"<<endl;
		return;
	}
	//ofstream out;
	////	out.open("Torsion_test_log.txt");
	//	out<<"Steric Force"<<"F_1 "<<"F_2 "<<"F_3 "<<"F_4 "<<"    Torsion angle "<<"    Bond Length"<<"  L_1 "<<" L_2"<<"  L_3 "<<endl;
	// iterate all torsions and update forces
	for (vector<TorsionData>::size_type i = 0; i<torsion_data_holder_.size(); ++i)
	{
		vector3 force_1 = torsion_data_holder_[i].atom1->get_force(), force_2 = torsion_data_holder_[i].atom2->get_force(), force_3 = torsion_data_holder_[i].atom3->get_force(),force_4 = torsion_data_holder_[i].atom4->get_force();
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
		double dEdphi = torsion_data_holder_[i].value.V1 * sin(torsion) - 2. * torsion_data_holder_[i].value.V2 * sin(torsion * 2.) + 3. * torsion_data_holder_[i].value.V3 * sin(torsion * 3.);
		dEdphi *= (-0.5 * FORCE_FACTOR);
		if ((t%u) * cb > 0.0)
			dEdphi = -dEdphi;
		// calculate vectors between v3 and v1, v4 and v2
		vector3 ca = torsion_data_holder_[i].atom3->get_position() - torsion_data_holder_[i].atom1->get_position();
		vector3 db = torsion_data_holder_[i].atom4->get_position() - torsion_data_holder_[i].atom2->get_position();
		vector3 dEdt = (double)(dEdphi/(t.length_2()*length_cb)) * (t%cb);
		vector3 dEdu = -(double)(dEdphi/(u.length_2()*length_cb)) * (u%cb);
		// calculate the forces on atom 1-4
		
		vector3 for_1,for_2,for_3,for_4;
		for_1 = dEdt % cb;
		for_2 = ca % dEdt + dEdu % dc;
		for_3 = dEdt % ba + db %dEdu;
		for_4 = dEdu % cb;
		//debug
			
		//out<<"            "<<for_1<<"  "<<for_2<<"  "<<for_3<<"  "<<for_4<<"  "<<torsion<<"   "<<"   "<<length_ab<<" "<<length_cb<<" "<<length_dc<<endl;
       
		force_1 += for_1;
		force_2 += for_2;
		force_3 += for_3;
		force_4 += for_4;
		// calculate the vectors between atom1 and atom2, atom2 and atom3, atom3 and atom4
		/*vector3 ij=torsion_data_holder_[i].atom2->get_position()-torsion_data_holder_[i].atom1->get_position();
		vector3 jk=torsion_data_holder_[i].atom3->get_position()-torsion_data_holder_[i].atom2->get_position();
		vector3 kl=torsion_data_holder_[i].atom4->get_position()-torsion_data_holder_[i].atom3->get_position();
		double r_ij=ij.length();
		double r_jk=jk.length();
		double r_kl=kl.length();
		if(isNearZero(r_ij)|| isNearZero(r_jk)|| isNearZero(r_kl))
		{
			continue;
		}
		double angle_ijk = vectorAngle(ij,jk);
		double angle_jkl = vectorAngle(jk,kl);
		//normalize the bond vectors
		ij.normalize();
		jk.normalize();
		kl.normalize();
		double sin_j = sin(angle_ijk);
		double sin_k = sin(angle_jkl);

		double rsj = r_ij * sin_j;
		double rsk = r_kl * sin_k;

		double rs2j = 1.0 / (rsj * sin_j);
		double rs2k = 1.0 / (rsk * sin_k);

		double rrj = r_ij / r_jk;
		double rrk = r_kl / r_jk;

		double rrcj = rrj * (-cos(angle_ijk));
		double rrck = rrk * (-cos(angle_jkl));

		vector3 a = ij % jk;
		vector3 b = jk % kl;
		vector3 c = a % b;
		double d1 = c * jk;
		double d2 = a * b;
		double angle = atan2(d1, d2);

		vector3 di = -a * rs2j;
		vector3 dl = b * rs2k;
		vector3 dj = di * (rrcj - 1.) - dl * rrck;
		vector3 dk = -(di + dj + dl);

		double c1 = torsion_data_holder_[i].value.V1 * sin(angle) - 2. * torsion_data_holder_[i].value.V2 * sin(angle * 2.) + 3. * torsion_data_holder_[i].value.V3 * sin(angle * 3.);
		c1 *= 0.5 * FORCE_FACTOR;
		force_1 += di * c1;
		force_2 += dj * c1;
		force_3 += dk * c1;
		force_4 += dl * c1;*/
		torsion_data_holder_[i].atom1->set_force(force_1);
		torsion_data_holder_[i].atom2->set_force(force_2);
		torsion_data_holder_[i].atom3->set_force(force_3);
		torsion_data_holder_[i].atom4->set_force(force_4);		
	}
	//out.close();
}
