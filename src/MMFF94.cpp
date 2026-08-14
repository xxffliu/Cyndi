#include "../include/MMFF94.h"
#include "../include/MMFF94Stretch.h"
#include "../include/MMFF94StretchBend.h"
#include "../include/MMFF94Bend.h"
#include "../include/MMFF94OOP.h"
#include "../include/MMFF94Torsion.h"
#include "../include/MMFF94VDW.h"
#include "../include/MMFF94Electrostatic.h"
#include <sstream>
using namespace std;
// default constructor
MMFF94::MMFF94():
param_is_initialized_(false),
file_name_(DEFAULT_MMFF94_PARAM_FILE)
{ 

	// set force name
	set_ff_name("MMFF94 ["+file_name_+"]");
	//creat the component list
	insert_component(new MMFF94Stretch(*this));
	insert_component(new MMFF94Bend(*this));
	insert_component(new MMFF94Str_Bend(*this));
	insert_component(new MMFF94OOP(*this));
	insert_component(new MMFF94Torsion(*this));
	insert_component(new MMFF94VDW(*this));
	insert_component(new MMFF94Ele(*this));
}
// constructor with a MOL
MMFF94::MMFF94 (MOL& mol):ForceField(),
param_is_initialized_(false),
file_name_(DEFAULT_MMFF94_PARAM_FILE)
{
	// set force name
	set_ff_name("MMFF94 ["+file_name_+"]");
	if(!mol.is_initialized())
		return;
	//creat the component list
	insert_component(new MMFF94Stretch(*this));
	insert_component(new MMFF94Bend(*this));
	insert_component(new MMFF94Str_Bend(*this));
	insert_component(new MMFF94Torsion(*this));
	insert_component(new MMFF94OOP(*this));
	insert_component(new MMFF94VDW(*this));
	insert_component(new MMFF94Ele(*this));
	//
	bool result = setup(mol);
	if(!result){
		cout<<"force field setup failed!"<<endl;
		valid_ = false;
	}
}
// special constructor with a MOL and only specified energy terms are included
MMFF94 ::MMFF94 (MOL& mol, vector<MMFF94COMPONENT> energy_terms_list):ForceField(),
param_is_initialized_(false),
file_name_(DEFAULT_MMFF94_PARAM_FILE){
	// set force name
	set_ff_name("MMFF94  ["+file_name_+"]");
	// loop the energy terms list and creat corresponding terms
	for(vector<MMFF94COMPONENT>::iterator it = energy_terms_list.begin();
		it != energy_terms_list.end(); ++it)
	switch(*it)
	{
		case MMFF94_BOND_STRETCH:
			insert_component(new MMFF94Stretch(*this));
			break;
		case MMFF94_ANGLE_BEND:
			insert_component(new MMFF94Bend(*this));
			break;
		case MMFF94_DIHEDRAL_TORSION:
			insert_component(new MMFF94Torsion(*this));
			break;
		case MMFF94_OOP_BEND:
			insert_component(new MMFF94OOP(*this));
			break;
		case MMFF94_VDW:
			insert_component(new MMFF94VDW(*this));
			break;
		case MMFF94_ELE:
			insert_component(new MMFF94Ele(*this));
			break;
		case MMFF94_STRETCH_BEND:
			insert_component(new MMFF94Str_Bend(*this));
			break;

			
	}
	//
	bool result = setup(mol);
	if(!result)
	{
		cout<<"force field setup failed!"<<endl;
		valid_ = false;
	}
}
// copy constructor
MMFF94::MMFF94(const MMFF94& mmff94):
ForceField(mmff94),
param_is_initialized_(mmff94.param_is_initialized_),
file_name_(mmff94.file_name_){}
// destructor
MMFF94::~MMFF94(){
}

void MMFF94::clear(){
	ForceField::clear();
	file_name_ = DEFAULT_MMFF94_PARAM_FILE;
	param_is_initialized_ = false;
}

// assignment operator
const MMFF94& MMFF94::operator=(const MMFF94& mmff94)
{
	if (this != &mmff94)
	{
		ForceField::operator=(mmff94);
		file_name_ = mmff94.file_name_;
		param_is_initialized_ = mmff94.param_is_initialized_;
	}
	return *this;
}

 bool MMFF94::specific_setup()
 {
	// check whether the molecule is aasigned
	if (get_mol() == 0)
		return false;
	//debug
#ifdef DEBUG
	cout<<"Atom Typering..."<<endl
#endif
	// bond all ff params with current FFParameter object 
	// assign the atom numeric type according to symbolic type implemented in mol2 file from the FFParameter object bonded with the FF
	parameter_.read_parameter(DEFAULT_MMFF94_PARAM_FILE);
	if ((!parameter_.is_initialized()) && (!parameter_.is_valid()))
	{
		parameter_.clear();
		param_is_initialized_ = false;
		return false;
	}
	else
	{
		param_is_initialized_ = true;
	}
	if(SetTypes())
		return true;
	else
		return false;
}

void  MMFF94::PerceiveAromatic()
{
	//bool done=false;
	vector<RING>::iterator _ring ;
	for(_ring=mol_->get_ring_vector().begin();_ring!=mol_->get_ring_vector().end();_ring++)
	{      
		if((*_ring).is_aromatic)
		{
		    ATOMVec::iterator atom_in_ring;
			BONDVec::iterator bond_in_ring;
			for(atom_in_ring=(*_ring).vatom.begin();atom_in_ring!=(*_ring).vatom.end();atom_in_ring++)
			{
				if(!((*atom_in_ring)->isAromaticRingAtom()))
				{
					(*atom_in_ring)->setAromaticRingAtom();
					//done=true;
				}
			}
			for(bond_in_ring=(*_ring).vbond.begin();bond_in_ring!=(*_ring).vbond.end();bond_in_ring++)
			{
				if(!((*bond_in_ring)->isAromaticRingBond()))
					(*bond_in_ring)->setAromaticRingBond();
			}					
		}		
	
	}	 
	return ;
}
  // Symbolic atom typing is skipped
  // 
  // atom typing is based on:
  //   MMFF94 I - Table III
  //   MMFF94 V - Table I
  //
MMFF94AtomType MMFF94::GetType(ATOM *atom)
{
	MMFF94AtomType type;
    BOND * bond;
    int oxygenCount, nitrogenCount, sulphurCount;
	//string doubleBondTo;
    ////////////////////////////////
    // Aromatic Atoms
    ////////////////////////////////
	if (atom->is_aromatic()||atom->isAromaticRingAtom())
	{
		// 5-member ring
		// fixed by xfliu, 20060226
		if (mol_->IsInAromaticRingSize(atom,5))
		{
			bool IsAromatic = true;;
			ATOMVec alphaPos, betaPos;
			ATOMVec alphaAtoms, betaAtoms;
			vector<int> alpha, beta;
			// fixed by xfliu, 20090227
			//if (atom->is_sulfur() && atom->get_num_bonded_atom() == 2)
			if(atom->is_sulfur())
			{
				type.numeric = 44;
				type.symbolic = "STHI";
				return type; // Aromatic 5-ring sulfur with pi lone pair (STHI)
			}
			else if (atom->is_oxygen())
			{
				type.numeric = 59;
				type.symbolic = "OFUR";
				return type; // Aromatic 5-ring oxygen with pi lone pair (OFUR)
			}
			else if (atom->is_nitrogen())
			{
				for ( ATOMVec::iterator aiter=atom->get_atom_list().begin();aiter!=atom->get_atom_list().end();aiter++)
				{
					if ((*aiter)->is_oxygen() && ((*aiter)->get_num_neighbor_bond()== 1))
					{
						type.numeric = 82;
						type.symbolic = "N5OX";
						return type; // N-oxide nitrogen in 5-ring alpha position, 
              //N-oxide nitrogen in 5-ring beta position, 
              // N-oxide nitrogen in other 5-ring  position, 
              // (N5AX, N5BX, N5OX) 
					}
				}
			}
			for(ATOMVec::iterator aiter=atom->get_atom_list().begin();aiter!=atom->get_atom_list().end();aiter++)
			{
				//if ( /*!(mol_->GetBond(atom, *aiter))->is_aromatic()) ||*/( !(mol_->GetBond(atom, *aiter))->isAromaticRingBond())|| !(mol_->IsInRingSize(*aiter,5)))
				if(!mol_->IsInRingSize(*aiter,5))
					continue;
				if(In_the_sameRing(atom,*aiter))
				{
					alphaPos.push_back(*aiter);
				}
				for(ATOMVec::iterator aiter2=(*aiter)->get_atom_list().begin();aiter2!=(*aiter)->get_atom_list().end();aiter2++)
				{
					if ((*aiter2)->get_id() == atom->get_id())
						continue;
					//if (/*!((mol_->GetBond(*aiter, *aiter2))->is_aromatic()) ||*/!((mol_->GetBond(*aiter, *aiter2))->isAromaticRingBond())|| !(mol_->IsInRingSize(*aiter,5)))
					if(!mol_->IsInRingSize(*aiter,5))
						continue;
					if (In_the_sameRing(atom, *aiter2))
					{
						betaPos.push_back(*aiter2);
					}
				}
			}
			if (IsAromatic)
			{
				for (unsigned int i = 0; i < alphaPos.size(); i++)
				{
					if (alphaPos[i]->is_sulfur())
					{
						alphaAtoms.push_back(alphaPos[i]);
						alpha.push_back(alphaPos[i]->get_atomic_num());
					}
					else if (alphaPos[i]->is_oxygen())
					{
						alphaAtoms.push_back(alphaPos[i]);
						alpha.push_back(alphaPos[i]->get_atomic_num());
					}
					else if (alphaPos[i]->is_nitrogen() && (alphaPos[i]->get_num_neighbor_bond() == 3))
					{
						bool IsNOxide = false;
						for(ATOMVec::iterator aiter=alphaPos[i]->get_atom_list().begin();aiter!=alphaPos[i]->get_atom_list().end();aiter++)
						{
							if ((*aiter)->is_oxygen() && ((*aiter) ->get_num_neighbor_bond() == 1))
							{
								IsNOxide = true;
							}
						}
						if (!IsNOxide)
						{
							alphaAtoms.push_back(alphaPos[i]);
							alpha.push_back(alphaPos[i]->get_atomic_num());
						}
					}
				}
				for (unsigned int i = 0; i < betaPos.size(); i++)
				{
					if (betaPos[i]->is_sulfur())
					{
						betaAtoms.push_back(betaPos[i]);
						beta.push_back(betaPos[i]->get_atomic_num());
					}
					else if (betaPos[i]->is_oxygen())
					{
						betaAtoms.push_back(betaPos[i]);
						beta.push_back(betaPos[i]->get_atomic_num());
					}
					else if (betaPos[i]->is_nitrogen() && (betaPos[i]->get_num_neighbor_bond() == 3))
					{
						bool IsNOxide = false;
						for(ATOMVec::iterator aiter=betaPos[i]->get_atom_list().begin();aiter!=betaPos[i]->get_atom_list().end();aiter++)
						{
							if ((*aiter)->is_oxygen() && ((*aiter)->get_num_neighbor_bond() == 1))
							{
								IsNOxide = true;
							}
						}
						if (!IsNOxide)
						{
							betaAtoms.push_back(betaPos[i]);
							beta.push_back(betaPos[i]->get_atomic_num());
						}
					}
				}
				if (!betaAtoms.size())
				{
					nitrogenCount = 0;
					for(ATOMVec::iterator aiter=atom->get_atom_list().begin();aiter!=atom->get_atom_list().end();aiter++)
					{
						if ((*aiter)->is_nitrogen() && ((*aiter)->get_num_neighbor_bond() == 3))
						{
							if (((*aiter)->BOSum() == 4) && ((*aiter)->is_aromatic()|| (*aiter)->isAromaticRingAtom()))
							{
								nitrogenCount++;
							}
							else if (((*aiter)->BOSum() == 3) && (!((*aiter)->is_aromatic()||(*aiter)->isAromaticRingAtom())||(!(*aiter)->isAromaticRingAtom())))
							{
								nitrogenCount++;
							}  
						}
					}
					if (nitrogenCount >= 2)
					{
						type.numeric = 80;
						type.symbolic = "CIM+";
						return type; // Aromatic carbon between N's in imidazolium (CIM+)
					}
				}
				// fixed by xfliu, 20090227
				if(!alphaAtoms.empty() && !betaAtoms.empty())
				{
					sort(alpha.begin(), alpha.end());
					reverse(alpha.begin(), alpha.end());
					sort(beta.begin(), beta.end());
					reverse(beta.begin(), beta.end());
					if(alpha[0] > beta[0])
					{
						betaAtoms.clear();
					}
					else if(alpha[0] < beta[0])
					{
						alphaAtoms.clear();
					}
				}

				if (alphaAtoms.empty() && betaAtoms.empty())
				{
					if (atom->is_carbon())
					{
              // there is no S:, O:, or N:
              // this is the case for anions with only carbon and nitrogen in the ring
						type.numeric = 78;
						type.symbolic = "C5";
						return type; // General carbon in 5-membered aromatic ring (C5)
					}
					else if (atom->is_nitrogen())
					{
						if (atom->get_num_neighbor_bond() == 3)
						{
						// this is the N: atom
							type.numeric = 39;
							type.symbolic = "NPYL";
							return type; // Aromatic 5 ring nitrogen with pi lone pair (NPYL)
						}
						else
						{
                // again, no S:, O:, or N:
							type.numeric = 76;
							type.symbolic = "N5M";
							return type; // Nitrogen in 5-ring aromatic anion (N5M)
						}
					}
				}
				if (alphaAtoms.size() == 2)
				{
					if (atom->is_carbon() && In_the_sameRing(alphaAtoms[0], alphaAtoms[1]))
					{
						if (alphaAtoms[0]->is_nitrogen() && alphaAtoms[1]->is_nitrogen())
						{
							if ((alphaAtoms[0]->get_num_neighbor_bond() == 3) && (alphaAtoms[1]->get_num_neighbor_bond() == 3))
							{
								type.numeric = 80;
								type.symbolic = "CIM+";
								return type; // Aromatic carbon between N's in imidazolium (CIM+)
							}
						}
					}
				}
				if (!alphaAtoms.empty() && betaAtoms.empty())
				{
					if (atom->is_carbon())
					{
						type.numeric = 63;
						type.symbolic = "C5A";
						return type; // Aromatic 5-ring C, alpha to N:, O:, or S: (C5A)
					}
					else if (atom->is_nitrogen())
					{
						if (atom->get_num_neighbor_bond() == 3)
						{
							type.numeric = 81;
							type.symbolic = "N5+";
							return type; // Posivite nitrogen in 5-ring alpha position (N5A+)
						}
						else
						{
							type.numeric = 65;
							type.symbolic = "N5A";
							return type; // Aromatic 5-ring N, alpha to N:, O:, or S: (N5A)
						}
					}
				}
				if (alphaAtoms.empty() && !betaAtoms.empty())
				{
					if (atom->is_carbon())
					{
						type.numeric = 64;
						type.symbolic = "C5B";
						return type; // Aromatic 5-ring C, beta to N:, O:, or S: (C5B)
					}
					else if (atom->is_nitrogen())
					{
						if (atom->get_num_neighbor_bond() == 3)
						{
							type.numeric = 81;
							type.symbolic = "N5+";
							return type; // Posivite nitrogen in 5-ring beta position (N5B+)
						}
						else
						{
							type.numeric = 66;
							type.symbolic = "N5B";
							return type; // Aromatic 5-ring N, beta to N:, O:, or S: (N5B)
						}
					}
				}
				if (!alphaAtoms.empty() && !betaAtoms.empty())
				{
					for (unsigned int i = 0; i < alphaAtoms.size(); i++)
					{
						for (unsigned int j = 0; j < betaAtoms.size(); j++)
						{
							// fixed by xfliu, 20090227
							if (In_the_sameRing(alphaAtoms[i], betaAtoms[j]))
							{
								if (atom->is_carbon())
								{
									type.numeric = 78;
									type.symbolic = "C5";
									return type; // General carbon in 5-membered aromatic ring (C5)
									//C5 is used if an atom is both C5A and C5B
								}
								else if (atom->is_nitrogen())
								{
									type.numeric = 79;
									type.symbolic = "N5";
									return type; // General nitrogen in 5-membered aromatic ring (N5)
								}
							}
						}
					}
					for (unsigned int i = 0; i < alphaAtoms.size(); i++)
					{
						if (alphaAtoms[i]->is_sulfur() || alphaAtoms[i]->is_oxygen())
						{
							if (atom->is_carbon())
							{
								type.numeric = 63;
								type.symbolic = "C5A";
								return type; // Aromatic 5-ring C, alpha to N:, O:, or S: (C5A)
							}
							else if (atom->is_nitrogen())
							{
								type.numeric = 65;
								type.symbolic = "N5A";
								return type; // Aromatic 5-ring N, alpha to N:, O:, or S: (N5A)
							}
						}
					}
					for (unsigned int i = 0; i < betaAtoms.size(); i++)
					{
						if (betaAtoms[i]->is_sulfur() || betaAtoms[i]->is_oxygen())
						{
							if (atom->is_carbon())
							{
								type.numeric = 64;
								type.symbolic = "C5B";
								return type; // Aromatic 5-ring C, beta to N:, O:, or S: (C5B)
							}
							else if (atom->is_nitrogen())
							{
								type.numeric = 66;
								type.symbolic = "N5B";
								return type; // Aromatic 5-ring N, beta to N:, O:, or S: (N5B)
							}
						}
					}
					if (atom->is_carbon())
					{
						type.numeric = 78;
						type.symbolic = "C5";
						return type; // General carbon in 5-membered aromatic ring (C5)
					}
					else if (atom->is_nitrogen())
					{
						type.numeric = 79;
						type.symbolic = "N5";
						return type; // General nitrogen in 5-membered aromatic ring (N5)
					}
				}
			}
		}
	
	// 6-member ring
		if (mol_->IsInRingSize(atom,6))
		{
			if (atom->is_carbon())
			{
				type.numeric = 37;
				type.symbolic = "CB";
				//cout<<atom->get_id()<<endl;
				return type; // Aromatic carbon, e.g., in benzene (CB)
			}
			else if (atom->is_nitrogen())
			{
				for ( ATOMVec::iterator aiter=atom->get_atom_list().begin();aiter!=atom->get_atom_list().end();aiter++)
				{
					if ((*aiter)->is_oxygen() && ((*aiter)->get_num_neighbor_bond() == 1))
					{
						type.numeric = 69;
						type.symbolic = "NPOX";
						return type; // Pyridinium N-oxide nitrogen (NPOX)
					}
				}
				if (atom->get_num_neighbor_bond() == 3)
				{
					type.numeric = 58;
					type.symbolic = "NPD+";
					return type; // Aromatic nitrogen in pyridinium (NPD+)
				}
				else
				{
					type.numeric = 38;
					type.symbolic = "NPYD";
					return type; // Aromatic nitrogen with sigma lone pair (NPYD)
				}
			}
		}
	}
    
    ////////////////////////////////
    // Hydrogen
    ////////////////////////////////
	if (atom->is_hydrogen())
	{
		for ( ATOMVec::iterator aiter=atom->get_atom_list().begin();aiter!=atom->get_atom_list().end();aiter++)
		{
			if ((*aiter)->is_carbon())
			{
				type.numeric = 5;
				type.symbolic = "HC";
				return type; // Hydrogen attatched to carbon (HC)
			}
			if((*aiter)->get_atomic_num()==14)
			{
				type.numeric = 5;
				type.symbolic = "HSI";
				return type;
			}
			if ((*aiter)->is_oxygen())
			{
				if ((*aiter)->BOSum() == 3)
				{
					if ((*aiter)->get_num_neighbor_bond() == 3)
					{
						type.numeric = 50;
						type.symbolic = "HO+";
						return type; // Hydrogen on oxonium oxygen (HO+)
					}
					else
					{
						type.numeric = 52;
						type.symbolic = "HO=+";
						return type; // Hydrogen on oxenium oxygen (HO=+)
					}
				} 
				int hydrogenCount = 0;
				for(ATOMVec::iterator aiter2=(*aiter)->get_atom_list().begin();aiter2!=(*aiter)->get_atom_list().end();aiter2++)
				{
					if ((*aiter2)->is_hydrogen())
					{
						hydrogenCount++;
						continue;
					}
					if ((*aiter2)->is_carbon())
					{
						// fixed by xfliu, 20090227
						if ((*aiter2)->is_aromatic()||(*aiter2)->isAromaticRingAtom())
						{
							type.numeric = 29;
							type.symbolic = "HOCC";
							return type; // phenol
						}
						for(ATOMVec::iterator aiter3=(*aiter2)->get_atom_list().begin();aiter3!=(*aiter2)->get_atom_list().end();aiter3++)
						{
							if ((*aiter3)->get_id() == (*aiter)->get_id())
								continue;
							bond = mol_->GetBond(*aiter2, *aiter3);
							if (bond->is_double())
							{
								if ((*aiter3)->is_oxygen())
								{
									type.numeric = 24;
									type.symbolic = "HOCO";
									return type; // Hydroxyl hydrogen in carboxylic acids (HOCO)
								}
								if ((*aiter3)->is_carbon() || (*aiter)->is_nitrogen())
								{
									type.numeric = 29;
									type.symbolic = "HOCC";
									return type;
								}
                    // Enolic or phenolic hydroxyl hydrogen,
                    // Hydroxyl hydrogen in HO-C=N moiety (HOCC, HOCN)
							}
						}
					}
					if ((*aiter2)->is_phosphorus())
					{
						type.numeric = 24;
						type.symbolic = "HOP";
						return type; // Hydroxyl hydrogen in H-O-P moiety (HOP)
					}
					if ((*aiter2)->is_sulfur())
					{
						type.numeric = 33;
						type.symbolic = "HOS";
						return type; // Hydrogen on oxygen attached to sulfur (HOS)
					}
				}
				if (hydrogenCount == 2)
				{
					type.numeric = 31;
					type.symbolic = "HOH";
					return type; // Hydroxyl hydrogen in water (HOH)
				}
				// if not the above cases 
				type.numeric = 21;
				type.symbolic = "HO";
				return type; // Hydroxyl hydrogen in alcohols, Generic hydroxyl hydrogen (HOR, HO)
			}
			if ((*aiter)->is_nitrogen())
			{
				switch (GetType((*aiter)).numeric)
				{
				case 81:
					type.numeric = 36;
					type.symbolic = "HIM+";
					return type; // Hydrogen on imidazolium nitrogen (HIM+)
				case 68:
					type.numeric = 23;
					type.symbolic = "HNOX";
					return type; // Hydrogen on N in N-oxide (HNOX)
				case 67:
					type.numeric = 23;
					type.symbolic = "HNOX";
					return type; // Hydrogen on N in N-oxide (HNOX)
				case 62:
					type.numeric = 23;
					type.symbolic = "HNR";
					return type; // Generic hydrogen on sp3 nitrogen, e.g., in amines (HNR)
				case 56:
					type.numeric = 36;
					type.symbolic = "HGD+";
					return type; // Hydrogen on guanimdinium nitrogen (HGD+)
				case 55:
					type.numeric = 36;
					type.symbolic = "HNN+";
					return type; // Hydrogen on amidinium nitrogen (HNN+)
				// fixed by xfliu, 20090227
				case 58:
					type.numeric = 36;
					type.symbolic = "NH+";
					return type; // Hydrogen on pyridinium-type N (FC=1) (aromatic)'
				case 43:
					type.numeric = 28;
					type.symbolic = "HN2";// general H on sp2 nitrogen
					return type; // Hydrogen on NSO2, or NSO3 nitrogen, Hydrogen on N triply bonded to C (HNSO, HNC%)
				case 39:
					type.numeric = 23;
					type.symbolic = "HPYL";
					return type; // Hydrogen on nitrogen in pyrrole (HPYL)
				case 8:
					type.numeric = 23;
					type.symbolic = "HN";
					return type; // Generic hydrogen on sp3 nitrogen, e.g., in amines, Hydrogen on nitrogen in ammonia (HNR, H3N)
				}
				if ((*aiter)->BOSum() == 4)
				{
					if ((*aiter)->get_num_neighbor_bond() == 2)
					{
						type.numeric = 28;
						type.symbolic = "HN2";
						return type; // Hydrogen on N triply bonded to C (HNC%)
					}
					else
					{
						type.numeric = 36;
						type.symbolic = "HN+"; // general H on any N+
						return type; // Hydrogen on pyridinium nitrogen, Hydrogen on protonated imine nitrogen (HPD+, HNC+)
					}
				}
				if ((*aiter)->get_num_neighbor_bond() == 2)
				{
					for(ATOMVec::iterator aiter2=(*aiter)->get_atom_list().begin();aiter2!=(*aiter)->get_atom_list().end();aiter2++)
					{
						if ((*aiter2)->is_hydrogen())
							continue;
						bond = mol_->GetBond((*aiter), (*aiter2));
						if (bond->is_double())
						{
							if ((*aiter2)->is_carbon() || (*aiter2)->is_nitrogen())
							{
								type.numeric = 27;
								type.symbolic = "HN=C";
								return type;
							}// Hydrogen on imine nitrogen, Hydrogen on azo nitrogen (HN=C, HN=N) 
						type.numeric = 28;
						type.symbolic = "HN2";
						return type; // Generic hydrogen on sp2 nitrogen (HSP2)
						}
					}
				}
				for(ATOMVec::iterator aiter2=(*aiter)->get_atom_list().begin();aiter2!=(*aiter)->get_atom_list().end();aiter2++)
				{
					if ((*aiter2)->is_hydrogen())
						continue;
					if ((*aiter2)->is_carbon())
					{
						if ((*aiter2)->is_aromatic()||(*aiter)->isAromaticRingAtom())
						{
							type.numeric = 28;
							type.symbolic = "HN2";
							return type; // deloc. lp pair
						}
						for(ATOMVec::iterator aiter3=(*aiter2)->get_atom_list().begin();aiter3!=(*aiter2)->get_atom_list().end();aiter3++)
						{
							if ((*aiter3)->get_id() == (*aiter)->get_id())
								continue;
							bond = mol_->GetBond((*aiter2), (*aiter3));
							if (bond->is_double())
							{
								if ((*aiter3)->is_carbon() || (*aiter3)->is_nitrogen() || (*aiter3)->is_oxygen() || (*aiter3)->is_sulfur())
								{
									type.numeric = 28;
									type.symbolic = "HN2";
									return type; // Hydrogen on amide nitrogen, Hydrogen on thioamide nitrogen,
								// Hydrogen on enamine nitrogen, Hydrogen in H-N-C=N moiety (HNCO, HNCS, HNCC, HNCN)
								}
							}
						}
					}
					if ((*aiter2)->is_nitrogen())
					{
						for(ATOMVec::iterator aiter3=(*aiter2)->get_atom_list().begin();aiter3!=(*aiter2)->get_atom_list().end();aiter3++)
						{
							if ((*aiter3)->get_id() == (*aiter)->get_id())
								continue;
							bond = mol_->GetBond((*aiter2), (*aiter3));
							if (bond->is_double())
							{
								if ((*aiter3)->is_carbon() || (*aiter3)->is_nitrogen())
								{
									type.numeric = 28;
									type.symbolic = "HN2";
									return type; // Hydrogen in H-N-N=C moiety, Hydrogen in H-N-N=N moiety (HNNC, HNNN)
								}
							}
						}
					}
					if ((*aiter2)->is_sulfur())
					{
						for(ATOMVec::iterator aiter3=(*aiter2)->get_atom_list().begin();aiter3!=(*aiter2)->get_atom_list().end();aiter3++)
						{
							if ((*aiter3)->get_id() == (*aiter)->get_id())
								continue;
							if ((*aiter3)->is_oxygen() || ((*aiter3)->get_num_neighbor_bond() == 1))
							{
								type.numeric = 28;
								type.symbolic = "HN2";
								return type; // Hydrogen on NSO, NSO2 or NSO3 nitrogen (HNSO)
							}
						}
					}
				}
				type.numeric = 23;
				type.symbolic = "HN";
				return type; // Generic hydrogen on sp3 nitrogen e.g., in amines, Hydrogen on nitrogen in pyrrole, Hydrogen in ammonia,Hydrogen on N in N-oxide (HNR, HPYL, H3N, HNOX)
			}
			if ((*aiter)->is_sulfur() || (*aiter)->is_phosphorus())
			{
				type.numeric = 71;
				type.symbolic = "HP";
				return type; // Hydrogen attached to sulfur, Hydrogen attached to >S= sulfur doubly bonded to N,
          // Hydrogen attached to phosphorus (HS, HS=N, HP)
			}
		}
	}

    ////////////////////////////////
    // Carbon
    ////////////////////////////////
	if (atom->is_carbon())
	{
      // 4 neighbours
		if (atom->get_num_neighbor_bond() == 4)
		{
			if (mol_->IsInRingSize(atom,3))
			{
				type.numeric = 22;
				type.symbolic = "CR3R";
				// Aliphatic carbon in 3-membered ring (CR3R)
			} 
			else if (mol_->IsInRingSize(atom,4))
			{
				type.numeric = 20;
				type.symbolic = "CR4R";
				// Aliphatic carbon in 4-membered ring (CR4R)
			}
			else
			{
				type.numeric = 1;
				type.symbolic = "CR";
			}
			// Alkyl carbon (CR)
		}
      // 3 neighbours
 		else if (atom->get_num_neighbor_bond() == 3)
		{
			int N2count = 0;
			int N3count = 0;
			oxygenCount = sulphurCount =  0;
			string doubleBondTo;
			for ( ATOMVec::iterator aiter=atom->get_atom_list().begin();aiter!=atom->get_atom_list().end();aiter++)
			{
				bond = mol_->GetBond(*aiter, atom);
				if (bond->is_double())
				{
					doubleBondTo = (*aiter)->get_element();
				} 
				if ((*aiter)->get_num_neighbor_bond() == 1)
				{
					if ((*aiter)->is_oxygen())
					{
						oxygenCount++;
					}
					else if ((*aiter)->is_sulfur())
					{
						sulphurCount++;
					}
				}
				else if ((*aiter)->get_num_neighbor_bond() == 3)
				{
					if ((*aiter)->is_nitrogen())
					{
						N3count++;
					}
				}
				else if (((*aiter)->get_num_neighbor_bond() == 2) && bond->is_double())
				{
					if ((*aiter)->is_nitrogen())
					{
						N2count++;
					}
				}
			}
			if ((N3count >= 2) && (doubleBondTo == "N") && !N2count)
			{
				// N3==C--N3
				type.numeric = 57;
				type.symbolic = "CNN+";
				// Guanidinium carbon, Carbon in +N=C-N: resonance structures (CGD+, CNN+)
			}
			else if ((oxygenCount == 2) || (sulphurCount == 2))
			{
				// O1-?-C-?-O1 or S1-?-C-?-S1
				type.numeric = 41;
				type.symbolic = "CO2M";
				// Carbon in carboxylate anion, Carbon in thiocarboxylate anion (CO2M, CS2M)
			}
			else if (mol_->IsInRingSize(atom,4) && (doubleBondTo == "C"))
			{
				type.numeric = 30;
				type.symbolic = "CR4E";
				// Olefinic carbon in 4-membered ring (CR4E)
			}
			else if ((doubleBondTo == "N") || (doubleBondTo =="O") || 
            (doubleBondTo == "P") || (doubleBondTo == "S"))
			{
				// C==N, C==O, C==P, C==S
				type.numeric = 3;
				type.symbolic = "C=";
				// Generic carbonyl carbon, Imine-type carbon, Guanidine carbon,
				// Ketone or aldehyde carbonyl carbon, Amide carbonyl carbon,
				// Carboxylic acid or ester carbonyl carbon, Carbamate carbonyl carbon,
				// Carbonic acid or ester carbonyl carbon, Thioester carbonyl (double
				// bonded to O or S), Thioamide carbon (double bonded to S), Carbon
				// in >C=SO2, Sulfinyl carbon in >C=S=O, Thiocarboxylic acid or ester 
				// carbon, Carbon doubly bonded to P (C=O, C=N, CGD, C=OR, C=ON, COO,
				// COON, COOO, C=OS, C=S, C=SN, CSO2, CS=O, CSS, C=P)
			}
			else
			{
				type.numeric = 2;
				type.symbolic = "CSP2";
			}
			// Vinylic Carbon, Generic sp2 carbon (C=C, CSP2)
		}
      // 2 neighbours
		else if (atom->get_num_neighbor_bond() == 2)
		{
			type.numeric = 4;
			type.symbolic = "CSP";
			// Acetylenic carbon, Allenic caron (CSP, =C=)
		}
      // 1 neighbours
		else if (atom->get_num_neighbor_bond() == 1)
		{
			type.numeric = 60;
			type.symbolic = "C%";
			// Isonitrile carbon (C%-)
		}
		return type;
	}

    ////////////////////////////////
    // Nitrogen
    ////////////////////////////////
	if (atom->is_nitrogen())
	{
      // 4 neighbours
		if (atom->get_num_neighbor_bond() == 4)
		{
			for ( ATOMVec::iterator aiter=atom->get_atom_list().begin();aiter!=atom->get_atom_list().end();aiter++)
			{
				if ((*aiter)->is_oxygen() && ((*aiter)->get_num_neighbor_bond() == 1))
				{
					type.numeric = 68;
					type.symbolic = "N3OX";
					return type; // sp3-hybridized N-oxide nitrogen (N3OX)
				}
			}
			type.numeric = 34;
			type.symbolic = "NR+";
			return type; // Quaternary nitrogen (NR+)
      }
      // 3 neighbours
		else if (atom->get_num_neighbor_bond() == 3)
		{
			// fixed by xfliu, 20060226
			if (atom->BOSum() >= 4)
			{
				oxygenCount = nitrogenCount =0;
				int doubleBondTo=0;
				for ( ATOMVec::iterator aiter=atom->get_atom_list().begin();aiter!=atom->get_atom_list().end();aiter++)
				{
					if ((*aiter)->is_oxygen() && ((*aiter)->get_num_neighbor_bond() == 1))
					{
						oxygenCount++;
					}
					else if ((*aiter)->is_nitrogen())
					{
						bond = mol_->GetBond(*aiter, atom);
						if (bond->is_double())
						{
							doubleBondTo = 7;
						}
					}
					else if ((*aiter)->is_carbon())
					{
						bond = mol_->GetBond(*aiter, atom);
						if (bond->is_double())
						{
							for(ATOMVec::iterator aiter2=(*aiter)->get_atom_list().begin();aiter2!=(*aiter)->get_atom_list().end();aiter2++)
							{
								if ((*aiter2)->is_nitrogen() && ((*aiter2)->get_num_neighbor_bond() == 3))
								{
									nitrogenCount++;
								}
							}
						}
					}
				}
				if (oxygenCount == 1)
				{
					type.numeric = 67;
					type.symbolic = "N2OX";
					//return type; // sp2-hybridized N-oxide nitrogen (N2OX)
				}
				else if (oxygenCount >= 2)
				{
					type.numeric = 45;
					type.symbolic = "NO3";
					//return 45; // Nitrogen in nitro group, Nitrogen in nitrate group (NO2, NO3)
				}
				else if (nitrogenCount == 1)
				{
					type.numeric = 54;
					type.symbolic = "N+=C";
					//return 54; // Iminium nitrogen (N+=C)
				}
				else if (nitrogenCount == 2)
				{
					type.numeric = 55;
					type.symbolic = "NCN+";
					//return 55; // Either nitrogen in N+=C-N: (NCN+)
				}
				else if (nitrogenCount == 3)
				{
					type.numeric = 56;
					type.symbolic = "NGD+";
					//return 56; // Guanidinium nitrogen (NGD+)
				}
				else if (doubleBondTo == 7)
				{
					type.numeric = 54;
					type.symbolic = "H+=";
					//return 54; // Positivly charged nitrogen doubly bonded to nitrogen (N+=N)
				}
				return type;
			}
			else if (atom->BOSum() == 3)
			{
				bool IsAmide = false;
				bool IsSulfonAmide = false;
				bool IsNNNorNNC = false;
				string tripleBondTo ;
				string  doubleBondTo ;
				for ( ATOMVec::iterator aiter=atom->get_atom_list().begin();aiter!=atom->get_atom_list().end();aiter++)
				{
					if ((*aiter)->is_sulfur() || (*aiter)->is_phosphorus())
					{
						oxygenCount = 0;
						for(ATOMVec::iterator aiter2=(*aiter)->get_atom_list().begin();aiter2!=(*aiter)->get_atom_list().end();aiter2++)
						{
							if((*aiter2)->is_oxygen() && ((*aiter2)->get_num_neighbor_bond() == 1))
							{
								oxygenCount++;
							}
						}
						if (oxygenCount >= 2)
						{
							IsSulfonAmide = true;
                //return 43; // Sulfonamide nitrogen (NSO2, NSO3)
						}
					}
				}
				for ( ATOMVec::iterator aiter=atom->get_atom_list().begin();aiter!=atom->get_atom_list().end();aiter++)
				{
					if ((*aiter)->is_carbon())
					{
						for(ATOMVec::iterator aiter2=(*aiter)->get_atom_list().begin();aiter2!=(*aiter)->get_atom_list().end();aiter2++)
						{
							bond = mol_->GetBond((*aiter), (*aiter2));
							if (bond->is_double() && ((*aiter2)->is_oxygen() || (*aiter2)->is_sulfur()))
							{
								IsAmide = true;
                  //return 10; // Amide nitrogen, Thioamide nitrogen (NC=O, NC=S)
							}
						}
					}
				}
				for ( ATOMVec::iterator aiter=atom->get_atom_list().begin();aiter!=atom->get_atom_list().end();aiter++)
				{
					if ((*aiter)->is_carbon())
					{
						int N2count = 0;
						int N3count = 0;
						oxygenCount = sulphurCount = 0;for(ATOMVec::iterator aiter2=(*aiter)->get_atom_list().begin();aiter2!=(*aiter)->get_atom_list().end();aiter2++)
						{
							bond = mol_->GetBond(*aiter, *aiter2);
							if (bond->is_double())
							{
								doubleBondTo = (*aiter2)->get_element();
							}
							if (bond->is_aromatic()||bond->isAromaticRingBond())
							{
								if (((*aiter2)->get_element() == "N") || ((*aiter2)->get_element() == "C"))
								{
									doubleBondTo = (*aiter2)->get_element();
								}
							}
							if (bond->is_triple())
							{
								tripleBondTo = (*aiter2)->get_element();
							}
							if ((*aiter2)->is_nitrogen() && ((*aiter2)->get_num_neighbor_bond() == 3))
							{
								int aiterOxygen = 0;
								for(ATOMVec::iterator aiter3=(*aiter2)->get_atom_list().begin();aiter3!=(*aiter2)->get_atom_list().end();aiter3++)
								{
									if ((*aiter3)->is_oxygen())
									{
										aiterOxygen++;
									}
								}
								if (aiterOxygen < 2)
								{
									N3count++;
								}
							}
							if ((*aiter2)->is_nitrogen() && ((*aiter2)->get_num_neighbor_bond() == 2) && (bond->is_double() || bond->is_aromatic()||bond->isAromaticRingBond()))
							{
								N2count++;
							}
							if ((*aiter2)->is_aromatic()||(*aiter2)->isAromaticRingAtom())
							{
								if ((*aiter2)->is_oxygen())
								{
									oxygenCount++;
								}
								else if ((*aiter2)->is_sulfur())
								{
									sulphurCount++;
								}
							}
						}
						if (N3count == 3)
						{
							type.numeric = 56;
							type.symbolic = "NGD+";
							return type; // Guanidinium nitrogen (NGD+)
						}
						else if (!IsAmide && !IsSulfonAmide && !oxygenCount && !sulphurCount && ((*aiter)->is_aromatic()||(*aiter)->isAromaticRingAtom()))
						{
							type.numeric = 40;
							type.symbolic = "NC=C";
							return type;
						}
						else if ((N3count == 2) && (doubleBondTo == "N") && !N2count)
						{
							type.numeric = 55;
							type.symbolic = "NCN+";
							return type; // Either nitrogen in N+=C-N: (NCN+)
						}
					}
					if ((*aiter)->is_nitrogen())
					{
						nitrogenCount = 0;
						for(ATOMVec::iterator aiter2=(*aiter)->get_atom_list().begin();aiter2!=(*aiter)->get_atom_list().end();aiter2++)
						{
							bond = mol_->GetBond(*aiter, *aiter2);
							if (bond->is_double())
							{
								if ((*aiter2)->is_carbon())
								{
									oxygenCount = sulphurCount = 0;
									for(ATOMVec::iterator aiter3=(*aiter2)->get_atom_list().begin();aiter3!=(*aiter2)->get_atom_list().end();aiter3++)
									{
										if ((*aiter3)->is_oxygen())
										{
											oxygenCount++;
										}
										else if ((*aiter3)->is_sulfur())
										{
											sulphurCount++;
										}
										else if ((*aiter3)->is_sulfur())
										{
											nitrogenCount++;
										}
									}
									if (!oxygenCount && !sulphurCount && (nitrogenCount == 1))
									{
										bool bondToAromC = false;
										for ( ATOMVec::iterator aiter4=atom->get_atom_list().begin();aiter4!=atom->get_atom_list().end();aiter4++)
										{
											if ((( *aiter4)->is_aromatic()||( *aiter4)->isAromaticRingAtom()) && (*aiter4)->is_carbon() && mol_->IsInRingSize(*aiter4,6))
											{
												bondToAromC = true;
											}
										}
										if (!bondToAromC)
										{
											IsNNNorNNC = true;
										}
									}
								}
								else if ((*aiter2)->is_nitrogen())
								{
									bool bondToAromC = false;
									for ( ATOMVec::iterator aiter4=atom->get_atom_list().begin();aiter4!=atom->get_atom_list().end();aiter4++)
									{
										if (((*aiter4)->is_aromatic()||(*aiter4)->isAromaticRingAtom()) &&(*aiter4)->is_carbon() &&mol_->IsInRingSize(*aiter4,6))
										{
											bondToAromC = true;
										}
									}
									if (!bondToAromC)
									{
										IsNNNorNNC = true;
									}
								}
							}
						}
					}
				}
				if (IsSulfonAmide)
				{
					type.numeric = 43;
					type.symbolic = "NC#N";
					return type; // Sulfonamide nitrogen (NSO2, NSO3)
				}
				else if (IsAmide)
				{
					type.numeric = 10;
					type.symbolic = "NC=O";
					return type; // Amide nitrogen, Thioamide nitrogen (NC=O, NC=S)
				}
				else if ((doubleBondTo == "C") || (doubleBondTo == "N") ||(doubleBondTo == "P") || (tripleBondTo == "C"))
				{
					type.numeric = 40;
					type.symbolic = "NC=C";
					return type; // Enamine or aniline nitrogen (deloc. lp), Nitrogen in N-C=N with deloc. lp,
            // Nitrogen in N-C=N with deloc. lp, Nitrogen attached to C-C triple bond
            // (NC=C, NC=N, NC=P, NC%C)
				}
				else if (tripleBondTo == "N")
				{
					type.numeric = 43;
					type.symbolic = "NC#C";
					return type; // Nitrogen attached to cyano group (NC%N)
				}
				else if (IsNNNorNNC)
				{
					type.numeric = 10;
					type.symbolic = "NC=O";
					return type; // Nitrogen in N-N=C moiety with deloc. lp
            // Nitrogen in N-N=N moiety with deloc. lp (NN=C, NN=N)
				}
				else
				{
					type.numeric = 8;
					type.symbolic = "NR";
					return type; // Amine nitrogen (NR)
				}
			}
		}
      // 2 neighbours
		else if (atom->get_num_neighbor_bond() == 2)
		{
			if (atom->BOSum() == 4)
			{
				for ( ATOMVec::iterator aiter=atom->get_atom_list().begin();aiter!=atom->get_atom_list().end();aiter++)
				{
					bond = mol_->GetBond(*aiter, atom);
					if (bond->is_triple())
					{
						type.numeric = 61;
						type.symbolic = "NR%";
						return type; // Isonitrile nitrogen (NR%)
					}
				}
				type.numeric = 53;
				type.symbolic = "=N=";
				return type; // Central nitrogen in C=N=N or N=N=N (=N=)
			}
			else if (atom->BOSum() == 3)
			{
				//doubleBondTo = "/0";
				for ( ATOMVec::iterator aiter=atom->get_atom_list().begin();aiter!=atom->get_atom_list().end();aiter++)
				{
					bond = mol_->GetBond(*aiter, atom);
					if ((*aiter)->is_oxygen() && bond->is_double() && ((*aiter)->get_num_neighbor_bond() == 1))
					{
						type.numeric = 46;
						type.symbolic = "N=O";
						return type; // Nitrogen in nitroso group (N=O)
					}
					else if (((*aiter)->is_carbon() || (*aiter)->is_nitrogen()) && bond->is_double())
					{
						type.numeric = 9;
						type.symbolic = "N=C";
						return type; // Iminie nitrogen, Azo-group nitrogen (N=C, N=N)
					}
				}
				for ( ATOMVec::iterator aiter=atom->get_atom_list().begin();aiter!=atom->get_atom_list().end();aiter++)
				{
					if ((*aiter)->is_sulfur())
					{
						oxygenCount = 0;
						for(ATOMVec::iterator aiter2=(*aiter)->get_atom_list().begin();aiter2!=(*aiter)->get_atom_list().end();aiter2++)
						{
							if ((*aiter2)->is_oxygen() && ((*aiter2)->get_num_neighbor_bond() == 1))
							{
								oxygenCount++;
							}
						}
						if (oxygenCount >= 2)
						{
							type.numeric = 43;
							type.symbolic = "NSO2";
							// Sulfonamide nitrogen (NSO2, NSO3)
						}
						// fixed by xfliu, 20090227
						// it's weired....
						else if(atom->get_symbol_type() == "N.2" && oxygenCount == 0)
						{
							type.numeric = 62;
							type.symbolic = "NM";
							// Divalent nitrogen replacing monovalent O in SO2 group (NSO)
						}
						else if(atom->get_symbol_type() == "N.2" && oxygenCount == 1)
						{
							type.numeric = 48;
							type.symbolic = "NSO";
							// Divalent nitrogen replacing monovalent O in SO2 group (NSO)
						}
					}
				}
				return type; 
			}
			else if (atom->BOSum() == 2)
			{
				oxygenCount = sulphurCount = 0;
				for ( ATOMVec::iterator aiter=atom->get_atom_list().begin();aiter!=atom->get_atom_list().end();aiter++)
				{
					if ((*aiter)->is_sulfur())
					{
						for(ATOMVec::iterator aiter2=(*aiter)->get_atom_list().begin();aiter2!=(*aiter)->get_atom_list().end();aiter2++)
						{
							if ((*aiter2)->is_oxygen() && ((*aiter2)->get_num_neighbor_bond() == 1))
							{
								oxygenCount++;
							}
						}
						if (oxygenCount == 1)
						{
							type.numeric = 48;
							type.symbolic = "NSO";
							return type; // Divalent nitrogen replacing monovalent O in SO2 group (NSO)
						}
					}
				}
				type.numeric = 62;
				type.symbolic = "NM";
				return type; // Anionic divalent nitrogen (NM)
			}
		}
      // 1 neighbours
		else if (atom->get_num_neighbor_bond() == 1)
		{
			for ( ATOMVec::iterator aiter=atom->get_atom_list().begin();aiter!=atom->get_atom_list().end();aiter++)
			{
				bond = mol_->GetBond(*aiter, atom);
				if (bond->is_triple())
				{
					type.numeric = 42;
					type.symbolic = "NSP";
					return type; // Triply bonded nitrogen (NSP)
				}
				else if ((*aiter)->is_nitrogen() && ((*aiter)->get_num_neighbor_bond() == 2))
				{
					type.numeric = 47;
					type.symbolic = "NAZT";
					return type; // Terminal nitrogen in azido or diazo group (NAZT)
				}
			}
		}
	}
    ////////////////////////////////
    // Oxygen
    ////////////////////////////////
	if (atom->is_oxygen())
	{
      // 3 neighbours
		if (atom->get_num_neighbor_bond() == 3)
		{
			type.numeric = 49;
			type.symbolic = "O+";
			return type; // Oxonium oxygen (O+)
		}
      // 2 neighbours
		else if (atom->get_num_neighbor_bond() == 2)
		{
			int hydrogenCount = 0;
			for ( ATOMVec::iterator aiter=atom->get_atom_list().begin();aiter!=atom->get_atom_list().end();aiter++)
			{
				if ((*aiter)->is_hydrogen())
				{
					hydrogenCount++;
				}
			}
			if (hydrogenCount == 2)
			{
				// H--O--H
				type.numeric = 70;
				type.symbolic = "OH2";
				return type; // Oxygen in water (OH2)
			}
			else if (atom->BOSum() == 3)
			{
				type.numeric = 51;
				type.symbolic = "O=+";
				return type; // Oxenium oxygen (O=+)
			}
			else
			{
				type.numeric = 6;
				type.symbolic = "O";
				return type; // Generic divalent oxygen, Ether oxygen, Carboxylic acid or ester oxygen,
        // Enolic or phenolic oxygen, Oxygen in -O-C=N- moiety, Divalent oxygen in
        // thioacid or ester, Divalent nitrate "ether" oxygen, Divalent oxygen in
        // sulfate group, Divalent oxygen in sulfite group, One of two divalent
        // oxygens attached to sulfur, Divalent oxygen in R(RO)S=O, Other divalent
        // oxygen attached to sulfur, Divalent oxygen in phosphate group, Divalent
        // oxygen in phosphite group, Divalent oxygen (one of two oxygens attached
        // to P), Other divalent oxygen (-O-, OR, OC=O, OC=C, OC=N, OC=S, ONO2, 
        // ON=O, OSO3, OSO2, OSO, OS=O, -OS, OPO3, OPO2, OPO, -OP)

        // 59 ar
			}
		}
      // 1 neighbour
		else if (atom->get_num_neighbor_bond() == 1)
		{
			oxygenCount = sulphurCount = 0;
			for ( ATOMVec::iterator aiter=atom->get_atom_list().begin();aiter!=atom->get_atom_list().end();aiter++)
			{
				bond = mol_->GetBond(*aiter, atom);
				//fixed by fbai,20090228
				if ((*aiter)->is_carbon() || (*aiter)->is_nitrogen()||(*aiter)->is_sulfur())
				{
					for(ATOMVec::iterator aiter2=(*aiter)->get_atom_list().begin();aiter2!=(*aiter)->get_atom_list().end();aiter2++)
					{
						if ((*aiter2)->is_oxygen() && ((*aiter2)->get_num_neighbor_bond() == 1))
						{
							oxygenCount++;
						}
						else if ((*aiter2)->is_sulfur() && ((*aiter2)->get_num_neighbor_bond() == 1))
						{
							sulphurCount++;
						}
					}
				}
				// O---H

				if ((*aiter)->is_hydrogen())
				{
					type.numeric = 35;
					type.symbolic = "OM";
					return type;
				}
				// O-?-C
				else if ((*aiter)->is_carbon())
				{
					if (oxygenCount == 2)
					{
						// O-?-C-?-O
						type.numeric = 32;
						type.symbolic = "O2CM";
						//return 32; // Oxygen in carboxylate group (O2CM)
					}
					else if (bond->is_single())
					{ 
						// O--C
						type.numeric = 35;
						type.symbolic = "OM";
						//return 35; // Oxide oxygen on sp3 carbon, Oxide oxygen on sp2 carbon (OM, OM2)
					}
					else
					{ 
						// O==C
						type.numeric = 7;
						type.symbolic = "O=";
						//return 7; // Generic carbonyl oxygen, Carbonyl oxygen in amides,
              // Carbonyl oxygen in aldehydes and ketones, Carbonyl
              // oxygen in acids or esters (O=C, O=CN, O=CR, O=CO)
					}
					return type;
				}
				// O-?-N
				else if ((*aiter)->is_nitrogen())
				{
					if (oxygenCount >= 2)
					{
						// O-?-N-?-O
						type.numeric = 32;
						type.symbolic = "OX";
						//return 32; // Oxygen in nitro group, Nitro-group oxygen in nitrate,
              // Nitrate anion oxygen (O2N, O2NO, O3N)
					}
					// fixed by xfliu, 20090227
					else if (bond->is_single())
					{
						if((*aiter)->is_ring && ((*aiter)->get_symbol_type() == "N.2" || (*aiter)->get_symbol_type() == "N.ar"))
						{
							// O--N
							type.numeric = 32;
							type.symbolic = "OX"; // Oxygen in N-oxides (ONX)
						}
						// fixed by xfliu, 20090227
						else if((*aiter)->get_mmff94_symbol_type() == "N3OX" || (*aiter)->get_mmff94_symbol_type() == "N2OX")
						{
							type.numeric = 32;
							type.symbolic = "OX"; // Oxygen in N-oxides (ONX)
						}
						else
						{
							type.numeric = 35;
							type.symbolic = "OM"; // terminal O, negatively charged
						}
					}
					else
					{ 
						// O==N
						type.numeric = 7;
						type.symbolic = "O=";
						//return 7; // Nitroso oxygen (O=N)
					}
					return type;
				}
				// O-?-S
				else if ((*aiter)->is_sulfur())
				{
					if (sulphurCount == 1)
					{ 
						// O1-?-S-?-S1
						type.numeric = 32;
						//fixed by fbai ,20090228
						type.symbolic = "OSMS";
						//return 32; // Terminal oxygen in thiosulfinate anion (OSMS)
					}
					if (bond->is_single())
					{ 
						// O--S
						type.numeric = 32;
						type.symbolic = "OX";
						//return 32; // Single terminal oxygen on sulfur, One of 2 terminal O's on sulfur, 
              // One of 3 terminal O's on sulfur, Terminal O in sulfate anion, 
              // (O-S, O2S, O3S, O4S)
					}
					// fixed by xfliu, 20090227
					else if((bond->is_double() || bond->is_aromatic()) && (*aiter)->get_num_bonded_atom() >= 3)
					{
						// O=S=O
						type.numeric = 32;
						type.symbolic = "OX";
						// oxygen on sulphone, the sulphur must connecte to 4 atoms and the bond connecting oxygen is double.
					}
					//fixed by fbai,20090228
					else if((bond->is_double()))
					{ 
						// O==S
						type.numeric = 7;
						type.symbolic = "O=";
						//return 7; // Doubly bonded sulfoxide oxygen, O=S on sulfur doubly bonded 
              // to, e.g., C (O=S, O=S=)
					}
					return type;
				}
				// fixed by xfliu, 20090227
				// O-?-P
				//else if((*aiter)->is_phosphate())
				//{

				else
				{
					type.numeric = 32;
					type.symbolic = "OX";
					return type; // Oxygen in phosphine oxide, One of 2 terminal O's on sulfur, 
          // One of 3 terminal O's on sulfur, One of 4 terminal O's on sulfur, 
          // Oxygen in perchlorate anion (OP, O2P, O3P, O4P, O4Cl)
				}
			}
		}
	}
    
    ////////////////////////////////
    // Flourine
    ////////////////////////////////
	if (atom->get_element() == "F")
	{
      // 1 neighbour
		if (atom->get_num_neighbor_bond() == 1)
		{
			type.numeric = 11;
			type.symbolic = "F";
			//return 11; // Fluorine (F)
		}
      // 0 neighbours
		else if (atom->get_num_neighbor_bond() == 0)
		{
			type.numeric = 89;
			type.symbolic = "F-";
			//return 89; // Fluoride anion (F-)
		}
		return type;
    }
    
 
    ////////////////////////////////
    // Phosphorus
    ////////////////////////////////
	if (atom->is_phosphorus())
	{
		if (atom->get_num_neighbor_bond() == 4)
		{
			type.numeric = 25;
			type.symbolic = "PTET";
			//return 25; // Phosphate group phosphorus, Phosphorus with 3 attached oxygens,
        // Phosphorus with 2 attached oxygens, Phosphine oxide phosphorus,
        // General tetracoordinate phosphorus (PO4, PO3, PO2, PO, PTET)
		}
		else if (atom->get_num_neighbor_bond() == 3)
		{
			type.numeric = 26;
			type.symbolic = "P";
			//return 26; // Phosphorus in phosphines (P)
		}
		else if (atom->get_num_neighbor_bond() == 2)
		{
			type.numeric = 75;
			type.symbolic = "-P=C";
			//return 75; // Phosphorus doubly bonded to C (-P=C)
		}
		return type;
	}
    
    ////////////////////////////////
    // Sulfur
    ////////////////////////////////
	if (atom->is_sulfur())
	{
		string doubleBondTo("None");
      // 4 neighbours
      if (atom->get_num_neighbor_bond() == 4)
	  {
		  type.numeric = 18;
		  type.symbolic = "SO2X";
		  //return 18; // Sulfone sulfur, Sulfonamide sulfur, Sulfonate group sulfur,
        // Sulfate group sulfur, Sulfur in nitrogen analog of sulfone 
        // (SO2, SO2N, SO3, SO4, SNO)
	  }
      // 3 neighbours
      else if (atom->get_num_neighbor_bond() == 3)
	  {
		  oxygenCount = sulphurCount =  0;
		  for ( ATOMVec::iterator aiter=atom->get_atom_list().begin();aiter!=atom->get_atom_list().end();aiter++)
		  {
			  bond = mol_->GetBond(*aiter, atom);
			  if (bond->is_double())
			  {
				  doubleBondTo = (*aiter)->get_element();
			  }
			  if ((*aiter)->get_num_neighbor_bond() == 1)
			  {
				  if ((*aiter)->is_oxygen())
				  {
					  oxygenCount++;
				  }
				  else if ((*aiter)->is_sulfur())
				  {
					  sulphurCount++;
				  }
			  } 
		  }
		  if (oxygenCount == 2)
		  {
			  if (doubleBondTo == "C")
			  {
				  type.numeric = 18;
				  type.symbolic = "SO2";
				  //return 18; // Sulfone sulfur, doubly bonded to carbon (=SO2)
			  }
			  else
			  {
				  type.numeric = 73;
				  type.symbolic = "SO2M";
				  //return 73; // Sulfur in anionic sulfinate group (SO2M)
			  }
			  //return type;
		  }
		  else if (oxygenCount && sulphurCount)
		  {
			  type.numeric = 73;
			  type.symbolic = "SO2M";
			  //return 73; // Tricoordinate sulfur in anionic thiosulfinate group (SSOM)
		  }
  
        //if ((doubleBondTo == 6) || (doubleBondTo == 8))
		  else
		  {
			  type.numeric = 17;
			  type.symbolic = "SN";
			  //return 17; // Sulfur doubly bonded to carbon, Sulfoxide sulfur (S=C, S=O)
		  }
		  //return type;
	  }
      // 2 neighbours
      else if (atom->get_num_neighbor_bond() == 2)
	  {
		  int doubleBondTo = 0;
		  for ( ATOMVec::iterator aiter=atom->get_atom_list().begin();aiter!=atom->get_atom_list().end();aiter++)
		  {
			  if ((*aiter)->is_oxygen())
			  {
				  bond = mol_->GetBond(*aiter, atom);
				  if (bond->is_double())
				  {
					  doubleBondTo = 8;
				  }
			  }
		  }
		  if (doubleBondTo == 8)
		  {
			  type.numeric = 74;
			  type.symbolic = "=S=O";
			  //return 74; // Sulfinyl sulfur, e.g., in C=S=O (=S=O)
		  }
		  else
		  {
			  type.numeric = 15;
			  type.symbolic = "S";
			  //return 15; // Thiol, sulfide, or disulfide sulfor (S)
		  }
		  //return type;
      }
      // 1 neighbour
      else if (atom->get_num_neighbor_bond() == 1)
	  {
		  sulphurCount = 0;
		  // fixed by xfliu, 20090226
		  string doubleBondTo = "None";
		  for ( ATOMVec::iterator aiter=atom->get_atom_list().begin();aiter!=atom->get_atom_list().end();aiter++)
		  {
			  for ( ATOMVec::iterator aiter2=(*aiter)->get_atom_list().begin();aiter2!=(*aiter)->get_atom_list().end();aiter2++)
			  {
				  if ((*aiter2)->is_sulfur() && ((*aiter2)->get_num_neighbor_bond() == 1))
				  {
					  sulphurCount++;
				  }
			  }
			  bond = mol_->GetBond(*aiter, atom);
			  if (bond->is_double())
			  {
				  doubleBondTo = (*aiter)->get_element();
			  }
		  }
		  if ((doubleBondTo == "C") && (sulphurCount != 2))
		  {
			  type.numeric = 16;
			  type.symbolic = "S=C";
			  //return 16; // Sulfur doubly bonded to carbon (S=C)
		  }
		  else
		  {
			  type.numeric = 72;
			  type.symbolic = "SX";
			  //return 72; // Terminal sulfur bonded to P, Anionic terminal sulfur,
        // Terminal sulfur in thiosulfinate group (S-P, SM, SSMO)
		  }
		  //return type;
	  }
	  return type;

      // 44 ar
	}
    
    ////////////////////////////////
    // Clorine
    ////////////////////////////////
	if (atom->is_clorine())
	{
      // 4 neighbour
      if (atom->get_num_neighbor_bond() == 4)
	  {
		  oxygenCount = 0;
		  for ( ATOMVec::iterator aiter=atom->get_atom_list().begin();aiter!=atom->get_atom_list().end();aiter++)
		  {
			  if ((*aiter)->is_oxygen())
			  {
				  oxygenCount++;
			  }
		  }
		  if (oxygenCount == 4)
		  {
			  type.numeric = 77;
			  type.symbolic = "CLO4";
			  //return 77; // Perchlorate anion chlorine (CLO4)
		  }
	  }
      // 1 neighbour
      else if (atom->get_num_neighbor_bond() == 1)
	  {
		  type.numeric = 12;
		  type.symbolic = "CL";
		  //return 12; // Chlorine (CL)
      }
      // 0 neighbours
      else if (atom->get_num_neighbor_bond() == 0)
	  {
		  type.numeric = 90;
		  type.symbolic = "CL-";
		  //return 90; // Chloride anion (CL-)
      }
	  return type;
	}
    
   
    ////////////////////////////////
    // Bromine
    ////////////////////////////////
	if (atom->is_bromine())
	{
      // 1 neighbour
		if (atom->get_num_neighbor_bond() == 1)
		{
			type.numeric = 13;
			type.symbolic = "BR";
			//return 13; // Bromine (BR)
		}
      // 0 neighbours
		else if (atom->get_num_neighbor_bond() == 0)
		{
			type.numeric = 91;
			type.symbolic = "BR-";
			//return 91; // Bromide anion (BR-)
		}
		return type;
	}
 
    ////////////////////////////////
    // Iodine
    ////////////////////////////////
	if (atom->is_iodine())
	{
      // 1 neighbour
		if (atom->get_num_neighbor_bond() == 1)
		{
			type.numeric = 14;
			type.symbolic = "I";
			return type; // Iodine (I)
		}
	}

	// Metal atom and other uncommon atoms in small molecules

	////////////////////////////////
    // Silicon
    ////////////////////////////////
	if (atom->is_silicon())
	{
		type.numeric = 19;
		type.symbolic = "SI";
		return type; // Iodine (I)
	}
	////////////////////////////////
    // Lithium
    ////////////////////////////////
	if (atom->get_element() == "LI" || atom->get_element() == "Li")
	{
      // 0 neighbours
		if (atom->get_num_neighbor_bond() == 0)
		{
			type.numeric = 92;
			type.symbolic = "LI+";
			return type; // Lithium cation (LI+)
		}
    }
	////////////////////////////////
    // Sodium
    ////////////////////////////////
	if (atom->get_element() == "NA" || atom->get_element() == "Na")
	{
		type.numeric = 93;
		type.symbolic = "NA+";
		return type; // Sodium cation (NA+)
    }
    
    ////////////////////////////////
    // Magnesium
    ////////////////////////////////
	if (atom->get_element() == "MG" || atom->get_element() == "Mg")
	{
		type.numeric = 99;
		type.symbolic = "MG+2";
		return type; // Dipositive magnesium cation (MG+2)
    }
	////////////////////////////////
    // Potasium
    ////////////////////////////////
    if (atom->get_element() == "K")
	{
		type.numeric = 94;
		type.symbolic = "K+";
		return type; // Potasium cation (K+)
    }
    
    ////////////////////////////////
    // Calcium
    ////////////////////////////////
	if (atom->get_element() == "CA" || atom->get_element() == "Ca")
	{
		// 0 neighbours
		if (atom->get_num_neighbor_bond() == 0)
		{
			type.numeric = 96;
			type.symbolic = "CA+2";
			return type; // Dipositive calcium cation (CA+2)
		}
    }
 
    ////////////////////////////////
    // Iron
    ////////////////////////////////
	if (atom->get_element() == "FE" || atom->get_element() == "Fe")
	{
		type.numeric = 87;
		type.symbolic = "FE+2";
		return type;
		// it's very hard to identify the FE+2 and FE+3 only from the topological information, so the Fe is assigned the default type as FE+2;
      //return 87; // Dipositive iron (FE+2)
      //return 88; // Tripositive iron (FE+3)
    }
    
    ////////////////////////////////
    // Copper
    ////////////////////////////////
	if (atom->get_element() == "CU" || atom->get_element() == "Cu")
	{
		type.numeric = 98;
		type.symbolic = "CU+2";
		return type;
		// it's very hard to identify the CU+1 and CU+2 only from the topological information, so the Cu is assigned the default type as CU+2;
      //return 97; // Monopositive copper cation (CU+1)
      //return 98; // Dipositive copper cation (CU+2)
    }
    
    ////////////////////////////////
    // Zinc
    ////////////////////////////////
	if (atom->get_element() == "ZN" || atom->get_element() == "Zn")
	{
		type.numeric = 95;
		type.symbolic = "ZN+2";
		return type; // Dipositive zinc cation (ZN+2)
    }
    return type;
	// end of mmff94 types perceiving
}

  bool MMFF94::SetTypes()
  {
	 
    MMFF94AtomType type; 
    // It might be needed to run this function more than once...
    bool done = false;
    /*while (!done) {
      done = PerceiveAromatic();
    }*/
	//PerceiveAromatic();
    
	for(ATOMVec::iterator aiter=mol_->get_atom_vector().begin();aiter!=mol_->get_atom_vector().end();aiter++)
	{
		if((*aiter)->get_mmff94_type() != 0)
			continue;
		type = GetType((*aiter));
		if(type.numeric == 0)
		{
#ifdef DEBUG
			cout<<"Warning: MMFF94::SetTypes(): No proper atom type assigned for :"<<(*aiter)->get_symbol_type()<<" "<<(*aiter)->get_id()<<endl<<"Default types will be assigned according to element type"<<endl;
#endif
			if((*aiter)->is_carbon())
			{
				type.numeric = 1;
				type.symbolic = "C";
				// sp3 carbon
			}
			else if((*aiter)->is_hydrogen())
			{
				type.numeric = 5;
				type.symbolic = "HC";
			}
			else if((*aiter)->is_oxygen())
			{
				type.numeric = 6;
				type.symbolic = "O";
			}
			else if((*aiter)->is_nitrogen())
			{
				type.numeric = 8;
				type.symbolic = "N";
			}
			else if((*aiter)->is_phosphorus())
			{
				type.numeric = 26;
				type.symbolic = "P";
			}
			else if((*aiter)->is_sulfur())
			{
				type.numeric = 15;
				type.symbolic = "S";
			}
			else
			{
				cout<<"Warning: MMFF94::SetTypes(): Unknown element for atom "<<(*aiter)->get_id()<<" "<<(*aiter)->get_element()<<endl;
				add_unassigned_atom(*aiter);
				return false;
			}

		}
		(*aiter)->set_mmff94_type(type.numeric);
		(*aiter)->set_mmff94_symbol_type(type.symbolic);
		// fixed by xfliu, 20090227
		// just make a patch for oxygen bonded to nitrogen
		if(type.symbolic == "N2OX")
		{
			for(ATOMVec::iterator it = (*aiter)->get_atom_list().begin(); it != (*aiter)->get_atom_list().end(); ++it)
			{
				if((*it)->is_oxygen() && (*it)->get_mmff94_symbol_type() == "OM")
				{
					(*it)->set_mmff94_type(32);
					(*it)->set_mmff94_symbol_type("OX");
					break;
				}
			}
		}
		//debug
		//cout<<(*aiter)->get_id()<<" "<<(*aiter)->get_symbol_type()<<" "<<type.symbolic<<endl;
    }
    return true;
  }


  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  //
  //  Calculate bond type, angle type, stretch-bend type, torsion type
  //
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
 
  //
  // MMFF part V - page 620
  //
  // BTij is 1 when:
  // a) single bond between atoms i and j, both i and j are not aromatic and both types have sbmb set in mmffprop.par, or
  // b) bewtween two aromatic atoms, but the bond is not aromatic (e.g. connecting bond in biphenyl)
  //
  int MMFF94::GetBondType(ATOM* a, ATOM* b)
  {
	  // fixed by xfliu, 20090227
	  if (!mol_->GetBond(a,b)->is_single())
		  return 0;
	  /*else if(mol_->GetBond(a,b)->is_amide())
	  {
		  if(HasAromSet(a->get_mmff94_type()) && HasAromSet(b->get_mmff94_type()))
			  return 0;
		  else
			  return 1;
	  }*/
	  //if (!mol_->GetBond(a, b)->is_aromatic() || !mol_->GetBond(a,b)->isAromaticRingBond())
	  else if (!mol_->GetBond(a,b)->isAromaticRingBond())
	  {
		  if (HasAromSet(a->get_mmff94_type()) && HasAromSet(b->get_mmff94_type()))
			  return 1;
		  // Extremely weired!!!! modified by xfliu, 20090227
		  /*else
		  {
			  if(a->get_mmff94_type() == 3)
			  {
				  if(b->get_mmff94_type() == 63 || b->get_mmff94_type() == 64)
					  return 1;
			  }
			  else if(b->get_mmff94_type() == 3)
			  {
				  if(a->get_mmff94_type() == 63 || a->get_mmff94_type() == 64)
					  return 1;
			  }
		  }*/
	  }
	  // fixed by xfliu, 20090227
	  if (HasSbmbSet(a->get_mmff94_type()) && HasSbmbSet(b->get_mmff94_type()))
	  {
		  return 1;
	  }
    return 0;
  }
  
  int MMFF94::GetAngleType(ATOM* a, ATOM* b, ATOM *c)
  {
    int sumbondtypes;

    sumbondtypes = GetBondType(a,b) + GetBondType(b, c);

	if (mol_->IsInRingSize(a,3) && mol_->IsInRingSize(b,3) && mol_->IsInRingSize(c,3) &&In_the_sameRing(a,c))
      switch (sumbondtypes) {
      case 0:
        return 3; 
      case 1:
        return 5; 
      case 2:
        return 6; 
      }
    
    if (mol_->IsInRingSize(a,4) && mol_->IsInRingSize(b,4) && mol_->IsInRingSize(c,4) && In_the_sameRing(a,c))
      switch (sumbondtypes) {
      case 0:
        return 4; 
      case 1:
        return 7; 
      case 2:
        return 8; 
      }
    
    return sumbondtypes;
  }
  
  int MMFF94::GetStrBndType(ATOM* a, ATOM* b, ATOM *c)
  {
    int btab, btbc, atabc;
    bool inverse;

    btab = GetBondType(a, b);
    btbc = GetBondType(b, c);
    atabc = GetAngleType(a, b, c);
	// modified by FBai,according to the Ball project.
	// there is some erro in CHARMM doc.
	if (a->get_mmff94_type() <= c->get_mmff94_type())
      inverse = false;
    else
      inverse = true;

    switch (atabc) {
    case 0:
      return 0;

    case 1:
      if (btab)
        if (!inverse)
          return 1;
        else
          return 2;
      if (btbc)
        if (!inverse)
          return 2;
        else
          return 1;

    case 2:
      return 3;

    case 3:
      return 5;

    case 4:
      return 4;

    case 5:
      if (btab)
        if (!inverse)
          return 6;
        else
          return 7;
      if (btbc)
        if (!inverse)
          return 7;
        else
          return 6;
      
    case 6:
      return 8;
      
    case 7:
      if (btab)
        if (!inverse)
          return 9;
        else
          return 10;
      if (btbc)
        if (!inverse)
          return 10;
        else
          return 9;
      
    case 8:
      return 11;
    }
    return -1; //???

  }
  
  //
  // MMFF part IV - page 609
  //
  // TTijkl = 1 when BTjk = 1
  // TTijkl = 2 when BTjk = 0 but BTij and/or BTkl = 1
  // TTijkl = 4 when i, j, k and l are all members of the same four-membered ring
  // TTijkl = 5 when i, j, k and l are members of a five-membered ring and at least one is a sp3-hybridized carbon (MMFF atom type 1)
  //
  int MMFF94::GetTorsionType(ATOM* a, ATOM* b, ATOM *c, ATOM *d)
  {
    int btab, btbc, btcd;

    btab = GetBondType(a, b);
    btbc = GetBondType(b, c);
    btcd = GetBondType(c, d);
    
    if (btbc == 1)
      return 1;
    
    if (mol_->IsInRingSize(a,4) && mol_->IsInRingSize(b,4) && mol_->IsInRingSize(c,4) && mol_->IsInRingSize(d,4))
      if (In_the_sameRing(a,b) && In_the_sameRing(b,c) && In_the_sameRing(c,d))
        return 4;
   
	if (mol_->GetBond(b,c)->is_single()) {
      if (btab || btcd)
        return 2;
      /*
        unsigned int order1 = GetCXT(0, atoi(d->GetType()), atoi(c->GetType()), atoi(b->GetType()), atoi(a->GetType())); 
        unsigned int order2 = GetCXT(0, atoi(a->GetType()), atoi(b->GetType()), atoi(c->GetType()), atoi(d->GetType()));
    
        cout << "GetTorsionType(" << a->GetType() << ", " << b->GetType() << ", " << c->GetType() << ", " << d->GetType() << ")" << endl;
        cout << "    order1 = " << order1 << endl;
        cout << "    order2 = " << order2 << endl;
        cout << "    btab = " << btab << endl;
        cout << "    btbc = " << btbc << endl;
        cout << "    btcd = " << btcd << endl;
      */
    }
    
    if (mol_->IsInRingSize(a,5) && mol_->IsInRingSize(b,5) && mol_->IsInRingSize(c,5) && mol_->IsInRingSize(d,5))
		if( a->get_mmff94_type() == 1 || b->get_mmff94_type() == 1 || c->get_mmff94_type() == 1 || d->get_mmff94_type() == 1)
			return 5;

      /*vector<RING>::iterator ri;
      vector<int>::iterator rj;
	  for (ri = mol_->get_ring_vector().begin();ri!=mol_->get_ring_vector().end();ri++) { // for each ring
		  if ((*ri).is_aromatic)
          continue;
	
        if ((*ri).size != 5)
          continue;
		
		if(!(mol_->IsInRingSize(a,5)))
        
        if (!(*ri).IsMember(a) || !(*ri).IsMember(b) || !(*ri).IsMember(c) || !(*ri).IsMember(d))
          continue;
	
        return 5;
      }
    }*/
	return 0;
  }

  // CXB = MC * (I * MA + J) + BTij
  unsigned int MMFF94::GetCXB(int type, int a, int b)
  {
    unsigned int cxb;
    cxb = 2 * (a * 136 + b) + type;
    return cxb;
  }
  
  // CXA = MC * (J * MA^2 + I * MA + K) + ATijk
  unsigned int MMFF94::GetCXA(int type, int a, int b, int c)
  {
    unsigned int cxa;
    cxa = 9 * (b * 18496 + a * 136 + c) + type;
    return cxa;
  }
  
  // CXS = MC * (J * MA^2 + I * MA + K) + STijk
  unsigned int MMFF94::GetCXS(int type, int a, int b, int c)
  {
    unsigned int cxs;
    cxs = 12 * (b * 18496 + a * 136 + c) + type;
    return cxs;
  }
  
  // CXO = J * MA^3 + I * MA^2 + K * MA + L
  unsigned int MMFF94::GetCXO(int a, int b, int c, int d)
  {
    unsigned int cxo;
    cxo = b * 2515456 + a * 18496 + c * 136 + d;
    return cxo;
  }
  
  // CXT = MC * (J * MA^3 + K * MA^2 + I * MA + L) + TTijkl
  unsigned int MMFF94::GetCXT(int type, int a, int b, int c, int d)
  {
    unsigned int cxt;
    cxt = 6 * (b * 2515456 + c * 18496 + a * 136 + d) + type;
    return cxt;
  }
  
  // CXQ = MC * (I * MA + J) + BTij
  unsigned int MMFF94::GetCXQ(int type, int a, int b)
  {
    unsigned int cxq;
    cxq = 2 * (a * 136 + b) + type;
    return cxq;
  }

    ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  //
  //  Various tables & misc. functions
  //
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
 
  // MMFF part V - TABLE I
  bool MMFF94::HasLinSet(int atomtype)
  {
	  bool defined=false;
	  if(MMFF94_LIN[atomtype-1])defined=true;
	  return defined;  
  }
 
  // MMFF part V - TABLE I
  bool MMFF94::HasPilpSet(int atomtype)
  {
	  bool defined=false;
	  if( MMFF94_PILP[atomtype-1])defined=true;
	  return defined;  
  }
  
  // MMFF part V - TABLE I
  bool MMFF94::HasAromSet(int atomtype)
  {
      bool defined=false;
	  if( MMFF94_AROM[atomtype-1])defined=true;
	  return defined;  
  }
 
  // MMFF part V - TABLE I
  bool MMFF94::HasSbmbSet(int atomtype)
  {
      bool defined = false;
	  if( MMFF94_SBMB[atomtype-1])
	  {
		  defined = true;
	  }
	  return defined;  
  }

  // MMFF part V - TABLE I
  int MMFF94::GetCrd(int atomtype)
  {
	 int crd=0;
	 crd= MMFF94_CRD[atomtype-1];
	 return crd;
  }

  // MMFF part V - TABLE I
  int MMFF94::GetVal(int atomtype)
  {
    int crd=0;
	 crd= MMFF94_VAL[atomtype-1];
	 return crd;
  }

  // MMFF part V - TABLE I
  int MMFF94::GetMltb(int atomtype)
  {
    int crd=0;
	 crd= MMFF94_MLTB[atomtype-1];
	 return crd;
  }
 
  // MMFF part V - TABLE VI
  double MMFF94::GetZParam(ATOM* atom)
  {
	  if (atom->is_hydrogen())
      return 1.395;
	  if (atom->is_carbon())
      return 2.494;
    if (atom->is_nitrogen())
      return 2.711;
    if (atom->is_oxygen())
      return 3.045;
    if (atom->get_atomic_num() == 9) // F
      return 2.847;
    if (atom->get_atomic_num() == 14) // Si
      return 2.350;
    if (atom->is_phosphorus())
      return 2.350;
    if (atom->is_sulfur())
      return 2.980;
    if (atom->get_atomic_num() == 17) // Cl
      return 2.909;
    if (atom->get_atomic_num() == 35) // Br
      return 3.017;
    if (atom->get_atomic_num() == 53) // I
      return 3.086;

    return 0.0;
  }
 
  // MMFF part V - TABLE VI
  double MMFF94::GetCParam(ATOM* atom)
  {
    if (atom->get_atomic_num() == 5) // B
      return 0.704;
	if (atom->is_carbon())
      return 1.016;
	if (atom->is_nitrogen())
      return 1.113;
	if (atom->is_oxygen())
      return 1.337;
   /* if (atom->GetAtomicNum() == 14) // Si
      return 0.811;*/
	if (atom->is_phosphorus())
      return 1.068;
	if (atom->is_sulfur())
      return 1.249;
	if (atom->get_atomic_num() == 17) // Cl
      return 1.078;
    /*if (atom->GetAtomicNum() == 33) // As
      return 0.825;*/

    return 0.0;
  }
 
  // MMFF part V - TABLE X
  double MMFF94::GetUParam(ATOM* atom)
  {
	  if (atom->is_carbon())
      return 2.0;
    if (atom->is_nitrogen())
      return 2.0;
	if (atom->is_oxygen())
      return 2.0;
    /*if (atom->GetAtomicNum() == 14) // Si
      return 1.25;*/
	if (atom->is_phosphorus())
      return 1.25;
	if (atom->is_sulfur())
      return 1.25;
    
    return 0.0;
  }
  
  // MMFF part V - TABLE X
  double MMFF94::GetVParam(ATOM* atom)
  {
	  if (atom->is_carbon())
      return 2.12;
	  if (atom->is_nitrogen())
      return 1.5;
	  if (atom->is_oxygen())
      return 0.2;
   /* if (atom->GetAtomicNum() == 14) // Si
      return 1.22;*/
    if (atom->is_phosphorus())
      return 2.4;
    if (atom->is_phosphorus())
      return 0.49;
    
    return 0.0;
  }
    
  // R Blom and A Haaland, J. Mol. Struct., 128, 21-27 (1985)
 double MMFF94::GetCovalentRadius(ATOM* a)
 {

    switch (a->get_atomic_num())
	{
    case 1:
      return 0.33; // corrected value from MMFF part V
    case 5:
      return 0.81;
    case 6:
      return 0.77; // corrected value from MMFF part V
    case 7:
      return 0.73;
    case 8:
      return 0.72;
    case 9:
      return 0.74;
    case 13:
      return 1.22;
    case 14:
      return 1.15;
    case 15:
      return 1.09;
    case 16:
      return 1.03;
    case 17:
      return 1.01;
    case 31:
      return 1.19;
    case 32:
      return 1.20;
    case 33:
      return 1.20;
    case 34:
      return 1.16;
    case 35:
      return 1.15;
    case 44:
      return 1.46;
    case 50:
      return 1.40;
    case 51:
      return 1.41;
    case 52:
      return 1.35;
    case 53:
      return 1.33;
    case 81:
      return 1.51;
    case 82:
      return 1.53;
    case 83:
      return 1.55;
    default:
      return 0.0;
    }
  }
  
 
  // MMFF part V - page 625
  double MMFF94::GetRuleBondLength(ATOM* a, ATOM* b)
  {
    double r0ab, r0a, r0b, c, Xa, Xb;
    int Ha, Hb, BOab;
	int type1,type2;
    r0a = GetCovalentRadius(a);
    r0b = GetCovalentRadius(b);
	Xa = GetAllredRochowElectroNeg(a);
	Xb = GetAllredRochowElectroNeg(b);
	type1=GetType(a).numeric;
	type2=GetType(b).numeric;
    
   
	if (a->is_hydrogen())
      r0a = 0.33;
	if (b->is_hydrogen())
      r0b = 0.33;
    
	if (a->is_hydrogen() || b->is_hydrogen())
      c = 0.050;
    else
      c = 0.085;

    if (GetMltb(type1) == 3)
      Ha = 1;
    else if ((GetMltb(type1) == 1) || (GetMltb(type2) == 2))
      Ha = 2;
    else
      Ha = 3;

    if (GetMltb(type2) == 3)
      Hb = 1;
    else if ((GetMltb(type2) == 1) || (GetMltb(type2) == 2))
      Hb = 2;
    else
      Hb = 3;

	BOab = mol_->GetBond(a,b)->get_bond_order();
    if (((GetMltb(type1) == 1) && (GetMltb(type2)) == 1))
      BOab = 4;
    if (((GetMltb(type1) == 1) && (GetMltb(type2)) == 2))
      BOab = 5;
    if (((GetMltb(type1) == 2) && (GetMltb(type2)) == 1))
      BOab = 5;
	if (mol_->GetBond(a,b)->is_aromatic())
      if (!HasPilpSet(type1) && !HasPilpSet(type2))
        BOab = 4;
      else
        BOab = 5;
     
    switch (BOab)
	{
    case 5:
      r0a -= 0.04;
      r0b -= 0.04;
      break;
    case 4:
      r0a -= 0.075;
      r0b -= 0.075;
      break;
    case 3:
      r0a -= 0.17;
      r0b -= 0.17;
      break;
    case 2:
      r0a -= 0.10;
      r0b -= 0.10;
      break;
    case 1:
      if (Ha == 1)
        r0a -= 0.08;
      if (Ha == 2)
        r0a -= 0.03;
      if (Hb == 1)
        r0b -= 0.08;
      if (Hb == 2)
        r0b -= 0.03;
    }
    
    /*
      cout << "Ha=" << Ha << "  Hb=" << Hb << "  BOab=" << BOab << endl;
      cout << "r0a=" << r0a << "  Xa=" << Xa << endl;
      cout << "r0b=" << r0b << "  Xb=" << Xb << endl;
      cout << "r0a + r0b=" << r0a +r0b << endl;
      cout << "c=" << c << "  |Xa-Xb|=" << fabs(Xa-Xb) << "  |Xa-Xb|^1.4=" << pow(fabs(Xa-Xb), 1.4) << endl;
    */
    r0ab = r0a + r0b - c * pow(fabs(Xa - Xb), 1.4) - 0.008; 

    return r0ab;
  }
  double MMFF94::GetAllredRochowElectroNeg(ATOM*a)
  {
	 if(a->get_element()=="H")
		 return 2.20;
	 else if(a->get_element()=="C")
		 return 1.56;
	 else if(a->get_element()=="O")
		 return 3.50;
	 else if(a->get_element()=="N")
		 return 3.07;
	 else if (a->get_element()=="P")
		 return 2.06;
	 else if(a->get_element()=="S")
		 return 2.44;
	 else if(a->get_element()=="F")
		 return 4.10;
	 else if(a->get_element()=="Cl")
		 return 2.83;
	 else if(a->get_element()=="Br")
		 return 2.74;
	 else if(a->get_element()=="I")
		 return 2.21;
	 else return 0.0;

  }
                           
int MMFF94::get_update_frequency() const
{
	return DEFAULT_UPDATEFREQUENCY;
}

double  MMFF94::get_stretch_energy() const
{
	FFComponent* component = get_component("MMFF94 Stretch");
	if (component != 0)
	{
		//debug
		//cout<<"stretch flag"<<endl;
		return component->get_energy();
	}
	else
		return 0.0;
}

double  MMFF94::get_bend_energy() const
{
	FFComponent* component = get_component("MMFF94 Bend");
	
	if (component != 0)
	{
		//debug
      //cout<<" Bend flag"<<endl;
		return component->get_energy();
	}
	else
		return 0.0;
}

double  MMFF94::get_torsion_energy() const
{
	FFComponent* component = get_component("MMFF94 Torsion");
	if (component != 0)
	{
		//debug
		//cout<<"torsion flag"<<endl;
		return component->get_energy();
	}
	else
		return 0.0;
}

double MMFF94::get_oop_energy() const
{
	FFComponent* component = get_component("MMFF94 OOP");
	if (component != 0)
	{
		//debug
		//cout<<"oop flag"<<endl;
		return component->get_energy();
	}
	else
		return 0.0;
}

double MMFF94:: get_stretch_bend_energy() const
{
	FFComponent* component = get_component("MMFF94 Str_Bend");
	if(component !=0)
	{
	//debug
	//cout<<" Str_Bend flag"<<endl;
		return component->get_energy();
	}
	return 0.0;
}
double  MMFF94::get_vdw_energy() const
{
	FFComponent* component = get_component("MMFF94 VDW");
	if (component != 0)
	{
		//debug
		//cout<<"vdw flag"<<endl;
		return component->get_energy();
	}
	return 0.0;
}

double  MMFF94::get_ele_energy() const
{
	FFComponent* component = get_component("MMFF94 Ele");
	if (component != 0)
	{
		//debug
		//cout<<"ele flag"<<endl;
		return component->get_energy();
	}
	return 0.0;
}

bool  MMFF94::has_initialized_param() const
{
	return param_is_initialized_;
}


string  MMFF94::get_results() const
{
	ostringstream os;

	os<<"\n"
		<<" Details of MMFF94 Energy:\n"
		<<" - Electrostatic       :     " <<get_ele_energy()<<endl 
		<<" - Van der Waals       :     " <<get_vdw_energy()<<endl
		<<" - Bond Stretching     :     " <<get_stretch_energy()<<endl
		<<" - Angle Bending       :     " <<get_bend_energy()<<endl 
		<<" - Torsion             :     " <<get_torsion_energy()<<endl
		<<" - OOP Bend            :     " <<get_oop_energy()<<endl 
		<<" - Str-Bend Interaction:     "<<get_stretch_bend_energy()<<endl
		<<"-------------------------------------------------\n" 
		<<"  Total Energy         :"<<get_energy()<<" Kcal/mol"<<endl;

	return os.str();
}

 //MMFF94 mmff94;                 