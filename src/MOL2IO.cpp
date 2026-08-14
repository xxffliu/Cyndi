#include "../include/MOL2IO.h"
#include <iomanip>
#include <sstream>
#include <ctime>

//using namespace std;
// default constructor, using individual default constructor of each member

MOL2IO::MOL2IO():
atoms_(),
bonds_(),
substructures_(),
sets_(),
has_sets_info_(false),
has_substr_info_(false),
first_mol_flag_(true),
mark_(fstream::beg)
{
	input_.close();
	input_.clear();
	output_.close();
	output_.clear();
}
// constructor
MOL2IO::MOL2IO(const string& filename, const string& openmode):
atoms_(),bonds_(),
substructures_(),
sets_(),
num_of_lines_(0),
has_sets_info_(false),
has_substr_info_(false),
first_mol_flag_(true),
mark_(0)
{

	bool SuccessOpen = open(filename, openmode);
	if(!SuccessOpen)
	{
		cout<<"MOL2IO::MOL2IO(): Failed to Open file "<<filename<<" Please check if the file exists"<<endl;
		exit(1);
	}
}

bool MOL2IO::open(const string& filename, const string& openmode)
{
	clear();
	stream_clear_();
	if (openmode == "in")
	{
		// establish the atomic radius and weight map
		AtomRadius_["C"] = 1.70, AtomRadius_["N"] = 1.55, AtomRadius_["O"] = 1.52,
		AtomRadius_["P"] = 1.80, AtomRadius_["S"] = 1.80, AtomRadius_["Cl"] = 1.75,
		AtomRadius_["Br"] = 1.85, AtomRadius_["I"] = 1.98, AtomRadius_["F"] = 1.47,
		AtomRadius_["H"] = 1.20, AtomRadius_["X"] = 0.0;

		AtomWeight_["C"] = 12.01, AtomWeight_["O"] = 16.00, AtomWeight_["P"] = 30.97,
		AtomWeight_["S"] = 32.07, AtomWeight_["Cl"] = 35.45, AtomWeight_["Br"] = 79.90,
		AtomWeight_["I"] = 126.9, AtomWeight_["F"] = 19.0, AtomWeight_["H"] = 1.008,
		AtomWeight_["N"] = 14.01, AtomWeight_["X"] = 0.0;
		input_.open(filename.c_str(), ifstream::binary);
		if(input_)
			return true;
		else
			return false;
	}
	else if(openmode == "out")
	{
		output_.open(filename.c_str());
		if(output_)
			return true;
		else
			return false;
	}
	else if(openmode == "app")
	{
		output_.open(filename.c_str(), ofstream::app);
		if(output_)
			return true;
		else
			return false;
	}
	else
	{
		cout<<"MOL2IO::MOL2IO(): Error: No appropriate open mode with the name of "<<openmode<<endl;
		exit(1);
	}
}

void MOL2IO::clear()
{

	// clear the structure for the molecule section
	molecule_.name = "";
	molecule_.num_of_atoms = 0;
	molecule_.num_of_bonds = 0;
	molecule_.num_of_substructures = 0;
	molecule_.num_of_features = 0;
	molecule_.num_of_sets = 0;

	// clear the vectors for the other sections
	atoms_.clear();
	bonds_.clear();
	substructures_.clear();
	sets_.clear();
	num_of_lines_ = 0;
	has_sets_info_ = false;
	has_substr_info_ = false;
}
void MOL2IO::stream_clear_()
{
	mark_ = fstream::beg;
	first_mol_flag_ = true;
	input_.close();
	input_.clear();
	output_.close();
	output_.clear();
}
MOL2IO::~MOL2IO()
{
	clear();
	stream_clear_();
}

bool MOL2IO::readConfs(deque<MOL>& conf)
{
	conf.clear();
	MOL tmpmol;
	bool first_flag = true;
	string previous_name_prefix, current_name_prefix;
	ifstream::pos_type tmp_pos;
	while(read(&tmpmol))
	{
		if(first_flag)
		{
			previous_name_prefix = tmpmol.get_name();
			// take each molecules are independent and return the current molecule
			conf.push_back(tmpmol);
			first_flag = false;
			tmp_pos = mark_;
		}
		else
		{
			current_name_prefix = tmpmol.get_name();
			if(current_name_prefix == previous_name_prefix)
			{
				conf.push_back(tmpmol);
				tmp_pos = mark_;
			}
			else
			{
				mark_ = tmp_pos;
				break;
			}
		}
	}
	if(!conf.empty())
		return true;
	else
		return false;
}
	
//void MOL2IO::bond_to_mol(MOL::MOL& mol){;}
bool MOL2IO::read(MOL* mol)
{
	mol->clear();
	clear();
	bool isFileEnd = false;
	input_.seekg(mark_);
	if(input_.eof())
		return false;
	if(input_)
	{
		first_mol_flag_ = true;
		string buffer;
		ifstream::pos_type tmp_pos = mark_;
		while(true)
		{
			if(input_.eof())
			{
				break;
			}
			else
			{
				mark_ = input_.tellg();
				getline(input_, buffer);
				if(buffer.find("@<TRIPOS>MOLECULE") != string::npos)
					if(first_mol_flag_)
					{
						first_mol_flag_ = false;
						continue;
					}
					else
					{
						break;
					}
			}
		}
		input_.clear();
		string line ;
		input_.seekg(tmp_pos);
		while(input_.tellg() != mark_ && !input_.eof())
		{
			getline(input_, line);
			num_of_lines_ ++;
			if(line.find("#") == 0)
			{
				continue;
			}
			else if(line == "")
				continue;
			else if(line.find("@<TRIPOS>ATOM") == 0)
			{
				read_atom_section_();
			}
			else if(line.find("@<TRIPOS>BOND") == 0)
			{
				read_bond_section_();
			}
			else if(line.find("@<TRIPOS>MOLECULE") == 0)
			{
				read_molecule_section_();

			}
			/*else if(line.find("@<TRIPOS>SET") == 0){
			has_sets_info_ = true;
			read_set_section_();

			}*/
			else if(line.find("@<TRIPOS>SUBSTRUCTURE") == 0)
			{
				has_substr_info_ = true;
				read_substructure_section_();

			}
			else if(line.find("@<TRIPOS>COMMENT") == 0)
				continue;
			else
			{
				//cout<<"Warning: MOL2IO::read(): section ignored in line"<<num_of_lines_<<endl;
				continue;
			}                                         
		}
		return build_all_(mol);
	}
	else
	{
		return false;
	}
}

void MOL2IO::read_atom_section_(){
	string buffer;
	while((input_.peek() != '@') && getline(input_, buffer)){
		num_of_lines_++;
		if(buffer.find("#") == 0)
			continue;
		vector<string> fields;
		istringstream sin(buffer);
		string token;
		while(sin>>token)
			fields.push_back(token);
		if(fields.size() < 6)
			continue;
		AtomInfo atom;
		atom.id = str2int(fields[0]);
		atom.name = fields[1];
		atom.position.SetX(str2double(fields[2]));
		atom.position.SetY(str2double(fields[3]));
		atom.position.SetZ(str2double(fields[4]));
		// here we patch the nomination problems for S.O2 and S.O
		if(fields[5] == "S.O2")
			fields[5] = "S.o2";
		else if(fields[5] == "S.O")
			fields[5] = "S.o";
		atom.symbolic_type = fields[5];
		// fixed by xfliu on 20090909
		// not all of the mol2 files have the columns of substructures information and charge
		if(fields.size() >= 7)
			atom.substructure_id = str2int(fields[6]);
		if(fields.size() >= 8)
			atom.substructure_name = str2int(fields[7]);
		if(fields.size() >= 9)
			atom.charge = str2double(fields[8]);
		atoms_.push_back(atom);
	}
}
void MOL2IO::read_bond_section_(){
	string buffer;
	while((input_.peek() != '@') && getline(input_, buffer)){
		num_of_lines_++;
		if(buffer.find("#") == 0)
			continue;
		vector<string> fields;

		istringstream sin(buffer);
		string token;
		while(sin>>token){
			fields.push_back(token);
		}
		if(fields.size() < 4)
			continue;
		BondInfo bond;
		bond.id = str2int(fields[0]);
		bond.atom1 = str2int(fields[1]);
		bond.atom2 = str2int(fields[2]);
		bond.type = fields[3];
		bonds_.push_back(bond);
	}
}
void MOL2IO::read_molecule_section_(){
	string buffer;
	int line_num = 0;
	while((input_.peek() != '@') && getline(input_,buffer)){
		//while(getline(input_,buffer)){
		//cout << buffer << endl;
		//if (buffer.find("@") != string::npos)
		//break;
		num_of_lines_++;
		if(buffer.find("#") == 0)
			continue;
		line_num ++;
		vector<string> fields;
		istringstream sin(buffer);
		string token;
		while(sin>>token)
			fields.push_back(token);
		switch(line_num){
			// we read five lines describing molecular information
							   case 1:
								   if(fields.size() == 1)
										molecule_.name = fields[0];
								   else
								   {
									   string temp;
									   for(vector<string>::iterator it = fields.begin(); it != fields.end(); ++it)
										   temp += *it;
									   molecule_.name = temp;
								   }

								   if(molecule_.name == "****")
									   molecule_.name = "MOLECULE";
								   break;
							   case 2:
								   molecule_.num_of_atoms = str2int(fields[0]);
								   if (fields.size() > 1)
									   molecule_.num_of_bonds = str2int(fields[1]);
								   if (fields.size() > 2)
									   molecule_.num_of_substructures = str2int(fields[2]);
								   if (fields.size() > 3)
									   molecule_.num_of_features = str2int(fields[3]);
								   if (fields.size() > 4)
									   molecule_.num_of_sets = str2int(fields[4]);

								   break;
							   case 3:
								   molecule_.type = fields[0];
								   break;
							   case 4:
								   molecule_.charge_type = fields[0];
								   break;
							   case 5:

								   if(fields.size() == 0)
									   break;
								   else{
									   molecule_.annotate = fields[0];
									   break;
								   }
		}
	}
}
// we ignore set, and substructure section at beta version
void MOL2IO::read_set_section_(){
	string buffer;
	int line_num = 0;
	while((input_.peek() != '@') && getline(input_,buffer)){
		num_of_lines_++;
		if(buffer.find("#") == 0)
			continue;
		line_num ++;
		vector<string> fields;
		istringstream sin(buffer);
		string token;
		while(sin>>token)
			fields.push_back(token);
		SetInfo set;
		switch(line_num){
			// we read 2 lines describing set information
			// currently we ignore "DYNAMIC" type and only "STATIC" type is read
							   case 1:
								   set.name = fields[0];
								   set.type = fields[1];
								   set.sub_type = fields[2];
								   set.annotate = fields[3];
								   break;
							   case 2:
								   if(set.type != "STATIC")
									   continue;
								   for(int i = 0; i!=fields.size();i++)
									   set.members[i] = str2int(fields[i]);
								   break;
		}
		sets_.push_back(set);
	}
}

void MOL2IO::read_substructure_section_(){/*????*/}

bool MOL2IO::build_all_(MOL* mol){
	bool read_anything = false;
	// construct atoms
	vector<ATOM*> atom_ptr(atoms_.size());
	for (vector<AtomInfo>::size_type i = 0; i<atoms_.size(); i++){
		read_anything = true;
		ATOM* atom = new ATOM;
		atom->set_id(atoms_[i].id);
		atom->set_name(atoms_[i].name);
		atom->set_symbol_type(atoms_[i].symbolic_type);
		if (atoms_[i].symbolic_type != "*"){
			string type = atoms_[i].symbolic_type;
			if(type.size()>1 && type[1] == '.')
				atom->set_element(type.substr(0,1));
			else
				atom->set_element(atoms_[i].symbolic_type);
		}
		else
			atom->set_element("X");
		string element = atom->get_element();
		// set the atomic weight and vdw radius
		atom->radius_ = AtomRadius_[element];
		atom->weight_ = AtomWeight_[element];
		mol->MW_ += AtomWeight_[element];
		//
		atom->set_position(atoms_[i].position);
		atom->set_charge(atoms_[i].charge);
		atom_ptr[i] = atom;
	}
	//cout<<"atom constructed"<<endl; //debug
	// construct bonds
	vector<BOND*> bond_ptr(bonds_.size());
	for (vector<BondInfo>::size_type i = 0; i<bonds_.size(); i++){
		if ((bonds_[i].atom1 > atom_ptr.size()) || (bonds_[i].atom2 > atom_ptr.size()))
			cout << "Error: MOL2File::read(): cannot build bond between atoms " 
			<< bonds_[i].atom1 << " and " << bonds_[i].atom2 << endl;
		BOND* bond = new BOND;
		bond->set_id(bonds_[i].id);
		bond->set_type(bonds_[i].type);
		bond->set_first_atom(atom_ptr[bonds_[i].atom1-1]);
		bond->set_second_atom(atom_ptr[bonds_[i].atom2-1]);
		float bond_order = 0;
		if(bonds_[i].type == "1" || bonds_[i].type == "am")
			bond_order = 1;
		else if(bonds_[i].type == "2")
			bond_order = 2;
		else if(bonds_[i].type == "3")
			bond_order = 3;
		else if(bonds_[i].type == "ar")
			bond_order = 1.5;
		bond->set_bond_order(bond_order);
		bond_ptr[i] = bond;
	}
	//cout<<"bond constructed"<<endl;
	// constuct fragment
	vector<FRAGMENT*> frag_ptr;
	if(has_sets_info_){
		frag_ptr.resize(sets_.size());
		for(vector<SetInfo>::size_type i = 0; i < sets_.size(); i++){
			FRAGMENT* frag = new FRAGMENT;
			if(sets_[i].sub_type == "ATOMS"){
				vector<ATOM*> atoms;
				for(vector<int>::size_type j = 0; j < sets_[i].members.size(); j++)
					atoms.push_back(atom_ptr[sets_[i].members[j]-1]);
				frag->set_member_atoms(atoms);
				frag->set_num_mem_atoms(int(atoms.size()));
			}
			else if(sets_[i].sub_type == "BONDS"){
				vector<BOND*> bonds;
				for(vector<int>::size_type j = 0; j< sets_[i].members.size(); j++)
					bonds.push_back(bond_ptr[sets_[i].members[j]-1]);
				frag->set_member_bonds(bonds);
				frag->set_num_mem_bonds(int(bonds.size()));
			}
			// set the id beginning with 1
			frag->set_id(int(i+1));
			frag->set_name(sets_[i].annotate);
			frag_ptr[i] = frag;
		}
	}
	// construct molecule
	//cout<<"Reading molecule "<<molecule_.name<<endl;
	mol->set_atom_vector(atom_ptr);
	mol->set_bond_vector(bond_ptr);
	if(has_sets_info_)
		mol->set_fragment_vector(frag_ptr);
	mol->set_name(molecule_.name);
	// fixed by xfliu, 20090316
	if(molecule_.num_of_atoms > atom_ptr.size())
		molecule_.num_of_atoms = atom_ptr.size();
	else if(molecule_.num_of_bonds > bond_ptr.size())
		molecule_.num_of_bonds = bond_ptr.size();
	mol->set_num_atom(molecule_.num_of_atoms);
	mol->set_num_bond(molecule_.num_of_bonds);
	mol->set_num_fragment(molecule_.num_of_sets);
	//mol->initialize();
	return read_anything;
}

bool MOL2IO::write(const MOL& mol0){
	MOL& mol = const_cast<MOL&>(mol0);
	if(output_){
		// output the file name, pharmacophore model name and creation time in comment part
		time_t it = time(NULL);
		string time(ctime(&it));		 
		output_<<"# Name:		"<<mol.get_name()<<endl;
		output_<<"# Pharmacophore name:		"<<mol.get_pharmacophore_name()<<endl;
		output_<<"# Created and modified by GAPS on "<<time<<endl;
		// if any, write the comment section, here including fitness value, energy value, etc.
		output_<<"# "<<mol.get_comment()<<endl;
		// output the fit value against the pharmacophore model and energy value;
		output_<<"# Fit value:		"<<mol.get_fitness()<<endl<<"# Energy:		"<<mol.get_energy()<<"kJ/mol"<<endl;
		output_<<endl;
		///////////////////////////////
		// write the molecule header
		output_<<"@<TRIPOS>MOLECULE"<<endl;
		// write the molecule name
		string name = mol.get_name();
		if(name == "")
			name = "****";
		output_<<name<<endl;
		// write the num of atoms, bonds and substructures
		output_<<setw(5)<<mol.get_num_atom()<<" "<<
			setw(5)<<mol.get_num_bond()<<setw(5)<<" "<<
			setw(5)<<mol.get_num_fragment()<<endl;
		// write molecule and charge type
		string charge_type = mol.get_charge_type();
		if(charge_type != "")
			output_<<"SMALL"<<endl<<charge_type<<endl;
		else
			output_<<"SMALL"<<endl<<"USER_CHARGE"<<endl;
		// write a blank line as required by mol2 format
		output_<<endl;
		// done with the molecule header;
		// now write the atom and bond section
		// atom section header first
		if(mol.get_num_atom() == mol.get_atom_vector().size()){
			output_<<"@<TRIPOS>ATOM"<<endl;
			//Format: atom_id atom_name x y z atom_type subst_id subst_name charge
			for (ATOMVec::iterator aiter = mol.get_atom_vector().begin();
				aiter != mol.get_atom_vector().end(); ++aiter){
					int id = (*aiter)->get_id();
					string name = (*aiter)->get_name();
					if(name == "")
						name = "****";
					double x = (*aiter)->get_position().x();
					double y = (*aiter)->get_position().y();
					double z = (*aiter)->get_position().z();
					string type = (*aiter)->get_symbol_type();
					if(type == "")
						type = "****";
					// method for writing substructure information will be append later
					string sub_id = "1 ";
					string sub_name = "**** ";
					double charge = (*aiter)->get_charge();
					// write all the information as the formation
					output_.setf(ios_base::right,ios_base::adjustfield);
					output_<<setw(4)<<id<<" ";
					output_.setf(ios_base::left,ios_base::adjustfield);
					output_<<setw(6)<<name<<" ";
					output_.setf(ios_base::right,ios_base::adjustfield);
					output_<<setiosflags(ios_base::fixed)<<setprecision(4)
						<<setw(8)<<x<<" "<<setw(8)<<y<<" "<<setw(8)<<z<<" ";
					output_.setf(ios_base::left,ios_base::adjustfield);
					output_<<setw(6)<<type<<" ";
					output_.setf(ios_base::right,ios_base::adjustfield);
					output_<<setw(3)<<sub_id;
					output_.setf(ios_base::left,ios_base::adjustfield);
					output_<<setw(8)<<sub_name;
					output_<<setiosflags(ios_base::showpoint)<<setw(9)<<charge<<endl;
			}
		}
		// done with the atom section

		// write the bond section
		if(mol.get_num_bond() == mol.get_bond_vector().size()){
			output_<<"@<TRIPOS>BOND"<<endl;
			for(BONDVec::iterator biter = mol.get_bond_vector().begin();
				biter != mol.get_bond_vector().end(); ++biter){
					int id = (*biter)->get_id();
					int atom1 = (*biter)->get_first_atom()->get_id();
					int atom2 = (*biter)->get_second_atom()->get_id();
					string type = (*biter)->get_type();
					//write the formated bond information 
					output_.setf(ios_base::right,ios_base::adjustfield);
					output_<<setw(4)<<id<<" "<<setw(4)<<atom1<<" "<<setw(4)<<atom2<<" ";
					output_.setf(ios_base::left,ios_base::adjustfield);
					output_<<setw(4)<<type<<endl;
			}
		}
		// done with the bond section
		// ignore the substructure information temperally
		output_<<"@<TRIPOS>SUBSTRUCTURE"<<endl<<"1 **** 2 TEMP 0 **** **** 0 ROOT"<<endl;
		// done with writing
		return true;
	}
	else
		return false;
}





