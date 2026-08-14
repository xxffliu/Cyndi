#include "../include/MMFF94VDW.h"
#include "../include/MMFF94.h"
#include "../include/utility.h"

using namespace std;

// constructors
MMFF94VDW::MMFF94VDW():FFComponent(),fast_access(),vdw_data_holder_(),vdw_energy_(0.0)
{
	//set component name
    set_name("MMFF94 VDW");
}
MMFF94VDW::MMFF94VDW(ForceField& ff):FFComponent(ff),fast_access(),vdw_data_holder_(),vdw_energy_(0.0)
{
	//set component name
    set_name("MMFF94 VDW");
}
//copy constructor
MMFF94VDW::MMFF94VDW(const MMFF94VDW& to_copy):FFComponent(to_copy)
{
	fast_access = to_copy.fast_access;
    vdw_data_holder_ = to_copy.vdw_data_holder_;
    vdw_energy_ = to_copy.vdw_energy_;
}
// destructor
MMFF94VDW::~MMFF94VDW()
{
	vdw_data_holder_.clear();
}
// extract vdw parameters from the non bonded FFParamter object and build a hashtable for fast access
bool MMFF94VDW::extract_VDW_parameters(FFParameter& ffp)
{
	if(!ffp.is_valid())
         return false;
    // 
    // build a 2 dim array of atom types and loop variables
    //FFParameter::AtomTypes atom_types = ffp.get_atomtypes();
    int num_types = ffp.get_num_types();
    int num_entries = num_types * num_types;
    num_of_atom_types_ = num_types;
    // clear internal data structures
	fast_access.is_defined.clear();
	fast_access.is_defined.resize(num_entries);
	fast_access.R_ij.clear();
	fast_access.R_ij.resize(num_entries);
	fast_access.R_ij_7.clear();
	fast_access.R_ij_7.resize(num_entries);
	fast_access.Epsilon_ij.clear();
	fast_access.Epsilon_ij.resize(num_entries);
	int index_1=0;


    for (int i = 0; i < num_types; i++)
	{
		for(int j=0;j<num_types;j++)
		{
			index_1= i +j * num_types;
			fast_access.is_defined[index_1] = false;
		}
	}
     // start pack the parameters into the vector fast_access;
    int type;
    string name_type;
    vector<double> alpha(num_types);
    vector<double> N(num_types);
	vector<double> A(num_types);
	vector<double> G(num_types);
	vector<int> DA(num_types);
	vector<bool>atom_defined(num_types);
    //debug
    //cout<<get_force_field()->get_parameters().params_in_each_section["vdw"].size()<<endl;
    for (vector<vector<string> >::iterator it = force_field_->get_parameters().params_in_each_section["vdw"].begin();it != force_field_->get_parameters().params_in_each_section["vdw"].end(); ++it)
	{
		type =str2int((*it)[0]);
		if (type == -1)
		{
			cout<<"MMFF94VDW::extract_VDW_parameters(): error! no numeric atom type defined for atom "<<name_type<<endl;
            return false;
        }
		alpha[type-1]=str2double((*it)[1]);
		N[type-1]=str2double((*it)[2]);
		A[type-1]=str2double((*it)[3]);
		G[type-1]=str2double((*it)[4]);
		DA[type-1]=str2int((*it)[5]);
		atom_defined[type-1]=true;
	}

    // now assemble all vdw parameters for all known atom type pairs
    for (int i = 0; i< num_types; i++)
	{
		for (int j = 0; j<num_types; j++)
		{
			// calculate index for R_IJ,R_ij,epsion_IJ for the field
            int index = i +j * num_types;
            int sym_index = i * num_types + j;
            //debug
            //cout<<atom_defined[i]<<" "<<atom_defined[j]<<endl;
			if (atom_defined[i] && atom_defined[j])
			{
				// calculate R_IJ,R_ij and epsilon_IJ
				// average method is arithmatic average
                //these calculations only need to be done once for each pair, 
				//we do them now and save them for later use
				double R_AA, R_BB, R_AB6, g_AB, g_AB2;
                double R_AB2, R_AB4, sqrt_a, sqrt_b;
 
                R_AA =A[i]*pow(alpha[i], 0.25);
                R_BB =A[j]*pow(alpha[j], 0.25);
                sqrt_a = sqrt(alpha[i]/N[i]);
                sqrt_b = sqrt(alpha[j] /N[j]);
      
                if (DA[i] == 1)
				{ // hydrogen bond donor
					fast_access.R_ij[index] = 0.5 * (R_AA + R_BB);
					fast_access.R_ij[sym_index]=fast_access.R_ij[index];
					R_AB2 = fast_access.R_ij[index] *  fast_access.R_ij[index];
					R_AB4 = R_AB2 * R_AB2;
                    R_AB6 = R_AB4 * R_AB2;
					   
					if (DA[j] == 2)
					{ // hydrogen bond acceptor
						fast_access.Epsilon_ij[index] = 0.5 * (181.16 * G[i]*G[j] *alpha[i] * alpha[j]) / (sqrt_a + sqrt_b) * (1.0 / R_AB6);
						fast_access.Epsilon_ij[sym_index]=fast_access.Epsilon_ij[index];
						// R_IJ is scaled to 0.8 for D-A interactions. The value used in the calculation of epsilon is not scaled. 
						fast_access.R_ij[index] = 0.8 *fast_access.R_ij[index];
						fast_access.R_ij[sym_index]=fast_access.R_ij[index];
					}
					else
						fast_access.Epsilon_ij[index]= (181.16 * G[i] * G[j] * alpha[i] * alpha[j]) / (sqrt_a + sqrt_b) * (1.0 / R_AB6);
					fast_access.Epsilon_ij[sym_index]=fast_access.Epsilon_ij[index];
					R_AB2 = fast_access.R_ij[index] *fast_access.R_ij[index];
					R_AB4 = R_AB2 * R_AB2;
					R_AB6 = R_AB4 * R_AB2;
					fast_access.R_ij_7[index]=R_AB6 *fast_access.R_ij[index];
					fast_access.R_ij_7[sym_index]=fast_access.R_ij_7[index];
				}
				else if(DA[j]== 1)
				{ // hydrogen bond donor
					fast_access.R_ij[index]=0.5*(R_AA + R_BB);
					fast_access.R_ij[sym_index]=fast_access.R_ij[index];
					R_AB2=fast_access.R_ij[index]*fast_access.R_ij[index];
					R_AB4 = R_AB2 * R_AB2;
					R_AB6 = R_AB4 * R_AB2;					   
     
					if ( DA[i]== 2)
					{ // hydrogen bond acceptor
						fast_access.Epsilon_ij[index]=0.5*(181.16*G[i]*G[j]*alpha[i]*alpha[j]) / (sqrt_a + sqrt_b) * (1.0 / R_AB6);
						fast_access.Epsilon_ij[sym_index]=fast_access.Epsilon_ij[index];
						// R_IJ is scaled to 0.8 for D-A interactions. The value used in the calculation of epsilon is not scaled. 
						fast_access.R_ij[index] = 0.8 *fast_access.R_ij[index];
						fast_access.R_ij[sym_index]=fast_access.R_ij[index];
					}
					else
						fast_access.Epsilon_ij[index] = (181.16 * G[i]* G[j] * alpha[i] * alpha[j]) / (sqrt_a + sqrt_b) * (1.0 / R_AB6);
					fast_access.Epsilon_ij[sym_index] =fast_access.Epsilon_ij[index];
					R_AB2 = fast_access.R_ij[index] *fast_access.R_ij[index];
					R_AB4 = R_AB2 * R_AB2;
					R_AB6 = R_AB4 * R_AB2;
					fast_access.R_ij_7[index]=fast_access.R_ij[index]*R_AB6;
					fast_access.R_ij_7[sym_index]=fast_access.R_ij_7[index];
				}
				else
				{
					g_AB = (R_AA - R_BB) / ( R_AA + R_BB);
                    g_AB2 = g_AB * g_AB;
					fast_access.R_ij[index]= 0.5 * (R_AA + R_BB) * (1.0 + 0.2 * (1.0 - exp(-12.0 * g_AB2)));
					fast_access.R_ij[sym_index]=fast_access.R_ij[index];
					R_AB2 = fast_access.R_ij[index] * fast_access.R_ij[index];
					R_AB4 = R_AB2 * R_AB2;
					R_AB6 = R_AB4 * R_AB2;
					fast_access.R_ij_7[index]=fast_access.R_ij[index]*R_AB6;
					fast_access.R_ij_7[sym_index]=fast_access.R_ij_7[index];
					fast_access.Epsilon_ij[index]=(181.16 *G[i]* G[j] *alpha[i] * alpha[j]) / (sqrt_a + sqrt_b) * (1.0 / R_AB6);
					fast_access.Epsilon_ij[sym_index]=fast_access.Epsilon_ij[index];
				} 
				fast_access.is_defined[index]=true;
				fast_access.is_defined[sym_index]=true;
			}
		}
	}
     return true;
}

// query a set of vdw parameters has defined for a given combination of atom types
bool MMFF94VDW::has_params(int i, int j) const
{
     if ((i > 0 && i <= num_of_atom_types_) && (j>0 && j <= num_of_atom_types_))
		 return (fast_access.is_defined[(i-1) +(j-1) * num_of_atom_types_]);
     else
         return false;
}
// assign the parameters for a given combination of atom types;
bool MMFF94VDW::assign_params(MMFF94VDW::VDWForceValues& param, int i, int j) const
{
	if (has_params(i, j))
	{
		int index=(i-1) + (j-1)*num_of_atom_types_;
		param.Epsilon_ij=fast_access.Epsilon_ij[index];
		param.R=fast_access.R_ij[index];
		param.is_defined=fast_access.is_defined[index];
		param.R_7=fast_access.R_ij_7[index];
        return true;
	}
    else
		return false;
}

// set up method
bool MMFF94VDW::setup()
{
	if (force_field_ == NULL)
	{
		cout<<"MMFF94VDW::setup(): force field bound component can not be found"<<endl;
        return false;
	}
    // clear the vdw parameter container
    vdw_data_holder_.clear();
    // tempararily set this component enabled;
    setenabled(true);
    // extract the L-J vdw parameters
    MMFF94* MMFF94_force_field = dynamic_cast<MMFF94*>(force_field_);
    if (MMFF94_force_field == NULL || MMFF94_force_field->has_initialized_param())
	{
		bool result = extract_VDW_parameters(force_field_->get_parameters());
        if (!result)
		{
			cout << "MMFF94VDW::setup(): can not access L-J vdw section"<<endl;
            return false;
		}
	}
    // now iterate all atom pairs and remove 1-2, 1-3 pairs
    VDWForceValues value;
    for (ATOMVec::iterator aiter1 = force_field_->get_atoms().begin();aiter1 != force_field_->get_atoms().end(); ++aiter1)
	{
		for (ATOMVec::iterator aiter2 = aiter1 + 1;aiter2 != force_field_->get_atoms().end(); ++aiter2)
		{
			// second check if is bonded
            if ((*aiter2)->is_bonded_to(*aiter1))
				continue;
            // then check if is geminal
            else if ((*aiter2)->is_geminal_to(*aiter1))
				continue;
			else
			{
				int atom_type_A = (*aiter1)->get_mmff94_type();
                int atom_type_B = (*aiter2)->get_mmff94_type();
                //debug
                //cout<<(*aiter1)->get_id()<<" "<<(*aiter2)->get_id()<<endl;
                vdw_data_holder_.push_back(MMFF94VDW::VDWData());
                vdw_data_holder_.back().atom1 = *aiter1;
                vdw_data_holder_.back().atom2 = *aiter2;                       
                // check if the given atom pair has predefined vdw parameters
                if(!assign_params(value, atom_type_A, atom_type_B))
                {
					cout<<"Warning: MMFF94VDW::setup(): cannot assign proper L-J vdw parameters for atom paris "<<(*aiter1)->get_id()<<"-"<<(*aiter2)->get_id()<<endl;
				}
                vdw_data_holder_.back().value = value;
			}
		}
	}
	return true;
}

double MMFF94VDW::update_energy()
{
       energy_ = 0.0;
       vdw_energy_ = 0.0;
       // loop all non bond atom pairs and calculate vdw and ele energies
       for (vector<VDWData>::size_type i = 0; i != vdw_data_holder_.size();++i)
	   {
		   const double r_ij = vdw_data_holder_[i].atom1->get_position().dist(vdw_data_holder_[i].atom2->get_position());
		   // fast integer powers: r^7 and (c)^7, avoiding expensive pow()
		   const double r_ij_2 = r_ij * r_ij;
		   const double r_ij_4 = r_ij_2 * r_ij_2;
		   const double r_ij_7 = r_ij_4 * r_ij_2 * r_ij;
		   const double c = 1.07*vdw_data_holder_[i].value.R/(r_ij+0.07*vdw_data_holder_[i].value.R);
		   const double c_2 = c * c;
		   const double c_4 = c_2 * c_2;
		   const double c_7 = c_4 * c_2 * c;
		   double first=vdw_data_holder_[i].value.Epsilon_ij*c_7;
		   double second=1.12*vdw_data_holder_[i].value.R_7/(r_ij_7+0.12*vdw_data_holder_[i].value.R_7)-2;
		   vdw_energy_ += first*second;
	   }
       energy_ = vdw_energy_;
       //energy_ = ele_energy_;
	   //cout<<"vdw energy:"<<energy_<<endl;
       return energy_;            
}

// calculate the forces imposed on each atoms by vdw and ele interaction
void MMFF94VDW::update_forces()
{
	if (force_field_ == NULL)
	{
		cout<<"MMFF94VDW::update_force(): error! this component doesn't bond to any force field"<<endl;
		return;
	}
    // iterate all non bond pairs and update forces
	double dbuf = 0.07;
    double gbuf = 0.12;
    for (vector<VDWData>::size_type i = 0; i != vdw_data_holder_.size();++i)
	{
		vector3 force_1 = vdw_data_holder_[i].atom1->get_force(), force_2 = vdw_data_holder_[i].atom2->get_force();
        vector3 direction(vdw_data_holder_[i].atom1->get_position() - vdw_data_holder_[i].atom2->get_position());
        double distance = direction.length();
        //double distance_2 = direction.length_2();
		direction.normalize();
        //direction = direction.normalize();
        if (!isNearZero(distance))
		{
			const double q = distance / vdw_data_holder_[i].value.R;
	    	const double q_2 = q * q;
	    	const double q6 = q_2 * q_2 * q_2;
			const double q7 = q6 * q;
			const double rpe = 1. / (q + dbuf);
			const double pe = (1. + dbuf) * rpe;
			const double rp7g = 1. / (q7 + gbuf);
			const double h = (1. + gbuf) * rp7g;
			const double gh = h - 2.;
			const double dgdp = -7. * q6 * h * rp7g;
			const double p2 = pe * pe;
			const double p4 = p2 * p2;
			const double f = vdw_data_holder_[i].value.Epsilon_ij * p4 * p2 * pe;
			const double dfdp = -7. * f * rpe;
			const double dgedp = f * dgdp + gh * dfdp;
			const double vdw_factor = dgedp / vdw_data_holder_[i].value.R;																	
			vector3 force = direction * vdw_factor * FORCE_FACTOR;
			force_1 -= force;
			force_2 += force;
			vdw_data_holder_[i].atom1->set_force(force_1);
            vdw_data_holder_[i].atom2->set_force(force_2);
		}
	}
}


                      
void MMFF94VDW::update(){;}                      


       
                                                                       
                                 
                                 
                            
                
