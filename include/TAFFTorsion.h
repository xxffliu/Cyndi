#ifndef TAFFTORSION_H
#define TAFFTORSION_H

#ifndef FFCOMPONENT_H
#include "FFComponent.h"
#endif

#ifndef FORCEFIELD_H
#include "ForceField.h"
#endif

/*#ifndef AVOGADRO
#define AVOGADRO  6.0221367E+23L;
#endif*/

//static double AVOGADRO = 6.0221367e+23;

class TAFFTorsion : public FFComponent{
      public:
             // define a structure holding TAFF torsion force parameters
             struct ForceValues{
                    double V;
                    double n;
                    };
             // define a structure holding TAFF torsion atom types and force parameters
             struct TorsionData{
                    ATOM* atom1;
                    ATOM* atom2;
                    ATOM* atom3;
                    ATOM* atom4;
                    int bond_order;
                    ForceValues value;
                    };
             // a structure holding all TAFF torsion parameters for fast access
             struct TorsionHashData{
                    //int num_of_atom_types;
                    //double* V;
                    vector<double> V;
                    //double* n;
                    vector<double> n;
                    //bool* is_defined;
                    vector<double> is_defined;
                    };
             // a holder for fast access through atoms type index and bond order
             std::map<int,TorsionHashData> fast_access;
             // extract torsion ff parameters from the TAFF parameter file
             // and build some data structure for fast access these data 
             bool extract_TOR_parameters(FFParameter& ffp);
             // query a set of parameters has defined for a given combination of atoms;
             bool has_params(int i, int j, int k, int l, int bond_order);
             // return the parameters for a given atom type combination
             TAFFTorsion::ForceValues get_params(int i, int j, int k, int l, int bond_order);
             // assign the parameters to a given atom type combination
             // if no parameters are define for this 4 combination, return false and nothing is changed
             bool assign_params(ForceValues& param, int i, int j, int k, int l, int bond_order);
             
             // default constructor
             TAFFTorsion();
             // constructor
             TAFFTorsion(ForceField& ff);
             // copy constructor
             TAFFTorsion(const TAFFTorsion& to_copy);
             // destructor
             virtual ~TAFFTorsion();
             
             // set up method
             virtual bool setup();
             
             // access methds
             // update current bend energy
             virtual double update_energy();
             // update forces imposed on each atoms by bend energy
             virtual void update_forces();
             // a method convert from symbolic bond type to numeric type;
		     inline int get_bond_order(const string& type) const{
                    if (type == "1") return 1;
                    else if (type == "2") return 2;
                    else if (type == "3") return 3;
                    else if (type == "ar") return 4;
                    else if (type == "am") return 5;
                    else
					{
						cout<<"TAFFTorsion::get_bond_order: No bond order returned with bond type: "<<type<<endl;
						return 0;
					}
             }
                    
       private:
               // a container holding all bond ff parameters
               std::vector <TorsionData> torsion_data_holder_; 
               int num_of_atom_types_;
};

#endif
