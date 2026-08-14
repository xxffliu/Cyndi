#ifndef TAFFSTRETCH_H
#define TAFFSTRETCH_H

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

class TAFFStretch : public FFComponent{
      public:
             // a structure holding bond stretch spring constant and equilium bond length
             struct ForceValues{
             double k;
             double r0;
             };
             // a structure holding whole bond stretch ff parameters
             struct StretchData{
             ATOM* atom1;
             ATOM* atom2;
             int bond_order;
             ForceValues value;
             };
             // a structure holding all bond stretch parameters for fast access
             struct StretchHashData{
                    
                    //double* k;
                    vector<double> k;
                    //double* r0;
                    vector<double> r0;
                    //bool* is_defined;
                    vector<bool> is_defined;
                    };
             // fast access through bond order and index of the two bond atoms
             map<int,StretchHashData> fast_access;
             // extract bond stretch ff parameters from the TAFF parameter file
             // and build some data structure for fast access these data 
             bool extract_BS_parameters(FFParameter& ffp);
             // query a set of parameters has defined for a given combination of atoms;
             bool has_params(int i, int j, int bond_order);
             // return the parameters for a given atom type combination
             ForceValues get_params(int i, int j, int bond_order);
             // assign the parameters to a given atom type combination
             // if no parameters are define for this 2 combination, return false and nothing is changed
             bool assign_params(ForceValues& param, int i, int j, int bond_order);
             
             
             // default constructor
             TAFFStretch();
             // constructor
             TAFFStretch(ForceField& ff);
             // copy constructor
             TAFFStretch(const TAFFStretch& to_copy);
             // destructor
             virtual ~TAFFStretch();
             
             // setup methods
             virtual bool setup();
             
             // assessors
             // Calculates and returns the component's energy.
		     virtual double update_energy();
		     // Calculates and returns the component's forces.
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
						cout<<"TAFFStretch::get_bond_order: No bond order returned with bond type: "<<type<<endl;
						return 0;
					}
             }
	     
      private:
              // a container holding all bond ff parameters
             std::vector <StretchData> stretch_data_holder_; 
             int num_of_atom_types_;
};

#endif
