#ifndef MMFF94STRETCH_H
#define MMFF94STRETCH_H

#ifndef FFCOMPONENT_H
#include "FFComponent.h"
#endif

#ifndef FORCEFIELD_H
#include "ForceField.h"
#endif

#include "MMFF94.h"

class MMFF94Stretch : public FFComponent{
      public:
             // a structure holding bond stretch spring constant and equilium bond length
             struct ForceValues{
             double k;
             double r0;
             };
             // a structure holding whole bond stretch ff parameters
             struct StretchData{
				 ATOM*atom1;
				 ATOM*atom2;
				 int bond_type;
				 ForceValues value;
             };
             // a structure holding all bond stretch parameters for fast access
			 struct StretchHashData{
                    //int bond class;
				   int bond_type;
					//int atom1 type;
				   int type1;
					//int atom2 type;
				   int type2;
				    //double* k;
				   double k;
                    //double* r0;
                   double r0;
                    //bool* is_defined;
                   bool is_defined;
                    };
			 struct StretchHashData_emprical{			
					//int atom1 type;
				   int type1;
					//int atom2 type;
				   int type2;
				    //double* k;
				   double k;
                    //double* r0;
                   double r0;
                    //bool* is_defined;
                   bool is_defined;

			 };
             // fast access through CXB and index of the two bond atoms
             map<int,StretchHashData> fast_access;
			 map<int,StretchHashData_emprical> fast_access_emperical;
             // extract bond stretch ff parameters from the TAFF parameter file
             // and build some data structure for fast access these data 
             bool extract_BS_parameters(FFParameter& ffp);
             // query a set of parameters has defined for a given combination of atoms;
             bool has_params(int i, int j, int bond_type);
			 bool has_params_emprical(int i,int j);
             // assign the parameters to a given atom type combination
             // if no parameters are define for this 2 combination, return false and nothing is changed
             bool assign_params(ForceValues& param, int i, int j, int bond_type);
			 bool assign_params_emperical(ForceValues& param, int i, int j);
	         double GetBondLength(ATOM* a, ATOM* b);
             
             
             // default constructor
             MMFF94Stretch();
             // constructor
             MMFF94Stretch(ForceField& ff);
             // copy constructor
             MMFF94Stretch(const MMFF94Stretch& to_copy);
			 // assignment operator
			 MMFF94Stretch& operator=(const MMFF94Stretch& str);
             // destructor
             virtual ~MMFF94Stretch();
             
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
						cout<<"MMFF94Stretch::get_bond_order: No bond order returned with bond type: "<<type<<endl;
						return 0;
					}
             }
	     
      private:
              // a container holding all bond ff parameters
             std::vector <StretchData> stretch_data_holder_; 
             MMFF94* mmff94_force_field;
};
extern MMFF94Stretch mmff94stretch;
// fast access through CXB and index of the two bond atoms
//extern map<int,MMFF94Stretch::StretchHashData> stretch_fast_access;
//extern map<int,MMFF94Stretch::StretchHashData_emprical> stretch_fast_access_emperical;
#endif   

