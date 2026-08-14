#ifndef MMFF94TORSION_H
#define MMFF94TORSION_H

#ifndef FFCOMPONENT_H
#include "FFComponent.h"
#endif

#ifndef FORCEFIELD_H
#include "ForceField.h"
#endif

#include "MMFF94.h"

//static double AVOGADRO = 6.0221367e+23;

class MMFF94Torsion : public FFComponent{
      public:
             // define a structure holding TAFF torsion force parameters
             struct ForceValues{			
                    double V1;
                    double V2;
					double V3;
					bool hasDefined;
					ForceValues():V1(0.0),V2(0.),V3(0.), hasDefined(false){}
                    };
             // define a structure holding TAFF torsion atom types and force parameters
             struct TorsionData{
                    ATOM* atom1;
                    ATOM* atom2;
                    ATOM* atom3;
                    ATOM* atom4;
					int torsion_type;
					ForceValues value;
                    };
             // a structure holding all MMFF94 torsion parameters for fast access
             struct TorsionHashData{
                    //int num_of_atom_types;
				 double V1;
                 double V2;
	             double V3;
                 bool  is_defined;
				 int torsion_type;
			 };
             // a holder for fast access through atoms type index and bond order
             map<int,TorsionHashData> fast_access_tor;
             // extract torsion ff parameters from the TAFF parameter file
             // and build some data structure for fast access these data 
             bool extract_TOR_parameters(FFParameter& ffp);
             // query a set of parameters has defined for a given combination of atoms;
             bool has_params(int i, int j, int k, int l, int torsion_type);
             // assign the parameters to a given atom type combination
             // if no parameters are define for this 4 combination, return false and nothing is changed
             bool assign_params(ForceValues& param, int i, int j, int k, int l, int torsion_type);
             
             // default constructor
             MMFF94Torsion();
             // constructor
             MMFF94Torsion(ForceField& ff);
             // copy constructor
             MMFF94Torsion(const MMFF94Torsion& to_copy);
             // destructor
             virtual ~MMFF94Torsion();
             
             // set up method
             virtual bool setup();
             
             // access methds
             // update current bend energy
             virtual double update_energy();
             // update forces imposed on each atoms by bend energy
             virtual void update_forces();
		   
                    
       private:
               // a container holding all bond ff parameters
               std::vector <TorsionData> torsion_data_holder_; 
               MMFF94* mmff94_force_field;
};

#endif
