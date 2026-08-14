#ifndef TAFFNONBOND_H
#define TAFFNONBOND_H

#ifndef FFCOMPONENT_H
#include "FFComponent.h"
#endif

#ifndef FORCEFIELD_H
#include "ForceField.h"
#endif

/*#ifndef AVOGADRO
#define AVOGADRO  6.0221367E+23L;
#endif

#ifndef VACUUM_PERMITTIVITY
#define VACUUM_PERMITTIVITY    	= 8.85419E-12L;
#endif

#ifndef PI
#define PI 3.1415926535897932384626433
#endif

#ifndef e0
#define e0 1.60217738E-19L
#endif
*/
//static double AVOGADRO = 6.0221367e+23;


class TAFFVDW : public FFComponent{
      public:
             // define a datastructure holding TAFF LJ vdw parameters
             struct VDWForceValues{
                    double A;
                    double B;
                    };
             // define a datastructure holding vdw parameters
             struct VDWData{
                    ATOM* atom1;
                    ATOM* atom2;
                    VDWForceValues value;
                    bool is_14_interaction;
                    };
             // define a hash datastructure for fast access
             struct HashVDWData{
                    std::vector<double> Aij;
                    std::vector<double> Bij;
                    std::vector<bool> is_defined;
                    };
             // a holder for fast access through index
             HashVDWData fast_access;
             // extract VDW ff parameters from the TAFF parameter file
             // and build some data structure for fast access these data 
             bool extract_VDW_parameters(FFParameter& ffp);
             // query a set of parameters has defined for a given combination of atoms;
             bool has_params(int i, int j) const;
             // return the parameters for a given atom type combination
             TAFFVDW::VDWForceValues get_params(int i, int j) const;
             // assign the parameters to a given atom type combination
             // if no parameters are define for this 2 combination, return false and nothing is changed
             bool assign_params(VDWForceValues& param, int i, int j) const;
             
             // default constructor
             TAFFVDW();
             // constructor
             TAFFVDW(ForceField& ff);
             // copy constructor
             TAFFVDW(const TAFFVDW& to_copy);
             // destructor
             virtual ~TAFFVDW();
             
             // set up method
             virtual bool setup();
             
             // access methds
             // update current non-bond energy
             virtual double update_energy();
             // return current vdw energy
             double get_vdw_energy()const;
             // update forces imposed on each atoms by non-bond energy
             virtual void update_forces();
             //Update the pair list.This method is called by the force field whenever ForceField::update is called.
             // It is used to recalculate the nonbonded pair list.
		     virtual void update();
		     
       protected:
                 double vdw_energy_;
                 
       private:
               std::vector<VDWData> vdw_data_holder_;
               int num_of_atom_types_;
};

#endif		     
