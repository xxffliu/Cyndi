// Molecular Mechanics: TRIPOS Force Field
#ifndef TAFF_H
#define TAFF_H

#ifndef FORCEFIELD_H
#include "ForceField.h"
#endif

#ifndef FFCOMPONENT_H
#include "FFComponent.h"
#endif

#ifndef MOL_H
#include "MOL.h"
#endif

#ifndef UTILITY_H
#include "utility.h"
#endif

#ifndef INPUT_H
	#include "ParamInput.h"
#endif

// TRIPOS Force Field class

// a enum type for the tripos force field componeents
enum TAFFCOMPONENT
{
	TAFF_BOND_STRETCH = 1,
	TAFF_ANGLE_BEND = 2,
	TAFF_DIHEDRAL_TORSION = 3,
	TAFF_OOP_BEND = 4,
	TAFF_VDW = 5,
	TAFF_ELE = 6
};
//static const string DEFAULT_TAFF_PARAM_FILE = MOGAParam_.FFParam_;
static const string DEFAULT_TAFF_PARAM_FILE = "TAFF.parm";
class TAFF : public ForceField{
      public:
             // default constructor
             TAFF();
             // constructor
             TAFF(MOL& mol);
             // special constructor with a MOL and only specified energy terms are included
             TAFF(MOL& mol, vector<TAFFCOMPONENT> energy_terms_list);
             // copy constructor
             TAFF(const TAFF& taff);
             // destructor
             virtual ~TAFF();
             // assessment operator
             const TAFF& operator=(const TAFF& taff);
             // clear method
             virtual void clear();
             
             // setup method
             virtual bool specific_setup();
             
             // Accessors attatched to TAFF
             // bond stretch
             double get_stretch_energy() const;
             // angle bend
             double get_bend_energy() const;
             // proper torsion
             double get_torsion_energy() const;
             // oop
             double get_oop_energy() const;
             // VDW
             double get_vdw_energy() const;
             // electrostatic
             double get_ele_energy() const;
             
             // return true if the parameters have already been initialized;
             bool has_initialized_param() const;
             // return update frequency, default is 20
             int get_update_frequency() const;
             
             //get final result in string form
             virtual string get_results() const;
             
      protected:
                string file_name_;
                bool param_is_initialized_;
};
extern TAFF taff;             
#endif             
