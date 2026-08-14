// Molecular Mechanics:  Force Field 94
#ifndef MMFF94_H
#define MMFF94_H

#ifndef FORCEFIELD_H
#include "ForceField.h"
#endif


#ifndef FFCOMPONENT_H
#include "FFComponent.h"
#endif

#ifndef VECTOR3_H
#include "vector3.h"
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

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// MMFF94 Force Field class

// a enum type for the mmff94 force field componeents
// a enum type for the mmff94 force field componeents
enum MMFF94COMPONENT
{
	MMFF94_BOND_STRETCH = 1,
	MMFF94_ANGLE_BEND = 2,
	MMFF94_DIHEDRAL_TORSION = 3,
	MMFF94_OOP_BEND = 4,
	MMFF94_STRETCH_BEND = 5,
	MMFF94_VDW = 6,
	MMFF94_ELE = 7
};
//static const string DEFAULT_MMFF94_PARAM_FILE = DEFAULT_MMFF94_PARAM_FILE;
static const string DEFAULT_MMFF94_PARAM_FILE = "MMFF94.parm";

// a structure depositing MMFF94 numeric and symbolic atom type
struct MMFF94AtomType
{
	int numeric;
	string symbolic;
	MMFF94AtomType():numeric(0),symbolic("UNK"){};
};

class MMFF94 : public ForceField
{     
   
      public:
             // default constructor
             MMFF94();
             // constructor
             MMFF94(MOL& mol);
             // special constructor with a MOL and only specified energy terms are included
             MMFF94(MOL& mol, vector<MMFF94COMPONENT> energy_terms_list);
             // copy constructor
             MMFF94(const MMFF94& mmff94);
             // destructor
             virtual ~MMFF94();
             // assessment operator
             const MMFF94& operator=(const MMFF94& mmff94);
             // clear method
             virtual void clear();
			 //bool extract_BS_parameters(FFParameter& ffp);
              //! \return Get the MMFF94 atom type for atom
             MMFF94AtomType GetType(ATOM *atom) ;
			 //perceiving  the atom or the bond whether is in an aromatic ring
             void PerceiveAromatic(); 
			 bool SetTypes();
			 //! \return The bond type (BTIJ)
			 int GetBondType(ATOM* a, ATOM* b);
			 //! \return The angle type (ATIJK)
			 int GetAngleType(ATOM* a, ATOM* b, ATOM *c);
			 //! \return The stretch-bend type (SBTIJK)
			 int GetStrBndType(ATOM* a, ATOM* b, ATOM *c);
			 //! \return The torsion type (TTIJKL)
			 int GetTorsionType(ATOM* a, ATOM* b, ATOM *c, ATOM *d);
			 //!we have assingned the content of mmffprop.par's in the head file of the FFComponent
			 //! \return true if atomtype has sbmb set in mmffprop.par
			 bool HasSbmbSet(int atomtype);
			 //! \return true if atomtype has pilp set in mmffprop.par
			 bool HasPilpSet(int atomtype);
			 //! \return true if atomtype has arom set in mmffprop.par
			 bool HasAromSet(int atomtype);
			 //! \return true if atomtype has lin set in mmffprop.par
			 bool HasLinSet(int atomtype);
			 //! \return the crd value for the atomtype in mmffprop.par
			 int GetCrd(int atomtype);
			 //! \return the val value for the atomtype in mmffprop.par
			 int GetVal(int atomtype); 
			 //! \return the mltb value for the atomtype in mmffprop.par
			 int GetMltb(int atomtype);
			//return the canonical bond index
			 unsigned int GetCXB(int type, int a, int b) ;
			 //! \return the canonical angle index
			 unsigned int GetCXA(int type, int a, int b, int c) ;
			 //! \return the canonical stretch-bend index
			 unsigned int GetCXS(int type, int a, int b, int c) ;
			 //! \return the canonical out-of-plane index
			 unsigned int GetCXO(int a, int b, int c, int d) ;
			 //! \return the canonical torsion index
			 unsigned int GetCXT(int type, int a, int b, int c, int d) ;
			 //! \return the canonical bond-charge-increment index
			 unsigned int GetCXQ(int type, int a, int b) ;
			 //! \return the U value for the atom from table X page 631
			 double GetUParam(ATOM* atom);
			 //! \return the Z value for the atom from table VI page 628
			 double GetZParam(ATOM* atom);
			 //! \return the C value for the atom from table VI page 628
			 double GetCParam(ATOM* atom);
			 //! \return the V value for the atom from table X page 631
			 double GetVParam(ATOM* atom);
			 //! return the covalent radius from Blom and Haaland, value from etab if not available
			 double GetCovalentRadius(ATOM* a);
		     //! return the bond length calculated with a modified version of the Schomaker-Stevenson rule
             double GetRuleBondLength(ATOM* a,ATOM* b);
			 //get the value of the Allred Rochow ElectroNegativity
			 double GetAllredRochowElectroNeg(ATOM*a); 
			 
             // setup method
             virtual bool specific_setup();
             
             // Accessors attatched to MMFF94
             // bond stretch
             double get_stretch_energy() const;
             // angle bend
             double get_bend_energy() const;
             // proper torsion
             double get_torsion_energy() const;
             // oop
             double get_oop_energy() const;
			 // stretch-bend
			 double get_stretch_bend_energy() const;
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
//extern MMFF94 mmff94;             
#endif             
