#ifndef MMFF94VDW_H
#define MMFF94VDW_H

#ifndef FFCOMPONENT_H
#include "FFComponent.h"
#endif

#ifndef FORCEFIELD_H
#include "ForceField.h"
#endif

class MMFF94VDW : public FFComponent{
      public:
             // define a datastructure holding MMFF94 LJ vdw parameters
             struct VDWForceValues{
				    double R;                 
					double Epsilon_ij;
                    bool  is_defined;
					double  R_7;
				 
                    };
             // define a datastructure holding vdw parameters
             struct VDWData{
                    ATOM* atom1;
                    ATOM* atom2;
                    VDWForceValues value;					
                    };
             // define a hash datastructure for fast access
             struct HashVDWData{
			     std::vector<bool> is_defined;
			     std::vector<double>R_ij;
			     std::vector<double> Epsilon_ij;
			     std::vector<double>R_ij_7;
				   
                    };
             // a holder for fast access through index
             HashVDWData fast_access;
             // extract VDW ff parameters from the TAFF parameter file
             // and build some data structure for fast access these data 
             bool extract_VDW_parameters(FFParameter& ffp);
             // query a set of parameters has defined for a given combination of atoms;
             bool has_params(int i, int j) const;
             // assign the parameters to a given atom type combination
             // if no parameters are define for this 2 combination, return false and nothing is changed
             bool assign_params(VDWForceValues& param, int i, int j) const;
             
             // default constructor
             MMFF94VDW();
             // constructor
             MMFF94VDW(ForceField& ff);
             // copy constructor
             MMFF94VDW(const MMFF94VDW& to_copy);
             // destructor
             virtual ~MMFF94VDW();
             
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
               //int num_of_atom_types_;
};

#endif		     
