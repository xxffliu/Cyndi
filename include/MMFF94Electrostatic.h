#ifndef MMFF94Electrostatic_H
#define MMFF94Electrostatic_H

#ifndef FFCOMPONENT_H
#include "FFComponent.h"
#endif

#ifndef FORCEFIELD_H
#include "ForceField.h"
#endif

#include "MMFF94.h"

class MMFF94Ele: public FFComponent{
      public:
             // define a datastructure holding MMFF94 LJ Ele.parameters
             struct EleForceValues{                 
					double q_i;
					double q_j;
                    bool is_defined;
					bool is_14_interaction;
					double r_ij;
									 
                    };
             // define a datastructure holding Ele.parameters
             struct EleData{
                    ATOM* atom1;
                    ATOM* atom2;
                    //EleForceValues value;
					bool is_14_interaction;
                    };
             // define a hash datastructure of the chg file for fast access
             struct HashEleData_chg{
				 int a;
				 int b;
				 double bci;
				 bool is_defined;

			 };
			 //define a hash datastructure of the pbci file for fast access
			 struct HashEleData_pbci{
				 int a ;
				 double pbci;
				 double fcadj;
				 bool is_defined;
			 };
             // a map for storing the chg  content ,key is cxq 
             map<int,HashEleData_chg> fast_access_chg;
			 //a map for storing the  pbci   content ,key is the mmff94 atom type
			 map<int,HashEleData_pbci> fast_access_pbci;
			 // extract Eleff parameters from the TAFF parameter file
             // and build some data structure for fast access these data 
             bool extract_Ele_parameters(FFParameter& ffp);
             // query a set of parameters has defined for a given combination of atoms;
             bool has_params_chg(int i, int j, int type ) ;
			 bool has_params_pbci(int i);
             // assign the parameters to a given atom type combination
             // if no parameters are define for this 2 combination, return false and nothing is changed
             bool assign_params(EleForceValues &value,int i, int j) ;
             
             // default constructor
             MMFF94Ele();
             // constructor
             MMFF94Ele(ForceField& ff);
             // copy constructor
             MMFF94Ele(const MMFF94Ele& to_copy);
             // destructor
             virtual ~MMFF94Ele();
             
             // set up method
             virtual bool setup();
             
             // access methds
             // update current non-bond energy
             virtual double update_energy();
             // return current Eleenergy
             double get_Ele_energy()const;
             // update forces imposed on each atoms by non-bond energy
             virtual void update_forces();
             //Update the pair list.This method is called by the force field whenever ForceField::update is called.
             // It is used to recalculate the nonbonded pair list.
		     virtual void update();
			  //  Sets formal charges
			 bool SetFormalCharges();
			 //  Sets partial charges
			 bool SetPartialCharges();
			 
		     
       protected:
                 double Ele_energy_;
                 
       private:
               std::vector<EleData> Ele_data_holder_;
               MMFF94* mmff94_force_field;
};

#endif		     
