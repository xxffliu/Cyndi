#ifndef MMFF94OOP_H
#define MMFF94OOP_H

#ifndef FFCOMPONENT_H
#include "FFComponent.h"
#endif

#ifndef FORCEFIELD_H
#include "ForceField.h"
#endif

#include "MMFF94.h"
//#include "MMFF94Bend.h"

//static double AVOGADRO = 6.0221367e+23;

class MMFF94OOP : public FFComponent{
      public:

             // a structure holding OOP force field parameters
             struct OOPData{
                    double k;
					//double willson_angle;
                    ATOM* atom1;
                    ATOM* atom2;
                    ATOM* atom3;
					ATOM* center;
                    };
             // a structure holding all OOP parameters for fast access
             struct OOPHashData{
				 double k0;
				 //double* k0;
				 bool is_defined;
				 //bool* is_defined;
				// bool willson_angle;
				 //bool*willson angle;
				 int type1;
				 //atom 1;
				 int type2;
				 //atom 2;
				 int type3;
				 //atom 3;
				 int type_center;
				 //center atom ;
                    };
			 struct ForceValues{
				 double k0;
			 };
			 
			 map<int ,OOPHashData> fast_access;
             // extract angle oop ff parameters from the MMFF94 parameter file
             // and build some data structure for fast access these data 
             bool extract_OOP_parameters(FFParameter& ffp);
             // query a set of parameters has defined for a given atom;
             bool has_params(int i,int j,int k,int l) ;
             // assign the parameters to a given atom type
             // if no parameters are define for this root atom type, return false and nothing is changed
             bool assign_params(ForceValues &param, int i,int j,int k,int l) ;
             
             // default constructor
             MMFF94OOP();
             // constructor
             MMFF94OOP(ForceField& ff);
             // copy constructor
             MMFF94OOP(const MMFF94OOP& to_copy);
             // destructor
             virtual ~MMFF94OOP();
             
             // set up method
             virtual bool setup();
             
             // access methds
             // update current bend energy
             virtual double update_energy();
             // update forces imposed on each atoms by bend energy
             virtual void update_forces();
      private:
              vector<OOPData> oop_data_holder_;
              MMFF94* mmff94_force_field;
};

#endif
