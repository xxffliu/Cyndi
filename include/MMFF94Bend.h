#ifndef MMFFBEND_H
#define MMFF94BEND_H

#ifndef FFCOMPONENT_H
#include "FFComponent.h"
#endif

#ifndef FORCEFIELD_H
#include "ForceField.h"
#endif

#ifndef VECTOR3_H
#include "vector3.h"
#endif

#include "MMFF94.h"

//static double AVOGADRO = 6.0221367e+23;

using namespace std;


class MMFF94Bend : public FFComponent{
      public:
             // a structure holding angle bend force parameters
             struct ForceValues{
                    double theta0;
                    double k;
                    };
             // a structure holding all the atoms and force parameters belong to an angle
             struct BendData{
                    ATOM* atom1;
                    ATOM* atom2;
                    ATOM* atom3;
					unsigned int angle_type;
                    ForceValues value;
                    };
             // a structure holding all TAFF angle bend parameters for fast access
             struct BendHashData{
                    //int angle class;
				    int angle_type;
					//int atom1 type;
					int type1;
					//int atom2 type;
					int type2;
					//int atom3 type
					int type3;
				    //double* k;
                    double  k;
                    //double* r0;
                    double  theta0;
                    //bool* is_defined;
                    bool  is_defined;
                    };
			 struct EqualValues{
				    int equl2;
					int equl3;
					int equl4;
					int equl5;			    
			        };
			 map<int,EqualValues>Equal_access;

             // fast access through CXA and index of the angles
             map<int,BendHashData> fast_access;
             // extract angle bend ff parameters from the MMFF94 parameter file
             // and build some data structure for fast access these data 
             bool extract_AB_parameters(FFParameter& ffp);
             // query a set of parameters has defined for a given combination of atoms;
             bool has_params(int i, int j, int k,int angle_type) ;
             // assign the parameters to a given atom type combination
             // if no parameters are define for this 3 combination, return false and nothing is changed
             bool assign_params(ForceValues& param, int i, int j, int k,int angle_type);
			 // return the reference bond angle for given atom types combination
			 double GetBondAngle(ATOM* a, ATOM* b, ATOM* c);
             // default constructor
            MMFF94Bend();
            MMFF94Bend(ForceField& ff);
             // copy constructor
            MMFF94Bend(const MMFF94Bend& to_copy);
			// overloading assignment operator
			MMFF94Bend& operator=(const MMFF94Bend& bend);
           
             // destructor
             virtual ~MMFF94Bend();
             
             // set up method
             virtual bool setup();
             
             // access methds
             // update current bend energy
             virtual double update_energy();
             // update forces imposed on each atoms by bend energy
             virtual void update_forces();
			 //! \return the level 2 equivalent atom type for type (mmffdef.par)
			 int EqLvl2(int type);
			 //! \return the level 3 equivalent atom type for type (mmffdef.par)
			 int EqLvl3(int type);
			 //! \return the level 4 equivalent atom type for type (mmffdef.par)
			 int EqLvl4(int type);
			 //! \return the level 5 equivalent atom type for type (mmffdef.par)
			 int EqLvl5(int type);
           
			 				
      private:
              std::vector<BendData> bend_data_holder_;
              MMFF94* mmff94_force_field;
};
extern MMFF94Bend mmff94bend;
//extern map<int,MMFF94Bend::BendHashData> bend_fast_access;
#endif
