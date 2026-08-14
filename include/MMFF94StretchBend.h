#ifndef MMFFStretchBend_H
#define MMFFStretchBend_H

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


class MMFF94Str_Bend : public FFComponent{
      public:
             // a structure holding angle Str_Str_Bend force parameters
             struct ForceValues{
				  
                    double k_ijk;
					double k_kji;
                    };
             // a structure holding all the atoms and force parameters belong to an angle
             struct Str_BendData{
                    ATOM* atom1;
                    ATOM* atom2;
                    ATOM* atom3;
					unsigned int strbnd_type;
                    ForceValues value;
					 double r_ij;
					double r_kj;
                    double theta0;
                    };
             // a structure holding all MMFF94 angle Str_Str_Bend parameters for fast access
             struct Str_BendHashData{
                    //int stretch bend class;
				    int strbnd__type;
					//int atom1 type;
					int type1;
					//int atom2 type;
					int type2;
					//int atom3 type
					int type3;
				    //double* k_ijk;
                    double k_ijk;
					//double*k_kji;
					double k_kji;				
                    bool  is_defined;
			 };
			
             // fast access through CXS and index of the stretch bend 
             map<int,Str_BendHashData> fast_access_strbnd;
			 //fast access through CXS and index  of the streth bond by emperical rule
			 map<int,Str_BendHashData>  emp_fast_acess_sb;
             // extract angle Str_Str_Bend ff parameters from the MMFF94 parameter file
             // and build some data structure for fast access these data 
             bool extract_AB_parameters(FFParameter& ffp);
             // query a set of parameters has defined for a given combination of atoms;
             bool has_params(int i, int j, int k,int strbnd_type) ;
			 bool has_params_emperical(int i,int j,int k);
             // assign the parameters to a given atom type combination
             // if no parameters are define for this 3 combination, return false and nothing is changed
             bool assign_params(ForceValues& param, int i, int j, int k,int strbnd_type);
			 bool assign_params_emprical(ForceValues& param, int i, int j, int k);
			 //! \return The row of the element atom in the periodic table
			 int GetElementRow(ATOM *atom);
             // default constructor
            MMFF94Str_Bend();
            MMFF94Str_Bend(ForceField& ff);
             // copy constructor
            MMFF94Str_Bend(const MMFF94Str_Bend& to_copy);
           
             // destructor
             virtual ~MMFF94Str_Bend();
             
             // set up method
             virtual bool setup();
             
             // access methds
             // update current Str_Str_Bend energy
             virtual double update_energy();
             // update forces imposed on each atoms by Str_Str_Bend energy
             virtual void update_forces();
			 
             
      private:
              std::vector<Str_BendData> Str_Bend_data_holder_;
              MMFF94* mmff94_force_field;
};
extern  MMFF94Str_Bend mmff94stretchbend;
#endif
