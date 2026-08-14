#ifndef TAFFBEND_H
#define TAFFBEND_H

#ifndef FFCOMPONENT_H
#include "FFComponent.h"
#endif

#ifndef FORCEFIELD_H
#include "ForceField.h"
#endif

#ifndef VECTOR3_H
#include "vector3.h"
#endif

/*#ifndef AVOGADRO
#define AVOGADRO  6.0221367E+23L;
#endif
*/

//static double AVOGADRO = 6.0221367e+23;

using namespace std;


class TAFFBend : public FFComponent{
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
                    ForceValues value;
                    };
             // a structure holding all TAFF angle bend parameters for fast access
             struct BendHashData{
                    //int num_of_atom_types;
                    //double* theta0;
                    vector<double> theta0;
                    //double* k;
                    vector<double> k;
                    //bool* is_defined;
                    vector<bool> is_defined;
                    void clear();
                    };
             // a holder for fast access through atoms type index
             BendHashData fast_access;
             
             // extract angle bend ff parameters from the TAFF parameter file
             // and build some data structure for fast access these data 
             bool extract_AB_parameters(FFParameter& ffp);
             // query a set of parameters has defined for a given combination of atoms;
             bool has_params(int i, int j, int k) const;
             // return the parameters for a given atom type combination
             ForceValues get_params(int i, int j, int k) const;
             // assign the parameters to a given atom type combination
             // if no parameters are define for this 3 combination, return false and nothing is changed
             bool assign_params(ForceValues& param, int i, int j, int k) const;
             
             // default constructor
             TAFFBend();
             // constructor
             TAFFBend(ForceField& ff);
             // copy constructor
             TAFFBend(const TAFFBend& to_copy);
             // destructor
             virtual ~TAFFBend();
             
             // set up method
             virtual bool setup();
             
             // access methds
             // update current bend energy
             virtual double update_energy();
             // update forces imposed on each atoms by bend energy
             virtual void update_forces();
             
      private:
              std::vector<BendData> bend_data_holder_;
              int num_of_atom_types_;
};

#endif
