#ifndef TAFFOOP_H
#define TAFFOOP_H

#ifndef FFCOMPONENT_H
#include "FFComponent.h"
#endif

#ifndef FORCEFIELD_H
#include "ForceField.h"
#endif

/*#ifndef AVOGADRO
#define AVOGADRO  6.0221367E+23L;
#endif
*/

//static double AVOGADRO = 6.0221367e+23;

class TAFFOOP : public FFComponent{
      public:
             enum{
                  NUM_OOP_TYPES = 6
                  };
             // a structure holding OOP force field parameters
             struct OOPData{
                    double k;
                    ATOM* root;
                    ATOM* atom1;
                    ATOM* atom2;
                    ATOM* atom3;
                    };
             // a structure holding all OOP parameters for fast access
             struct OOPHashData{
                    std::vector<double> k0;
                    //double* k0;
                    std::vector<bool> is_defined;
                    //bool* is_defined;
                    };
             OOPHashData fast_access;
             // extract angle oop ff parameters from the TAFF parameter file
             // and build some data structure for fast access these data 
             bool extract_OOP_parameters(FFParameter& ffp);
             // query a set of parameters has defined for a given atom;
             bool has_params(int i) const;
             // return the parameters for a given atom type
             double get_params(int i) const;
             // assign the parameters to a given atom type
             // if no parameters are define for this root atom type, return false and nothing is changed
             bool assign_params(double& k, int i) const;
             
             // default constructor
             TAFFOOP();
             // constructor
             TAFFOOP(ForceField& ff);
             // copy constructor
             TAFFOOP(const TAFFOOP& to_copy);
             // destructor
             virtual ~TAFFOOP();
             
             // set up method
             virtual bool setup();
             
             // access methds
             // update current bend energy
             virtual double update_energy();
             // update forces imposed on each atoms by bend energy
             virtual void update_forces();
      private:
              vector<OOPData> oop_data_holder_;
              int num_of_atom_types_;
};

#endif
