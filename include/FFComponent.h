#ifndef FFCOMPONENT_H
#define FFCOMPONENT_H

/*#ifndef ATOM_H
#define ATOM_H
#endif
*/
#include <string>

using namespace std;
class ForceField;
// const vectors storing the properties for each mmff94 atom type
const  unsigned int MMFF94_ASPEC[] = {	6 ,	6 ,	6 ,	6 ,	1,	8 ,	8 ,	7 ,	7 ,	7 ,	9 ,	17,	35,	53,	16,	16,	16,	16,	14,	6 ,	1 ,	6 ,	1 ,	1 ,	15,	15,	1 ,	1 ,	1 ,	6 ,	1 ,	8 ,	1 ,	7 ,	8 ,	1 ,	6 ,	7 ,	7 ,	7 ,	6 ,	7 ,	7 ,	16,	7 ,	7 ,	7 ,	7 ,	8 ,	1 ,	8 ,	1 ,	7 ,	7 ,	7 ,	7 ,	6 ,	7 ,	8 ,	6 ,	7 ,	7 ,	6 ,	6 ,	7 ,	7 ,	7 ,	7 ,	7 ,	8 ,	1 ,	16,	16,	16,	15,	7 ,	17,	6 ,	7 ,	6 ,	7 ,	7 ,	0 ,	0 ,	0 , 0 ,	26,	26,	9 ,	17,	35,	3 ,	11,	19,	30,	20,	29,	29,	12
};
const unsigned int MMFF94_CRD[]=     {  4,	3,	3,	2,	1,	2,	1,	3,	2,	3,	1,	1,	1,	1,	2,	1,	3,	4,	4,	4,	1,	4,	1,	1,	4,	3,	1,	1,	1,	3,	1,	1,	1,	4,	1,	1,	3,	2,	3,	3,	3,	1,	3,	2,	3,	2,	1,	2,	3,	1,	2,	1,	2,	3,	3,	3,	3,	3,	2,	1,	2,	2,	3,	3,	2,	2,	3,	4,	3,	2,	1,	1,	3,	2,	2,	2,	4,	3,	2,	3,	3,	3,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0
};
const unsigned int MMFF94_VAL[]=     {  4 ,	4 ,	4 ,	4 ,	1 ,	2 ,	2 ,	3 ,	3 ,	3 ,	1 ,	1 ,	1 ,	1 ,	2 ,	2 ,	4 ,	4 ,	4 ,	4 ,	1 ,	4 ,	1 ,	1 ,	4 ,	3 ,	1 ,	1 ,	1 ,	4 ,	1 ,	12,	1 ,	4 ,	1 ,	1 ,	4 ,	3 ,	3 ,	3 ,	4 ,	3 ,	3 ,	2 ,	4 ,	3 ,	2 ,	2 ,	3 ,	1 ,	3 ,	1 ,	4 ,	4 ,	34,	34,	4 ,	4 ,	2 ,	3 ,	4 ,	2 ,	4 ,	4 ,	3 ,	3 ,	4 ,	4 ,	4 ,	2 ,	1 ,	1 ,	3 ,	4 ,	3 ,	2 ,	4 ,	4 ,	3 ,	4 ,	4 ,	4 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 ,	0 
};
const unsigned int MMFF94_PILP[]=    {  0,	0,	0,	0,	0,	1,	0,	1,	0,	1,	1,	1,	1,	1,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	1,	0,	0,	1,	0,	0,	0,	1,	1,	0,	0,	1,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	0,	0,	1,	0,	0,	0,	0,	0,	0,	0,	1,	0,	1,	0,	0,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0
};
const unsigned int MMFF94_MLTB[]=    {	0,	2,	2,	3,	0,	0,	2,	0,	2,	1,	0,	0,	0,	0,	0,	2,	2,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	2,	0,	1,	0,	0,	1,	0,	2,	2,	1,	0,	1,	3,	0,	1,	2,	2,	2,	0,	0,	0,	2,	0,	2,	2,	1,	1,	2,	1,	1,	3,	3,	0,	2,	2,	2,	2,	2,	0,	1,	0,	0,	1,	0,	2,	2,	0,	0,	2,	2,	2,	1,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0
};
const unsigned int MMFF94_AROM[]=    {	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	1,	1,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	1,	0,	0,	0,	1,	1,	1,	1,	0,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	1,	1,	0,	1,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0
};
const unsigned int MMFF94_LIN[]=     {  0,	0,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0
};
const unsigned int MMFF94_SBMB[]=    {	0,	1,	1,	1,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	0,	0,	1,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	0,	0,	1,	1,	0,	0,	0,	0,	1,	1,	0,	0,	1,	0,	0,	0,	0,	0,	0,	0,	1,	0,	0,	1,	0,	1,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0
};


class FFComponent{
      friend class ForceField;
      public:
      public:

             // constructors
             FFComponent();
             FFComponent(ForceField& ff);
             // copy constructor;
             FFComponent(const FFComponent& ffc);
             // destructor
             virtual ~FFComponent();
             
             //setup
             virtual bool setup();
             
             // return the FF this component registered in
             ForceField* get_force_field() const;

             // set the force field
             void set_force_field(ForceField& ff);
             
             // set the component name
             void set_name(const string& name);

             // get component name
             string get_name() const;
             
             inline bool isenabled() const {return enabled_;}
             inline void setenabled(bool flag){ enabled_ = flag;}
             
             // FF calculation. these mothods are inherited by specific FF component;
             virtual double get_energy() const;
             virtual double update_energy();
             virtual void update_forces();
             virtual void update();
			 // set the number of atom types
			 inline void set_num_atom_types(int num)
			 {
				 num_of_atom_types_ = num;
				 return;
			 }
			 // return the number of atom types defined in the FF
			 inline int get_num_atom_types()
			 {
				 return num_of_atom_types_;
			 }
      protected:
                // the force field this component registered in
                ForceField* force_field_;
                // the energy of the component
                double energy_;
				// number of atom types from the force field
				int num_of_atom_types_;
      private:
              std::string name_;
              bool enabled_;

};
#endif              
