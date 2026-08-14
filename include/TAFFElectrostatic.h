#ifndef TAFFELE_H
#define TAFFELE_H

#ifndef FFCOMPONENT_H
#include "FFComponent.h"
#endif

#ifndef FORCEFIELD_H
#include "ForceField.h"
#endif

static double VACUUM_PERMITTIVITY = 8.85419e-12;
static double e0 = 1.60217738e-19;

class TAFFEle : public FFComponent
{
	public:
		// define a datastructure holding vdw parameters
        struct ElePair
		{
			ATOM* atom1;
            ATOM* atom2;
            bool is_14_interaction;
        };
		TAFFEle();
		TAFFEle(ForceField& ff);
		TAFFEle(const TAFFEle& to_copy);
		virtual ~TAFFEle();
		// set up method
		virtual bool setup();

		// access methds
		// update current electrostatic energy
		virtual double update_energy();
		// return current electrostatic energy
		double get_ele_energy()const; 
		// update forces imposed on each atoms by non-bond energy
		virtual void update_forces();
		//Update the pair list.This method is called by the force field whenever ForceField::update is called.
		// It is used to recalculate the nonbonded pair list.
		virtual void update();

	protected:
		vector<ElePair> ele_pair_list_;
		double ele_energy_;
};
#endif