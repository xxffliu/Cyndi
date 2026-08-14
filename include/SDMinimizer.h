// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: steepestDescent.h,v 1.24.20.2 2007-05-16 15:56:13 aleru Exp $
// Line Search Minimizer: A special class for the line search minimization algorithm

#ifndef SDMINIMIZER_H
#define SDMINIMIZER_H

#ifndef ENERGYMINIMIZER_H
# include "EnergyMinimizer.h"
#endif

#ifndef CONICLINESEARCH_H
# include "ConicLineSearch.h"
#endif

//#ifndef BALL_MOLMEC_MINIMIZATION_LINESEARCH_H
//# include <BALL/MOLMEC/MINIMIZATION/lineSearch.h>
//#endif


	/** A minimizer for geometry optimization based on steepest descent steps.
	 *  \ingroup  MolmecEnergyMinimizer
	 */
	class  SDMinimizer 
		: public EnergyMinimizer
	{
		
		public:
			
			/** @name Constructors and Destructors
			*/
			//@{
			
			
			/** Default constructor.
			*/
			SDMinimizer();
			
			/** Constructor.
			*/
			SDMinimizer(ForceField& force_field);
			
			/** Constructor.
			*/
			SDMinimizer(ForceField& force_field, const MIN_OPTION& options);
			
			/** Copy constructor
			*/
			SDMinimizer(const SDMinimizer& minimizer);
			
			/** Destructor.
			*/
			virtual ~SDMinimizer();
			
			//@}
			/** @name Assignments
			*/
			//@{
			
			/** Assignment operator
			*/
			const SDMinimizer&	operator = (const SDMinimizer& minimizer);
			
			//@}
			/** @name Setup methods
			*/
			//@{
			
			/** Specific setup
			*/
			virtual bool specific_setup();
			
			//@}
			/** @name Accessors
			 */
			//@{
			
			/** Minimize the energy of the system using scaled steepest descent steps.
			*/
			virtual bool minimize(int steps = 0, bool resume = false);
			
			/** Find the next step using a line search.
			*/
			virtual double find_step();
			
			/** Update the search direction.
			 *  (Scaled) steepest descent searches along the current gradient only.
			 *  Therefore, updateDirection only assigns direction to the (scaled) last gradient
			 *  computed (current_gradient_).
			 */
			virtual void update_direction();
			
		protected:
			
			//@}
			/** @name Protected Attributes
			 */
			//@{
			
			/*_ The line search minimizer.
					This member is used to perform the line search in findStep
			*/
			ConicLineSearch line_search_;
			//LineSearch line_search_;
			//@}
			
	};


#endif // BALL_MOLMEC_MINIMIZATION_STEEPESTDESCENT_H
