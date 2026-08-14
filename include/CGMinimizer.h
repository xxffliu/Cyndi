
#ifndef CGMinimizer_H 
#define CGMinimizer_H 

#ifndef	ENERGYMINIMIZER_H
#	include "EnergyMinimizer.h"
#endif

#ifndef CONICLINESEARCH_H
# include "ConicLineSearch.h"
#endif


/*#ifndef LINESEARCH_H
# include "LineSearch.h"
#endif*/

	/** A minimizer for geometry optimization based on different
	 *  conjugate gradient (CG) algorithms, see
	 *
	 *  [1] Fletcher, R. and Reeves, C.M.: Function Minimization by 
	 *  Conjugate Gradients, Comp. J., 7, 149-154, 1964
	 *
	 *  [2] David Shanno "Conjugate Gradient Methods With Inexact Searches"
	 *  Mathematics of Operations Research, Vol. 3, No. 3, August 1978, pp. 244-256
	 *
	 *  [3] S. Watowich et. al "A Stable Rapidly Converging Conjugate Gradient Method
	 *  for Energy Minimization", Journal of Computational Chemistry, Vol. 9,
	 *  No. 6, pp. 650-661 (1988)
	 *
	 *  [4] M.J.D. Powell "Convergence properties of algorithms for nonlinear
	 *  optimization, Report No. DAMPT 1985/NA1, University of Cambridge, 1985
	 *
	 *  [5] D:F: Shanno "Globally convergent conjugate gradient algorithms",
	 *  Mathematical Programming 33 (1985), pp. 61-67
	 *
	 *  \ingroup  MolmecEnergyMinimizer
	 */
	class CGMinimizer 
		: public EnergyMinimizer
	{
		public:
			
			/** @name Enums
			*/
			//@{
			
			/** The different conjugate gradient methods implemented.
			 *  @see updateDirection for details on the implementation and references
			 */
			enum UpdateMethod
			{
				/** Polak-Ribiere method
				*/
				POLAK_RIBIERE = 1,
				
				/** Fletcher-Reeves method
				*/
				FLETCHER_REEVES = 2,
				
				/** Shanno
				*/
				SHANNO = 3
			};
			
			//@}
			/** @name Options and Defaults
			 */
			//@{
			
			/** Options names
			*/
			/*struct MIN_OPTION
			{
				///The initial step length used in the line search 
				int UPDATE_METHOD; 
			};
			
			/**	Defaults for all options
			struct DEFAULT_OPTION
			{
				/** The initial step length used in the line search 
				static const int UPDATE_METHOD;
			};
			*/
			//@}
			/** @name Constructors and Destructors	
			*/
			//@{
			
			
			/** Default constructor.
			*/
			CGMinimizer();
			
			/** Constructor expecting a valid force field
			*/
			CGMinimizer(ForceField& force_field);
			
			
			/** Constructor expecting a valid force field and options 
			*/
			CGMinimizer(ForceField& force_field, const MIN_OPTION& options);
			
			/** Copy constructor
			*/
			CGMinimizer(const CGMinimizer& rhs);
			
			//assignment operator
			const CGMinimizer& operator=(const CGMinimizer& rhs);
			
			/** Destructor.
			*/
			virtual ~CGMinimizer();
			
			//@}
			/** @name Assignments
			*/
			//@{
			
			/** Assignment operator
			*/
			//const CGMinimizer& operator = (const CGMinimizer& rhs);
			
			//@}
			/** @name Setup methods. They do all necessary preparations.
			*/
			//@{
			
			/** Specific setup
			*/
			virtual bool specific_setup();
			
			//@}
			/** @name Accessors
			*/
			//@{
			
			/** Set explicitly the criterion used for updateDirection. If the
			 *  method is not recognized, this function fails.
			 */
			void set_update_method(UpdateMethod method);
			
			/** Returns the current method for updateDirection.
			 */
			UpdateMethod get_update_method() const;
			
			/** Calculate the next step.
			 *  First this method updates the model and performs a line search
			 *  into the calculated direction afterwards.
			 *  @return double <b>\geq 0.0</b> if the line search found an acceptable solution, otherwise -1.0.
			 *  @see EnergyMinimizer::findStep
			 */
			virtual double find_step();
			
			/** Update the search direction.
			 *  This method updates the search direction.
			 *  It uses the different conjugate gradient caculations
			 *  dependend on the options.
			 */
			virtual void update_direction();
			
			/** Minimize the energy of the system.
			 *  This method executes at most <tt>iterations</tt> minimization steps.
			 *  If the number of iterations is not given, the number specified in the
			 *  options is taken.
			 *  @param iterations the maximum number of iterations
			 *  @param resume <b>true</b> to resume a previous run
			 *  @see EnergyMinimizer::minimize
			 */
			virtual bool minimize(int iterations = 0, bool resume = false); 
			
		protected:
		
			//@}
			/** @name Protected Attributes
			*/
			//@{
			
			/*_ The line search
			*/
			ConicLineSearch line_search_;
			//LineSearch line_search_;
			/*_ The last step size (in respect to the length of the computed direction vector),
				  so the length of the last step was step_*||direction_||.
			*/
			double step_;
			
			/*_ The unscaled last search direction
			*/
			Gradient unscaled_direction_;
			
			/*_ The number of movable atoms.
			*/
			int number_of_atoms_;
			
			/*_ The update method used for the CG
			*/
			int update_method_;
			
			/*_ Is this the (initial) first iteration?
			*/
			bool first_iter_;
			
			/*_ g^t*g where g is the gradient of the last iteration
			*/
			double old_gtg_;
			
			/*_ Some variables that are needed for the Shanno direction calculation
			*/
			vector<vector3> a_i_, b_i_, p_t_, y_t_, p_i_, y_i_;
			double D_1_, D_4_;
			
			/*_ Frequency for restarts.
			*/
			int restart_frequency_;
			
			/*_ We count the iterations since the last restart.
				  There are other possibilities for a restart so the
				  iteration counter may not coincide (by modulo) with the
				  restart frequency.
			*/
			int last_restart_iter_;
			
			//@}
		
	};

#endif 
