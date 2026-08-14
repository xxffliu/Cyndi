#ifndef ENERGY_MINIMIZER_H
#define ENERGY_MINIMIZER_H

#ifndef GRADIENT_H
#include "Gradient.h"
#endif

#ifndef VECTOR3_H
#include "Vector3.h"
#endif
#include "utility.h"
class ForceField;

class EnergyMinimizer{
      
      public:
             struct MIN_OPTION{
				int MAXIMAL_NUMBER_OF_ITERATIONS; //Max number of iterations
				int ENERGY_OUTPUT_FREQUENCY; // Energy ouput frequency
				int NUMBER_OF_ITERATION; // Number of iteration
				float ENERGY_DIFFERENCE_BOUND; // The energy difference needed for assuming 'equal energy'
				int MAX_SAME_ENERGY;  // The number of iterations without any change in energy.
				float MAX_GRADIENT; // Default value for the maximum RMS gradient.
				float MAXIMUM_DISPLACEMENT; // Max shift of single atom.
				int UPDATE_METHOD; // Specify update method using in CG minimizer (default Shannon method)
				MIN_OPTION();
            };
            // default constructor
            EnergyMinimizer();
            // constructor
            EnergyMinimizer(ForceField& ff);
            // constructor
            EnergyMinimizer(ForceField& ff, const MIN_OPTION & opt);
            // copy constructor
            EnergyMinimizer(const EnergyMinimizer& minimizer);
            // assignment operator
            const EnergyMinimizer& operator=(const EnergyMinimizer& rhs);
            // destructor
            virtual ~EnergyMinimizer();
            
            // check the minimizer is valid
            bool isValid() const;
            
            //set up methods
            bool setup(ForceField& ff);
            bool setup(ForceField& ff, MIN_OPTION option);
            //specific setup
            virtual bool specific_setup();
            
            // implements the convergence criterion. if any of following 2 criterions is fullfilled,
            // the method return "true":
            //  1) RMS gradient is below max_rms_gradient_;
            //  2) same_energy_counter_ is above max_same_energy_;
            // this method is inherited by other minimizer classes;
            virtual bool is_converged() const;
            
            // Calculate the next step;
            // Inherited by other minizer classes and tries to determine the next step to be taken;
            // It performs a line search and returns the step length with respect to the current direction.
            // If failure raised returns -1;
            virtual double find_step();
            
            // inherited by other minimizer classes to update the search direction
            virtual void update_direction();
            
            // update current energy;
            // this method calls ForceField->update_energy() and stores the result in current_energy_;
            virtual double update_energy();
            // update current force;
            // this method calls ForceField->update_force();
            virtual void update_force();
            
            // Store the current energy and gradient. The current gradient and current energy is copied into
            // initial energy and initial gradient. This is usually done at the start of an iteration.
			void store_gradient_energy();
			
			// Print the energy.
		    // This method is called by finishIteration after every energy_output_frequency_ steps.
			//  It prints the current RMS gradient and the current energy to std::ostream or log.
			virtual void print_energy() const;
			
			// Finishing step for this iteration.
			//  This method should be called at the end of the main iteration
			//  loop implemented in minimize . It takes over some administrative stuff:
			//  - increment the iteration counter number_of_iterations_
			//  - call printEnergy if necessary
			//  - call ForceField::update ForceField::update if necessary (to rebuild the pair lists!)
			//  - update the same_energy_counter_ tested in isConverged
			//  This method should be overwritten only in rare cases. Even then, the programmer
			//  should make sure to call EnergyMinimizer::finishIteration or
			//  has to take care of the above items himself.
			//  All derived classes should call this method at the end of the minimize main loop.
			//  Otherwise strange things might happen.

			virtual void finish_iteration();
			/** Return the number of iterations performed.
			 */
			int get_number_of_iterations() const;
			
			/** Return a reference to the current search direction
			 */
			Gradient& get_direction();
			
			/** Return a reference to the current gradient
			 */
			Gradient& get_gradient();
			
			/** Return a reference to the initial gradient
			 */
			Gradient& get_initial_gradient();
			
			/** Return the current energy
			 */
			double get_energy() const;
			
			/** Return a reference to the current energy
			 */
			double& get_energy();
			
			/** Return the initial energy
			 */
			double get_initial_energy() const;
			
			/** Return a mutable reference to the initial energy
			 */
			double& get_initial_energy();
			
			/** Set the number of iterations performed so far.
			 */
			void set_number_of_iterations(int number_of_iterations);
			
			/** Get the maximum number of iterations
			 */
			int get_max_number_of_iterations() const;
			
			/** Set the maximum number of iterations
			 */
			void set_max_number_of_iterations(int number_of_iterations);
			
			/** Set the maximum number of iterations allowed with equal energy
			 *  (second convergence criterion)
			 */
			void  set_max_same_energy(int number);
			
			/** Get the maximum number of iterations allowed with equal energy
			 *  (second convergence criterion)
			 */
			int get_max_same_energy() const;
			
			/** Set the energy output frequency
			 */
			void set_energy_output_frequency(int energy_output_frequency);
			
			/** Get the energy ouput frequency
			 */
			int get_energy_output_frequency() const;
			/** Set the energy difference bound for convergence
			*/
			void set_energy_difference_bound(float energy_difference_bound);
			
			/** Get the energy difference bound
			 */
			float get_energy_difference_bound() const;
			
			/** Set the maximum RMS gradient (first convergence criterion).
			 *  The gradient unit of the gradient is <b>kJ/(mol \AA)</b>.
			 */
			void set_max_gradient(float max_gradient);
			
			/** Get the maximum RMS gradient (first convergence criterion).
			 *  The gradient unit of the gradient is <b>kJ/(mol \AA)</b>.
			 */
			float get_max_gradient() const;
			
			/** Set the maximum displacement value.
			 *  This is the maximum distance an atom may be moved by the minimizer in one iteration.
			 */
			void  set_maximum_displacement(float maximum_displacement);
			
			/** Get the maximum displacement value
			*/
			float get_maximum_displacement() const;
			/** Return the force field of the energy minimizer
			*/
			ForceField* get_force_field();
			
			/** Return the number of force updates since the start of the minimization.
			*/
			int get_force_update_counter() const;
			
			/** Return the number of energy updates since the start of the minimization.
			*/
			int get_energy_update_counter() const;
			/** Minimize the energy of the system bound to the force field.
			 *  If a number of steps is given, the minimization is aborted after
			 *  that number of steps, regardless of the number of steps given in
			 *  the options (<tt>MAX_STEPS</tt>). Together with the <tt>resume</tt> option
			 *  this feature is used to extract properties or visualize the results
			 *  in the course of the minimization. If <tt>resume</tt> is set to <b>true</b>,
			 *  the minimization continues with the former step width and settings.
			 */
			virtual bool minimize(int steps = 0, bool resume = false);
			
			//toggle  the state of abort condition
			void enable_energy_abort_condition(bool state);
			
			// check if abort condition is enabled
			bool energy_abort_condition_enabled() const;
			
			// set the energy criteria to abort
			void set_energy_to_abort(float value);
			
			// return the energy criteria to abort
			float get_energy_to_abort() const;
			
			/** Return true, if the minimization was aborted, e.g. because of strange
			 *  energies or gradient.
			 */
			bool was_aborted() const;
			
			//Store the current atom positions in saved_position_
            //considers these coordinates as start coordinates.
            virtual void save_position();
            
            // Resets the atom positions to the saved positions.
			// If coordinates weres stored using  savePositions() , the atoms
			// coordinates are reset to the saved positions.
			// If no savedPositions exist the coordinates remain unchanged.
		    virtual void reset_position();
		    
		    //  Move all atoms along a direction vector.
			//	The method translates all atoms a long a given direction.
			//	The direction vector is multiplied with a step length <tt>step</tt>.
			//	If a saved position exists ( savePositions() ), it is used as a start
			//	position (i.e. the final positions are vec{start} + {step} \ vec{direction}).
			//	Otherwise, the current atom positions are used.
			//	If the gradient's size differs from the number of atoms, nothing is done.
		    virtual void move_to(const Gradient& direction, float step = 1.0);
			
			// public attributes
            //const MIN_OPTION default_option;
			MIN_OPTION option;
			
      protected:
             /** The gradient at the beginning of the current minimization step.
			*/
			Gradient initial_grad_;
			
			/** The current gradient.
			*/
			Gradient current_grad_;
			
			/** The energy at the beginning of the current minimization step.
			*/
			double initial_energy_;
			
			/** The current energy.
			*/
			double current_energy_;
			
			/** The gradient from the last step
			*/
			Gradient old_grad_;
			
			/** The energy from the last step
			*/
			double old_energy_;
			
			/** The current search direction
			*/
			Gradient direction_;
		 
			/** The boolean variable indicates if the setup of the energy minimizer was successful
			*/
			bool valid_;
			 
			
			/** The force field bound to the energy minimizer.
			 *  Among other data the force field contains the molecular system
			 *  whose energy will be minimized by the energy minimizer.
			 */
			ForceField* force_field_;
			
			/** The current iteration number
			*/
			int number_of_iterations_;
			
			/** Maximum number of iterations 
			*/
			int maximal_number_of_iterations_;
			
			/** Frequency of energy output
			*/
			int energy_output_frequency_;
			
			
			/** If the energy difference (before and after an iteration)
			 * is smaller than this bound, the minimization procedure stops.
			 */
			double energy_difference_bound_;
			
			/** The maximum RMS gradient tolerated (first convergence criterion)
			*/
			double max_gradient_;
			
			/** The maximum number of iterations with same energy.
			 *  When this number is reached, we assume the system to have converged
			 *  (second convergence criterion)
			 */
			int max_same_energy_;
			
			/** A counter for the number of steps with a similar energy.
			*/
			int same_energy_counter_;
			
			/** The maximal shift of an atom per iteration step (in Angstrom).
			*/
			float maximum_displacement_;
			
			/** Internal counter: how often is a force update done.
			 *  Measure for the speed of minimization.
			 */
			int force_update_counter_;
			
			/** Internal counter: how often is an energy update done.
			 *  Measure for the speed of minimization.
			 */
			int energy_update_counter_;
			
			/** Numerical lower bound: we don't want to compute the reciprocal 
			 *  of a number which is lower than 'cutlo_'.
			 */
			float cutlo_;
			
			/** The last step size (in respect to the length of the computed direction vector),
			 *  so the length of the last step was step_*||direction_||.
			 */
			double step_;
			
			//_ 
			bool abort_by_energy_enabled_;
			
			//_ 
			float abort_energy_;
			
			//_
			bool aborted_;
			
			// a temp vector holding for position of each atom due to changes by minimization
			vector<vector3> saved_position_;
};

#endif
			 
			
              
            
            
             
