
#include "../include/SDMinimizer.h"
#include "../include/ForceField.h"
#include <math.h>
#include <limits>
//#define BALL_DEBUG
#undef BALL_DEBUG

using namespace std;


	
	// Default constructor
	 SDMinimizer::SDMinimizer()
		:	EnergyMinimizer(),
			line_search_(*this)
	{
	}
	
	// Copy constructor 
	SDMinimizer::SDMinimizer 
		(const SDMinimizer & minimizer)
		: EnergyMinimizer(minimizer),
			line_search_(minimizer.line_search_)
	{
		line_search_.set_minimizer(*this);
	}
	
	// Constructor initialized with a force field
	SDMinimizer::SDMinimizer(ForceField& force_field)
		:	EnergyMinimizer(),
			line_search_(*this)
	{
		valid_ = setup(force_field);
		
		if (!valid_)
		{
			cout << "SDMinimizer: setup failed! " << std::endl;
		}
	}

	// Constructor initialized with a force field and a set of options
	SDMinimizer::SDMinimizer
		(ForceField& force_field, const MIN_OPTION& new_options)
		:	EnergyMinimizer(),
			line_search_(*this)
	{
		option = new_options;
		valid_  = setup (force_field, new_options);
		
		if (!valid_)
		{
			cout << "SDMinimizer: setup failed! " << endl;
		}
	}

	
	// Destructor
	SDMinimizer::~SDMinimizer()
	{
	}
	
	// Assignment operator
	const SDMinimizer& SDMinimizer::operator=
			(const SDMinimizer& minimizer)
	{
		if (&minimizer != this)
		{
			EnergyMinimizer::operator=(minimizer);
			line_search_      = minimizer.line_search_;
			line_search_.set_minimizer(*this);
		}
		return *this;
	}
	
	// virtual function for the specific setup of derived classes
	bool SDMinimizer::specific_setup()
	{
		// Make sure the force field is assigned and valid!
		if (force_field_ == 0 || !force_field_->isValid())
		{
			return false;
		}
		
		// Invalidate the initial gradient in order to ensure
		// its re-evaluation at the start of minimize().
		initial_grad_.invalidate();
		
		return true;
	}
	
	/*  The minimizer optimizes the energy of the system
	*/
	bool SDMinimizer::minimize(int iterations, bool resume)
	{
		aborted_ = false;
		
		// Check for validity of minimizer and force field
		if (!isValid() || get_force_field() == 0 || !get_force_field()->isValid())
		{
			cout << "SDMinimizer: minimizer is not initialized correctly!" << std::endl;
			aborted_ = true;
			return false;
		}
		
		// Make sure we have something worth moving.
		if (get_force_field()->get_num_atoms() == 0)
		{
			return true;
		}
		
		// Some aliases
		ATOMVec& atoms(get_force_field()->get_atoms());
		
		// If the run is to be continued, don't reset the iteration counter and the initial step size
		if (!resume)
		{
			// reset the number of iterations for a restart
			set_number_of_iterations(0);
			same_energy_counter_ = 0;
			initial_grad_.invalidate();
			current_grad_.invalidate();
			
			// Obviously, we don't have "old" energies yet, so we initialize it a with 
			// sensible value. We don't need "old" gradients.
			old_energy_ = (std::numeric_limits<float>::max)();
		}
		int max_iterations = std::min(get_number_of_iterations() + iterations, get_max_number_of_iterations());
		
		// Save the current atom positions
		save_position();
		bool converged = false;	
		while (!converged && (get_number_of_iterations() < max_iterations))
		{
			// Try to take a new step
			double stp = find_step();
			
			// Check whether we were successful.
			if (stp > 0.0)
			{
				// Use this step as new reference step if findStep was successful
				save_position();
			}
			
			// Store the energy, there's no need to store the old gradient
			old_energy_ = initial_energy_;
			
			// Store the current gradient and energy
			store_gradient_energy();
			
			#ifdef BALL_DEBUG
				cout << "SDM::minimize: end of main: current grad RMS = " << current_grad_.rms << std::endl;
			#endif
			
			// Check for convergence.
			converged = is_converged();
			
			// Increment iteration counter, take snapshots, print energy,
			// update pair lists, and check the same-energy counter
			finish_iteration();
			
			if ((!converged) && (stp < 0.))
			{
				// Nasty case: No convergence and the step computation failed.
				// We must give up:-(
				aborted_ = true;
				return false;
			}
			/*if (_isnanf(force_field_->get_energy()))
			{
				aborted_ = true;
				return false;
			}
			
			if (isnanf(get_gradient().rbegin()->x()) ||
		        isnanf(get_gradient().rbegin()->y()) ||
			    isnanf(get_gradient().rbegin()->z())) 
			{
				aborted_ = true;
				return false;
			}*/
			
			if (abort_by_energy_enabled_)
			{
				if (force_field_->get_energy() > abort_energy_) 
				{
					aborted_ = true;
					return false;
				}
			}
		}
		
		return converged;
	}
	
	double SDMinimizer::find_step()
	{
		// Compute the new direction
		update_direction();
		
		bool result = false;
		int iter = 0;
		while ((!result) && (iter < 10))
		{
			double step;
			
			// No need to assure the maximum displacement here since our
			// line search pays attention to this constraint.
			
			result = line_search_.minimize(step);
			
			if (!result)
			{
				// Some aliases
				//ATOMVec& atoms(get_force_field()->get_atoms());
				int n = get_force_field()->get_num_atoms();
				for(int i = 0; i < n; ++i)
				{
					direction_[i] *= 0.5;
				}
				direction_.norm     *= 0.5;
				direction_.rms      *= 0.5;
				direction_.inv_norm *= 2.;
				reset_position();
			}
			else
			{
				return step;
			}
			++iter;
		}
		
		// If we are here something went wrong
		return -1.0;
	}
	
	void SDMinimizer::update_direction()
	{
		// If we do not have a valid gradient, recalculate the gradient, the energy,
		// and the search direction
		if (!initial_grad_.isValid())
		{
			// Compute the initial energy and the initial forces
			update_energy();
			update_force();
			
			// Store them for later use
			store_gradient_energy();
		}
		
		// The direction is the normalized negative gradient
		direction_ = initial_grad_;
		direction_.minus();
		direction_.normalize();
	}
	

