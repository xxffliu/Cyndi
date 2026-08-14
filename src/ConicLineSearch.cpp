
#include "../include/ConicLineSearch.h"
#include "../include/EnergyMinimizer.h"
//#include <BALL/MOLMEC/COMMON/atomVector.h>
#include "../include/ForceField.h"

// Parameter alpha for 'sufficient energy decrease'.
#define CONICLINESEARCH_DEFAULT_ALPHA 1.e-4

// Parameter beta for 'sufficient gradient reduction'.
#define CONICLINESEARCH_DEFAULT_BETA 0.9

// Lower bound for energy values.
#define CONICLINESEARCH_DEFAULT_MIN_ENERGY -1.e+10

// Maximum number of interpolation steps for a line search.
#define CONICLINESEARCH_DEFAULT_MAX_STEPS 50



using namespace std; 
	// Default constructor
	ConicLineSearch::ConicLineSearch()
		:	alpha_(CONICLINESEARCH_DEFAULT_ALPHA),
			beta_(CONICLINESEARCH_DEFAULT_BETA),
			max_steps_(CONICLINESEARCH_DEFAULT_MAX_STEPS),
			lower_energy_bound_(CONICLINESEARCH_DEFAULT_MIN_ENERGY),
			is_bracketed_(false),
			ctrl_stp_(0.),
			ctrl_val_(1.e+100),
			minimizer_(0)
	{
	}
	
	
	// Copy constructor 
	ConicLineSearch::ConicLineSearch(const ConicLineSearch& line_search)
		:	alpha_(line_search.alpha_),
			beta_(line_search.beta_),
			max_steps_(line_search.max_steps_),
			lower_energy_bound_(line_search.lower_energy_bound_),
			is_bracketed_(line_search.is_bracketed_),
			ctrl_stp_(line_search.ctrl_stp_),
			ctrl_val_(line_search.ctrl_val_),
			minimizer_(line_search.minimizer_)
	{
	}
	
	// Assignment operator
	/*const ConicLineSearch& ConicLineSearch::operator = (const ConicLineSearch& line_search)
	{
		alpha_              = line_search.alpha_;
		beta_               = line_search.beta_;
		max_steps_          = line_search.max_steps_;
		lower_energy_bound_ = line_search.lower_energy_bound_;
		is_bracketed_       = line_search.is_bracketed_;
		ctrl_stp_           = line_search.ctrl_stp_;
		ctrl_val_           = line_search.ctrl_val_;
		minimizer_          = line_search.minimizer_;
		
		return *this;
	}*/
	
	// Detailed constructor
	ConicLineSearch::ConicLineSearch(EnergyMinimizer& minimizer)
		:	alpha_(CONICLINESEARCH_DEFAULT_ALPHA),
			beta_(CONICLINESEARCH_DEFAULT_BETA),
			max_steps_(CONICLINESEARCH_DEFAULT_MAX_STEPS),
			lower_energy_bound_(CONICLINESEARCH_DEFAULT_MIN_ENERGY),
			is_bracketed_(false),
			ctrl_stp_(0.),
			ctrl_val_(1.e+100),
			minimizer_(&minimizer)
	{
	}
	
	// Destructor
	ConicLineSearch::~ConicLineSearch()
		throw()
	{
	}
	
	// Set the parameter alpha_.
	void ConicLineSearch::set_alpha(double alpha)
	{
		alpha_ = alpha;
	}
	
	// Get the parameter alpha_.
	double ConicLineSearch::get_alpha() const
	{
		return alpha_;
	}
	
	// Set the parameter beta_.
	void ConicLineSearch::set_beta(double beta)
	{
		beta_ = beta;
	}
	
	// Get the parameter beta_.
	double ConicLineSearch::get_beta() const
	{
		return beta_;
	}
	
	// Set the parameter max_steps_.
	void ConicLineSearch::set_max_steps(int max_steps)
	{
		max_steps_	= max_steps;
	}
	
	// Get the parameter max_steps_.
	int ConicLineSearch::get_max_steps() const
	{
		return max_steps_;
	}
	
	// Set the lower bound on energy values.
	void ConicLineSearch::set_lower_bound(double lbound)
	{
		lower_energy_bound_ = lbound;
	}
	
	// Get the lower bound on energy values.
	double ConicLineSearch::get_lower_bound() const
	{
		return lower_energy_bound_;
	}
	
	// Set the flag is_bracketed_.
	void ConicLineSearch::set_bracketed_flag(bool bracktd)
	{
		is_bracketed_ = bracktd;
	}
	
	// Return whether a minimizer has already been bracketed.
	bool ConicLineSearch::is_bracketed() const
	{
		return is_bracketed_;
	}
	
	// Set the calling minimizer class.
	void ConicLineSearch::set_minimizer(EnergyMinimizer& minimizer)
	{
		minimizer_ = &minimizer;
	}
	
	/*	The minimizer optimizes the energy of the system using a two stage line 
			search algorithm similar to that proposd by More and Thuente.
	*/
	bool ConicLineSearch::minimize(double& stp)
	{
		//debug
			//cout << "oops! LS:minimize(" << stp << ")" << std::endl;

		// Check whether a valid minimizer and a valid force field exist.
		if ((minimizer_ == 0) || (minimizer_->get_force_field() == 0))
		{
			return false;
		}
		
		// Define some aliases for convenience.
		//AtomVector&				atoms(const_cast<AtomVector&>(minimizer_->getForceField()->getAtoms()));
		const Gradient&		direction(minimizer_->get_direction());
		//EnergyMinimizer&	minimizer(*minimizer_);
		Gradient&					gradient(minimizer_->get_gradient());
		Gradient&					initial_gradient(minimizer_->get_initial_gradient());
		
		// If we do not have a valid gradient for the first step,
		// calculate it.
		if (!minimizer_->isValid())
		{   
            // debug
            //cout<<"The minimizer is not valid"<<endl;
			// Reset the atoms to the start position (stp = 0)
			minimizer_->reset_position();
			
			// Calculate the initial energy and forces
			minimizer_->update_force();
			minimizer_->update_energy();
			initial_gradient = gradient;
			
			// Force a recalculation of the current gradient
			// as well since updateForces overwrote everything!
			gradient.invalidate();
		}
		
		// Initial energy value.
		double f_init = minimizer_->get_initial_energy();
		
		// Initial directional derivative.
		double g_init = (initial_gradient * direction);
		
		// Minimum and maximum stepsizes.
		
		// Compute the maximum step size by the minimizers 'maximum displacement'
		maxstp_ = minimizer_->get_maximum_displacement();
		double minstp = 1.e-10;
		if (maxstp_ < 0.)
		{   //debug
		    //cout<<"No maximum displacement given, estimate the maximum stepsize"<<endl;
			// No maximum displacement given, estimate the maximum stepsize
			double tmp = gradient.norm/direction.norm;
			maxstp_ = 1.e+5*tmp;
			
			// Useful safeguard for minimum stepsize
			minstp = 1.e-10*tmp;
			//debug
			//cout<<"minimum stepsize is: "<<minstp<<endl;
		}
		else
		{   
            //debug
            //cout<<"Find the maximum translation"<<endl;
			// Find the maximum translation
			vector<vector3>::const_iterator git(direction.begin());
			double max = 0.;
			double cur = 0.;
			for (; git != direction.end(); ++git)
			{
				cur = git->length_2();
				if (cur > max)
				{
					max = cur;
				}
			}
			max = sqrt(max);
			if (max > 1.e-16)
			{
				maxstp_ /= max;
			}
			else
			{   
                //debug
                //cout<<"Something went wrong, we estimate the maximum stepsize"<<endl;
				// Something went wrong, we estimate the maximum stepsize
				double tmp = gradient.norm/direction.norm;
				maxstp_ = 1.e+5*tmp;
				
				// Useful safeguard for minimum stepsize
				minstp = 1.e-10*tmp;
			}
		}
		
		// A minimum has not been bracketed so far.
		is_bracketed_ = false;
		
		// Compute initial stepsize
		//debug
		//cout<<"g_init is: "<<g_init<<endl;
		stp = min(1., (lower_energy_bound_ - f_init)*4./g_init);
		stp = max(stp, minstp);
		stp = min(stp, maxstp_);
		// Objective value at stp
		double f;
		
		// Directional derivative at stp
		double g = 0.;
		
		// Interval data
		double stp_up = 0.;
		double f_up = f_init;
		double g_up = g_init;
		double stp_lo = 0.;
		double f_lo = f_init;
		double g_lo = g_init;
		
		// Main loop

		for(int iteration = 0; iteration < max_steps_; ++iteration)
		{
			// Get energy and directional derivative at 'stp'
			
			// Move the atoms to the new position.
			minimizer_->move_to(direction, stp);
			// Update energy and gradient
			f = minimizer_->update_energy();
			minimizer_->update_force();
			g = (gradient * direction);
			
			// We check a few things
			//debug
			//cout<<f<<" "<<f_init<<" "<<alpha_<<" "<<g_init<<endl;
			bool con1  = f - f_init <= alpha_*stp*g_init;
			//cout<<g<<" "<<beta_<<" "<<g_init<<endl;
			bool con2  = g >= beta_*g_init;
			
			// Check whether we cannot proceed because the minimum
			// tolerable stepsize is achieved.
			if ((stp <= minstp) && !con1)
			{   //debug
			    //cout<<"minimum tolerable stepsize is achieved. "<<stp<<" "<<minstp<<" "<<con1<<endl;
				return false;
			}
			
			// Check whether we have the maximum stepsize achieved
			if ((stp >= maxstp_) && con1 && !con2)
			{   
                //debug
                //cout<<"maximum stepsize achieved"<<endl;
				return false;
			}
			
			// Convergence test (weak Wolfe conditions)
			if (con1 && con2)
			{   
				return true;
			}
			
			if (!is_bracketed_)
			{
				// Bracketing stage
				ctrl_stp_ = stp_lo;
				ctrl_val_ = f_lo;
				
				stp_lo = stp_up;
				f_lo   = f_up;
				g_lo   = g_up;
				stp_up = stp;
				f_up   = f;
				g_up   = g;
				
				if (!con1)
				{
					is_bracketed_ = true;
				}
			}
			else
			{
				// Interpolation stage
				if (!con1)
				{
					ctrl_stp_ = stp_up;
					ctrl_val_ = f_up;
					
					stp_up = stp;
					f_up   = f;
					g_up   = g;
				}
				else
				{
					ctrl_stp_ = stp_lo;
					ctrl_val_ = f_lo;
					
					stp_lo = stp;
					f_lo   = f;
					g_lo   = g;
				}
			}
			
			if (stp >= maxstp_)
			{
				// New stepsize
				stp = take_step(stp_lo, f_lo, g_lo, stp_up, f_up, g_up);
				if (stp >= maxstp_)
				{   
                    //debug
                    //cout<<"stp >= maxstp_"<<endl;
					return false;
				}
			}
			else
			{
				// New stepsize
				stp = take_step(stp_lo, f_lo, g_lo, stp_up, f_up, g_up);
			}
		}
		
		// If we're here something went wrong
		return false;
	}
	
	// Computes a step for a search procedure by case differentiation 
	// dependend on whether a minimum could already be bracketed or not.
	double ConicLineSearch::take_step(double stp_lo, double f_lo, double g_lo, 
																	 double stp_up, double f_up, double g_up) const
	{
		double stp;
	
		// Compute needed values
		double q = (f_up - f_lo) / (g_lo*(stp_up - stp_lo));
		double p = g_up / g_lo;
		
		double conic_stp = -1.;
		double cubic_stp = -1.;
		double quad_stp_val = -1.;
		
		// Compute the conic step
		double tmp = q*q - p;
		if (tmp >= 0.)
		{
			tmp = q + sqrt(tmp);
			if (tmp > 0.)
			{
				tmp = 1. - p/tmp/tmp/tmp;
				if (!compAndSafeguardStep_(stp_lo, stp_up, tmp, conic_stp))
				{
					conic_stp = -1.;
				}
			}
		}
		
		// Compute the cubic step
		tmp = p - q*3. + 2.;
		double tmp2 = tmp*tmp - (p - q*2. + 1.)*3.;
		if (tmp2 >= 0.)
		{
			tmp += sqrt(tmp2);
			if (!compAndSafeguardStep_(stp_lo, stp_up, tmp, cubic_stp))
			{
				cubic_stp = -1.;
			}
		}
		
		// Compute the quadratic step based on one directional derivative
		tmp = (1. - q)*2.;
		if (!compAndSafeguardStep_(stp_lo, stp_up, tmp, quad_stp_val))
		{
			quad_stp_val = -1.;
		}
		
		// Check in which stage we are
		if (!is_bracketed_)
		{
			// In the (extrapolation) bracketing stage
			// we choose the step which is closest to stp_up
			stp = conic_stp;
			if (fabs(cubic_stp - stp_up) < fabs(stp - stp_up))
			{
				stp = cubic_stp;
			}
			if (fabs(quad_stp_val - stp_up) < fabs(stp - stp_up))
			{
				stp = quad_stp_val;
			}
		}
		else
		{
			// We are in the interpolation stage
			
			double val_diff = 1.e+100;
			stp = conic_stp;
			
			// Check whether the conic step is possible
			if (conic_stp != -1.)
			{
				// Check the value of the conic interpolating function at the
				// control step
				double tmp = (q + sqrt(q*q-p))/p - 1.;
				double s   = stp_up-stp_lo;
				
				double d   = 2.*tmp/s;
				double x   = ctrl_stp_ - stp_lo;
				double den = d*x + 1.;
				den       *= den;
				
				if (den < 1.e-100)
				{
					// We cannot check the value since the control
					// step is too close to a singularity
					conic_stp = -1.;
					stp       = cubic_stp;
				}
				else
				{
					double a   = -(p - ((4.*tmp+2.)*q + 1.)*q)/p/s/g_lo;
					double b   = g_lo;
					double nom = a*x + b;
					nom       *= x;
					
					double val = nom/den + f_lo;
					
					val_diff   = fabs(val - ctrl_val_);
					stp        = conic_stp;
				}
			}
			
			// Check whether the cubic step is possible
			if (cubic_stp != -1.)
			{
				// Check the value of the cubic interpolating function
				// at the control step
				double s = stp_up - stp_lo;
				
				double a = -g_lo*(-1. + 2.*q - p)/s/s;
				double b = (-2. + 3.*q - p)/p;
				double c = g_lo;
				
				double x   = ctrl_stp_ - stp_lo;
				double val = a*x + b;
				val        = val*x + c;
				val        = val*x + f_lo;
				
				double tmp = fabs(val - ctrl_val_);
				if (tmp < val_diff)
				{
					val_diff = tmp;
					stp = cubic_stp;
				}
			}
			
			// Check whether the quadratic step is possible
			if (quad_stp_val != -1.)
			{
				double s   = stp_up - stp_lo;
			
				double a   = -g_lo*(1.-q)/s;
				double b   = g_lo;
				double x   = ctrl_stp_ - stp_lo;
				double val = a*x + b;
				val        = val*x + f_lo;
				
				double tmp = fabs(val - ctrl_val_);
				if (tmp < val_diff)
				{
					val_diff = tmp;
					stp = quad_stp_val;
				}
			}
		}
		
		// Check whether all interpolations have failed
		if (stp == -1.)
		{
			
			//We try a quadratic step with two directional derivatives
			double tmp = 1. - p;
			if (!compAndSafeguardStep_(stp_lo, stp_up, tmp, stp))
			{
				return stp;
			}
			
			// All other cases failed, we do a bisection
			if (!is_bracketed_)
			{
				stp = stp_up*2.;
			}
			else
			{
				stp = (stp_lo + stp_up)*0.5;
			}
			stp = min(stp, maxstp_);
		}
		return stp;
	}
	
	
	// Compute the safeguarded step and decide whether it can be accepted
	bool ConicLineSearch::compAndSafeguardStep_(double stp_lo, double stp_up, double den, double& stp) const
	{
		if (!is_bracketed_)
		{
			if ((den > 0.) && (den < 1.))
			{
				// We are in the extrapolation stage and accept 
				// this step
				stp = stp_lo + (stp_up - stp_lo)/den;
				stp = max(stp, stp_up*1.1);
				stp = min(stp, stp_up*2.);
				
				stp = min(stp, maxstp_);
				return true;
			}
		}
		else
		{
			if (den > 1.)
			{
				// We are in the interpolation stage and accept
				// this step
				stp = stp_lo + (stp_up - stp_lo)/den;
				stp = min(stp, stp_lo + (stp_up - stp_lo)*0.9999);
				return true;
			}
		}
		return false;
	}
	
	

