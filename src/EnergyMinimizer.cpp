#include "../include/EnergyMinimizer.h"
#include "../include/ForceField.h"
#include <limits>
using namespace std;
// set the initial values for DEFAULT
/*int EnergyMinimizer::DEFAULT::MAXIMAL_NUMBER_OF_ITERATIONS = 1000;
int EnergyMinimizer::DEFAULT::ENERGY_OUTPUT_FREQUENCY = 50;
int EnergyMinimizer::DEFAULT::MAX_SAME_ENERGY = 20; 
int EnergyMinimizer::DEFAULT::NUMBER_OF_ITERATION = 0;              // start number 
float EnergyMinimizer::DEFAULT::ENERGY_DIFFERENCE_BOUND = 1e-2;      // in kJ/mol
float EnergyMinimizer::DEFAULT::MAX_GRADIENT = 0.01;                 // in kJ/(mol A) 
// maximum displacement == -1. means: let the line search the maximum stepsize estimate
float EnergyMinimizer::DEFAULT::MAXIMUM_DISPLACEMENT = -1.0;         // Angstrom
int EnergyMinimizer::DEFAULT::UPDATE_METHOD = 3;  //Shannon method;
*/
// set the initial values for OPTIONS (at test stage all are set the same to DEFAULT)
EnergyMinimizer::MIN_OPTION::MIN_OPTION():
          MAXIMAL_NUMBER_OF_ITERATIONS(1000),
          ENERGY_OUTPUT_FREQUENCY(50),
          MAX_SAME_ENERGY(20),
          NUMBER_OF_ITERATION(0),              // start number 
          ENERGY_DIFFERENCE_BOUND(1e-2),     // in kJ/mol
          MAX_GRADIENT(0.01),                 // in kJ/(mol A) 
          
// maximum displacement == -1. means: let the line search the maximum stepsize estimate
          MAXIMUM_DISPLACEMENT(-1.0),         // Angstrom
          UPDATE_METHOD(3){}

//constructors
EnergyMinimizer::EnergyMinimizer():
                                   initial_grad_(),
			                       current_grad_(),
			                       initial_energy_(0.0),
			                       current_energy_(0.0),
			                       old_grad_(),
			                       old_energy_(0.0),
			                       direction_(),
			                       valid_(false),
			                       force_field_(0),
			                       number_of_iterations_(0),
			                       maximal_number_of_iterations_(0),
			                       energy_output_frequency_(0),
			                       energy_difference_bound_(0.0),
			                       max_gradient_(0.0F),
			                       max_same_energy_(0),
			                       same_energy_counter_(0),
			                       maximum_displacement_(0.0F),
			                       force_update_counter_(0),
			                       energy_update_counter_(0),
			                       cutlo_(0.),
			                       step_(0.),
			                       abort_by_energy_enabled_(true),
			                       abort_energy_(1000000000.0),
			                       aborted_(false),
                                   saved_position_(){}

EnergyMinimizer::EnergyMinimizer(const EnergyMinimizer& to_copy):
			                           initial_grad_(to_copy.initial_grad_),
			                           current_grad_(to_copy.current_grad_),
			                           initial_energy_(to_copy.initial_energy_),
			                           current_energy_(to_copy.current_energy_),
			                           old_grad_(to_copy.old_grad_),
			                           old_energy_(to_copy.old_energy_),
			                           direction_(to_copy.direction_),
			                           valid_(to_copy.valid_),
			                           force_field_(to_copy.force_field_),
			                           number_of_iterations_(to_copy.number_of_iterations_),
			                           maximal_number_of_iterations_(to_copy.maximal_number_of_iterations_),
			                           energy_output_frequency_(to_copy.energy_output_frequency_),
			                           energy_difference_bound_(to_copy.energy_difference_bound_),
			                           max_gradient_(to_copy.max_gradient_),
			                           max_same_energy_(to_copy.max_same_energy_),
			                           same_energy_counter_(to_copy.same_energy_counter_),
			                           maximum_displacement_(to_copy.maximum_displacement_),
			                           force_update_counter_(to_copy.force_update_counter_),
			                           energy_update_counter_(to_copy.energy_update_counter_),
			                           cutlo_(to_copy.cutlo_),
			                           step_(to_copy.step_),
			                           abort_by_energy_enabled_(to_copy.abort_by_energy_enabled_),
			                           abort_energy_(to_copy.abort_energy_),
			                           aborted_(false),
                                       saved_position_(to_copy.saved_position_){}

// constructor with a force field
EnergyMinimizer::EnergyMinimizer(ForceField& ff){
                                            valid_ = setup(ff);
                                            if(!valid_)
                                                cout<<"EnergyMinimizer::EneryMinimizer(): Energy minimizer set up failed"<<endl;
                                            }

// constructor with a force field and a set of options
EnergyMinimizer::EnergyMinimizer(ForceField& ff, const MIN_OPTION & opt){
                            valid_ = setup(ff, opt);
                            if(!valid_)
                                cout<<"EnergyMinimizer::EneryMinimizer(): Energy minimizer set up failed"<<endl;
                            }
                            
// assignment operator
const EnergyMinimizer& EnergyMinimizer::operator = (const EnergyMinimizer& energy_minimizer)
	{
		// Guard against self assignment
		if (&energy_minimizer != this) 
		{
			// Copy the attributes
			option                       = energy_minimizer.option;
			valid_                        = energy_minimizer.valid_;
			force_field_                  = energy_minimizer.force_field_;
			number_of_iterations_         = energy_minimizer.number_of_iterations_;
			maximal_number_of_iterations_ = energy_minimizer.maximal_number_of_iterations_ ;
			energy_output_frequency_      = energy_minimizer.energy_output_frequency_;
			energy_difference_bound_      = energy_minimizer.energy_difference_bound_ ;
			max_gradient_                 = energy_minimizer.max_gradient_ ;
			max_same_energy_              = energy_minimizer.max_same_energy_; 
			same_energy_counter_          = energy_minimizer.same_energy_counter_;
			maximum_displacement_         = energy_minimizer.maximum_displacement_;
			force_update_counter_         = energy_minimizer.force_update_counter_;
			energy_update_counter_        = energy_minimizer.energy_update_counter_;
			cutlo_                        = energy_minimizer.cutlo_;
			step_                         = energy_minimizer.step_;
			abort_by_energy_enabled_      = energy_minimizer.abort_by_energy_enabled_;
			abort_energy_                 = energy_minimizer.abort_energy_;
			saved_position_               = energy_minimizer.saved_position_;
		}
		return (*this);
	}
                            
// destructor
EnergyMinimizer::~EnergyMinimizer(){}

//set up methods
bool EnergyMinimizer::setup(ForceField& ff){
    //store the force field
    force_field_ = &ff;
    valid_ = force_field_->isValid();
    if(!valid_)
        cout<<"EnergyMinimizer::setup():The force field bound to the minimizer is not valid"<<endl;
    //compute cutlo
    float epsilon = 1.;
    float eps = 1.;
    while (1. + eps > 1.)
		{
			epsilon = eps;
			eps /= 2.;
		}
    cutlo_ = sqrt((std::numeric_limits<float>::min)()/epsilon);
    //initilize default options;
    maximal_number_of_iterations_ =option.MAXIMAL_NUMBER_OF_ITERATIONS;
			
    energy_output_frequency_ = option.ENERGY_OUTPUT_FREQUENCY;
			
    number_of_iterations_ = option.NUMBER_OF_ITERATION;
			
    max_same_energy_ = option.MAX_SAME_ENERGY;
			
    energy_difference_bound_ = option.ENERGY_DIFFERENCE_BOUND;
			
    max_gradient_ = option.MAX_GRADIENT;
			
    maximum_displacement_ = option.MAXIMUM_DISPLACEMENT;
    
    energy_update_counter_ = 0;
    force_update_counter_ = 0;
    
    // specific setup method
    valid_ = specific_setup();
    if(!valid_)
         cout<<"EnergyMinimizer::setup(): EnergyMinimizer specific setup method failed"<<endl;

    return valid_;
}

// set up method with a force field and a set of option
bool EnergyMinimizer::setup(ForceField& ff, MIN_OPTION new_option){
                                   option = new_option;
                                   return setup(ff);
}

// specific setup
bool EnergyMinimizer::specific_setup(){
     return true;
     }
     
// minimize the bounded force field, virtual

bool EnergyMinimizer::minimize(int steps, bool resume){return false;}

//update the search direction
void EnergyMinimizer::update_direction(){return;}

// find the new step
double EnergyMinimizer::find_step(){return 0.0;}

// calculate current energy
double EnergyMinimizer::update_energy(){
       if(force_field_ != 0){
           current_energy_ = force_field_->update_energy();
           energy_update_counter_ ++;
       }
       //debug
       //cout<<"EnergyMinimizer::update_energy(): Current energy is "<<current_energy_<<endl;
       return current_energy_;
}

//calculate  current forces
void EnergyMinimizer::update_force(){
     force_field_->update_forces();
     // assign current gradient
     current_grad_.set(force_field_->get_atoms());
     force_update_counter_ ++;
     //debug
     //cout<<"EnergyMinimizer::update_force(): Current force RMS is "<<current_grad_.rms<<endl;
	 return;
     }

// check if is converged
bool EnergyMinimizer::is_converged() const{
     bool converged = ((current_grad_.rms <= max_gradient_)||(same_energy_counter_ >= max_same_energy_));
     //debug
     /*if(current_grad_.rms <= max_gradient_)
         cout<<"Is converged becaused of current_grad_.rms "<<current_grad_.rms<<" <= max_gradient_ "<<max_gradient_<<endl;
     else if(same_energy_counter_ >= max_same_energy_)
         cout<<"Is converged becaused of same_energy_counter_ >= max_same_energy_"<<endl;
         */
     return converged;
}
     
void EnergyMinimizer::print_energy() const{
          if(isValid()){
			cout << "iteration " << number_of_iterations_
								 << "  RMS gradient " << current_grad_.rms
								 << " kcal/(mol A)  total energy " << force_field_->get_energy() << " kcal/mol"
								 << endl;
                              }
}

void EnergyMinimizer::finish_iteration(){
     //perform a force field update in regular interval
     float max = 0.0;
     for(int i = 0; i<direction_.size(); ++i){
             float tmp = direction_[i].length();
             if(max<tmp)
                 max = tmp;
             }
             max = step_*sqrt(max);
             if (((force_field_->get_update_frequency() != 0)
				&& (number_of_iterations_ % force_field_->get_update_frequency() == 0))
						|| (max > 8.))
		                          {
			                          force_field_->update();
			                          initial_grad_.invalidate();
		                           }
             
		    // print the energy every energy_output_frequency_ iterations
#ifdef DEBUG
		     if ((energy_output_frequency_ != 0)
				&& (number_of_iterations_ % energy_output_frequency_ == 0))
		     {
                                          print_energy();
		      }
#endif
            // Check whether there the new energy and the old energy differ significantly
		     if (fabs(initial_energy_ - old_energy_) < energy_difference_bound_)
			// count if there is the same energy between last iteration and
			// this iteration
			   same_energy_counter_++;
		    else
		       same_energy_counter_ = 0;
		
		    // Increment the iteration counter
            ++number_of_iterations_;
          }
          


// Set the number of the current iteration
void EnergyMinimizer::set_number_of_iterations(int number_of_iterations){
		number_of_iterations_ = number_of_iterations;
	}
	
Gradient& EnergyMinimizer::get_gradient(){
		return current_grad_;
	}
	
Gradient& EnergyMinimizer::get_initial_gradient(){
		return initial_grad_;
	}
	
double EnergyMinimizer::get_energy() const{
		return current_energy_;
	}
	
double& EnergyMinimizer::get_energy(){
		return current_energy_;
	}
	
double EnergyMinimizer::get_initial_energy() const{
		return initial_energy_;
	}
	
double& EnergyMinimizer::get_initial_energy(){
		return initial_energy_;
	}
	
Gradient& EnergyMinimizer::get_direction(){
		return direction_;
	}
	
// Get the number of the current iteration
int EnergyMinimizer::get_number_of_iterations() const{
		return number_of_iterations_;
	}
	
void EnergyMinimizer::store_gradient_energy(){
		initial_energy_ = current_energy_;
		initial_grad_ = current_grad_;
	}
	
// Set the maximal number of iterations
void EnergyMinimizer::set_max_number_of_iterations(int maximal_number_of_iterations){
		maximal_number_of_iterations_ = maximal_number_of_iterations;
	}
	
// Get the maximal number of iterations
int EnergyMinimizer::get_max_number_of_iterations() const{
		return maximal_number_of_iterations_;
	}
	
// Is the energy minimizer valid: did the setup work?
bool EnergyMinimizer::isValid() const{
		return valid_;
	}
	
// Set the energy output frequency
void EnergyMinimizer::set_energy_output_frequency(int energy_output_frequency){
		energy_output_frequency_ = energy_output_frequency;
	}
	
// Get the energy ouput frequency
int EnergyMinimizer::get_energy_output_frequency() const{
		return energy_output_frequency_;
	}
	
// Set the energy difference bound
void EnergyMinimizer::set_energy_difference_bound(float energy_difference_bound){
		energy_difference_bound_ = energy_difference_bound;
	}
	
// Set explicitly the option max_gradient_
void  EnergyMinimizer::set_max_gradient(float max_gradient){
		max_gradient_ = max_gradient;
	}
	
// Get the current value of the maximum gradient bound
float EnergyMinimizer::get_max_gradient() const{
		return max_gradient_;
	}
	
// Set explicitly the number of iterations for detecting convergence due to invariant energy 
void  EnergyMinimizer::set_max_same_energy(int number){
		max_same_energy_ = number;
	}
	
// Get the value of max_same_energy, i.e. the number of iterations after which the algorithm is stopped when the
// energy remains constant
int EnergyMinimizer::get_max_same_energy() const{
		return max_same_energy_;
	}
	
// Get the energy difference bound
float EnergyMinimizer::get_energy_difference_bound() const{
		return energy_difference_bound_;
	}
	
// Set the maximal shift
void EnergyMinimizer::set_maximum_displacement(float displacement){
		maximum_displacement_ = displacement;
	}
	
// Get the maximal shift
float EnergyMinimizer::get_maximum_displacement() const{
		return maximum_displacement_;
	}

// Get the force field of the energy minimizer
ForceField*	EnergyMinimizer::get_force_field(){
		return force_field_;
	}
	
int EnergyMinimizer::get_force_update_counter() const{
		return force_update_counter_;
	}
	
int EnergyMinimizer::get_energy_update_counter() const{
		return energy_update_counter_;
	}
	
void EnergyMinimizer::enable_energy_abort_condition(bool state){
		abort_by_energy_enabled_ = state;
	}
	
bool EnergyMinimizer::energy_abort_condition_enabled() const{
		return abort_by_energy_enabled_;
	}
	
void EnergyMinimizer::set_energy_to_abort(float value){
		abort_energy_ = value;
	}
	
float EnergyMinimizer::get_energy_to_abort() const{
		return abort_energy_;
	}
	
bool EnergyMinimizer::was_aborted() const{
		return aborted_;
	}

void EnergyMinimizer::save_position(){
                 saved_position_.clear();
                 saved_position_.resize(force_field_->get_num_atoms());

		         // copy all positions
		         ATOMVec::const_iterator it(force_field_->get_atoms().begin());
		         vector<vector3>::iterator pos_it(saved_position_.begin());
		         for (; it != force_field_->get_atoms().end(); ++it, ++pos_it)
		  			      *pos_it = (*it)->get_position();
}
                                                       
void EnergyMinimizer::reset_position(){
     // move only if a saved position exists for every atom
		if (saved_position_.size() == force_field_->get_num_atoms())
		{
			ATOMVec::iterator it(force_field_->get_atoms().begin());
			vector<vector3>::const_iterator pos_it(saved_position_.begin());
			for (; it != force_field_->get_atoms().end(); ++it, ++pos_it)
			{
				(*it)->set_position(*pos_it);
			}
		}
}

void EnergyMinimizer::move_to(const Gradient& gradient, float step)
	{
		// move only if a saved position exists for every atom
		if (gradient.size() == force_field_->get_num_atoms())
		{   
			// use the saved positions
			if (saved_position_.size() == force_field_->get_num_atoms())
			{    

				ATOMVec::iterator it(force_field_->get_atoms().begin());
				vector<vector3>::const_iterator pos_it(saved_position_.begin());
				vector<vector3>::const_iterator grad_it(gradient.begin());
				for (; it != force_field_->get_atoms().end(); ++it, ++pos_it, ++grad_it)
				{
					(*it)->set_position(*pos_it + *grad_it * step);
					//debug
						//cout<< "   - atom " << (*it)->get_symbol_type() << " @ " << (*it)->get_position() << std::endl;
				}
			}
			// we don't have saved positions, use the current atom positions
			else 
			{
				ATOMVec::iterator it(force_field_->get_atoms().begin());
				vector<vector3>::const_iterator grad_it(gradient.begin());
				for (; it != force_field_->get_atoms().end(); ++it, ++grad_it)
				{    
                    vector3 temp = (*it)->get_position();
					temp += *grad_it * step;
					(*it)->set_position(temp);
					//debug
						//cout << "   - atom " << (*it)->get_symbol_type() << " @ " << (*it)->get_position() << std::endl;
				}				
			}
		}
	}
