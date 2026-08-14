#include "../include/Gradient.h"
#include "../include/Atom.h"

using namespace std;

static double NA = 6.0221367E+23L;   	 // 1 / mol

Gradient::Gradient():norm(0.0),inv_norm(0.0),rms(0.0),valid_(false){}

Gradient::Gradient(const ATOMVec& avec){
                         set(avec);
                         }

Gradient::Gradient(const Gradient& grad):
                         vector<vector3>(){
                         set(grad);
                         }
                         
const Gradient& Gradient::operator=(const Gradient& rhs){
      set(rhs);
      return *this;
      }


Gradient::~Gradient(){
                      valid_ = false;
                      }

void Gradient::set(const ATOMVec& avec){
     // resize to hold all vectors
     size_type max_index = (size_type)avec.size();
     resize(max_index);
     // copy all forces, the gradient is the negative force and is
     // stored in units of kJ/(mol A). The forces
     // are in units of Newton, so we have to use
	// a conversion factor of -1.0 / 1e3 (J->kJ) / 1e10 (m->A) * NA (1->mol)
	norm = 0.0;
	vector<vector3>::iterator it(begin()); 
	for (int i = 0; i < max_index; ++i, ++it){
        *it = avec[i]->get_force() * NA / -1.0e13;
        norm += (*it) * (*it);
		}

    // calculate the norm and its inverse
    norm = sqrt(norm);
    inv_norm = 1.0 / norm;
    if (max_index > 0)
        rms = norm / sqrt(3.0 * (double)max_index);
    else 
		rms = 0.0;

    // the gradient is now valid
		valid_ = true;
	}

void Gradient::set(const Gradient& grad){
     //copy the gradient
     resize(grad.size());
     copy(grad.begin(), grad.end(),begin());
     // copy the members
     norm = grad.norm;
     inv_norm = grad.inv_norm;
     rms = grad.rms;
     valid_ = grad.valid_;
     }
void Gradient::minus(){
     // iterate vector and switch the sign
     for (vector<vector3>::iterator giter = begin(); giter != end(); ++giter)
         *giter *= -1.0;
}

void Gradient:: normalize(){
     for (vector<vector3>::iterator it = begin(); it != end(); ++it)
			*it *= inv_norm;
    // reset the norm and its inverse and calculate the root mean square
     norm = 1.0;
     inv_norm = 1.0;
     if (size() > 0)
        rms = 1.0 / sqrt(3.0 * (double)size());
     else 
		rms = 0.0;
}

double Gradient::operator * (const Gradient& gradient) const{
    size_type max_index = size();
    if (gradient.size() != max_index){
			cout<<"Error:Gradient::operator *(): The sizes of 2 gradient are not the same"<<endl;
			return 0.0;
            }
	double result = 0.0;
	for (size_type i = 0; i < max_index; i++)
		result += operator[](i) * gradient[i];
	return result;
	}

bool Gradient::isValid() const{
     return valid_;
     }
void Gradient::invalidate(){
     valid_ = false;
     }
