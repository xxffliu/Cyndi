#ifndef GRADIENT_H
#define GRADIENT_H

#ifndef VECTOR3_H
#include "Vector3.h"
#endif

#ifndef MOL_H
#include "Mol.h"
#endif

#include <vector>

class ATOM;
class Gradient:public vector<vector3>{
      public:
             Gradient();
             Gradient(const ATOMVec& avec);
             Gradient(const Gradient& grad);
             virtual ~Gradient();
             const Gradient& operator=(const Gradient& rhs);
             //const Gradient& operator=(const ATOMVec& rhs);
             
             void set(const Gradient& grad);
             void set(const ATOMVec& avec);
             
             // Make the force value negative
             void minus();
             
             // Normalize the gradient, rescale to unity length;
             
             void normalize();
             
             // Dot product operator, invalid if the two gradients have different sizes
             
             double operator*(const Gradient& grad) const;
             
             const vector3& operator [] (int i) const { return std::vector<vector3>::operator [] (i); }
		     vector3& operator [] (int i) { return std::vector<vector3>::operator [] (i); }

		    // Invalidate the gradient.	
		     void invalidate();

		    //	Return the validity flag.
		     bool isValid() const;

		/**	The gradient norm.
		*/
		     double norm;

		/**	The inverse of the gradient norm.
		*/
		     double inv_norm;

		/**	The root mean square of the gradient.
		*/
		     double rms;

       protected:
			
		//	The validity flag.

             bool valid_;
		
};

#endif
