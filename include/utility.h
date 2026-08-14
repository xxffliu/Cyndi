#ifndef UTILITY_H
#define UTILITY_H


#include <vector>
#include <cmath>
#include "Vector3.h"
#include "Matrix3x3.h"

/*#ifndef DEBUG
#define DEBUG
#endif*/

/*#ifndef TMP_FILE
#define TMP_FILE
#endif*/

//define the upper infinite boundary
const float GAINFINITY = 1e7;
//define the lower psitive boundary
const double EPSILON = 1e-6;

inline bool isNearZero(double num)
{
	return fabs(num) < EPSILON;
}
/*inline ifstream& fopener(ifstream &in, const string& file){
          in.close();
          in.clear();
          in.open(file.c_str());
          return in;
          }
vector<vector<string> >& ParameterSpliter(vector<vector<string> >& lines, ifstream& fin){
                      string buffer;
                      vector<string> fields;
                      while (getline(fin, buffer)){
                          istringstream sin(buffer);
                          string token;
                          while (cin>>token)
                              fields.push_back(token);
                          lines.push_back(fields);
                          }
                      return lines;
                      }

enum features{
     HYDRO = 1,
     HBA = 2,
     HBD = 3,
     POS_ION = 4,
     NEG_ION = 5,
     POS_CHAR = 6,
     NEG_CHAR = 7
     };
 
*/ 
// convert string type to corresponding int type
inline int str2int(string str){
             return atoi(str.c_str());
             }
// convert string type to correspoding double type
inline double str2double(string str){
       return atof(str.c_str());
       }
       
// rotate this point around a specified axis
// take the axis in the world coordinates system
// return the final position after rotation
// we must provide an axis vector, a atom postion of the bond and the rotation angle       
inline vector3 myrotate(vector3& p, vector3& v, vector3& orig, double angle){
     matrix3x3 mat;
     mat.RotAboutAxisByAngle(v, angle);
     vector3 tmp = p - orig;
     return mat*tmp + orig;
     }
// the overloaded version of myrote, we provide the axis vector, the vector to be rotate and the rotate angle,
// return the rotated angle
inline vector3 myrotate(vector3& vec, vector3& axis, double angle)
{
	matrix3x3 mat;
	mat.RotAboutAxisByAngle(axis, angle);
	return mat * vec;
}

inline matrix3x3 KabschFit(const vector<vector3>& r, const vector<vector3>& f,int size)
{
    double xxyx(0.), xxyy(0.), xxyz(0.);
    double xyyx(0.), xyyy(0.), xyyz(0.);
    double xzyx(0.), xzyy(0.), xzyz(0.);
    double d[4] = {0, 0, 0, 0}, q[4] = {0,0,0,0};
    double c[16],v[16];

    /* generate the upper triangle of the quadratic form matrix */

    for (int i = 0; i < size; ++i)
	{
		double rx = r[i].x(), ry = r[i].y(), rz = r[i].z();
        double fx = f[i].x(), fy = f[i].y(), fz = f[i].z();
        xxyx += fx * rx;
        xxyy += fx * ry;
        xxyz += fx * rz;
        xyyx += fy * rx;
        xyyy += fy * ry;
        xyyz += fy * rz;
        xzyx += fz * rx;
        xzyy += fz * ry;
        xzyz += fz * rz;
	}

    c[4*0+0] = xxyx + xyyy + xzyz;

    c[4*0+1] = xzyy - xyyz;
    c[4*1+1] = xxyx - xyyy - xzyz;

    c[4*0+2] = xxyz - xzyx;
    c[4*1+2] = xxyy + xyyx;
    c[4*2+2] = xyyy - xzyz - xxyx;

    c[4*0+3] = xyyx - xxyy;
    c[4*1+3] = xzyx + xxyz;
    c[4*2+3] = xyyz + xzyy;
    c[4*3+3] = xzyz - xxyx - xyyy;

    /* diagonalize c */

    matrix3x3::jacobi(4, c, d, v);
    /* extract the desired quaternion */

    q[0] = v[4*0+3];
    q[1] = v[4*1+3];
    q[2] = v[4*2+3];
    q[3] = v[4*3+3];

    /* generate the rotation matrix */
	double u[3][3];
    u[0][0] = q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3];
    u[1][0] = 2.0 * (q[1] * q[2] - q[0] * q[3]);
    u[2][0] = 2.0 * (q[1] * q[3] + q[0] * q[2]);

    u[0][1] = 2.0 * (q[2] * q[1] + q[0] * q[3]);
    u[1][1] = q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3];
    u[2][1] = 2.0 * (q[2] * q[3] - q[0] * q[1]);

    u[0][2] = 2.0 * (q[3] * q[1] - q[0] * q[2]);
    u[1][2] = 2.0 * (q[3] * q[2] + q[0] * q[1]);
    u[2][2] = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3];
    matrix3x3 rot(u);
	return rot;
    //return d[3];
  }
  


inline bool redundantRing(vector<int> ring)
{
	return ring.size() >= 8 ? true : false;
}
/*void clear_all(molecule& all){
      all.type = 0; all.x1 = 0; all.x2 = 0;all.y1 = 0; all.y2 = 0;all.z1 = 0;all.z2 = 0;
      }
*/      
/*       
//typedef vector< vector<string> > 

// Some constants defination

// EPSILON (used fr comparisons)
		double EPSILON = 1e-6;

		// PI
		const double  pi = 3.14159265358979323846L;

		// Euler's number - base of the natural logarithm
		const double  E  = 2.718281828459045235L;
		
		//	Elementary charge.
		const double	ELEMENTARY_CHARGE = 1.60217738E-19L;  	 // C     
	
		/// Elementary charge (alias)
		const double	e0								=	ELEMENTARY_CHARGE;

		// Electron mass.
		const double	ELECTRON_MASS   	= 9.1093897E-31L;   	 // kg

		// Proton mass.
		const double	PROTON_MASS     	= 1.6726230E-27L;   	 // kg

		// Neutron mass.
		const double	NEUTRON_MASS    	= 1.6749286E-27L;   	 // kg

		// Avogadro constant.
		const double	AVOGADRO        	= 6.0221367E+23L;   	 // 1 / mol

		// Avogadro constant (alias)
		const double	NA								= AVOGADRO;

		// Avogadro constant (alias)
		const double	MOL             	= AVOGADRO;

		// Boltzmann constant.
		const double	BOLTZMANN       	= 1.380657E-23L;  	   // J / K

		// Boltzmann constant (alias)
		const double	k	        				= BOLTZMANN;
		
		// Planck constant.
		const double	PLANCK 	         	= 6.6260754E-34L;      // J * sec

		// Planck constant (alias)
		const double	h       	   			= PLANCK;

		// Gas constant (= NA * k)	
		const double	GAS_CONSTANT 	  	= NA * k;

		// Gas constant (alias)
		const double R 								= GAS_CONSTANT;

		// Faraday constant (= NA * e0)
		const double	FARADAY         	= NA * e0;

		// Faraday constant (alias)
		const double	F    							= FARADAY;

		// Bohr radius.
		const double	BOHR_RADIUS     	= 5.29177249E-11L;     // m

		// Bohr radius (alias)
		const double	a0     						= BOHR_RADIUS;

		//  the following values from: 
		//  P.W.Atkins: Physical Chemistry, 5th ed., Oxford University Press, 1995

		// Vacuum permittivity.
		const double	VACUUM_PERMITTIVITY    	= 8.85419E-12L;     // C^2 / (J * m)

		// Vacuum permeability.
		const double	VACUUM_PERMEABILITY     = (4 * pi * 1E-7L);	// J s^2 / (C^2 * m)

		// Speed of light.
		const double	SPEED_OF_LIGHT          = 2.99792458E+8L;	  // m / s

		// Speed of Light (alias)
		const double	c												= SPEED_OF_LIGHT;

		// Gravitational constant.
		const double	GRAVITATIONAL_CONSTANT  = 6.67259E-11L;    	// N m^2 / kg^2

		// Fine structure constant.
		const double	FINE_STRUCTURE_CONSTANT = 7.29735E-3L;   		// 1      
			
		// Degree per rad.
		const double	DEG_PER_RAD				= 57.2957795130823209L;

		// Rad per degree.
		const double	RAD_PER_DEG			 	= 0.0174532925199432957L;

		// mm per inch.
		const double	MM_PER_INCH 			= 25.4L;

		// m per foot.
		const double	M_PER_FOOT  			= 3.048L;

		// Joule per calorie
		const double	JOULE_PER_CAL     = 4.184;

		// Calories per Joule.
		const double	CAL_PER_JOULE     = (1 / 4.184);
*/
#endif
