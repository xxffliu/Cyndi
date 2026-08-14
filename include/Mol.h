#ifndef MOL_H
#define MOL_H
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <deque>

#ifndef ATOM_H
#include "Atom.h"
#endif

#ifndef BOND_H
#include "Bond.h"
#endif

#ifndef FRAGMENT_H
#include "Fragment.h"
#endif

#ifndef GRAPH_H
#include "Graph.h"
#endif

#ifndef MATRIX3X3_H
#include "Matrix3x3.h"
#endif
#include "utility.h"


#define GASTEIGER_DENOM  20.02
#define GASTEIGER_DAMP   0.5
#define GASTEIGER_ITERS  10
typedef vector<ATOM*> ATOMVec;
typedef vector<BOND*> BONDVec;
typedef vector<FRAGMENT*> FRAGVec;

// a structure data type for characterize ring
struct RING
{
       ATOMVec vatom;
       BONDVec vbond;
       int size;
       bool is_aromatic;
	   // get the center of this ring
       vector3 get_center();
	   //finding that the atom a whether is one of the member of the ring
	   bool IsMember(ATOM *a);
	   // return number of the PI electrons
	   int NumPIElectrons();
	   // constructors
	   RING();
	   RING(const RING& to_copy);
	   RING& operator=(const RING& rhs);
	   virtual ~RING();
       void clear();
};

// check if bridged or fused rings
inline bool bridgedRing(const RING& ring)
{
	return ring.size < ring.vbond.size() ? true : false;
}
inline bool largeRing(const vector<int>& ring)
{
	return ring.size() >= 8 ? true : false;
}

 class  GasteigerState
	{
	public:
		GasteigerState():a(0.), b(0.), c(0.),denom(0.),chi(0.),q(0.){}
				~GasteigerState() {}

				void SetValues(double _a,double _b,double _c,double _q)
				{
					a = _a;
					b = _b;
					c = _c;
					denom=a+b+c;
					q = _q;
				}

				double a, b, c;
				double denom;
				double chi;
				double q;
};

// a class 

class MOL{
      public:
             friend bool operator<(const MOL& lhs, const MOL& rhs);
			 friend bool GreaterMol(const MOL& lhs, const MOL& rhs);
			 friend class MOL2IO;	
             // member
             // mapping the atoms to their types.
             map<ATOM*, vector<string> > atom_type_map;
             // mapping the rotatable parts of the molecule to their rotor bonds
             map<int, ATOMVec> rotation_map;
             // mapping the rotor bonds id to a increcent order
             map<int, int> rotor_id_map;
             
             //member methods
             // constuctor list, copy constructor and assignment operator overloading
             MOL();
             MOL(const MOL& mol0);
             MOL& operator=(const MOL& mol0);
             virtual ~MOL();
             void clear();
             // a genereal method to initialize the mol, including neighbot atoms and bonds
             // lists, searching ring systems, perceive the rotatable bonds and fragmentization
             // as well as guanidino or carboxyl group atom types
             // patching (assign guanidino N to N.pl3, C to C.cat, carboxyl O to O.co2
             // C to C.2, all bonds type to ar)
             void initialize();
			 bool is_initialized();
             //backup current position
             void bk_position();
             // 'get' and 'set' methods
             string get_name();
             string get_charge_type();
             int get_num_atom();
             int get_num_bond();
             int get_num_fragment();
             int get_num_ring();
			 inline int get_num_arom_ring()
			 {
				 int counter(0);
				 for(vector<RING>::iterator it  = _vring.begin(); it != _vring.end(); ++it)
				 {
					 if(it->is_aromatic)
						 counter += 1;
				 }
				 return counter;
			 }
             ATOMVec& get_atom_vector();
             vector<vector3> get_coordinates();
			 vector<vector3> get_heavy_coordinates();
             ATOM* get_atom(int);
             ATOM* get_first_atom();
             BONDVec& get_bond_vector();
             BOND* get_bond(int);
             BOND* get_first_bond();
             vector<RING>& get_ring_vector();
             int get_num_HB_acceptor();//to do
             int get_num_HB_donor(); // to do
             inline float get_MW()
			 {
				 return MW_;
			 }
             int get_num_of_rot_bonds();
             //BOND* get_bond(ATOM&,ATOM&);
             FRAGVec& get_fragment_vector();
             FRAGMENT* get_fragment(int);
             string get_comment();
			 string get_pharmacophore_name();
             double get_rmsd();
             // calculate the rmsd between current position and the original one
             double updateRMSD();
             // calculate the rmsd value between two molecules
             double get_RMSD(MOL& mol);
             // return the barycenter
             vector3 get_barycenter();
             // return the fitness score
             float get_fitness();
             // return the current energy
             // Caution! the current energy only record the energy calculated in ForceField, and
             // can only be updated in ForceField
             float get_energy();
             // set up methods
			 void set_comment(string comment = "");
			 void set_pharmacophore_name(const string& name);
             void set_name(const string& name0);
             void set_num_atom(int n);
             void set_num_bond(int n);
             void set_num_fragment(int n);
             void set_atom_vector(ATOMVec&);
             void set_bond_vector(BONDVec&);
             void set_coordinates(vector<vector3> coord);
             void set_fitness(float fit);
             void set_energy(float energy);
             void set_rmsd(float rmsd);
             void set_charge_type(string type);
             // initialize neighbor list of each atom
             void init_neighbor_list();
             //////////////////////////////////////
             // fragilize the molecule
             bool fragmentize();
             void set_fragment_vector(FRAGVec&);
             void set_comment(const string&);
             friend istream& operator >> (istream& is, MOL& mol0);
             friend ostream& operator << (ostream& os, MOL& mol0);
             
             // ring finding related methods
             void find_rings();
             //////////////////////////////////////
             // rotatable bonds finding related methods
             void find_rotatable_bonds();

             
             // other methods
             void swap(MOL& mol1,MOL& mol2); 
             
             // reset the molecule to its original position
             void reset();
             // transform the coordinates according to a transform matrix
             void transform(matrix3x3& mat);
             void transform(double (*mat)[4]);
             // rotate part of the molecule around the specified rotor bond
			 // CAUTION: the angle parameter is DEG format; 
             void rotate(int rotor_id, double rotate_angle);
             void apply_rotor(int order, double angle);
             // find the rotatable parts by each rotors and packed them into rotation_map
             void set_rotation_mark();
             void clear_rotation_mark();
             
             //compute the geomitrical barycenter and the atom weight weighted barycenter of the molecue
             void compute_center();
             vector3 update_center();
             void compute_weight_center();
             // center the molecule according to its barycenter
             void center_pos();
             inline vector3 get_center(){
                    compute_center();
                    return barycenter_;
                    }
             void apply_trans_rot(const vector3& center, const vector3& target_center, 
                                  const vector3& rot = VZero, const vector3& trans = VZero);
			 // calculate the gyration radius of the molecule
			 // a gyration radius is defined as the mean square distance between all heavy atoms and the molecular cenroid
			 float ComputGyrationRadius();
             
             //RMSD minimizer
			 // This method computes the optimal transformation mapping
			 // one set of three-dimensional points onto another set of
			 // points. It implements the algorithm by Coutsalis et al.
			 //(J. Comput. Chem., 25(15), 1849 (2004)), which computes
			 //the RMSD-optimal transformation by solving an eigenvalue
			 // problem.
			 // the position of fixed_mol is fixed during the process
			 // return the minimized rmsd
			 double minimizeRMSD();
			 double minimizeRMSD(MOL& fixed_mol);
			 // calculate gasteiger masili charge for each atom
		     bool calculateGasteigerCharge();
              // judging  whether the atom  in the ring of the given size
			 bool IsInRingSize(ATOM * ,int ) ;
			 // fixed by xfliu, 20090226
			 bool IsInAromaticRingSize(ATOM* atom, int size);
			 // perceiving the aromatic rings
			 void PerceiveAromaticRings();
		     //return the bond between the atoms @p bgn and @p end or NULL if none exists
		      BOND  *GetBond(ATOM* bgn, ATOM* end) const;
			  //get the atom when we knowthe numeric atom type of  mmff94
			  ATOM * get_map(int type) ;
			  map<int,ATOM*>map_atom;
			
      private:
              string name;
              int num_atom;
              int num_bond;
              int num_fragment;
              int num_rotatable_bonds_;
              FRAGVec _vfragment;
              ATOMVec _vatom;
              BONDVec _vbond;
              vector3 barycenter_;
              vector<vector3> bk_pos_;
              // vertex vector -- initiate from atom information
              deque<VERTEX> vertex_;
              // edge vector -- initiate from bond information
              deque<EDGE> edge_;
              // rings vector
              vector<vector<int> > rings_;
              vector<RING> _vring;
              string comment_;
              bool is_fragmentized_;
			  bool is_initialized_;
              
			  // molecular weight
			  float MW_;
			  // the pharmacophore name to be aligned with
			  string pharmacophore_name_;
              // the fitness score of this molecule
              float fitness_;
              float energy_;
              // this rmsd refers to the rmsd value between current coordinates sys and the backup one
              float rmsd_;
              string charge_type_;
			  std::vector <GasteigerState> _gsv; //!< vector of internal GasteigerState (for each atom)

			  // initiate vertex and edge
              void InitiateGraph_();
              // prune the endpoint bonds attatched to the ring
              void PruneEndBonds_();
              // main method collapsing the ring
              void CollapsRing_(VERTEX& curr_vertex);
			  // combine fused rings into single ring
			  void CombineFusedRings_();
			  void InitialPartialCharges();
			  bool GasteigerSigmaChi(ATOM *atom,double &a,double &b,double &c);
             
              
};   
               
#endif