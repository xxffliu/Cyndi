/*
  Name: ConGen.h
  Copyright: DDDC
  Author: Xiaofeng Liu
  Date: 08-01-08 09:34
  Description: Using MOGA method to generate a series of low-energy diverse conformers       
*/

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <vector>
#include <iterator>

#include "TAFF.h"
#include "MMFF94.h"

#ifndef INPUT_H
	#include "ParamInput.h"
#endif
#include "random.h"

//const int MAX = 4;
//const double DELTA = 1e-6;
//const int MAX = 4; // The number of objective functions 
const int MAX_ARCHIVE_SIZE = 10000; // maximum size possible (the upper bound) of the archive
// epsilon matrix, this controls the num of solutions generated into the archive. The smaller of
// these value, the more solutions remained in the final archive
// Shoule be treated as a user-defined parameter
const double epsilon[4] = {15, 8, 0.5, 0.1};
const int supremumx = 71, infimumx = -72; // the upper and lower bounds for variables

const float torsionStep = 2.5;


struct Conformer
{
	friend bool operator<(const Conformer& lhs, const Conformer& rhs);
    vector<int> torsions;
    float TorsionEnergy;
	float VDWEnergy;
	float TotalEnergy;
    float rmsd;
	float GyrationRdius;
};

inline bool operator<(const Conformer& lhs, const Conformer& rhs)
{
	return lhs.VDWEnergy + lhs.TorsionEnergy < rhs.VDWEnergy + rhs.TorsionEnergy;
}
// Class definitions for an individual
struct individual
{
       vector<int> xreal;
       vector<double> f;        // the objective functions
       double cv;        // cv - parameter for constraint violation
	   float energy;
	   float gyration_radius;
       vector<double> box;        // box vectors for the individual
       individual();
	   individual& operator=(const individual& rhs);
       void clear();
};
              
class MOGA
{
      public:
             friend struct individual;
             // constructors and destructor
             MOGA();
             //MOGA(const MOGA& to_copy);
             MOGA(const MOL& mol);
             virtual ~MOGA();
             
             // clear method
             void clear();
             // setup method
             bool setup(const MOL& mol);
             //execuation method
             vector<Conformer> execuateMOGA();
             
      private:
              // member variables
			  float gyration_;
			  float MinEnergy_;
			  float Ecutoff_;
              vector<individual> population;
              vector<individual> archive;
              individual newchild, newchild1, newchild2;
              // Number of decision variables
              int N_of_x;
              MOL mol_;
			  // original position of the input mol
			  vector<vector3> input_pos_;
              TAFF taff_;
			  MMFF94 mmff94_;
              int archive_size;    //starting archive size = 0
              // member funtions
              // objective function
              void objective(individual* indv);
              // Create the random pop
              void create_random_pop ();
              // Strict Domination check function for population members indexed by m and n
              // returns 1 ( 0 ) meaning m dominates (does not dominate) n in minimization sense
              // only for population members _to_be_called_in_  create_archive()
              int strict_dom_check (int m, int n);
              // Create the initial archive
              // Consists of the non-dominated (strict domination) members
              // of the population
              void create_archive ();
              // Box funtion definitons, here done for minimization case
              void box_func (individual * ind1);
              // function to check for Box domination
              int box_dom (individual ind1, individual ind2);
              //to check whether 2 individuals are in the same box or not
              int same_box_check (individual ind1, individual ind2);
              // Check for domination (usual sense)
              int dom_check (individual ind1, individual ind2);
              // Finds the distance between corner[] and between the decision variable
              // coordinates of the individual
              double distance (individual ind1, vector<double> point);
              int con_update ();
              // The Update function, called to update the archive in every iteration
              int update ();
              void generate_replace ();
              int tournament (int p1, int p2);
              void compete ();
              // xrossover and mutation method
              void realcross (individual parent1, individual parent2);
              void real_mutate (individual * new_pop_ptr);
			  void remove_redundency_archive();

			  /************************* Output function for debug ****************/
			  // print rmsd and gyration radius for all individules
			  void print_rmsd_gyration_(int no_gen)
			  {
				  ofstream out("all_rmsd_gyration.out", ofstream::app);
				  for (int i = 0; i < MOGAParam_.PopSize_; i++)
				  {
					  out<<population[i].f[0]<<setw(16)<<population[i].f[1]<<setw(16)<<population[i].f[2]<<endl;
				  }
				  out.close();
			  }

			  // Prints the population members, whenever called
			  void print_pop_(int no_gen)
			  {
				  string name = "Cyndi_" + mol_.get_name() + "_all_population.out";
				  std::ofstream out(name.c_str(), ofstream::app);
				  out << "*************** Generation "<<no_gen<<" ***************"<<endl;
				  for (int i = 0; i < MOGAParam_.PopSize_; i++)
				  {
					  out<<"population "<<i<<": "<<setw(16);
					  /*for(int k = 0; k < population[i].xreal.size(); k++)
						  out<<population[i].xreal[k]<<setw(16);*/
					  for (int j = 0; j < MOGAParam_.NumObjects_; j++)
						  out << population[i].f[j] <<setw(16);
					  out <<endl;
				  }
				  out.close();
			  }


			  // Prints the function values of the archive members, when called
			  void print_archive_(int no_gen)
			  {
				  std::ofstream out(("Cyndi_" + mol_.get_name() + "_all_archives.out").c_str(), ofstream::app);
				  out << "*************** Generation "<<no_gen<<" ***************"<<endl;
				  //out << "\n Printing the function values of the archive members\n";
				  for (int i = 0; i < archive_size; i++)
				  {
					  out<<"archive "<<i<<": "<<setw(16);
					  /*for(int k = 0; k < archive[i].xreal.size(); k++)
						  out<<archive[i].xreal[k]<<setw(16);*/
					  for (int j = 0; j < MOGAParam_.NumObjects_; j++)
						  out << archive[i].f[j] <<setw(16);
					  //out <<archive[i].gyration_radius<<endl;
				  }
				  out.close();
			  }
			  // Prints into files the population and archive members
				// Also prints the Box vectors and Decision Variables of the archive members
				void print_func_values_()
				{
					FILE *fpop, *farch, *fvar, *fbox;
					string name = "Cyndi_" + mol_.get_name();
					fpop = fopen ((name + "_popltn.out").c_str(), "w");
					farch = fopen ((name + "_archive.out").c_str(), "w");
					fvar = fopen ((name + "_ar_variables.out").c_str(), "w");
					fbox = fopen ((name + "_box_vectors.out").c_str(), "w");
					cout << "\n\n The size of the archive is " << archive_size << "\n";
					//cout << " The final archive is in file archive.out\n\n";
					for (int i = 0; i < MOGAParam_.PopSize_; i++)
					{
						for (int j = 0; j < MOGAParam_.NumObjects_; j++)
						fprintf (fpop, "%f\t", population[i].f[j]);
						fprintf (fpop, "\n");
					}
					for (int i = 0; i < archive_size; i++)
					{
					for (int j = 0; j < MOGAParam_.NumObjects_; j++)
					{
						fprintf (farch, "%f\t", archive[i].f[j]);
						fprintf (fbox, "%f\t", archive[i].box[j] * MOGAParam_.EPSILON_[j]);
					}

						for (int j = 0; j < N_of_x; j++)
						fprintf (fvar, "%f\t", archive[i].xreal[j]);

						fprintf (farch, "\n");
						fprintf (fbox, "\n");
						fprintf (fvar, "\n");
					}
				}

};
