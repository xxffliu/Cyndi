#include "../include/ConGen.h"

individual::individual():cv(0), energy(0.0)
{
	xreal.clear();
	f.resize(4);
	box.resize(4);
}
individual& individual::operator =(const individual &rhs)
{
	xreal.assign(rhs.xreal.begin(),rhs.xreal.end());
	f = rhs.f;
	box = rhs.box;
	energy = rhs.energy;
	gyration_radius = rhs.gyration_radius;
	return *this;
}
MOGA::MOGA():archive_size(0), taff_()
{
	// v1 (2026): keep the electrostatic component enabled so the search
	// objective includes full MMFF94/TAFF energy (H-bonds & charge distribution).
	// Originally commented out: taff_.remove_component("TAFF ELE");
	// Originally commented out: mmff94_.remove_component("MMFF94 ELE");
	for(int i = 0; i < 4; i++)
	{
		if(isNearZero(MOGAParam_.EPSILON_[i]))
			MOGAParam_.EPSILON_[i] = epsilon[i];
	}
}
MOGA::MOGA(const MOL& mol):archive_size(0),taff_()
{
	// v1 (2026): electrostatics enabled during search (see note in ctor above)
	for(int i = 0; i < 4; i++)
	{
		if(isNearZero(MOGAParam_.EPSILON_[i]))
			MOGAParam_.EPSILON_[i] = epsilon[i];
	}
	setup(mol);
}

bool MOGA::setup(const MOL& mol)
{
	MinEnergy_ = 1000.0;
	mol_ = mol;
	mol_.initialize();
	mol_.bk_position();
	input_pos_ = mol_.get_coordinates();
	//MOGAParam_.EPSILON_[1] = mol_.get_num_of_rot_bonds() * 1.0;
	//debug
	//gyration_ = mol_.ComputGyrationRadius();
	bool FFSetup = false;
	if(MOGAParam_.FFType_ == FF_TAFF)
		FFSetup = taff_.setup(mol_);
	else
		FFSetup = mmff94_.setup(mol_);
	if(!FFSetup)
	{
		cout<<"Conformation generation failed due to FF assignment failed."<<endl;
		return false;
	}
	N_of_x = mol_.get_num_of_rot_bonds();
	population.clear();
	archive.clear();
	int popSize = MOGAParam_.PopSize_ + 20 * mol_.get_num_of_rot_bonds();
	population.resize(popSize);
	archive.resize(MAX_ARCHIVE_SIZE);
	newchild.xreal.resize(N_of_x);
	newchild1.xreal.resize(N_of_x);
	newchild2.xreal.resize(N_of_x);
	for(vector<individual>::iterator it = archive.begin(); it != archive.end(); ++it)
		it->xreal.resize(N_of_x);
	return true;
}

void individual::clear()
{
	xreal.clear();
	xreal.resize(4);
	f.clear();
	f.resize(4);
	box.clear();
	box.resize(4);
	energy = 0.;
}

void MOGA::clear()
{
	mol_.clear();
	input_pos_.clear();
	taff_.clear();
	N_of_x = 0;
	population.clear();
	archive.clear();
}

MOGA::~MOGA()
{
	clear();
}

// Main function starts here
vector<Conformer> MOGA::execuateMOGA()
{
	// check the number of objectives
	if(MOGAParam_.NumObjects_ <= 1 || MOGAParam_.NumObjects_ > 4)
	{
		cout<<"Error! Number of MOGA objectives is neither 2, 3 nor 4"<<endl;
		exit(1);
	}
	//define the returned vector containing the final archive variables
	vector<Conformer> arch;
	Conformer tmp;

	int child_flag1 = 0, flg = 0, child_flag2 = 0;
	if(!MOGAParam_.UseInputRandomSeed_)
		// generate the random number seed based on current time
		MOGAParam_.BasicSeed_ = TimeRandomSeed();

	warmup_random (MOGAParam_.BasicSeed_);
	create_random_pop ();
	create_archive();
#ifdef DEBUG
	cout <<endl<<"Starting GA loop........";
#endif
	for (int no_gen = 0; no_gen < MOGAParam_.MaxNumGen_; no_gen++)
	{

		generate_replace ();
		// cout <<"GEN = " << no_gen + 1<<endl;
		// cout << "Number of members in Archive at present :" << archive_size <<endl;
		newchild = newchild1;
		//update the archive
		child_flag1 = con_update ();
		// update the population
		compete ();
		newchild = newchild2;
		child_flag2 = con_update ();
		compete ();

		//print_pop_(no_gen);
		//print_archive_(no_gen);

	}
#ifdef DEBUG
	cout<<"Done"<<endl<<"Remove redundant conformers..."<<endl;
#endif
	//remove_redundency_archive();
	// temperary files
#ifdef TMP_FILE
	print_func_values_();
#endif
	//int OutNum = (archive_size <= MOGAParam_.MaxNumConformer_) ? archive_size : MOGAParam_.MaxNumConformer_;
	int OutNum = archive_size;
	for (int i = 0; i < OutNum; i++)
	{
		tmp.torsions = archive[i].xreal;
		tmp.VDWEnergy = archive[i].f[0];
		tmp.TorsionEnergy = archive[i].f[1];
		tmp.TotalEnergy = archive[i].energy;
		tmp.rmsd = archive[i].f[2];
		//tmp.GyrationRdius = 1/archive[i].f[3];
		tmp.GyrationRdius = -archive[i].f[3];
		arch.push_back(tmp);
	}
	mol_.reset();
	arch.resize(OutNum);
	return arch;
}

// objective functions
void MOGA::objective(individual* indv)
{
	//first check the val size is N_of_x
	if(indv->xreal.size() != N_of_x)
		return;
	// apply the torsion angles
	for(int i = 0; i < N_of_x; i++)
	{
		mol_.apply_rotor(i, torsionStep * indv->xreal[i]);
	}
	double vdwEnergy(0.), torsionEnergy(0.), eleEnergy(0.);
	if(MOGAParam_.FFType_ == FF_TAFF)
	{
		indv->energy = taff_.update_energy();
		vdwEnergy = taff_.get_vdw_energy();
		torsionEnergy = taff_.get_torsion_energy();
		eleEnergy = taff_.get_ele_energy();
	}
	else
	{
		indv->energy = mmff94_.update_energy();
		vdwEnergy = mmff94_.get_vdw_energy();
		torsionEnergy = mmff94_.get_torsion_energy();
		eleEnergy = mmff94_.get_ele_energy();
	}
	indv->gyration_radius = mol_.ComputGyrationRadius();

	// v1 (2026) objective design:
	//   f0 = conformation-dependent energy (torsion + vdw + ele). The full FF
	//       energy also contains stretch/bend/oop terms that are CONSTANT under
	//       dihedral rotation; including them flattens the Pareto front into a
	//       single epsilon box. Keep only terms that respond to torsion changes.
	//   f1 = aligned RMSD vs. input conformer -- diversity
	//   f2 = -gyration radius (optional) -- molecular extension
	//   f3 = -electrostatic energy (optional) -- polar/H-bond complement
	double confEnergy = torsionEnergy + vdwEnergy + eleEnergy;
	// 2 objectives
	if(MOGAParam_.NumObjects_ == 2)
	{
		indv->f[0] = confEnergy;
		indv->f[1] = mol_.minimizeRMSD();
	}
	// 3 objectives
	else if(MOGAParam_.NumObjects_ == 3)
	{
		indv->f[0] = confEnergy;
		indv->f[1] = mol_.minimizeRMSD();
		indv->f[2] = -mol_.ComputGyrationRadius();
	}
	// 4 objectives
	else if(MOGAParam_.NumObjects_ == 4)
	{
		indv->f[0] = confEnergy;
		indv->f[1] = mol_.minimizeRMSD();
		indv->f[2] = -mol_.ComputGyrationRadius();
		indv->f[3] = -eleEnergy;
	}
	mol_.reset();
	return ;
}

void MOGA::create_random_pop ()
{
	double ran(0.);
	for (int i = 0; i < MOGAParam_.PopSize_; i++)
	{
		population[i].xreal.resize(N_of_x);
		for (int j = 0; j < N_of_x; j++)
		{
			//ran = randomperc ();
			population[i].xreal[j] = rnd(infimumx, supremumx);
		}
		objective(&population[i]);
	}
}

/*This is the file used for crossover for Real Coded GA*/
// v1 (2026): torsions are periodic variables. A value maps to a dihedral angle
// via xreal*2.5 deg; the domain [infimumx, supremumx] spans exactly one full
// turn (144 steps x 2.5 deg = 360 deg). wrap_torsion maps any integer back into
// the canonical domain along the circle (shortest-arc equivalent).
inline int wrap_torsion(int v)
{
	const int period = supremumx - infimumx + 1;   // 144
	int r = (v - infimumx) % period;
	if (r < 0)
		r += period;
	return infimumx + r;
}

void MOGA::realcross(individual parent1, individual parent2)
{
	float betaq(0.), beta(0.), alpha(0.),expp(0.);
	int y1(0), y2(0), yu(0), yl(0), chld1(0), chld2(0);
	/*Loop over no of variables */
	for (int j = 0; j < N_of_x; j++)
	{
		/*Selected Two Parents */
		int par1 = parent1.xreal[j];
		int par2 = parent2.xreal[j];

		yl = infimumx;
		yu = supremumx;
		float rnd = randomperc ();

		/* Check whether variable is selected or not */
		if (rnd <= 0.5)
		{
			/*Variable selected */
			//ncross++;

			// v1 (2026): map par2 onto the shortest arc around par1 so that
			// crossover respects the periodicity of the dihedral variable
			// (e.g. 170 deg and -170 deg are only 20 deg apart on the circle).
#ifndef TORSION_LINEAR
			{
				int period = supremumx - infimumx + 1;   // 144
				int d = par2 - par1;
				if (d > period / 2)      d -= period;
				else if (d <= -period/2) d += period;
				// wrap back into [infimumx, supremumx] so the SBX bound terms
				// (y1-yl) and (yu-y2) stay non-negative
				par2 = wrap_torsion(par1 + d);
			}
#endif
			if (abs(par1 - par2) > 0)	// changed by Deb (31/10/01)
			{
				if (par2 > par1)
				{
					y2 = par2;
					y1 = par1;
				}
				else
				{
					y2 = par1;
					y1 = par2;
				}

				/*Find beta value */
				if ((y1 - yl) > (yu - y2))
				{
					beta = 1 + (2 * (yu - y2) / (y2 - y1));
					//printf("beta = %f\n",beta);
				}
				else
				{
					beta = 1 + (2 * (y1 - yl) / (y2 - y1));
					//printf("beta = %f\n",beta);
				}

				/*Find alpha */
				expp = MOGAParam_.n_distribution_c + 1.0;

				beta = 1.0 / beta;

				alpha = 2.0 - pow (beta, expp);

				if (alpha < 0.0)
				{
					cout<<"Warning: MOGA::realcross: Crossover Faild for individual "<<par1<<" "<<par2<<endl;
					continue;
				}

				rnd = randomperc ();

				if (rnd <= 1.0 / alpha)
				{
					alpha = alpha * rnd;
					expp = 1.0 / (MOGAParam_.n_distribution_c + 1.0);
					betaq = pow (alpha, expp);
				}
				else
				{
					alpha = alpha * rnd;
					alpha = 1.0 / (2.0 - alpha);
					expp = 1.0 / (MOGAParam_.n_distribution_c + 1.0);
					if (alpha < 0.0)
					{
					cout<<"Warning: MOGA::realcross: Crossover Faild for individual "<<par1<<" "<<par2<<endl;
					continue;
					}
					betaq = pow (alpha, expp);
				}

				/*Generating two children */
				chld1 = 0.5 * ((y1 + y2) - betaq * (y2 - y1));
				chld2 = 0.5 * ((y1 + y2) + betaq * (y2 - y1));

			}
			else
			{

				betaq = 1.0;
				y1 = par1;
				y2 = par2;

				/*Generation two children */
				chld1 = 0.5 * ((y1 + y2) - betaq * (y2 - y1));
				chld2 = 0.5 * ((y1 + y2) + betaq * (y2 - y1));

			}
			// v1 (2026): wrap children onto the circle instead of clamping,
			// so crossover can move across the +/-180 deg seam.
#ifdef TORSION_LINEAR
			if (chld1 < yl) chld1 = yl;
			if (chld1 > yu) chld1 = yu;
			if (chld2 < yl) chld2 = yl;
			if (chld2 > yu) chld2 = yu;
#else
			chld1 = wrap_torsion(chld1);
			chld2 = wrap_torsion(chld2);
#endif
		}
		else
		{

			/*Copying the children to parents */
			chld1 = par1;
			chld2 = par2;
		}
		newchild1.xreal[j] = chld1;
		newchild2.xreal[j] = chld2;
	}
	return;
}

void MOGA::real_mutate (individual * new_pop_ptr)
{
	for (int i = 0; i < N_of_x; i++)
	{
		double rnd = randomperc();
		int y(0);
		/*For each variable find whether to do mutation or not */
		if (rnd <= MOGAParam_.PMutation_)
		{
			y = new_pop_ptr->xreal[i];
			int yl = infimumx;
			int yu = supremumx;

			if (y > yl)
			{
				/*Calculate delta */
				double delta(0.);
				if ((y - yl) < (yu - y))
					delta = (y - yl) / (yu - yl);
				else
					delta = (yu - y) / (yu - yl);
				rnd = randomperc();

				double indi = 1.0 / (MOGAParam_.n_distribution_m + 1.0);
				double deltaq(0.);
				if (rnd <= 0.5)
				{
					double xy = 1.0 - delta;
					// v1 (2026): was 'int val = ...' which truncated the polynomial
					// mutation to 0/1 and destroyed the mutation distribution.
					double val = 2 * rnd + (1 - 2 * rnd) * (std::pow(xy, double(MOGAParam_.n_distribution_m + 1)));
					deltaq = pow(val, indi) - 1.0;
				}
				else
				{
					double xy = 1.0 - delta;
					double val = 2.0 * (1.0 - rnd) + 2.0 * (rnd - 0.5) * (pow (xy, double(MOGAParam_.n_distribution_m + 1)));
					deltaq = 1.0 - (pow (val, indi));
				}

				/*Change the value for the parent */
				//  *ptr  = *ptr + deltaq*(yu-yl);
				// Added by Deb (31/10/01)
				y = y + deltaq * (yu - yl);
				// v1 (2026): wrap onto the circle instead of clamping, so mutation
				// can cross the +/-180 deg seam.
#ifdef TORSION_LINEAR
				if (y < yl) y = yl;
				if (y > yu) y = yu;
				new_pop_ptr->xreal[i] = y;
#else
				new_pop_ptr->xreal[i] = wrap_torsion(y);
#endif
			}
			else
			{
				// v1 (2026): randomperc() returns [0,1); was assigned to int 'xy'
				// giving 0 always -- a real bug. Use a proper integer draw.
				int xy = (int)(randomperc() * (yu - yl + 1)) + yl;
				new_pop_ptr->xreal[i] = wrap_torsion(xy);
			}
		}
	}

	return;
}

void MOGA::compete ()
{
	int flag = 0, lg = 0, pp = 0;
	double rndp(0.);

	for (int i = 0; i < MOGAParam_.PopSize_; i++)
	{
		if (dom_check (newchild, population[i]) == 1)
		{
			population[i] = newchild;
			flag = 1;
			return;
		}
	}
	if (flag == 1)
		return;
	else
	{
		for (int i = 0; i < MOGAParam_.PopSize_; i++)
		{
			if (dom_check (population[i], newchild) == 1)
			{
				lg = -1;
				return;
			}
		}
		if (lg == -1)
			return;
		else
		{

			rndp = randomperc ();
			pp = (rndp < 1.0) ? (int) (rndp * MOGAParam_.PopSize_) : MOGAParam_.PopSize_ - 1;
			population[pp] = newchild;
			return;
		}
	}
	return;
}

int MOGA::tournament (int p1, int p2)
{
	if (dom_check (population[p1], population[p2]) == 1)
		return p1;
	else if (dom_check (population[p2], population[p1]) == 1)
		return p2;
	else
	{
		if (randomperc () > 0.5)
			return p1;
		else
			return p2;
	}
}

void MOGA::generate_replace ()
{
	int ar_parent, pop_parent, ch1_flg = 0, ch2_flg = 0;
	int pop_parent1, pop_parent2;
	double u = 0.0;
	double rnds = randomperc ();
	pop_parent1 = (rnds < 1) ? (int) (rnds * MOGAParam_.PopSize_) : MOGAParam_.PopSize_ - 1;
	rnds = randomperc ();
	pop_parent2 = (rnds < 1) ? (int) (rnds * MOGAParam_.PopSize_) : MOGAParam_.PopSize_ - 1;
	pop_parent = tournament (pop_parent1, pop_parent2);

	rnds = randomperc ();
	ar_parent = (rnds < 1) ? (int) (rnds * archive_size) : archive_size - 1;

	if (randomperc () <= MOGAParam_.PXover_)
	{
		realcross (population[pop_parent], archive[ar_parent]);
	}
	else
	{
		newchild1 = population[pop_parent];
		newchild2 = archive[ar_parent];
	}

	real_mutate (&(newchild1));
	real_mutate (&(newchild2));
	objective(&newchild1);
	objective(&newchild2);
}

// The Update function, called to update the archive in every iteration
int MOGA::update()
{
	int d(0), ch_flag(0), accept_flag(0), no_dom = 0, same_box = 0;
	double da = 0, dc = 0;
	double child_dist(0.), archive_dist(0.);
	vector<double> middle(MOGAParam_.NumObjects_, 0.); 

	box_func (&newchild);
	/*  added by xfliu on 20080417*/
	// modified by xfliu on 20080428
	// modefied by xfliu on 20080508
	for (int i = 0; i < archive_size; i++)
	{
		if(archive[i].energy < MinEnergy_)
			MinEnergy_ = archive[i].energy;
	}
	//float Ew = 1000. + 100.0 * (N_of_x - 1);
	float Ew = MOGAParam_.EnergyCutoff_ + 0.5 * (N_of_x - 1);
	Ecutoff_ = MinEnergy_ + Ew;
	for (int i = 0; i < archive_size; i++)
	{

		if (archive[i].energy > Ecutoff_)
		{
			for (int j = i; j < archive_size - 1; j++)
			{
				archive[j] = archive[j + 1];
			}
			archive_size--;
		}
	}
	//modified by xfliu on 20080507
	if(newchild.energy > Ecutoff_)
		return -1;

	for (int i = 0; i < archive_size; i++)
	{
		box_func (&archive[i]);

		if (box_dom (archive[i], newchild) == 1)
		{
			accept_flag = -1;
		}
	}

	if (accept_flag == -1)
		return -1;
	else
	{
		for (int i = 0; i < archive_size; i++)
		{
			no_dom = 0;
			if (box_dom (newchild, archive[i]) == 1)
			{
				no_dom++;
				for (int j = i; j < archive_size - 1; j++)
				{
					archive[j] = archive[j + 1];
				}
				archive_size--;
				d = 1;
				i = 0;
			}
		}
		if (d == 1)
		{
			archive_size++;
			archive[archive_size - 1] = newchild;
			accept_flag = 1;
		}
		else
		{
			for (int i = 0; i < archive_size; i++)
			{
				if (same_box_check (newchild, archive[i]) == 1)
				{
					same_box++;
					if (dom_check (newchild, archive[i]) == 1)
					{
						archive[i] = newchild;
						accept_flag = 1;
					}
					else
					{
						if (dom_check (newchild, archive[i]) == 0)	//changed
							accept_flag = -1;
						else
						{
							for (int p = 0; p < MOGAParam_.NumObjects_; p++)
							{
								middle[p] = newchild.box[p] * MOGAParam_.EPSILON_[p];	//changed
							}
							child_dist = distance (newchild, middle);
							archive_dist = distance (archive[i], middle);

							if (child_dist < archive_dist)
							{
								archive[i] = newchild;
								accept_flag = 1;
							}
							else
								accept_flag = -1;
						}
					}
				}
			}
		}
		if (d == 0 && accept_flag == 0 && same_box == 0)
		{
			archive_size++;
			archive[archive_size - 1] = newchild;

		}

	}

	return accept_flag;
}

// returns 0 if not accepted and 1 if the child is accepted
int MOGA::con_update ()
{
	int flg = 0, d = 0;

	// First Case : the new Child is a feasible solution,
	// i.e. it has zero constraint violation
	if (newchild.cv <= EPSILON)
	{
		// removing all infeasible solutions from the archive
		for (int i = 0; i < archive_size; i++)
		{
			if (archive[i].cv > EPSILON)
			{
				for (int j = i; j < archive_size - 1; j++)
				{
					archive[j] = archive[j + 1];
				}
				archive_size--;
				d = 1; // set flag to indicate that some archive members have been deleted
			}
		}
		// insert the new child in the usual manner into the archive after
		// removing the infeasible points.
		flg = update ();
	}
	else
	{
		// Second Case : New child is infeasible
		for (int k = 0; k < archive_size; k++)
		{
			if (newchild.cv > archive[k].cv)
			{
				flg = -1;  // reject the child, as there is a better solution in the archive
			}
		}
		if (flg == 0)
		{
			// removing all those members whose constraint violation
			// is higher than that of the new child
			for (int k = 0; k < archive_size; k++)
			{
				if (newchild.cv < archive[k].cv)
				{
					for (int j = k; j < archive_size - 1; j++)
					{
						archive[j] = archive[j + 1];
					}
					archive_size--;
					d = 1; // set flag to indicate that some archive members have been deleted
				}
			}
		}

		if (flg == 0)
		{
			archive_size++;
			archive[archive_size - 1] = newchild;
			// new child inserted into archive
			flg=1; // flag to indicate this
		}
	}
	return flg;
}

// Finds the distance between corner[] and between the decision variable
// coordinates of the individual
double MOGA::distance (individual ind1, vector<double>point)
{
	double dist = 0.0;
	for (int i = 0; i < MOGAParam_.NumObjects_; i++)
	{
		dist += pow ((ind1.f[i] - point[i]), 2);
	}
	dist = sqrt (dist);
	return (dist);
}

// Check for domination (usual sense)
int MOGA::dom_check (individual ind1, individual ind2)
{
	int flag = 1;
	int elag = 0;
	for (int j = 0; j < MOGAParam_.NumObjects_; j++)
	{
		if ((flag == 1) && (ind1.f[j] <= ind2.f[j]))
			flag = 1;
		else
			flag = 0;

		if ((ind1.f[j] < ind2.f[j]) || (elag == 1))
			elag = 1;
	}
	if ((flag == 1) && (elag == 1))
		return 1;
	else
		return 0;
}

//to check whether 2 individuals are in the same box or not
int MOGA::same_box_check (individual ind1, individual ind2)
{
	int flag = 1;
	for (int i = 0; i < MOGAParam_.NumObjects_; i++)
	{
		if ((flag == 1) && (ind1.box[i] == ind2.box[i]))
		{
			flag = 1;
		}
		else
		{
			flag = 0;
			break;
		}
	}
	return flag;
}

// function to check for Box domination
int MOGA::box_dom (individual ind1, individual ind2)
{
	int flag = 1;
	int elag = 0;
	for (int j = 0; j < MOGAParam_.NumObjects_; j++)
	{
		if ((flag == 1) && (ind1.box[j] < ind2.box[j]))
			flag = 1;
		else
			flag = 0;

		if ((ind1.box[j] < ind2.box[j]) || (elag == 1))
			elag = 1;
	}
	if ((flag == 1) && (elag == 1))
		return 1;
	else
		return 0;
}

int MOGA::strict_dom_check (int m, int n)
{
	int flag = 1;
	for (int p = 0; p < MOGAParam_.NumObjects_; p++)
	{
		if ((flag == 1) && (population[m].f[p] < population[n].f[p]))
			flag = 1;
		else
		{
			flag = 0;
			break;
		}
	}
	return flag;
}    

// Create the initial archive
// Consists of the non-dominated (strict domination) members
// of the population
void MOGA::create_archive ()
{
	int flag = 0;
	for (int i = 0; i < MOGAParam_.PopSize_; i++)
	{
		flag = 0;
		for (int j = 0; j < MOGAParam_.PopSize_; j++)
		{
			if (strict_dom_check (j, i) == 1)
				flag = -1;
		}
		if (flag == 0)
		{
			archive[archive_size] = population[i];
			archive_size++;
		}
	}
	return;
}

// Box funtion definitons, here done for minimization case
void MOGA::box_func (individual * ind1)
{
	for (int i = 0; i < MOGAParam_.NumObjects_; i++)
	{
		ind1->box[i] = floor ((ind1->f[i] / MOGAParam_.EPSILON_[i]));
	}
}

void MOGA::remove_redundency_archive()
{
	float RMSD_TOL = MOGAParam_.rmsdScaleFactor_;
	//float RMSD_TOL = MOGAParam_.rmsdScaleFactor_ * sqrt(float(1 + N_of_x));
	archive.resize(archive_size);
	vector<bool> del_archive(archive_size, false);
	for (int ar = 0; ar < archive_size; ar++)
	{
		if(del_archive[ar])
			continue;
		for(int i = 0; i < N_of_x; i++)
		{
			mol_.apply_rotor(i, torsionStep * archive[ar].xreal[i]);
		}
		mol_.bk_position();
		for (int ar_1 = ar + 1; ar_1 < archive_size; ar_1++)
		{
			if(del_archive[ar_1])
				continue;
			for(int i = 0; i < N_of_x; i++)
			{
				mol_.apply_rotor(i, torsionStep * archive[ar_1].xreal[i]);
			}
			float rmsd = mol_.minimizeRMSD();
			if(rmsd < RMSD_TOL)
				if((archive[ar].f[0] + archive[ar].f[1]) <= (archive[ar_1].f[0] + archive[ar_1].f[1]))
				{
					del_archive[ar_1] = true;
				}
				else
				{
					del_archive[ar] = true;
					break;
				}
		}
	}
	int i = 0;
	for(vector<individual>::iterator it = archive.begin(); it != archive.end(); ++i)
	{
		if(del_archive[i])
			it = archive.erase(it);
		else
			++it;
	}
	archive_size = archive.size();
}
