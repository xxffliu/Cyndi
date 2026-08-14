#include "../include/Cyndi.h"

int main(int argc, char* argv[])
{
	// the wall clock
	clock_t start, termination;
	// starting timing
	start = clock();
	CyndiOptionParser(argc, argv);
	execuateCONGEN();
	termination = clock();
	float ElapseTime = (double)(termination - start) / CLOCKS_PER_SEC;
	cout<<"Elapsed Time: "<<ElapseTime<<endl;
	return EXIT_SUCCESS;
}

int RemoveRedundantConf(deque<MOL>& confvec, int num_rot)
{
	float RMSD_TOL = MOGAParam_.rmsdScaleFactor_;
	int num_conf = confvec.size();
	vector<bool> del_conf(num_conf, false);
	MOL mol_1, mol_2;
	for (int ar = 0; ar < num_conf; ar++)
	{
		if(del_conf[ar])
			continue;
		mol_1 = confvec[ar];
		for (int ar_1 = ar + 1; ar_1 < num_conf; ar_1++)
		{
			if(del_conf[ar_1])
				continue;
			mol_2 = confvec[ar_1];
			float rmsd = mol_1.minimizeRMSD(mol_2);
			//float rmsd = mol_1.get_RMSD(mol_2);
			if(rmsd < RMSD_TOL)
				if( mol_1.get_energy() <= mol_2.get_energy())
				{
					del_conf[ar_1] = true;
				}
				else
				{
					del_conf[ar] = true;
					break;
				}
		}
	}
	int i = 0;
	for(deque<MOL>::iterator it = confvec.begin(); it != confvec.end(); ++i)
	{
		if(del_conf[i])
			it = confvec.erase(it);
		else
			++it;
	}
	return confvec.size();
}

void execuateCONGEN()
{
	cout<<"Name"<<setw(8)<<"Rot Bonds"<<setw(24)<<"Confs(before filtered)"<<setw(24)<<"Confs(after filtered)"<<setw(8)<<"time(s)"<<endl;
	MOL mol_, tmpmol;
	TAFF MOGATaff;
	MMFF94 MOGAmmff94;
	CGMinimizer cgm;
	string MolFile = MOGAParam_.InputFile_;
	vector<MOL> FailedVector;
	// herer we want to discard the delimiter such as "/" and ".", so we only take the substring between the last
	// ".'(extension mark) and "/"(the directory mark)
	// we use reversed iterator to locate
	string::reverse_iterator rpoint = find(MolFile.rbegin(), MolFile.rend(), '.');
	string::reverse_iterator rsplash = find(MolFile.rbegin(), MolFile.rend(), '/');
	// here to tackle the MS format directory (using '\' as delimiter)
	if(MolFile.find('\\') != string::npos)
		rsplash = find(MolFile.rbegin(), MolFile.rend(), '\\');
	string MolName = string(rsplash.base(), rpoint.base()-1);
	MOL2IO infile(MolFile, "in");
	MOL2IO outfile;
	outfile.open(MOGAParam_.OutputFile_, "out");
	int Counter(0), failCounter(0);
	int GlobalIdx = 0;   // v3 (2026): global molecule index (0-based)
	while (infile.read(&mol_))
	{
		// v3 (2026): multiprocess batch support -- skip molecules before
		// StartIndex_ and stop after MaxMolecules_ (0 = to EOF).
		if(GlobalIdx < MOGAParam_.StartIndex_)
		{
			GlobalIdx += 1;
			mol_.clear();
			continue;
		}
		if(MOGAParam_.MaxMolecules_ > 0 && Counter >= MOGAParam_.MaxMolecules_)
			break;
		Counter += 1;
		mol_.initialize();
		string name = mol_.get_name();
		int NumConfBelowCutoff = 0;
		float CurrGlobalMinEnergy = 1000.0;
		priority_queue<MOL> conf;
		double ElapseTime(0.);

		if(MOGAParam_.KeepInputConformer_)
			//conf.push_back(mol_);
			conf.push(mol_);
		if(MOGAParam_.OptimizeInputConformer_)
		{
			if(MOGAParam_.FFType_ == FF_TAFF)
			{
				MOGATaff.setup(mol_);
				cgm.setup(MOGATaff);
				cgm.minimize(MOGAParam_.MaxNumIteration_);
				MOGATaff.update_energy();
			}
			else
			{
				MOGAmmff94.setup(mol_);
				//MOGAmmff94.remove_component("MMFF94 Str_Bend");
				cgm.setup(MOGAmmff94);
				cgm.minimize(MOGAParam_.MaxNumIteration_);
				MOGAmmff94.update_energy();
			}
		}
		mol_.bk_position();
		if(mol_.get_num_of_rot_bonds() < 1)
		{
			cout<<name<<" "<<"Warning: Cyndi:excecuateCONGEN(): The number of rotational bonds of input molecule is ZERO! Skipping MOGA conformer generation"<<endl;
			outfile.write(mol_);
			continue;
		}
		// here we skip some huge molecules with over 30 rotatable bonds
		else if(mol_.get_num_of_rot_bonds() > 30)
		{
			cout<<name<<" "<<"Warning: Cyndi:excecuateCONGEN(): Too many rotatable bonds! Skipping MOGA conformer generation"<<endl;
			outfile.write(mol_);
			continue;
		}
		else
		{
			vector<MOL> archive;
			vector<Conformer> TorVec,tmp_tor_vec;
			TorVec.reserve(1000);
			// the wall clock
			clock_t start, termination;
			// starting timing
			start = clock();
			MOGA congen;
		// v3 (2026): tell the MOGA which global molecule index this is, so its
		// RNG seed = BasicSeed_ + 0.001 * GlobalIdx -- deterministic per
		// molecule and identical whether run solo or as a batch chunk.
		congen.set_global_idx(GlobalIdx);
			bool success = congen.setup(mol_);
			if(!success)
			{
				FailedVector.push_back(mol_);
				failCounter += 1;
				continue;
			}
			for(int i = 0; i < MOGAParam_.MaxNumRun_; i++)
			{
				tmp_tor_vec = congen.execuateMOGA();
				copy(tmp_tor_vec.begin(), tmp_tor_vec.end(), back_inserter(TorVec));
			}
			termination = clock();
			ElapseTime = (double)(termination - start) / CLOCKS_PER_SEC;
			int num_conf = TorVec.size();

			if(num_conf == 0)
			{
				cout<<"MOGA conformation generation failed for "<<name<<endl;
				FailedVector.push_back(mol_);
				failCounter += 1;
				continue;
			}

			// Setup the force field and minimizer ONCE per molecule. The FF binds to
			// atom pointers and reads coordinates on the fly during minimize(), so
			// repeated setup() per conformer only re-reads parameter files and
			// re-assigns atom types (wasted work). Optimization by Hermes, 2026-08.
			if(MOGAParam_.OptimizeConformer_)
			{
				if(MOGAParam_.FFType_ == FF_TAFF)
				{
					MOGATaff.setup(mol_);
					cgm.setup(MOGATaff);
				}
				else
				{
					MOGAmmff94.setup(mol_);
					cgm.setup(MOGAmmff94);
				}
			}

			for(int i = 0 ;i < num_conf; i++)
			{
				for(int j = 0; j < mol_.get_num_of_rot_bonds(); j++)
					mol_.apply_rotor(j, TorVec[i].torsions[j] * 2.5);

				// minimize the confomers if required
				if(MOGAParam_.OptimizeConformer_)
				{
					if(MOGAParam_.FFType_ == FF_TAFF)
					{
						cgm.minimize(MOGAParam_.MaxNumIteration_);
						MOGATaff.update_energy();
					}
					else
					{
						cgm.minimize(MOGAParam_.MaxNumIteration_);
						MOGAmmff94.update_energy();
					}
					mol_.minimizeRMSD();
					//debug
					//cout<<MOGATaff.get_vdw_energy()<<setw(16)<<MOGATaff.get_torsion_energy()<<setw(16)<<mol_.get_rmsd()<<setw(16)<<mol_.ComputGyrationRadius()<<endl;
				}
				else
				{
					mol_.set_energy(TorVec[i].TotalEnergy);
					mol_.set_rmsd(TorVec[i].rmsd);
					//cout<<TorVec[i].VDWEnergy<<setw(16)<<TorVec[i].TorsionEnergy<<setw(16)<<TorVec[i].TotalEnergy<<setw(16)<<TorVec[i].rmsd<<setw(16)<<TorVec[i].GyrationRdius<<endl;
				}
				if(CurrGlobalMinEnergy > mol_.get_energy())
				{
					CurrGlobalMinEnergy = mol_.get_energy();
				}
				archive.push_back(mol_);
				mol_.reset();
			}
			// remove conformers according to energy cutoff
			for(vector<MOL>::iterator mit = archive.begin(); mit != archive.end(); ++mit)
			{
				if(mit->get_energy() <= (CurrGlobalMinEnergy + MOGAParam_.EnergyCutoff_ + 0.5 * mit->get_num_of_rot_bonds()))
				{
					mit->set_fitness(-mit->get_energy());
					// rename the name of the molecule according to increment number
					//mol_.set_name(name + "_" + string(itoa(i, buffer, 10)));
					// push the conformers into the stack
					conf.push(*mit);
					NumConfBelowCutoff += 1;
				}
			}
			//sort(conf.begin(), conf.end());
			LeaderClustering lc(MOGAParam_.rmsdScaleFactor_);
			lc.cluster(conf);
			NumConfBelowCutoff = conf.size() > MOGAParam_.MaxNumConformer_ ? MOGAParam_.MaxNumConformer_ : conf.size();
			//NumConfBelowCutoff = RemoveRedundantConf(conf, mol_.get_num_of_rot_bonds());
			/*if(conf.size() > MOGAParam_.MaxNumConformer_)
			{
				conf.resize(MOGAParam_.MaxNumConformer_);
				NumConfBelowCutoff = MOGAParam_.MaxNumConformer_;
			}*/
			//cout<<"Done."<<endl<<num_conf<<" conformers generated."<<endl;
			//cout<<"Cyndi reserved "<<NumConfBelowCutoff<<" conformers because of energy cutoff"<<endl;
			cout<<name<<setw(8)<<mol_.get_num_of_rot_bonds()<<setw(24)<<num_conf<<setw(24)<<NumConfBelowCutoff<<setw(8)<<ElapseTime<<endl;
		}
		//debug
		//cout<<"VDW Energy"<<setw(16)<<"Torsion Energy"<<setw(16)<<"Total Energy"<<setw(16)<<"Rmsd"<<setw(16)<<"GyrationRadius"<<endl;
		// now dumping the molecules
		if(conf.empty())
			FailedVector.push_back(mol_);
		else
		{
			int counter = 1;
			while(!conf.empty())
			{
				//outfile.write(conf.front());
				//conf.pop_front();
				if(counter <= NumConfBelowCutoff)
				{
					outfile.write(conf.top());
					counter += 1;
				}
				conf.pop();
			}
		}
		GlobalIdx += 1;   // v3 (2026): advance the global index for the next molecule
	}
	if(!FailedVector.empty())
	{
		MOL2IO failfile("failed.mol2", "out");
		for(vector<MOL>::iterator it = FailedVector.begin(); it != FailedVector.end(); ++it)
			failfile.write(*it);
	}
	cout<<"Number of Molecules Processed: "<<Counter<<endl;
	cout<<"Number of Molecules Failed: "<<failCounter<<endl;
	mol_.clear();
	tmpmol.clear();
	return;
}
