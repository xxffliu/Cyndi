#include "../include/Cluster.h"

LeaderClustering::LeaderClustering(double dist) :_leaders(0), _nbrClusters(0), _threshold(dist){}

void LeaderClustering::cluster(priority_queue<MOL>& confqueue)
{
	// first pop all the conformers into a vector
	vector<MOL> points;
	while(!confqueue.empty())
	{
		points.push_back(confqueue.top());
		confqueue.pop();
	}
	
	// initialisation by giving each data point a random cluster number;
	unsigned int n = points.size();	
	std::vector<unsigned int> labels(n);
	std::vector<unsigned int> index(n);
	for ( unsigned int i=0; i<n; ++i )
	{	
		labels[i] = 0;
		index[i] = i;
	}

	// permutate the indices
	//srand(rand());
	//random_shuffle(index.begin(),index.end());
	
	// make the first element the first leader
	Leader newCenter;
	newCenter.center = points[index[0]];
	newCenter.index = index[0];
	newCenter.numMembers = 1;
	newCenter.members.push_back(points[index[0]]);
	// first clear current set of leaders
	if ( !_leaders.empty() )
		_leaders.clear();

	// add first one
	_leaders.push_back(newCenter);
	_nbrClusters = 1; 
	labels[index[0]] = 0; // set label of first data point
	
	double mseDistance(0.0);	
	for ( unsigned int i=1; i<n; ++i )
	{
		// std::cerr << " " << index[i];
		//vector3 r(points[index[i]]);
		
		double dClosest(10000.0);
		unsigned int c(_leaders.size());
		for ( unsigned int j = 0; j<_leaders.size(); ++j )
		{
			double d = _leaders[j].center.minimizeRMSD(points[index[i]]);
			if ( d < dClosest )
			{
				dClosest = d; 
				c = j; //store current cluster number
			}
		}
		
		if ( dClosest < _threshold ) // add to selected cluster
		{
			labels[index[i]] = c;
			_leaders[c].numMembers += 1;
			_leaders[c].members.push_back(points[index[i]]);
			mseDistance += (dClosest * dClosest);
		}
		else // create a new cluster
		{
			Leader newCenter;
			newCenter.center = points[index[i]];
			newCenter.index = index[i];
			newCenter.numMembers = 1;
			newCenter.members.push_back(points[index[i]]);
			_leaders.push_back(newCenter);
			_nbrClusters++;
			labels[index[i]] = _leaders.size() - 1;
			//std::cerr << _nbrClusters << "\r";
		}
	}
	
	// now we have clustered leaders, we pop all the cluster centers back into the confqueue
	for(int i = 0; i < _nbrClusters; ++i)
	{
		confqueue.push(_leaders[i].center);
	}
	//std::cerr << "MSE = " << mseDistance/n << std::endl;
	//std::cerr << "Nbr. Clusters = " << _nbrClusters << std::endl;
	//return _leaders;
	return;
}