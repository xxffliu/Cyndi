/*
 *  LeaderClustering.h
 *
 *  Created by Gert Thijs on 09/08/06.
 *  Copyright 2006 Silicos NV. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     + Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     + Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY SILICOS NV AND CONTRIBUTORS ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL SILICOS NV AND CONTRIBUTORS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
		* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
 *
 */

#ifndef  LEADERCLUSTERING__
#define  LEADERCLUSTERING__

#include <queue>
#include "Mol.h"
#include <vector>
#include <algorithm>
using namespace std;
/**
	\class LeaderClustering
	\brief Implenentation of leader clustering
 \ingroup _clustering_
 
	This class implements the leader clustering algorithm, which is a one pass 
	clustering algorithm. It starts by randomly selecting a data point and makes 
	it a cluster center (or leader). In the next steps a new data point is selected
	and the distance to the closest cluster is computed. If the distance to this 
	closest leader is smaller than a predefined threshold the point is added to this 
	cluster otherwise a new cluster is started of which this point will be the leader.
	This procedure is repeated until al points have been processed.
 
	The main advantage of this method is its simplicity and its speed of computation. 
	The major drawbacks of this approach are the dependance on the order of the data points
	and the need for the cutoff value to be defined.
 */
class LeaderClustering
{
public:
	struct Leader
	{
		MOL center;
		unsigned int index;
		unsigned int numMembers;
		vector<MOL> members;
			
		Leader() : center(), index(0), numMembers(0), members() {};
		void clear()
		{
			center.clear();
			index = 0;
			numMembers  = 0;
			members.clear();
		}
	};

private:
		
	/**
		\struct Leader 
		\brief Local representation of a cluster leader
		\param coord Vector with the coordinates of the leader
		\param index Cluster number of the leader
		\param members Number of data points assigned to this cluster
	 */
	vector<Leader> _leaders;   ///< Storage of cluster centers or leaders
	
	unsigned int _nbrClusters;      ///< Number of clusters found
	unsigned int _nbrVar;           ///< dimension of the data matrix
	double _threshold;              ///< threshold on the distance to start a new cluster 

public:
	
	/**
		\brief Constructor
		\param p Reference to a predefined set of parameters
	 */
	LeaderClustering(double dist);
	
	/**
		\brief Destructor
	 */
	~LeaderClustering()
	{
		_leaders.clear();
	}
	
	/**
		\brief Main clustering step
	
		\param data Reference to a Matrix that holds the data points
		\return Vector with the cluster labels of all data points
		 */
	void cluster(priority_queue<MOL>& confqueue);
	
	/**
		\brief Clustering step that starts from the currently defined leaders
	 
		\param data Reference to a Matrix that holds the data points
		\return A Vector with the cluster labels of all data points
	 */
	vector<unsigned int> update(const vector<vector3>& points);
	
	
	/**
		\brief Get the current number of clusters
	 */
	inline unsigned int getNbrOfClusters() { return _nbrClusters; };
	
	/**
		\brief Get the i-th leader as a SiMath::Vector 
	 */
	inline vector3 getLeader(const unsigned int i);

	// get the geometrical center of the i-th cluster
	inline vector3 getClusterCenter(const unsigned int i);
	
	/**
		\brief Reset the clustering to an empty set of leaders.
	 */
	void reset();

	/**
		\brief Set the threshold on the distance below which a data point is added to a cluster
	 
		\param d Threshold
		\throws LeaderClusteringError if the threshold is smaller than 0.
	 */
	inline void setThreshold(double d)
	{
		_threshold = d;
	}
};


#endif  LEADERCLUSTERING__
