/*
 * CiGRAM - Circular Gaussian Random grAph Model
 * Copyright (C) 2017  James Gilbert
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 3
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details:
 * http://www.gnu.org/licenses/gpl.txt
 */
#include "generate_graph.h"

/*
 * Rewire procedure
 */

bool comparator (const std::pair<double, int>  &l, const std::pair<double, int>  &r) {
		return l.first < r.first;
}

void rewire(Graph &G,
			std::vector<double> &theta,
			std::vector<double> &beta,
			std::vector<int> &component,
			int &num_components,
			boost::mt19937 &generator,
			double a,
			const std::vector< std::vector<int> > &nodeCommunities = std::vector< std::vector<int> >())
{
	/*
	* Sort edges in connected component by probability of existence
	* Wire up components
	* Remove least probable edges that don't create disconnected components
	*/
	std::vector< std::vector<int> > componentNodes (num_components);
	std::vector<int> componentSize (num_components);
	std::vector< std::vector<double> > compBeta (num_components);


	for (int i=0; i < num_components; ++i) {
		componentSize[i] = 0;
		componentNodes[i] = std::vector<int> ();
	}

	for (size_t i=0; i < component.size(); ++i) {
		componentNodes[component[i]].push_back(i);
		++componentSize[component[i]];
	}

	for (size_t i=0; i < componentSize.size(); ++i) {
		compBeta[i] = std::vector<double> (componentSize[i]);
	}

	// index of the giant component
	int giantC = distance(componentSize.begin(), max_element(componentSize.begin(), componentSize.end()));
	int n_i = 0;
	int n_j = 0;

	// store edges in sortable format
	std::vector< std::pair<int, int> > edgeSet (num_edges(G), std::make_pair(-1, -1) );
	std::vector< std::pair<double, int> > edgeScore (num_edges(G), std::make_pair(0.0, -1));
	int i = 0;

	graph_traits<Graph>::edge_iterator ei, ei_end;
	for (tie(ei, ei_end) = edges(G); ei != ei_end; ++ei) {
		edgeSet[i] = std::make_pair(source(*ei, G), target(*ei, G));
		edgeScore[i] = std::make_pair(beta[edgeSet[i].first] + beta[edgeSet[i].second], i);
		++i;
	}

	for (int i=0; i < num_components; ++i) {
		// for each connected component that is not maximal, select node
		// connect to node in largest component with selection matrix
		if (i == giantC)
			continue;

		for (size_t n=0; n < componentNodes[i].size(); ++n) {

			compBeta[i][n] = beta[n];
		}

		boost::random::discrete_distribution<> dist(compBeta[i]);
		n_i = componentNodes[i][dist(generator)];


		std::vector<double> selectionScoresj = selectionScores(theta[n_i], componentNodes[giantC], beta, theta, a);

		boost::random::discrete_distribution<> dist_j(selectionScoresj);
		n_j = componentNodes[giantC][dist_j(generator)];

		add_edge(n_i, n_j, G);

	}

	// Remove num_components - 1 edges
	// rank edges by probability and iterate through in this order
	// Use shortest paths method to check if edge can be removed

	int removals = 0;
	std::sort(edgeScore.begin(), edgeScore.end(), comparator);
	IndexMap index = get(boost::vertex_index, G);

	// This rewiring procedure is extremely costly, worst case a graph contains no cycles in which case every edge is tested
	for (uint i=0; i < edgeScore.size(); ++i) {

		if (removals == num_components - 1)
			break;

		n_i = edgeSet[edgeScore[i].second].first;
		n_j = edgeSet[edgeScore[i].second].second;

		// We don't want to remove edges between communities
		if (nodeCommunities.size()) {
			std::vector<int> intersection;
			std::set_intersection (nodeCommunities[n_i].begin(), nodeCommunities[n_i].end(),
											nodeCommunities[n_j].begin(), nodeCommunities[n_j].end(),
											std::back_inserter(intersection));

			if (intersection.size()){
				continue;
			}
		}

		// check to see if there is an alternative path
		remove_edge(n_i, n_j, G);

		// dfs until we find the edge pair or remove
		Graph::adjacency_iterator neighbourIt, neighbourEnd;

		std::set<int> visited;
		std::stack<int> toVisit;
		bool pathExists = false;
		int j, top;

		visited.insert(n_i);

		tie(neighbourIt, neighbourEnd) = adjacent_vertices(n_i, G);
		for (; neighbourIt != neighbourEnd; ++neighbourIt)
			toVisit.push(index[*neighbourIt]);

		// either we find our node and the search ends or we traverse every traverseable vertex
		while (!toVisit.empty()) {
			top = toVisit.top();
			tie(neighbourIt, neighbourEnd) = adjacent_vertices(top, G);
			for (; neighbourIt != neighbourEnd; ++neighbourIt) {

				j = index[*neighbourIt];
				if (j == n_j) {
					pathExists = true;
					++removals;
					break;
				}

				if (visited.find(j) == visited.end()) {
					visited.insert(j);
					toVisit.push(j);
				}
			}

			if (pathExists)
				break;

			toVisit.pop();
		}

		if (!pathExists) {
			add_edge(n_i, n_j, G);
		}
	}
}

/*
 * number of edges a node will be part of, sampling with replacement
 * TODO: move to sample class.
 */
void setFirstFreqCounter(int num, std::vector<int> &freqCounter,  std::vector<int> &n_degree, std::vector<double> &beta, boost::mt19937 &generator)
{
	int k = 0;
	int N = beta.size();
	int freeCount = N;
	int i;
	boost::random::discrete_distribution<> dist(beta);
	
	boost::variate_generator< boost::mt19937&, boost::random::discrete_distribution<> >
        generateRandom(generator, dist);

	while (k < num && freeCount > 0) {
		i = generateRandom();
		
		if (freqCounter[i] + n_degree[i] < N - 1) {
			++freqCounter[i];
			++k;
		} else{
			
			if (beta[i] > 0){
				--freeCount;
				beta[i] = 0;
				boost::random::discrete_distribution<> dist(beta);
				boost::variate_generator< boost::mt19937&, boost::random::discrete_distribution<> >
					generateRandom(generator, dist);
			}
			else{
				// This point means that the random generator is fubar
				break;
			}
		}
	}
}


int getNodeDegree(int i, Graph &G){
	int deg_count = 0;
	Graph::adjacency_iterator vi, vi_end;
	for (tie(vi, vi_end) = adjacent_vertices(i, G); vi != vi_end; ++vi){
		++deg_count;
	}

	return deg_count;
}

/*
 * create a subgraph with desired params
 *
 * nodes should be a subset of nodes to fill
 * theta should be of the same length as nodes, each index corresponding to the position of n
 */
Graph fillGraph(int N,
				int m,
				std::vector<double> &theta,
				std::vector<double> &beta,
				std::vector<int> &n_degree,
				double sigmaR,
				double a,
				boost::mt19937 &generator,
				int minDegree=1,
				int formConnected=1
			   )
{

	Graph subGraph(N);
	
	// No need to sample if we're building a fully dense graph, just add every edge
	if ( (2.0 * m)/(N * (N - 1)) == 1.0) {
		for (int i = 0; i < N; ++i) {
			for (int j = i; j < N; ++j) {
				if (j != i) {
					add_edge(i, j, subGraph);
					n_degree[i]++;
					n_degree[j]++;
				}
			}
		}
		return subGraph;
	}

	std::vector<int> nodes (N);
	
	for (int i=0; i < N; ++i) {
		n_degree[i] = 0;
		nodes[i] = i;
	}

	if (m > N - 1) {
		// ensure all nodes have at least one edge
		for (int i=0; i < N; ++i) {
			
			if (getNodeDegree(i, subGraph) < minDegree) {

				std::vector<int> available_population;
				// Add nodoes that aren't adjacent to the avialble population
				for (size_t t=0; t < nodes.size(); t++){
					int j = nodes[t];
					if (i != j)
						available_population.push_back(t);
				}

				if (available_population.size() == 0)
					break;
				
				std::vector<double> selectionScoresj = selectionScores(theta[i], available_population, beta, theta, a) ;

				boost::random::discrete_distribution<> dist_j(selectionScoresj);
				SampleARes sampler(available_population, selectionScoresj);
				std::vector<int> sample = sampler.getSamplesWRS(minDegree - getNodeDegree(i, subGraph), generator);
				for (size_t t=0; t < sample.size(); ++t) {
					int j = sample[t];
					add_edge(i, j, subGraph);
				}
			}
		}
	}

	for (int i=0; i < N; ++i) {
		n_degree[i] = getNodeDegree(i, subGraph);
	}
	
	uint fsum = 0;
	int assign = m - num_edges(subGraph);
	std::vector<int> freqCounter (N);
	for (int i=0; i < N; ++i) {
		 freqCounter[i] = 0;
	}
	
	/* Pain is caused by heterogentiy */
	while(fsum < m - num_edges(subGraph)){
		setFirstFreqCounter(assign, freqCounter, n_degree, beta, generator);
		
		fsum = 0;
		double bsum = 0.0;
		int freeCount = N;
		for (int i=0; i < N; ++i){
			if(n_degree[i] + freqCounter[i] >= N - 1){
				beta[i] = 0;
				--freeCount;
			}
			bsum += beta[i];
			fsum += freqCounter[i];
		}
		
		if(bsum < 1/N) {
			for (int i=0; i < N; ++i){
				if(n_degree[i] + freqCounter[i] < N - 1){
					beta[i] = 1/freeCount;
				}
			}
			bsum = 1.0;
		}
		
		
		for (int i=0; i < N; ++i){
			if(n_degree[i] + freqCounter[i] < N - 1){
				beta[i] /= bsum;
			}
		}

	}
		
	int reassign;
	int assigned = 0;
	int fail_count = 0;
	// fill graph with edges
	for (int i=0; i < N; ++i) {
		
		if (freqCounter[i] == 0)
			continue;

		reassign = 0;
		assigned = 0;
		std::vector<int> available_population;
		// Add nodoes that aren't adjacent to the avialble population
		for (size_t t=0; t < nodes.size(); t++){
			int j = nodes[t];
			if (i != j && !edge(i, j, subGraph).second )
				available_population.push_back(j);
		}

		std::vector<double> selectionScoresj = selectionScores(theta[i], available_population, beta, theta, a);

		SampleARes sampler(available_population, selectionScoresj);
		std::vector<int> sample = sampler.getSamplesWRS(freqCounter[i], generator);
		for (size_t t=0; t < sample.size(); ++t) {
			int j = sample[t];

			add_edge(i, j, subGraph);
			assigned += 1;
		}

		if (assigned < freqCounter[i]) {
			reassign += freqCounter[i] - int(sample.size());
			beta[i] = 0;
		}

		freqCounter[i] = 0;
		
		if (reassign > 0 && fail_count < N) {
			++fail_count;
			//  correction for first selection, probabilities for filled nodes will now be 0
			for (int ti=0; ti < N; ++ti) {
				n_degree[ti] = getNodeDegree(ti, subGraph);
			}
			setFirstFreqCounter(reassign, freqCounter, n_degree, beta, generator);
			// Must restart loop
			i = -1;
		}
	}

	std::vector<int> component(num_vertices(subGraph));
	int num_components = connected_components(subGraph, &component[0]);

	if (num_components > 1 && formConnected && m >= N - 1) {
		rewire(subGraph, theta, beta, component, num_components, generator, a);
	}
	
	for (int i=0; i < N; ++i) {
		n_degree[i] = getNodeDegree(i, subGraph);
	}
	
	return subGraph;
}


void assignCommunityNodes (uint K,
						  std::vector<int> &nodes,
						  std::vector< std::vector<int> > &communities,
						  std::vector< std::vector<int> > &nodeCommunities,
						  double overlapProb,
						  std::vector<double> &beta_k,
						  std::vector<double> &theta_k,
						  double a,
						  boost::mt19937 &generator,
						  double sigmaF,
						  double sigmaR,
						  int minCommunitySize)
{

	// Allows community positions to be fixed
	if (theta_k.size() < K) {
		theta_k = std::vector<double> (K);
		setTheta(theta_k, sigmaF, generator);
	}

	setBeta(theta_k, beta_k, sigmaR);

	boost::random::discrete_distribution<> dist(beta_k);
	boost::uniform_real<> rand_trial(0, 1) ;
	int k;
	
	for (size_t i=0; i < nodes.size(); ++i) {
		
		k = dist(generator);
		
		if (minCommunitySize) {
			for (size_t k_l=0; k_l < communities.size(); ++k_l) {
				if (int(communities[k_l].size()) < minCommunitySize) {
					k = int(k_l);
					break;
				}
			}
		}
		
		communities[k].push_back(nodes[i]);
		nodeCommunities[nodes[i]].push_back(k);
		// Now assign overlap
		if (overlapProb > 0.0) {

			double sumProbs = 0.0;
			for (uint j=0; j < K; ++j) {
				sumProbs += beta_k[j] * (-a * fabs(fabs(theta_k[k]) - fabs(theta_k[j])));
			}

			double prob;
			for (uint j=0; j < K; ++j) {
				if (int(j) == k)
					continue;

				prob = (beta_k[j] * (-a * fabs(fabs(theta_k[k]) - fabs(theta_k[j]))))/sumProbs;
				// chance of overlap is community selection probability, given first community multiplied by p_o parameter
				if (overlapProb * prob >= rand_trial(generator) ){
					communities[j].push_back(nodes[i]);
					nodeCommunities[nodes[i]].push_back(j);
				}
			}
		}
	}
}

int assignCommunityEdges (uint K,
						  int m,
						  std::vector< std::vector<int> > &communities,
						  std::vector<int> &communityEdgeCount,
						  std::vector<double> &beta_k,
						  boost::mt19937 &generator)
{

	std::set<int> filled;
	boost::random::discrete_distribution<> dist(beta_k);
	int interCommunityEdges = 0;
	uint k;
	double sum_of_elems=0;
	for(std::vector<double>::iterator j=beta_k.begin();j!=beta_k.end();++j)
		sum_of_elems += *j;
		
	while (interCommunityEdges < m && filled.size() < K && sum_of_elems > 0.0001) {
		k = dist(generator);
		// check maximal density is not reached, if it is, remove from selectable
		if ( double(communityEdgeCount[k]) < ((communities[k].size() * ( communities[k].size() - 1))/ 2.0)) {
			++communityEdgeCount[k];
			++interCommunityEdges;
		} else {
			if (filled.find(k) == filled.end()) {
				filled.insert(k);
				beta_k[k] = 0;
				boost::random::discrete_distribution<> dist(beta_k);
			}
		}
		
		sum_of_elems = 0;
		for(std::vector<double>::iterator j=beta_k.begin();j!=beta_k.end();++j)
			sum_of_elems += *j;
	}
	
	return interCommunityEdges;
}

int overlapReassignement (Graph &G,
						   int K,
						   int reassignEdges,
						   std::vector<double> &beta,
						   std::vector<double> &theta,
						   double a,
						   std::vector<double> &beta_communities,
						   std::vector< std::vector<int> > &communities,
						   std::vector<int> &communityEdgeCount,
						   std::vector< std::vector<double> > &beta_k,
						   std::vector< std::vector<int> > &n_degree_k,
						   std::vector<int> &n_degree,
						   boost::mt19937 &generator
  						)
{

	std::vector<int> communityRefill (K);
	std::vector<int> availableEdges (K);
	int k, n_i, N, EC;
	int assignable = 0;
	double betaNorm = 0;


	for (size_t k=0; k < communityEdgeCount.size(); ++k) {

		EC = 0;
		for (size_t i=0; i < communities[k].size(); ++i) {
			EC += n_degree_k[k][i];
		}

		N = communities[k].size();
		availableEdges[k] = ((N * (N -1) )/ 2) - (EC/2);

		if(availableEdges[k] == 0) {
			beta_communities[k] = 0;
		}

		assignable += availableEdges[k];
		betaNorm += beta_communities[k];

	}

	for (int k=0; k < K; ++k) {
		beta_communities[k] /= betaNorm;
	}

	// ensure we aren't assigning more edges than is possible
	if (assignable < reassignEdges)
		reassignEdges = assignable;

	boost::random::discrete_distribution<> dist(beta_communities);
	// When overlap occurs we need to reassign edges
	for (int i=0; i < reassignEdges; ++i) {

		do {
			k = dist(generator);
		} while(availableEdges[k] <= 0);

		++communityRefill[k];
		++communityEdgeCount[k];
		--availableEdges[k];
		// check maximal density is not reached, if it is, remove from selectable
		if ( availableEdges[k] == 0) {
			beta_communities[k] = 0;
			betaNorm = 0;

			for (int k=0; k < K; ++k)
				betaNorm += beta_communities[k];

			for (int k=0; k < K; ++k)
				beta_communities[k] /= betaNorm;

			boost::random::discrete_distribution<> dist(beta_communities);
		}

	}
	EC = 0;
	// Reassignment of edges inside communities
	for (size_t k=0; k < communityRefill.size(); ++k) {

		if (communityRefill[k] == 0)
			continue;

		std::vector<int> nodes (communities[k].size());
		std::vector<int> firstCounter (communities[k].size());

		for(size_t i=0; i < firstCounter.size(); i++){
			firstCounter[i] = 0;
			nodes[i] = i;
		}


		int count = 0;
		int N = nodes.size();
		int exit_early = 0;
		
		boost::random::discrete_distribution<> dist(beta_k[k]);
		while (count < communityRefill[k] && exit_early < N * 10) {
			int i = dist(generator);
			if (firstCounter[i] + n_degree_k[k][i] < N - 1) {
				++firstCounter[i];
				++count;
			}
			++exit_early;
		}


		for(size_t i=0; i < firstCounter.size(); i++){

			if (firstCounter[i] == 0)
				continue;
			int reassign = 0;

			n_i = communities[k][i];

			std::vector<int> available_population;
			// Add nodoes that aren't adjacent to the avialble population
			for (size_t t=0; t < communities[k].size(); t++){
				int j = communities[k][t];
				if (n_i != j && !edge(n_i, j, G).second )
					available_population.push_back(j);
			}

			std::vector<double> selectionScoresj = selectionScores(double(theta[n_i]), available_population, beta, theta, a) ;

			SampleARes sampler(available_population, selectionScoresj);
			std::vector<int> sample = sampler.getSamplesWRS(firstCounter[i], generator);

			for (size_t t=0; t < sample.size(); ++t) {
				int n_j = sample[t];
				add_edge(n_i, n_j, G);
				++EC;
				++n_degree[n_i];
				++n_degree[n_j];

			}

			if (int(sample.size()) < firstCounter[i]) {
				reassign = firstCounter[i] - int(sample.size());
				beta_k[k][i] = 0;
			}

			firstCounter[i] = 0;

			if (reassign > 0) {
				//  correction for first selection, probabilities for filled nodes will now be 0
				setFirstFreqCounter(reassign, firstCounter, n_degree_k[k], beta_k[k], generator);
				// Must restart loop
				i = -1;
			}

		}
	}
	return EC;
}

// Generate community graph
Graph generateGraph(uint N,
					uint K,
					double Density,
					double sigmaF,
					double sigmaR,
					double a,
					double communitySigmaF,
					double communitySigmaR,
					double ek_per,
					double overlapProb,
					std::vector<double> &theta,
					std::vector< std::vector<int> > &communities,
					std::vector<double> &theta_communities,
					int formConnected,
					int minDegree,
					int minCommunitySize,
					unsigned long int seed
				   )
{

	Graph G(N);

	boost::mt19937 generator;
	generator.seed(seed);

	std::vector<double> beta (N);

	std::vector<int> nodes (N);
	std::vector<int> n_degree (N);

	for (size_t i=0; i < nodes.size(); i++) {
		nodes[i] = i;
		n_degree[i] = 0;
	}

	int E = int (round((Density * N * (N-1))/2));

	// Allows theta to be a fixed
	if (theta.size() < N) {
		// Assign node positions
		theta = std::vector<double> (N);
		setTheta(theta, sigmaF, generator);
	}

	// normalised first selection probs
	setBeta(theta, beta, sigmaR);

	// No need to do all the other junk if we only generate a single community
	if (K == 1) {
		G = fillGraph(N, E, theta, beta, n_degree, sigmaR, a, generator, minDegree, formConnected);
		return G;
	}

	// vectors of inderemined size
	communities = std::vector< std::vector<int> >  (K);
	for (uint i=0; i < K; ++i) {
		communities[i] = std::vector<int> ();
	}

	int betweenEdges = int(round(E * ek_per));
	std::vector<double> beta_communities (K);

	std::vector<int> communityEdgeCount (K);
	std::vector< std::vector<int> > nodeCommunities (N);
	
	if (minCommunitySize > int(N/K)){
		minCommunitySize = int(N/K);
	}
	
	assignCommunityNodes(K, nodes,  communities, nodeCommunities, overlapProb, beta_communities, theta_communities, a, generator, communitySigmaF, communitySigmaR, minCommunitySize);

	std::vector<double> max_community (N);
	for (size_t n=0; n < nodes.size(); ++n) {
		double maxCom = 0;
		for (size_t t=0; t < nodeCommunities[n].size(); ++t) {
			int k = nodeCommunities[n][t];
			if (theta_communities[k] > maxCom) {
				maxCom = theta_communities[k];
			}
		}
		max_community[n] = maxCom;
	}

	// calculate maximal possible overlap  and number of available edges between communities
	int maxBE = 0;
	std::vector<int>::iterator it;

	// For all pairs of nodes
	for (size_t i=0; i < nodes.size(); ++i) {

		std::sort(nodeCommunities[i].begin(), nodeCommunities[i].end());
		for (size_t j = i + 1; j < nodes.size(); ++j) {
			if (i == 0)
				std::sort(nodeCommunities[j].begin(), nodeCommunities[j].end());

			std::vector<int> intersection;
			std::set_intersection (nodeCommunities[i].begin(), nodeCommunities[i].end(),
										nodeCommunities[j].begin(), nodeCommunities[j].end(),
										std::back_inserter(intersection));

			// If communities do not intersect, maxBE + 1
			if (intersection.size() == 0) {
				maxBE++;
			}
		}
	}

	// We cannot place more edges between communities than maxBE
	if (betweenEdges > maxBE) {
		betweenEdges = maxBE;
	}

	// Communities can end up being fully dense, we store the actual number of edges assigned
	betweenEdges = assignCommunityEdges(K, E - betweenEdges, communities, communityEdgeCount, beta_communities, generator);

	graph_traits<Graph>::edge_iterator ei, ei_end;
	int reassignEdges = 0;

	// selection for each community
	std::vector< std::vector<double> > theta_k (K);
	std::vector< std::vector<double> > beta_k (K);
	std::vector< std::vector<int> > n_degree_k (K);


	for (uint k=0; k < K; ++k) {

		theta_k[k] = std::vector<double> (communities[k].size());
		beta_k[k] = std::vector<double> (communities[k].size());
		n_degree_k[k] = std::vector<int> (communities[k].size());

		for (size_t i = 0; i < communities[k].size(); ++i) {
			theta_k[k][i] = theta[communities[k][i]];
			n_degree_k[k][i] = 0;
		}

		setBeta(theta_k[k], beta_k[k], sigmaR);

		Graph subGraph = fillGraph(communities[k].size(), communityEdgeCount[k], theta_k[k], beta_k[k], n_degree_k[k], sigmaR, a, generator, minDegree, formConnected);

		graph_traits<Graph>::edge_iterator ei, ei_end;
		// Merge the two graphs
		for (tie(ei, ei_end) = edges(subGraph); ei != ei_end; ++ei) {

			std::vector<int> epair (2);
			//track edges that occur multiple times
			epair[0] = communities[k][source(*ei, subGraph)];
			epair[1] = communities[k][target(*ei, subGraph)];

			std::sort(epair.begin(), epair.end());

			if (edge(epair[0], epair[1], G).second) {
				++reassignEdges;
			} else {
				add_edge(epair[0], epair[1], G);
				++n_degree[epair[0]];
				++n_degree[epair[1]];
			}
		}
	}

	if (reassignEdges){
		// additional edges as a result of overlap
		overlapReassignement(G, K, reassignEdges, beta, theta, a, beta_communities, communities, communityEdgeCount, beta_k, n_degree_k, n_degree, generator);
	}
	
	// selection for edges between communities
	setBeta(theta_communities, beta_communities, communitySigmaR);


	std::vector<int> freqCounter (K);
	for (size_t k=0; k < freqCounter.size(); ++k) {
		setBeta(theta_k[k], beta_k[k], sigmaR);
		freqCounter[k] = 0;
	}

	std::vector<int> nodeFreqCounter (N);

	int exit_early = 0;
	int interCount = 0;

	// until target edges achieved or possible edges are fill
	while (int(num_edges(G)) < int(E) and interCount < maxBE and exit_early < 300) {
		// select first community
		boost::random::discrete_distribution<> first_k_dist(beta_communities);
		int num = 0;

		for (size_t k=0; k < freqCounter.size(); ++k) {
			freqCounter[k] = 0;
		}

		for (size_t i=0; i < nodeFreqCounter.size(); ++i){
			nodeFreqCounter[i] = 0;
			
			if(n_degree[i] < minDegree){
				nodeFreqCounter[i] = minDegree - n_degree[i];
				num += minDegree - n_degree[i];
			}
		}
		
		while (num < int(E - num_edges(G))) {
			int k = first_k_dist(generator);
			if (communities[k].size() > 0) {
				++freqCounter[k];
				++num;
			}
		}

		for (size_t k=0; k < freqCounter.size(); ++k) {
			if (freqCounter[k] == 0)
				continue;

			boost::random::discrete_distribution<> dist(beta_k[k]);
			for (int l=0; l < freqCounter[k]; ++l) {
				// select node from community
				int n_i = dist(generator);
				int i = communities[k][n_i];
				++nodeFreqCounter[i];
			}
		}

		// select second nodes
		for (uint i=0; i < nodeFreqCounter.size(); ++i) {

			if (nodeFreqCounter[i] == 0)
				continue;

			std::vector<int> available_population;
			for (size_t t=0; t < nodes.size(); ++t) {

				uint j = nodes[t];

				if (edge(i, j, G).second || i == j)
					continue;

				std::vector<int> intersection;
				std::set_intersection (nodeCommunities[i].begin(), nodeCommunities[i].end(),
											nodeCommunities[j].begin(), nodeCommunities[j].end(),
											std::back_inserter(intersection));


				if (intersection.size() == 0) {
					available_population.push_back(t);
				}
			}

			// get correct weights
			std::vector<double>  selectionScoresj = interCommunitySelectionScores (i, theta[i], available_population, beta, theta, nodeCommunities, beta_communities, max_community, a);

			SampleARes sampler(available_population, selectionScoresj);
			std::vector<int> sample = sampler.getSamplesWRS(nodeFreqCounter[i], generator);

			for (size_t t=0; t < sample.size(); ++t) {
				int j = sample[t];

				add_edge(i, j, G);

				++n_degree[i];
				++n_degree[j];

				interCount += 1;
			}

		}

		// early termination in case of bugs
		++exit_early;
	}

	std::vector<int> component(num_vertices(G));
	int num_components = connected_components(G, &component[0]);


	if (num_components > 1 && formConnected && uint(ek_per * E)  >= K - 1) {
		rewire(G, theta, beta, component, num_components, generator, a, nodeCommunities);
	}
	return G;
}


/*
 * Converts graph to an adjacency list, not used but useful for testing
 */
std::vector< std::pair <int, int> > createEdgeList(Graph &G)
{
	std::vector< std::pair <int, int> > edgeList;
	graph_traits<Graph>::edge_iterator ei, ei_end;

	for (tie(ei, ei_end) = edges(G); ei != ei_end; ++ei) {
		// we want a sorted tuple to ensure the edge isn't already in our set
		std::vector<int> epair (2);
		epair[0] = source(*ei, G);
		epair[1] = target(*ei, G);
		std::sort(epair.begin(), epair.end());

		std::pair<int, int> p = std::make_pair(epair[0], epair[1]);

		bool found = false;

		for (size_t i=0; i < edgeList.size(); ++i) {
			if (edgeList[i].first == p.first && edgeList[i].second == p.second) {
				found = true;
				break;
			}
		}

		if (!found) {
			edgeList.push_back(p);
		}
	}

	return edgeList;
}
