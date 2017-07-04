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
#include <math.h>
#include <iostream>
#include <fstream>
#include <stack>
#include <set>

#include <boost/config.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/subgraph.hpp>

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/math/constants/constants.hpp>

#include "distributions.h"
#include "sample.h"

using namespace boost;
typedef boost::adjacency_matrix<boost::undirectedS> Graph;
typedef boost::property_map<Graph, boost::vertex_index_t>::type IndexMap;
typedef boost::graph_traits<Graph>::vertex_descriptor VertexDesc;

/*
 * Written in horrible procedural manner
 * TODO: rewrite in OOP
 */

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
				int formConnected=1,
				int minDegree=1,
				int minCommunitySize=0,
				unsigned long int seed=time(0));
