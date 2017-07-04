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
#include "sample.h"

SampleARes::SampleARes(std::vector<int> &pop, std::vector<double> &w) {

	if (pop.size() != w.size())
		throw std::invalid_argument( "weight and population must be the same size" );

	population = pop;
	weights = w;
	weightsSet = false;
	// Normalise the weights
	double wSum = 0;

	for (size_t i = 0; i < weights.size(); ++i)
		wSum += w[i];

	if (wSum > 0) {
		for (size_t i = 0; i < weights.size(); ++i){
			weights[i] = w[i]/wSum;
		}
		weightsSet = true;
	}
}


/*
 * Weighted sampling without replacement.
 * See "Weighted random sampling with a reservoir" by Efraimidis and Spirakis
 * doi: 10.1016/j.ipl.2005.11.003
 */
std::vector<int> SampleARes::getSamplesWRS(int sampleSize, boost::mt19937 &generator) {
	if (population.size() < uint(sampleSize))
		sampleSize = int(population.size());

	// If weights sum to zero, return vector of length 0
	if (!weightsSet) {
		std::vector<int> sample;
		return sample;
	}

	double key, randV, t_w;
	double weightsSum = 0;
	uint i = 0;
	boost::uniform_real<> rand_double(0, 1) ;
	std::vector< std::pair<int, double> > resivoir;
	// assign first sampleSize keys to resivoir
	while (resivoir.size() < uint(sampleSize) &&  i < weights.size()) {
		weightsSum += weights[i];
		if (weights[i] > 0 ) {
			randV = rand_double(generator);
			key = pow(randV, 1/weights[i]);
			resivoir.push_back(std::make_pair(int(i), key));
		}
		++i;
	}

	uint sampleEndPoint = i;
	// sort the resivoir by the key
	std::make_heap(resivoir.begin(), resivoir.end(), resivoirComparitor);
	std::pop_heap(resivoir.begin(), resivoir.end(), resivoirComparitor);
	std::pair<int, double> thresh = resivoir.back();

	// comparing key check allows exponential jumps reduing the number of ranom checks from O(n) to O(s log(n/s)
	double keyCheck = log(rand_double(generator))/log(thresh.second);

	//  Go through each population member and assign a key, if key is larger than lowest value in resivoir, insert it.
	for (uint i = sampleEndPoint; i < weights.size(); ++i) {

		weightsSum += weights[i];
		if (weights[i] > 0 && weightsSum > keyCheck) {

			resivoir.pop_back();
			// create new key must be greater t
			t_w = pow(thresh.second, weights[thresh.first]);
			boost::uniform_real<> rand_key(t_w, 1);
			randV = rand_key(generator);
			key = pow(randV, 1/weights[i]);
			// insert new, remove lowest
			resivoir.push_back(std::make_pair(int(i), key));
			std::make_heap(resivoir.begin(), resivoir.end(), resivoirComparitor);

			// threshold is the lowest key in the list
			std::pop_heap(resivoir.begin(), resivoir.end(), resivoirComparitor);
			thresh = resivoir.back();
			keyCheck = log(rand_double(generator))/log(thresh.second);
		}
	}

	// return the selected members of the population
	std::vector<int> sample (resivoir.size());
	for (size_t i=0; i < resivoir.size(); ++i)  {
		sample[i] = population[resivoir[i].first];
	}

	return sample;
}
