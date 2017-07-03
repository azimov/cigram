/*
 * CiGRAM - Circular Gaussian Random grAph Model
 * Copyright (C) 2014  James Gilbert
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
#include <boost/random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <algorithm>

/*
 * Weighted sampling without replacement using a resivoir method
 */
class SampleARes {

	std::vector<int> population;
	std::vector<double> weights;
	bool weightsSet;
	/*
	 * TODO : weighted sample with replacement
	 * TODO : unweighted sample with and without replacement (constructor creates unform weights)
	 * TODO : population as a generic template
	 */
	public:
		// weighted sample constructor
		SampleARes (std::vector<int> &pop, std::vector<double> &w);
		// sampling without replacement
		std::vector<int> getSamplesWRS(int sampleSize, boost::mt19937 &generator);

	private:
		// sort function for heap
		static bool resivoirComparitor (const std::pair<int, double>  &l, const std::pair<int, double>  &r) { return l.second > r.second; }
};
