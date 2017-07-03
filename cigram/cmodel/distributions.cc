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

#include "distributions.h"

/*
 * Wrapped gaussian pdf
 */
double alphaValue(double theta, double sigma, double mu)
{
	double alpha = 0;

	for (int k=-1; k <= 1; ++k)
		alpha += exp (-pow( theta - mu + (2 * pi * k), 2) / ( 2 * pow(sigma, 2) ));

	alpha = (1/(sigma * sqrt(2 * pi) )) * alpha;

	return alpha;
}

/*
 * Wrapped gaussian positions
 */
void setTheta(std::vector<double> &theta, double sigmaF, boost::mt19937 &generator, double mu) {
	// generator object, Is it possible to use a single generator for normal distribution and weighted stuff??
	boost::normal_distribution<> normal (mu, sigmaF);

	// Assign node positions
	for (size_t i=0; i < theta.size(); ++i) {
		theta[i] = normal(generator);

		if (theta[i] >= 0)
			theta[i] = fmod(theta[i], pi);
		else
			theta[i] = -fmod(fabs(theta[i]), pi);

		if (theta[i] > pi)
			theta[i] = -pi + (theta[i] - pi);
		if (theta[i] < -pi)
			theta[i] = pi - fabs(theta[i] - -pi);

	}

}

void setBeta(std::vector<double> &theta, std::vector<double> &beta, double sigmaR) {

	double betaSum = 0;
	// Generate our alpha coordiantes, wrapped gaussians
	for (size_t i=0; i < theta.size(); ++i) {
		beta[i] = alphaValue(theta[i], sigmaR);
		betaSum += beta[i];
	}

	// Normalise betas, fill selection matrix
	for (size_t i=0; i < beta.size(); ++i)
		beta[i] /= betaSum;
}


std::vector<double> selectionScores(double position, std::vector<int> &population, std::vector<double> &beta, std::vector<double> &theta, double a) {
	std::vector<double> scores (population.size());
	double sum = 0;
	for (size_t i=0; i < population.size(); ++i) {
		int n_j = population[i];
		scores[i] = beta[n_j] * exp(-a * fabs(fabs(position) - fabs(theta[n_j])));
		sum += scores[i];
	}

	if (sum == 0){
		/* Where values extremely close to 0, set to uniform probability */
		for (size_t i=0; i < population.size(); ++i) {
			scores[i] = 1/population.size();
		}	
		
		return scores;
	}
		

	for (size_t i=0; i < scores.size(); ++i)
		scores[i] /= sum;

	return scores;
}


std::vector<double> interCommunitySelectionScores (int n_i,
												   double position,
												   std::vector<int> &population,
												   std::vector<double> &beta,
												   std::vector<double> &theta,
												    std::vector< std::vector<int> > &nodeCommunities,
												   std::vector<double> &beta_community,
												   std::vector<double> &max_community,
												   double a)
{
	std::vector<double> scores (population.size());
	double sum = 0;


	for (size_t i=0; i < population.size(); ++i) {
		int n_j = population[i];

		double commProbs = 0;
		double dist;
		for (uint k=0; k < nodeCommunities[n_i].size(); ++k) {
			int k1 = nodeCommunities[n_i][k];
			for (uint l=0; l < nodeCommunities[n_j].size(); ++l) {
				int k2 = nodeCommunities[n_j][l];

				dist = fabs(fabs(beta_community[k1]) - fabs(beta_community[k2]));
				commProbs += beta_community[k1] * beta_community[k2] * exp(-a * dist);;

			}
		}
		scores[i] = commProbs * beta[n_j] * exp(-a * (fabs(fabs(position) - fabs(theta[n_j])) + fabs(fabs(max_community[n_i]) - fabs(max_community[n_j]))));
		sum += scores[i];
	}

	if (sum == 0)
		return scores;

	for (size_t i=0; i < scores.size(); ++i)
		scores[i] /= sum;

	return scores;
}
