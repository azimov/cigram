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
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/math/constants/constants.hpp>

const double pi = boost::math::constants::pi<double>();

double alphaValue(double theta, double sigma, double mu=pi);

void setTheta(std::vector<double> &theta, double sigmaF, boost::mt19937 &generator, double mu=0);

void setBeta(std::vector<double> &theta, std::vector<double> &beta, double sigmaR);

std::vector<double> selectionScores(double position, std::vector<int> &population, std::vector<double> &beta, std::vector<double> &theta, double a);

std::vector<double> interCommunitySelectionScores (int n_i,
												   double position,
												   std::vector<int> &population,
												   std::vector<double> &beta,
												   std::vector<double> &theta,
												    std::vector< std::vector<int> > &nodeCommunities,
												   std::vector<double> &beta_community,
												   std::vector<double> &max_community,
												   double a);
