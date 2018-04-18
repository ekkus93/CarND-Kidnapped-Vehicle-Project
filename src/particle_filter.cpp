/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <sstream>
#include <string>
#include <iterator>
#include "Eigen/Dense"
#include "particle_filter.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// *TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	num_particles = 100;

	normal_distribution<double> N_x(x, std[0]);
	normal_distribution<double> N_y(y, std[1]);
	normal_distribution<double> N_theta(theta, std[2]);

	for (int i=0; i < num_particles; i++)
	{
		Particle particle;
		particle.id = i;
		particle.x = N_x(gen);
		particle.y = N_y(gen);
		particle.theta = N_theta(gen);
		particle.weight = 1;

		particleCollection_.particles_.push_back(particle);
		weights.push_back(1);
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// *TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	for (int i=0; i<num_particles; i++)
	{
		double new_x;
		double new_y;
		double new_theta;

		double curr_x = particleCollection_.particles_[i].x;
		double curr_y = particleCollection_.particles_[i].y;
		double curr_theta = particleCollection_.particles_[i].theta;

		if (fabs(yaw_rate) > 0.001)
		{
			new_x = curr_x + velocity/yaw_rate*(sin(curr_theta+yaw_rate*delta_t) - sin(curr_theta));
			new_y = curr_y + velocity/yaw_rate*(cos(curr_theta)-cos(curr_theta+yaw_rate*delta_t));
			new_theta = curr_theta + yaw_rate*delta_t;
		}
		else
		{
			new_x = curr_x + velocity*delta_t*cos(curr_theta);
			new_y = curr_y + velocity*delta_t*sin(curr_theta);
			new_theta = curr_theta;
		}

		normal_distribution<double> N_x(new_x, std_pos[0]);
		normal_distribution<double> N_y(new_y, std_pos[1]);
		normal_distribution<double> N_theta(new_theta, std_pos[2]);

		particleCollection_.particles_[i].x = N_x(gen);
    particleCollection_.particles_[i].y = N_y(gen);
    particleCollection_.particles_[i].theta = N_theta(gen);
	}
}

void ParticleFilter::dataAssociation(LandmarkObsCollection &predictedCollection, 
																			LandmarkObsCollection &observationCollection) {
	// *TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

  for (int i=0; i<observationCollection.landmarkObsVect_.size(); i++) {
    double min_dist = numeric_limits<double>::max();

    for (int j=0; j<predictedCollection.landmarkObsVect_.size(); j++) {
			LandmarkObs currPred = predictedCollection.landmarkObsVect_[j];
      double d = dist(observationCollection.landmarkObsVect_[i].x, 
											observationCollection.landmarkObsVect_[i].y, 
											currPred.x, currPred.y);
      if (d < min_dist) 
			{
        observationCollection.landmarkObsVect_[i].id	= currPred.id;
        min_dist = d;
      }
  	}
  }
}

LandmarkObs ParticleFilter::ConvertObsCoordsToMapCoords(const LandmarkObs &obs, const Particle &particle)
{
	LandmarkObs mapObs;
	mapObs.id = obs.id;
	mapObs.x = cos(particle.theta) * obs.x - sin(particle.theta) * obs.y + particle.x;
	mapObs.y = sin(particle.theta) * obs.x + cos(particle.theta) * obs.y + particle.y;

	return mapObs;
}

LandmarkObsCollection ParticleFilter::FilterMapLandmarks(const Map &map, double x, double y, double max_range)
{
	LandmarkObsCollection filteredMapLandmarkCollection;

	for(int i=0; i<map.landmark_list.size(); i++)
	{
		int curr_id = map.landmark_list[i].id_i;
		double curr_x = map.landmark_list[i].x_f;
		double curr_y = map.landmark_list[i].y_f;

		double curr_dist = dist(x, y, curr_x, curr_y);
		if (curr_dist < max_range)
		{
			LandmarkObs currMapLandmark;
			currMapLandmark.id = curr_id;
			currMapLandmark.x = curr_x;
			currMapLandmark.y = curr_y;

			filteredMapLandmarkCollection.landmarkObsVect_.push_back(currMapLandmark);
		}
	}

	return filteredMapLandmarkCollection;
}

double ParticleFilter::CalcWeight(const LandmarkObsCollection &transformedObsCollection, 
const LandmarkObsCollection &predictedLandmarkCollection, 
double std_x, double std_y)
{
  double particle_likelihood = 1.0;

  double mu_x, mu_y;
  for (int i=0; i<transformedObsCollection.landmarkObsVect_.size(); i++) {
		LandmarkObs tObs = transformedObsCollection.landmarkObsVect_[i];

    // Find corresponding landmark on map for centering gaussian distribution
    for (int j=0; j<predictedLandmarkCollection.landmarkObsVect_.size(); j++)
		{
			LandmarkObs pLandmark = predictedLandmarkCollection.landmarkObsVect_[j];
      if (tObs.id == pLandmark.id) {
        mu_x = pLandmark.x;
        mu_y = pLandmark.y;
        break;
      }
		}

		double norm_factor = 2 * M_PI * std_x * std_y;
    double prob = exp( -( pow(tObs.x - mu_x, 2) / (2 * std_x * std_x) + 
										pow(tObs.y - mu_y, 2) / (2 * std_y * std_y) ) );

    particle_likelihood *= prob / norm_factor;
	}	

	return particle_likelihood;
}

void ParticleFilter::NormalizeParticleWeights() 
{
	double denom = numeric_limits<double>::epsilon();
	for(int i=0; i<particleCollection_.particles_.size(); i++)
	{
		denom += particleCollection_.particles_[i].weight;
	}

	for(int i=0; i<particleCollection_.particles_.size(); i++)
	{
		particleCollection_.particles_[i].weight /= denom;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                      const std::vector<LandmarkObs> &observations, const Map &map_landmarks)																	
{
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	double std_x = std_landmark[0];
	double std_y = std_landmark[1];

	for(int i=0; i<particleCollection_.particles_.size(); i++)
	{
		double p_x = particleCollection_.particles_[i].x;
		double p_y = particleCollection_.particles_[i].y;
		double p_theta = particleCollection_.particles_[i].theta;

		LandmarkObsCollection predictedLandmarkCollection =
			FilterMapLandmarks(map_landmarks, p_x, p_y, sensor_range);

		LandmarkObsCollection transformedObsCollection;
		for(int j=0; j<observations.size(); j++)
		{
			LandmarkObs transformedObs = ConvertObsCoordsToMapCoords(observations[j], particleCollection_.particles_[i]);

			transformedObsCollection.landmarkObsVect_.push_back(transformedObs);
		}

		dataAssociation(predictedLandmarkCollection, transformedObsCollection);

		particleCollection_.particles_[i].weight = CalcWeight(transformedObsCollection, predictedLandmarkCollection, 
																			std_x, std_y);
	}

	NormalizeParticleWeights();
}

void ParticleFilter::resample() {
	// *TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	int N = particleCollection_.particles_.size();
	weights.clear();
	for(int i=0; i<particleCollection_.particles_.size(); i++)
	{
		weights.push_back(particleCollection_.particles_[i].weight);
	}

	discrete_distribution<int> distribution(weights.begin(),weights.end());

	vector<Particle> newParticles;

	int index = distribution(gen);

	double beta = 0.0;
	double weight_max = *(max_element(begin(weights), end(weights)));
	normal_distribution<double> N_w(0.0, 2.0*weight_max);

	for(int i=0; i<N; i++)
	{
		beta += N_w(gen);
		while (beta > weights[index])
		{
			beta -= weights[index];
			index = (index+1) % N;
		}
		newParticles.push_back(particleCollection_.particles_[index]);
	}

	particleCollection_.particles_ = newParticles;

	// reset the weights to 1.0
	for(int i=0; i<particleCollection_.particles_.size(); i++)
	{
		particleCollection_.particles_[i].weight = 1.0;
	}
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;

		return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
