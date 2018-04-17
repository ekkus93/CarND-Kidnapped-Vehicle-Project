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
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include "Eigen/Dense"

#include "particle_filter.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

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

		particles.push_back(particle);
		weights.push_back(1);

		assert(particles.size()>0);
		assert(weights.size()>0);
	}

	cout << "particles.size(): " << particles.size() << "\n";
	cout << "weights.size(): " << weights.size() << "\n";
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

		double curr_x = particles[i].x;
		double curr_y = particles[i].y;
		double curr_theta = particles[i].theta;

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

		particles[i].x = N_x(gen);
    particles[i].y = N_y(gen);
    particles[i].theta = N_theta(gen);
	}
}

LandmarkObs ParticleFilter::GetBestMapLandmarkForObservationLandmark(const vector<LandmarkObs> &mapLandmarks, const LandmarkObs &observation)
{
	LandmarkObs bestMapLandmark;
	bestMapLandmark.id = 0;

	double bestDistance = std::numeric_limits<double>::infinity();

	for(int i=0; i<mapLandmarks.size(); i++)
	{
		double currDistance = dist(mapLandmarks[i].x, 
																mapLandmarks[i].y,
																observation.x, observation.y);
		if (currDistance < bestDistance)
		{
			bestMapLandmark = mapLandmarks[i];
			bestDistance = currDistance;
		}
	}
			
	return bestMapLandmark;
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// *TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

  for (int i=0; i<observations.size(); i++) {
		LandmarkObs currObs = observations[i];
    double min_dist = numeric_limits<double>::max();

    for (int j=0; j<predicted.size(); j++) {
			LandmarkObs currPred = predicted[j];
      double d = dist(currObs.x, currObs.y, currPred.x, currPred.y);
      if (d < min_dist) 
			{
        currObs.id	= currPred.id;
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

double ParticleFilter::CalcGaussianNorm(double sig_x, double sig_y)
{
	return 1.0/(2.0 * M_PI * sig_x * sig_y);
}

vector<LandmarkObs> ParticleFilter::FilterMapLandmarks(const Map &map, double x, double y, double max_range)
{
	vector<LandmarkObs> filteredMapLandmarks;

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

			filteredMapLandmarks.push_back(currMapLandmark);
		}
	}

	return filteredMapLandmarks;
}

double ParticleFilter::CalcWeight(const vector<LandmarkObs> &transformedObservations, 
vector<LandmarkObs> &predictedLandmarks, 
double norm_factor, double std_x, double std_y)
{
  double particle_likelihood = 1.0;

  double mu_x, mu_y;
  for (int i=0; i<transformedObservations.size(); i++) {
		LandmarkObs tObs = transformedObservations[i];

    // Find corresponding landmark on map for centering gaussian distribution
    for (int j=0; j<predictedLandmarks.size(); j++)
		{
      if (tObs.id == predictedLandmarks[i].id) {
        mu_x = predictedLandmarks[i].x;
        mu_y = predictedLandmarks[i].y;
        break;
      }
		}

    double prob = exp( -( pow(tObs.x - mu_x, 2) / (2 * std_x * std_x) + pow(tObs.y - mu_y, 2) / (2 * std_y * std_y) ) );

    particle_likelihood *= prob / norm_factor;
	}	

	return particle_likelihood;
}

void ParticleFilter::NormalizeParticles(vector<Particle> &particles) 
{
	double denom = numeric_limits<double>::epsilon();
	for(int i=0; i<particles.size(); i++)
	{
		denom += particles[i].weight;
	}

	for(int i=0; i<particles.size(); i++)
	{
		particles[i].weight /= denom;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   std::vector<LandmarkObs> observations, Map map_landmarks) 
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
	double gausNorm = CalcGaussianNorm(std_x, std_y);

	for(int i=0; i<particles.size(); i++)
	{
		double p_x = particles[i].x;
		double p_y = particles[i].y;
		double p_theta = particles[i].theta;

		vector<LandmarkObs> predictedLandmarks =
			FilterMapLandmarks(map_landmarks, p_x, p_y, sensor_range);		

		vector<LandmarkObs> transformedObservations;

		for(int j=0; j<observations.size(); j++)
		{
			LandmarkObs transformedObs = ConvertObsCoordsToMapCoords(observations[j], particles[i]);

			transformedObservations.push_back(transformedObs);
		}

		dataAssociation(predictedLandmarks, transformedObservations);

		particles[i].weight = CalcWeight(transformedObservations,
																			predictedLandmarks, gausNorm, 
																			std_x, std_y);
	}

	NormalizeParticles(particles);
}



void ParticleFilter::resample() {
	// *TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	int N = particles.size();
	weights.clear();
	for(int i=0; i<particles.size(); i++)
	{
		weights.push_back(particles[i].weight);
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
		newParticles.push_back(particles[index]);
	}

	particles = newParticles;

	// reset the weights to 1.0
	for(int i=0; i<particles.size(); i++)
	{
		particles[i].weight = 1.0;
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
