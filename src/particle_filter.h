/*
 * particle_filter.h
 *
 * 2D particle filter class.
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#ifndef PARTICLE_FILTER_H_
#define PARTICLE_FILTER_H_

#include <random>
#include "helper_functions.h"
#include "Eigen/Dense"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

struct Particle {
	int id;
	double x;
	double y;
	double theta;
	double weight;
	std::vector<int> associations;
	std::vector<double> sense_x;
	std::vector<double> sense_y;
};

class ParticleCollection 
{
	public:
		vector<Particle> particles_;

		MatrixXd GetPositionMatrix();
		MatrixXd GetPositionMatrix_1col();
		MatrixXd ConvertToMapCoords(const MatrixXd &positionMatrix_1col, const Particle &particle);
		VectorXd GetDistances(const MatrixXd &positionMatrix, const VectorXd &landmarkMapPosVec);
};

class ParticleFilter {
	// Number of particles to draw
	int num_particles; 
	
	std::default_random_engine gen;
	
	
	// Flag, if filter is initialized
	bool is_initialized;
	
	// Vector of weights of all particles
	std::vector<double> weights;
	
public:
	
	// Set of current particles
	//std::vector<Particle> particles;
	ParticleCollection particleCollection_;

	// Constructor
	// @param num_particles Number of particles
	ParticleFilter() : num_particles(0), is_initialized(false) {}

	// Destructor
	~ParticleFilter() {}

	/**
	 * init Initializes particle filter by initializing particles to Gaussian
	 *   distribution around first position and all the weights to 1.
	 * @param x Initial x position [m] (simulated estimate from GPS)
	 * @param y Initial y position [m]
	 * @param theta Initial orientation [rad]
	 * @param std[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
	 *   standard deviation of yaw [rad]]
	 */
	void init(double x, double y, double theta, double std[]);

	/**
	 * prediction Predicts the state for the next time step
	 *   using the process model.
	 * @param delta_t Time between time step t and t+1 in measurements [s]
	 * @param std_pos[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m]
	 *   standard deviation of yaw [rad]]
	 * @param velocity Velocity of car from t to t+1 [m/s]
	 * @param yaw_rate Yaw rate of car from t to t+1 [rad/s]
	 */
	void prediction(double delta_t, double std_pos[], double velocity, double yaw_rate);
	
	/**
	 * dataAssociation Finds which observations correspond to which landmarks (likely by using
	 *   a nearest-neighbors data association).
	 * @param predicted Vector of predicted landmark observations
	 * @param observations Vector of landmark observations
	 */
	void dataAssociation(LandmarkObsCollection &predictedCollection, 
																			LandmarkObsCollection &observationCollection);
	
	/**
	 * ConvertObsCoordsToMapCoords converts the local coordinates of
	 * obs to map coords in reference to particle.
	 * @param obs - landmark observation with local coordinates 
	 * @param particle - reference particle
	 * @return obs with map coordinates
	 */
	LandmarkObs ConvertObsCoordsToMapCoords(const LandmarkObs &obs, const Particle &particle);

	/**
	 * FilterMapLandmarks filters out map landmark observations which are
	 * outside of max_range from x,y
	 * @param map - Map with landmark observations
	 * @param x - x reference
	 * @param y - y reference
	 * @param max_range - max range from x,y
	 * @return subset of map landmark observations
	 */
	LandmarkObsCollection FilterMapLandmarks(const Map &map, double x, double y, double max_range);

	/**
	 * CalcWeight calculates the weight of a particle
	 * @param transformedObservations - observations based on map coordinates
	 * @param predictedLandmarks - filtered map landmarks
	 * @param std_x - standard dev of x
	 * @param std_y - standard dev of y
	 * @return the weight for a particle 
	 */
	double CalcWeight(const LandmarkObsCollection &transformedObsCollection, 
										LandmarkObsCollection &predictedLandmarkCollection, 
										double std_x, double std_y);

	/**
	 * NormalizeParticleWeights normalizes all of the 
	 * weights for each of the particles.
	 */
	void NormalizeParticleWeights();

	/**
	 * updateWeights Updates the weights for each particle based on the likelihood of the 
	 *   observed measurements. 
	 * @param sensor_range Range [m] of sensor
	 * @param std_landmark[] Array of dimension 2 [Landmark measurement uncertainty [x [m], y [m]]]
	 * @param observations Vector of landmark observations
	 * @param map Map class containing map landmarks
	 */
	void updateWeights(double sensor_range, double std_landmark[],
                      std::vector<LandmarkObs> observations, Map map_landmarks);
	
	/**
	 * resample Resamples from the updated set of particles to form
	 *   the new set of particles.
	 */
	void resample();

	/*
	 * Set a particles list of associations, along with the associations calculated world x,y coordinates
	 * This can be a very useful debugging tool to make sure transformations are correct and assocations correctly connected
	 */
	Particle SetAssociations(Particle& particle, const std::vector<int>& associations,
		                     const std::vector<double>& sense_x, const std::vector<double>& sense_y);

	
	std::string getAssociations(Particle best);
	std::string getSenseX(Particle best);
	std::string getSenseY(Particle best);

	/**
	* initialized Returns whether particle filter is initialized yet or not.
	*/
	const bool initialized() const {
		return is_initialized;
	}
};



#endif /* PARTICLE_FILTER_H_ */
