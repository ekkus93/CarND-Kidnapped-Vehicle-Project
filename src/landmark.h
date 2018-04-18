#ifndef LANDMARK_H_
#define LANDMARK_H_

#include <vector>
#include "particle.h"
#include "Eigen/Dense"

/*
 * Struct representing one landmark observation measurement.
 */
typedef struct {
	int id;				// Id of matching landmark in the map.
	double x;			// Local (vehicle coordinates) x position of landmark observation [m]
	double y;			// Local (vehicle coordinates) y position of landmark observation [m]
} LandmarkObs;

class LandmarkObsCollection
{
	public:
		std::vector<LandmarkObs> landmarkObsVect_;
		
		Eigen::MatrixXd GetPositionMatrix();
		Eigen::MatrixXd GetPositionMatrix_1col();
		LandmarkObsCollection ConvertToMapCoords(Particle &particle);
};

#endif /* LANDMARK_H_*/