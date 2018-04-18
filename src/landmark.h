#ifndef LANDMARK_H_
#define LANDMARK_H_

#include "helper_functions.h"
#include "Eigen/Dense"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Struct representing one landmark observation measurement.
 */
struct LandmarkObs {
	
	int id;				// Id of matching landmark in the map.
	double x;			// Local (vehicle coordinates) x position of landmark observation [m]
	double y;			// Local (vehicle coordinates) y position of landmark observation [m]
};

class LandmarkObsCollection
{
	public:
		vector<LandmarkObs> landmarkObsVect_;
		
		MatrixXd GetPositionMatrix();
		MatrixXd GetPositionMatrix_1col();
		MatrixXd ConvertToMapCoords(const MatrixXd &positionMatrix_1col, const LandmarkObs &landmarkObj);
};

#endif /* LANDMARK_H_*/