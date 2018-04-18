#include "landmark.h"
#include "Eigen/Dense"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

MatrixXd LandmarkObsCollection::GetPositionMatrix()
{
	MatrixXd positionMatrix = MatrixXd(landmarkObsVect_.size(), 2);
	for(int i=0; i<landmarkObsVect_.size(); i++)
	{
		positionMatrix.row(i) << landmarkObsVect_[i].x, landmarkObsVect_[i].y;
	}

	return positionMatrix;
}

MatrixXd LandmarkObsCollection::GetPositionMatrix_1col()
{
	MatrixXd positionMatrix = MatrixXd(landmarkObsVect_.size(), 3);
	for(int i=0; i<landmarkObsVect_.size(); i++)
	{
		positionMatrix.row(i) << landmarkObsVect_[i].x, landmarkObsVect_[i].y, 1;
	}

	return positionMatrix;
}

 LandmarkObsCollection LandmarkObsCollection::ConvertToMapCoords(Particle &particle)
{
	MatrixXd positionMatrix_1col = GetPositionMatrix_1col();

	MatrixXd homogeneousTransform = MatrixXd(3, 3);
  homogeneousTransform << cos(particle.theta), -sin(particle.theta), particle.x,
                            sin(particle.theta), cos(particle.theta), particle.y,
                            0, 0, 1;

	MatrixXd mapCoordsMaxtrix = positionMatrix_1col * homogeneousTransform.transpose();

	LandmarkObsCollection mapCoordsLandmarkObsCollection;
	for(int i=0; i<mapCoordsMaxtrix.rows(); i++)
	{
		LandmarkObs currLandmarkObs;
		currLandmarkObs.id = landmarkObsVect_[i].id;
		currLandmarkObs.x = mapCoordsMaxtrix(i, 0);
		currLandmarkObs.y = mapCoordsMaxtrix(i, 1);

		mapCoordsLandmarkObsCollection.landmarkObsVect_.push_back(currLandmarkObs);
	}

	return mapCoordsLandmarkObsCollection;
}