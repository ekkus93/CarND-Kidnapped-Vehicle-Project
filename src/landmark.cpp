#include "landmark.h"

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
	MatrixXd positionMatrix = MatrixXd(landmarkObsVect_.size(), 2);
	for(int i=0; i<landmarkObsVect_.size(); i++)
	{
		positionMatrix.row(i) << landmarkObsVect_[i].x, landmarkObsVect_[i].y, 1;
	}

	return positionMatrix;
}

MatrixXd LandmarkObsCollection::ConvertToMapCoords(const MatrixXd &positionMatrix_1col, const Particle &particle)
{
	MatrixXd homogeneousTransform = MatrixXd(3, 3);
  homogeneousTransform << cos(particle.theta), -sin(particle.theta), particle.x,
                            sin(particle.theta), cos(particle.theta), particle.y,
                            0, 0, 1;

	return positionMatrix_1col * homogeneousTransform.transpose();
}