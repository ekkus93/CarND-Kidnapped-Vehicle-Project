#include <numeric>
#include <math.h> 
#include <sstream>
#include <string>
#include <iterator>
#include "Eigen/Dense"
#include "particle.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

MatrixXd ParticleCollection::GetPositionMatrix()
{
	MatrixXd positionMatrix = MatrixXd(particles_.size(), 2);
	for(int i=0; i<particles_.size(); i++)
	{
		positionMatrix.row(i) << particles_[i].x, particles_[i].y;
	}

	return positionMatrix;
}

MatrixXd ParticleCollection::GetPositionMatrix_1col()
{
	MatrixXd positionMatrix = MatrixXd(particles_.size(), 3);
	for(int i=0; i<particles_.size(); i++)
	{
		positionMatrix.row(i) << particles_[i].x, particles_[i].y, 1;
	}

	return positionMatrix;
}

MatrixXd ParticleCollection::ConvertToMapCoords(const MatrixXd &positionMatrix_1col, const Particle &particle)
{
	MatrixXd homogeneousTransform = MatrixXd(3, 3);
  homogeneousTransform << cos(particle.theta), -sin(particle.theta), particle.x,
                            sin(particle.theta), cos(particle.theta), particle.y,
                            0, 0, 1;

	return positionMatrix_1col * homogeneousTransform.transpose();
}

VectorXd ParticleCollection::GetDistances(const MatrixXd &positionMatrix, const VectorXd &landmarkMapPosVec)
{
	MatrixXd diffMatrix = positionMatrix.array().rowwise() - landmarkMapPosVec.transpose().array();
	
	return diffMatrix.rowwise().norm();
}