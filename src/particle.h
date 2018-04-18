#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <random>
#include "Eigen/Dense"

typedef struct {
	int id;
	double x;
	double y;
	double theta;
	double weight;
	std::vector<int> associations;
	std::vector<double> sense_x;
	std::vector<double> sense_y;
} Particle;

class ParticleCollection 
{
	public:
		std::vector<Particle> particles_;

		Eigen::MatrixXd GetPositionMatrix();
		Eigen::MatrixXd GetPositionMatrix_1col();
		Eigen::MatrixXd ConvertToMapCoords(const Eigen::MatrixXd &positionMatrix_1col, const Particle &particle);
		Eigen::VectorXd GetDistances(const Eigen::MatrixXd &positionMatrix, const Eigen::VectorXd &landmarkMapPosVec);
};

#endif /* PARTICLE_H_ */