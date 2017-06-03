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

#include "particle_filter.h"
using namespace std;
static default_random_engine gen;
void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	//std is the sigma_pos from main.cpp, these are the uncertainties

	num_particles = 200;
	//Create normal distributions
	normal_distribution<double> distX(0, std[0]);	
	normal_distribution<double> distY(0, std[1]);
	normal_distribution<double> distTheta(0, std[2]);

	//Initialize the particles
	for (int i = 0; i < num_particles; i++)
	{
		Particle p;
		p.id = i;
		p.x = x;
		p.y = y;
		p.theta = theta;
		p.weight = 1.0;

		//add the noise in.
		p.x += distX(gen);
		p.y += distY(gen);
		p.theta += distTheta(gen);

		particles.push_back(p);
	}
	is_initialized = true;
	

	

}
void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	//Create normal distributions
	normal_distribution<double> distX(0, std_pos[0]);
	normal_distribution<double> distY(0, std_pos[1]);
	normal_distribution<double> distTheta(0, std_pos[2]);

	//iterate through the particles
	for (int i = 0; i < num_particles; i++)
	{
		if (fabs(yaw_rate) < 0.00001)
		{
			particles[i].x += velocity * delta_t * cos(particles[i].theta);
			particles[i].y += velocity * delta_t * sin(particles[1].theta);
		}
		else
		{
			particles[i].x += velocity / yaw_rate *(sin(particles[i].theta + yaw_rate *delta_t) - sin(particles[i].theta));
			particles[i].y += velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
			particles[i].theta += yaw_rate*delta_t;
		}

		//factor in the noise
		particles[i].x += distX(gen);
		particles[i].y += distY(gen);
		particles[i].theta += distTheta(gen);
	}



}


void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	for (int i = 0; i < observations.size(); i++)
	{
		//current observation
		LandmarkObs obs = observations[i];
		
		//get min distance to max possible distance
		double minDistance = numeric_limits<double>::max();

		int mapID = -1;
		for (int j = 0; j<predicted.size(); j++)
		{
			LandmarkObs p = predicted[j];

			//current distance
			double currentDistance = dist(obs.x, obs.y, p.x, p.y);

			//get the predicted landmark closest to the current observed landmark

			if( currentDistance < minDistance)
			{
				minDistance = currentDistance;
				mapID = p.id;
			}
		}
		observations[i].id = mapID;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
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

	//iterate through the particles
	for (int i = 0; i < num_particles; i++)
	{
		double partX = particles[i].x;
		double partY = particles[i].y;
		double partTheta = particles[i].theta;

		//use a vector to hold the landmark positions
		vector<LandmarkObs> predictions;

		//iterate through landmarks
		for ( int j = 0; j < map_landmarks.landmark_list.size(); j++)
		{
			float landMarkX = map_landmarks.landmark_list[j].x_f;
			float landMarkY = map_landmarks.landmark_list[j].y_f;
			int landMarkID = map_landmarks.landmark_list[j].id_i;

			//Only look at landmarks within sensor range of the particle
			if (fabs(landMarkX - partX)<= sensor_range && fabs(landMarkY - partY) <= sensor_range)
			{
				predictions.push_back(LandmarkObs{ landMarkID, landMarkX, landMarkY });
			}
		}

		//get a copy of the list of observations, transformed from vehicle to map coords.
		vector<LandmarkObs> transformedObs;

		for (int j = 0; j < transformedObs.size(); j++)
		{
			double transX = cos(partTheta)*observations[j].x - sin(partTheta)*observations[j].y + partX;
			double transY = sin(partTheta)*observations[j].x + cos(partTheta)*observations[j].y + partY;
			transformedObs.push_back(LandmarkObs{ observations[j].id, transX, transY });
		}

		//data associateon for the predictions and transformed observations
		dataAssociation(predictions, transformedObs);

		particles[i].weight = 1.0;

		for (int j = 0; j < transformedObs.size(); j++)
		{
			double obsX, obsY, predX, predY;
			obsX = transformedObs[j].x;
			obsY = transformedObs[j].y;

			int predID = transformedObs[j].id;

			//get the X Y coords of the prediction in relation to the observation
			for (int k = 0; k < predictions.size(); k++)
			{
				if(predictions[k].id == predID)
				{
					predX = predictions[k].x;
					predY = predictions[k].y;
				}
			}

			//calculate the weights with Gaussian
			double stdX = std_landmark[0];
			double stdY = std_landmark[1];
			double observedWeight = (1 / (2 * M_PI*stdX*stdY)) * exp(-(pow(predX - obsX, 2)/(2*pow(stdX,2)) + (pow(predY-obsY,2)/(2*pow(stdY,2)))));

			//use the product with the total observations weights
			particles[i].weight *= observedWeight;
		}		
	}

}




void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	vector<Particle> newParticles;

	//retrieve all of the weights
	vector<double> weights;
	for (int i = 0; i < num_particles; i++)
	{
		weights.push_back(particles[i].weight);
	}

	//generate a random index for the sample wheel
	uniform_int_distribution<int> uniformIntDist(0, num_particles - 1);
	int index = uniformIntDist(gen);

	double maxWeight = *max_element(weights.begin(), weights.end());
	uniform_real_distribution<double> uniformRealDist(0.0, maxWeight);

	double beta = 0.0;

	//
	for (int i = 0; i<num_particles; i++)
	{
		beta += uniformRealDist(gen) * 2.0*maxWeight;
		while(beta > weights[index])
		{
			beta -= weights[index];
			index = (index + 1) % num_particles;
		}
		newParticles.push_back(particles[index]);
	}
	particles = newParticles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

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
