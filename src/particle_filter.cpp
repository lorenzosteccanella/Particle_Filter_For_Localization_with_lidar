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

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	num_particles=100;
	
	std::default_random_engine gen;
	
	std::normal_distribution<double> N_x(x, std[0]);
	std::normal_distribution<double> N_y(y, std[1]);
	std::normal_distribution<double> N_theta(theta, std[2]);
	
	for( int i = 0; i < num_particles; i++)
	{
		Particle particle;
		particle.id=i;
		particle.x=N_x(gen);
		particle.y=N_y(gen);
		particle.theta = N_theta(gen);
		particle.weight = 1;

		particles.push_back(particle);
		weights.push_back(1);
	}
	is_initialized = true;
	 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	std::default_random_engine gen;
	for (int i = 0; i < num_particles; i++)
	{
		double new_x;
		double new_y;
		double new_theta;
		
		if (yaw_rate == 0)
		{
			new_x = particles[i].x + velocity*cos(particles[i].theta)*delta_t;
			new_y = particles[i].y + velocity*sin(particles[i].theta)*delta_t;
			new_theta = particles[i].theta;
		}

		else
		{
			new_x = particles[i].x + velocity/yaw_rate*(sin(particles[i].theta+yaw_rate*delta_t)-sin(particles[i].theta));
			new_y = particles[i].y + velocity/yaw_rate*(cos(particles[i].theta)-cos(particles[i].theta + yaw_rate*delta_t));
			new_theta = particles[i].theta + yaw_rate*delta_t;
		}
		
		std::normal_distribution<double> N_x(new_x, std_pos[0]);
		std::normal_distribution<double> N_y(new_y, std_pos[1]);
		std::normal_distribution<double> N_theta(new_theta, std_pos[2]);

		particles[i].x = N_x(gen);
		particles[i].y = N_y(gen);
		particles[i].theta = N_theta(gen);

	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.


	for (int i = 0; i < observations.size(); i++)
	{
		double min_dist= 99999999999.0;
		int index=-1;
		for (int j = 0; j < predicted.size(); j++)
		{
			double distance = dist(observations[i].x,observations[i].y,predicted[j].x,predicted[j].y);
			if(distance<min_dist){
				min_dist=distance;
				index=j;
			}
		}
		
		observations[i].id = predicted[index].id;
		
	}

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], std::vector<LandmarkObs> observations, Map map_landmarks) {

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

	for (int i = 0; i < num_particles; i++)
	{

		vector<LandmarkObs> Observation_mapTransformed;
		for (int j = 0; j < observations.size(); j++)
		{
			double obs_mapt_x = particles[i].x + (cos(particles[i].theta) * observations[j].x - sin(particles[i].theta) * observations[j].y);
			double obs_mapt_y = particles[i].y + (sin(particles[i].theta) * observations[j].x + cos(particles[i].theta) * observations[j].y);
			Observation_mapTransformed.push_back(LandmarkObs{observations[j].id, obs_mapt_x, obs_mapt_y});

		}
		
		vector<LandmarkObs> predictions;
		for (int j = 0; j < map_landmarks.landmark_list.size(); j++)
		{
			predictions.push_back(LandmarkObs{map_landmarks.landmark_list[j].id_i, map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f});	
		}

		dataAssociation(predictions, Observation_mapTransformed);

		particles[i].weight = 1.0;

		double mu_x, mu_y, x, y, sigma_x, sigma_y;

		for (int j = 0; j < Observation_mapTransformed.size(); j++)
		{
			
			for (int k = 0; k < predictions.size(); k++)
			{
				if (predictions[k].id == Observation_mapTransformed[j].id)
				{
					x = predictions[k].x;
					y = predictions[k].y;
				}
			}
			mu_x = Observation_mapTransformed[j].x;
			mu_y = Observation_mapTransformed[j].y;
			sigma_x = std_landmark[0];
			sigma_y = std_landmark[1];

			double w = (1 / (2 * M_PI * sigma_x * sigma_y)) * exp(-(pow(x - mu_x, 2) / (2 * pow(sigma_x, 2)) + (pow(y - mu_y, 2) / (2 * pow(sigma_y, 2)))));
			particles[i].weight *= w;
			
			
		}
		
		
		
	}

	
}

void ParticleFilter::resample()
{
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	std::default_random_engine gen;

	vector<Particle> particles_updated;
	vector<double> weights;
	for (int i = 0; i < num_particles; i++)
	{
		weights.push_back(particles[i].weight);
	}

	std::discrete_distribution<> d(weights.begin(), weights.end());

	for (int i = 0; i < particles.size(); i++) {
        	particles_updated.push_back(particles[d(gen)]);
	}

	particles = particles_updated;
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
