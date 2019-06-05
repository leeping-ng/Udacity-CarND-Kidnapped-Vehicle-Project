/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using namespace std;
static default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) 
{
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 100;  // TODO: Set the number of particles
  
  // Create Gaussian distributions for x, y and theta for random noise
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);  
  
  // Initialize particles
  for (int i = 0; i < num_particles; ++i) 
  {
    Particle p;
    
    // Fill in details on Particle struct
    p.id = i;
    p.x = dist_x(gen);
    p.y = dist_y(gen);
    p.theta = dist_theta(gen);
    p.weight = 1.0;
    
    // Add to set of current particles
    particles.push_back(p);
  }
 
  is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) 
{
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

  // Create Gaussian distributions for x, y and theta for random noise
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);
  
  for (int i = 0; i < num_particles; ++i) 
  {
    // Equations for condition when yaw rate is not zero
    if (abs(yaw_rate) != 0) 
    {
    	particles[i].x += (velocity/yaw_rate) * (sin(particles[i].theta + (yaw_rate * delta_t)) - sin(particles[i].theta));
    	particles[i].y += (velocity/yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + (yaw_rate * delta_t)));
    	particles[i].theta += yaw_rate * delta_t;
    }
    
    // Equations for condition when yaw rate is zero
    else 
    {
    	particles[i].x += velocity * delta_t * cos(particles[i].theta);
    	particles[i].y += velocity * delta_t * sin(particles[i].theta);
    }
    
    // Add Gaussian noise to the particles
    particles[i].x += dist_x(gen);
    particles[i].y += dist_y(gen);
    particles[i].theta += dist_theta(gen);    
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) 
{
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  
  // "predicted" contains the prediction measurements between one particular particle and all map landmarks within sensor range
  // "observations" contains actual landmark measurements gathered from LiDAR
  
  // for each observation
  // good practice to declare variables as unsigned if they will be compared to sizes, else warning will be flagged
  for (unsigned int i=0; i < observations.size(); i++) 
  {
    // initialize to largest possible value for type double
    double min_dist = std::numeric_limits<double>::max();
    // initialize mapID
    int mapID = -1;
    
    // for each observation, compare with each prediction
    for (unsigned int j=0; j < predicted.size(); j++) 
    {
    	// Use the helper function dist() that returns sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1))
      	double distance = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
      
      	if (distance < min_dist) 
        {
        	min_dist = distance;
          	mapID = predicted[j].id;
        }
      	// assign observed measurement to this landmark
      	observations[i].id = mapID;
    }     
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  
  // Loop for one particle at a time
  for (int i=0; i < num_particles; i++) 
  {
  	
    double p_x = particles[i].x;
    double p_y = particles[i].y;
    double p_theta = particles[i].theta;
    
    // This section creates a vector to hold map landmark locations that are within sensor range of particle
    // In map coordinate system
    vector<LandmarkObs> predictions;
    for (unsigned int j=0; j < map_landmarks.landmark_list.size(); j++) 
    {
   		float lm_x = map_landmarks.landmark_list[j].x_f;
      	float lm_y = map_landmarks.landmark_list[j].y_f;
      	int lm_id = map_landmarks.landmark_list[j].id_i;
      
      	// Use helper function dist() to check if its within range
      	if (dist(lm_x, lm_y, p_x, p_y) <= sensor_range) 
        {
        	predictions.push_back(LandmarkObs{lm_id, lm_x, lm_y});
        }
    }
    
    // This section creates a vector to hold observations transformed from vehicle coordinates to map coordinates
    vector<LandmarkObs> transformed_obs;
    for (unsigned int j = 0; j < observations.size(); j++) 
    {
      double t_x = p_x + cos(p_theta)*observations[j].x - sin(p_theta)*observations[j].y;
      double t_y = p_y + sin(p_theta)*observations[j].x + cos(p_theta)*observations[j].y;
      transformed_obs.push_back(LandmarkObs{observations[j].id, t_x, t_y});
    }
    
    // This section finds the predicted measurement that is closest to each observed measurement 
    // and assigns the id of this particular landmark to the transformed observation
    dataAssociation(predictions, transformed_obs);
    
    // Initialize particle weight to 1
    particles[i].weight = 1.0;
    
    // This section updates the weight of each particle
    for (unsigned int j=0; j < transformed_obs.size(); j++)
    {
    	double o_x, o_y, pred_x, pred_y;
      	o_x = transformed_obs[j].x;
      	o_y = transformed_obs[j].y;
      	
      	// For each observation, compare with each prediction
        for (unsigned int k = 0; k < predictions.size(); k++) 
        {
          // Match the observation with prediction using id
          if (predictions[k].id == transformed_obs[j].id) 
          {
            pred_x = predictions[k].x;
            pred_y = predictions[k].y;
          }
        }
      	
      	// The whole point of the above was to find, for each observation, the x and y value
      	// and for each prediction, the x and y value, and with these, 
      	// calculate the observation weight with multivariate Gaussian
      	double s_x = std_landmark[0];
      	double s_y = std_landmark[1];
      	double obs_w = ( 1/(2*M_PI*s_x*s_y)) * exp( -( pow(pred_x-o_x,2)/(2*pow(s_x, 2)) + (pow(pred_y-o_y,2)/(2*pow(s_y, 2))) ) );
      
      	// Calculate final weight of observation by multiplying all weights
      	particles[i].weight *= obs_w;
    }
  }
}

void ParticleFilter::resample() 
{
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  vector<Particle> resampled_particles;
  
  // This section gets all current weights and max weight
  vector<double> weights;
  double max_weight = numeric_limits<double>::min();
  for (int i = 0; i < num_particles; i++) 
  {
    weights.push_back(particles[i].weight);
   	
    if(particles[i].weight > max_weight) 
    {
		max_weight = particles[i].weight;
    }
  }
  
  // generate random starting index for resampling wheel
  uniform_int_distribution<int> intdist(0, num_particles-1);
  auto index = intdist(gen);
  uniform_real_distribution<double> realdist(0.0, max_weight);

  double beta = 0.0;

  // This section moves along the resampling wheel
  // The particles with higher weights have a higher change of being resampled
  for (int i=0; i < num_particles; i++) 
  {
    beta += realdist(gen) * 2.0;
    while (beta > weights[index]) 
    {
      beta -= weights[index];
      index = (index + 1) % num_particles;
    }
    resampled_particles.push_back(particles[index]);
  }
  particles = resampled_particles;

}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}