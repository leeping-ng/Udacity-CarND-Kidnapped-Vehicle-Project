# Udacity-CarND-Kidnapped-Vehicle-Project
A particle filter for localizing a vehicle using LiDAR measurements, for the Udacity Self-Driving Car Engineer Nanodegree.


A particle filter is a technique for localising, or finding out a robot's current location. The image below from Udacity shows a particle filter at work, where the red dots are the "particles" guessing the location, randomly scattered around the map initially. Each particle has x, y coordinates, and angle orientation.
<img src='https://github.com/leeping-ng/Udacity-CarND-Kidnapped-Vehicle-Project/blob/master/images/Particle%20Filter%20Diagram.png'>

There are 4 main steps in this particle filter implementation:
1. Initialisation: The position is estimated from GPS, which is usually not as accurate as we want
2. Prediction: The control input (velocity and yaw rate) are added for all particles
3. Update: Particle weights are updated using map landmark positions and feature measurements
4. Resampling: A resampling wheel is used, whereby particles with higher weights have a higher chance of surviving

The image below shows the successful implementation of the particle filter using a simulator.

<img src='https://github.com/leeping-ng/Udacity-CarND-Kidnapped-Vehicle-Project/blob/master/images/Simulator.JPG'>
