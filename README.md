# Unscented Kalman Filter 
## Target Tracking using Lidar and Radar Measurements

In this project utilize an Unscented Kalman Filter to estimate the state of a moving object of interest with noisy lidar and radar measurements.  

## Simulator Details:

This project involves the Term 2 Simulator which can be downloaded [here](https://github.com/udacity/self-driving-car-sim/releases)

## Code Details:
Code for the project are src/ukf.cpp, src/ukf.h, tools.cpp, and tools.h
main.cpp uses for uWebSocketIO in communicating with the simulator.


1. mkdir build
2. cd build
3. cmake ..
4. make
5. ./UnscentedKF


## INPUT: 
values provided by the simulator to the c++ program

["sensor_measurement"] => the measurment that the simulator observed (either lidar or radar)


## OUTPUT: 
values provided by the c++ program to the simulator

["estimate_x"] <= kalman filter estimated position x
["estimate_y"] <= kalman filter estimated position y
["rmse_x"]
["rmse_y"]
["rmse_vx"]
["rmse_vy"]

### Dataset1:
![Dataset1](https://github.com/aranga81/Unscented-Kalman-Filter/tree/master/results/dataset1.png)

### Dataset2:
![Dataset2](https://github.com/aranga81/Unscented-Kalman-Filter/tree/master/results/dataset2.png)

---

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./UnscentedKF` The current state uses i/o from the simulator.

