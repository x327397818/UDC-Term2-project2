# Unscented Kalman Filter Project

[//]: # (Image References)

[Output_img]: ./output/Result.PNG "Output_img"
[NIS_img]: ./output/NIS.png "NIS_img"
 


Self-Driving Car Engineer Nanodegree Program

In this project a unscented kalman filter is utilized to estimate the state of a moving object of interest with noisy lidar and radar measurements. Passing the project requires obtaining RMSE values that are lower than the tolerance outlined in the project rubric. 

This project involves the Term 2 Simulator which can be downloaded [here](https://github.com/udacity/self-driving-car-sim/releases)


---

## Dependencies

* cmake >= 3.5
  * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1 (Linux, Mac), 3.81 (Windows)
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools](https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make` 
   * On windows, you may need to run: `cmake .. -G "Unix Makefiles" && make`
4. Run it: `./UnscentedKF `

## Editor Settings

We've purposefully kept editor configuration files out of this repo in order to
keep it as simple and environment agnostic as possible. However, we recommend
using the following settings:

* indent using spaces
* set tab width to 2 spaces (keeps the matrices in source code aligned)

## Code Style

Please (do your best to) stick to [Google's C++ style guide](https://google.github.io/styleguide/cppguide.html).

## Generating Additional Data

This is optional!

If you'd like to generate your own radar and lidar data, see the
[utilities repo](https://github.com/udacity/CarND-Mercedes-SF-Utilities) for
Matlab scripts that can generate additional data.

## Result

Here is the result of simulation based on given dataset `obj_pose-laser-radar-synthetic-input.txt`

![alt text][Output_img]

```
Red circles: lidar measurements.
Blue circles: meansurements.
Green marker: Kalman filter output.
```

## Accuracy

`RMSE: [px, py, vx, vy] = [0.0661, 0.0829, 0.3283, 0.2286]`

Meet project rubric `RMSE <= [.11, .11, 0.52, 0.52]` 

And it is also better than the Extended Kalman Filter result `[0.0973, 0.0855, 0.4513, 0.4399]`

## Variances choise and NIS
The values for the process noise 

`std_a_ = 1.0 `

`std_yawdd_ = 0.5`

They are choosen based on the NIS result. Here is the final result of NIS.

![alt text][NIS_img]