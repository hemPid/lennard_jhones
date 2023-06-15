#pragma once
#include <cstdint>
#include <string>
//mass of particle
const double mass = 1;
//const double sigma = 3.4e-10;
//const double epsilon = 1.66e-21;
//dt variable, which used in equations
const double dtime = 0.001;
//paths to log files
const std::string energy_path = "C:\\Users\\iynb2\\lenard_jhones_system_v2\\logs\\energy.txt";
const std::string velocity_path = "C:\\Users\\iynb2\\lenard_jhones_system_v2\\logs\\vels.txt";
const std::string pos_path = "C:\\Users\\iynb2\\lenard_jhones_system_v2\\logs\\pos.xyz";
const std::string distance_path = "C:\\Users\\iynb2\\lenard_jhones_system_v2\\logs\\min_distance.txt";
const std::string free_run_path = "C:\\Users\\iynb2\\lenard_jhones_system_v2\\logs\\free_run.txt";
const std::string diffusion_path = "C:\\Users\\iynb2\\lenard_jhones_system_v2\\logs\\diffusion.txt";