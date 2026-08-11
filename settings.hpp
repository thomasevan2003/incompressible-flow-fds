#ifndef SETTINGS_HPP
#define SETTINGS_HPP

#include <cmath>

typedef double real;

// physics parameters
constexpr real V_inlet = 1.0; // [m/s]
constexpr real domain_width = 1.0; // [m]
constexpr real Re = 5'000.0;
constexpr real rho = 1.0; // [kg/m^3]
constexpr real mu = rho*V_inlet*domain_width/Re; // [kg/m*s]
constexpr real nu = mu/rho; // [m^2/s]

// simulation parameters
constexpr real total_time = 10.0;
constexpr int grid_width = 900;
constexpr int grid_height = 600;
constexpr real dx = domain_width/(real)(grid_width-1);
constexpr real domain_height = (real)(grid_height-1)*dx;
constexpr real dt = 0.0001; // [s]
constexpr real diffusion_number = 4.0*(nu*dt/dx/dx); // stability of viscous term, should be less than 1.0
constexpr int num_pressure_iterations = 2; 

// visualization settings
constexpr int num_particles = 2000;
constexpr real particle_lifetime = 1.0;
constexpr real particle_speed_multiplier = 1.0;

// output settings
constexpr real frame_time = 0.01;
constexpr int total_frames = (int)std::round(total_time/frame_time);
constexpr int steps_per_frame = (int)std::round(frame_time/dt);
constexpr int image_width = 1500;
constexpr int image_height = 1000;
constexpr real video_length = total_frames/(real)30.0;
constexpr real video_time_per_step = video_length/(real)(total_frames*steps_per_frame);

#endif