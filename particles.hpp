#ifndef PARTICLES_HPP
#define PARTICLES_HPP

#include "settings.hpp"
#include <vector>
#include <random>

struct Particle {
    real x;
    real y;
    real age;
};

real getRand() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<real> dist(0.0, 1.0);
    return dist(gen);
}

void initialize_particles(std::vector<Particle>& particles) {
    for (int i = 0; i < num_particles; ++i) { 
        particles[i].x = getRand()*domain_width;
        particles[i].y = getRand()*domain_height;
        particles[i].age = getRand()*particle_lifetime;
    }
}

void step_particles(std::vector<Particle>& particles, std::vector<Cell>& grid) {
    for (int k = 0; k < num_particles; ++k) {
        particles[k].age += video_time_per_step;
        if (particles[k].age > particle_lifetime) {
            particles[k].age = 0.0;
            particles[k].x = getRand()*domain_width;
            particles[k].y = getRand()*domain_height;
        }
        int i = (particles[k].x/domain_width)*grid_width;
        int j = (particles[k].y/domain_height)*grid_height;
        if (i >= 0 && i < grid_width && j >= 0 && j < grid_height) {
            int idx = i + j*grid_width;
            particles[k].x += grid[idx].u*dt*particle_speed_multiplier;
            particles[k].y += grid[idx].v*dt*particle_speed_multiplier;
        }
    }
}

#endif