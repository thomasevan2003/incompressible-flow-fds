#ifndef GRAPHICS_HPP
#define GRAPHICS_HPP

#include "simulation.hpp"
#include "settings.hpp"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <string>
#include <vector>
#include <iostream>
#include "Particles.hpp"
#include <cstring>
#include "dye.hpp"

void clear_image(std::vector<unsigned char>& image) {
    std::memset(image.data(), 0, image.size()*sizeof(unsigned char));
}

void draw_walls(std::vector<unsigned char>& image, const std::vector<Cell>& grid) {
    #pragma omp parallel for
    for (int idx = 0; idx < image_width*image_height; ++idx) {
        int i = idx % image_width;
        int j = idx / image_width;
        int m = (int)(((float)i+0.5f)*(float)(grid_width-1)/(float)image_width);
        int n = (int)(((float)j+0.5f)*(float)(grid_height-1)/(float)image_height);
        int gdx = m + n*grid_width;
        if (grid[gdx].wall) {
            image[3*idx] = 128;
            image[3*idx+1] = 128;
            image[3*idx+2] = 128; 
        }
    }
}

void draw_particles(std::vector<unsigned char>& image, const std::vector<Particle>& particles) {
    #pragma omp parallel for
    for (int k = 0; k < num_particles; ++k) {
        int i = (particles[k].x/domain_width)*image_width;
        int j = (particles[k].y/domain_height)*image_height;
        if (i >= 0 && i < image_width && j >= 0 && j < image_height) {
            int idx = i + j*image_width;
            image[3*idx] = 255;
            image[3*idx+1] = 255;
            image[3*idx+2] = 255;
        }
    }
}

void draw_velocity(std::vector<unsigned char>& image, const std::vector<Cell>& grid) {
    float V_max = 0.0f;
    for (int idx = 0; idx < grid_width*grid_height; ++idx) {
        real V = sqrt(grid[idx].u*grid[idx].u + grid[idx].v*grid[idx].v);
        if (V > V_max) {
            V_max = V;
        }
    }
    #pragma omp parallel for
    for (int idx = 0; idx < image_width*image_height; ++idx) {
        int i = idx % image_width;
        int j = idx / image_width;
        int m = (int)(((float)i+0.5f)*(float)grid_width/(float)image_width);
        int n = (int)(((float)j+0.5f)*(float)grid_height/(float)image_height);
        int gdx = m + n*grid_width;
        real u = grid[gdx].u;
        real v = grid[gdx].v;
        real V = sqrt(u*u+v*v);
        if (V/V_max < 0.66) {
            image[3*idx] = 255*(V/V_max*1.5);
        } else {
            image[3*idx] = 255;
            image[3*idx+1] = 255*(V-0.66*V_max)/(0.66*V_max);
        }
    }
}

void draw_dye(std::vector<unsigned char>& image, std::vector<Dye>& dye) {
    for (int idx = 0; idx < image_width*image_height; ++idx) {
        int i = idx % image_width;
        int j = idx / image_width;
        int m = (int)(((float)i+0.5f)*(float)(grid_width-1)/(float)image_width);
        int n = (int)(((float)j+0.5f)*(float)(grid_height-1)/(float)image_height);
        int gdx = m + n*grid_width;
        real red = dye[gdx].red;
        real blue = dye[gdx].blue;
        if (red > 0.0) {
            if (red < 0.66) {
                image[3*idx] = 255*(red/0.66);
            } else if (red < 0.99) {
                image[3*idx] = 255;
                image[3*idx+1] = 255*(red-0.66)/0.66;
            } else {
                image[3*idx] = 255;
                image[3*idx+1] = 128;
            }
        }
        if (blue > 0.0) {
            if (blue < 0.66) {
                image[3*idx+2] = 255*(blue/0.66);
            } else if (blue < 0.99) {
                image[3*idx+2] = 255;
                image[3*idx+1] = std::max((unsigned char)(255*(blue-0.66)/0.66), image[3*idx+1]);
            } else {
                image[3*idx+2] = 255;
                image[3*idx+1] = 128;
            }
        }
    }
}

void draw_pressure(std::vector<unsigned char>& image, const std::vector<Cell>& grid) {
    real p_min = grid[grid_width/4 + grid_height*grid_width/2].p;
    real p_max = grid[grid_width/4 + grid_height*grid_width/2].p;
    for (int idx = 0; idx < grid_width*grid_height; ++idx) {
        int i = idx % grid_width;
        int j = idx / grid_width;
        if (grid[idx].wall || i == grid_width-1 || j == grid_height - 1) {
            continue;
        }
        if (grid[idx].p < p_min) {
            p_min = grid[idx].p;
        }
        if (grid[idx].p > p_max) {
            p_max = grid[idx].p;
        }
    }
    #pragma omp parallel for
    for (int idx = 0; idx < image_width*image_height; ++idx) {
        int i = idx % image_width;
        int j = idx / image_width;
        int m = (int)(((float)i+0.5f)*(float)(grid_width-1)/(float)image_width);
        int n = (int)(((float)j+0.5f)*(float)(grid_height-1)/(float)image_height);
        int gdx = m + n*grid_width;
        real p = grid[gdx].p;
        real frac = 2.0*(p - p_min)/(p_max - p_min) - 1.0;
        if (frac > 0.0) {
            image[3*idx] = (unsigned char)(255*frac);
        } else {
            image[3*idx+2] = (unsigned char)(-255*frac);
        }
    }
}

void save_image(std::vector<unsigned char>& image, const std::string& dir, int frame) {
    char filename[256];
    sprintf(filename, "output/%s/%d.jpg", dir.c_str(), frame);
    stbi_write_jpg(filename, image_width, image_height, 3, image.data(), 90);
}

#endif