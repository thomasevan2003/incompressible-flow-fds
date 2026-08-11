#ifndef DYE_HPP
#define DYE_HPP

#include "settings.hpp"
#include "simulation.hpp"
#include <omp.h>

struct Dye {
    real red;
    real blue;
    real red_next;
    real blue_next;
};

void apply_boundary_conditions(std::vector<Dye>& dye) {
    for (int j = 0.8*grid_height; j < grid_height; ++j) {
        float frac = ((float)j - 0.8f*grid_height)/(0.2f*grid_height);
        if (frac > 0.5f) {
            dye[0 + j*grid_width].red = 0.95;
        } else {
            dye[0 + j*grid_width].blue = 0.95;
        }
        dye[grid_width-1 + j*grid_width].red = dye[grid_width-2 + j*grid_width].red;
        dye[grid_width-1 + j*grid_width].blue = dye[grid_width-2 + j*grid_width].blue;
    }
}

void step_dye(std::vector<Dye>& dye, std::vector<Cell>& grid) {
    apply_boundary_conditions(dye);
    #pragma omp parallel for
    for (int i = 1; i < grid_width-1; ++i) {
        for (int j = 1; j < grid_height-1; ++j) {
            int idx = i + j*grid_width;
            real red_e;
            real red_w;
            real red_n;
            real red_s;
            real blue_e;
            real blue_w;
            real blue_n;
            real blue_s;
            if (grid[idx].u > 0.0) {
                red_w = dye[idx-1].red;
                blue_w = dye[idx-1].blue;
            } else {
                red_w = dye[idx].red;
                blue_w = dye[idx].blue;
            }
            if (grid[idx+1].u > 0.0) {
                red_e = dye[idx].red;
                blue_e = dye[idx].blue;
            } else {
                red_e = dye[idx+1].red;
                blue_e = dye[idx+1].blue;
            }
            if (grid[idx].v > 0.0) {
                red_s = dye[idx-grid_width].red;
                blue_s = dye[idx-grid_width].blue;
            } else {
                red_s = dye[idx].red;
                blue_s = dye[idx].blue;
            }
            if (grid[idx+grid_width].v > 0.0) {
                red_n = dye[idx].red;
                blue_n = dye[idx].blue;
            } else {
                red_n = dye[idx+grid_width].red;
                blue_n = dye[idx+grid_width].blue;
            }
            dye[idx].red_next = dye[idx].red - dt/dx*(grid[idx+1].u*red_e - grid[idx].u*red_w + 
                                              grid[idx+grid_width].v*red_n - grid[idx].v*red_s);
            dye[idx].blue_next = dye[idx].blue - dt/dx*(grid[idx+1].u*blue_e - grid[idx].u*blue_w + 
                                              grid[idx+grid_width].v*blue_n - grid[idx].v*blue_s);
        }
    }
    #pragma omp parallel for
    for (int i = 1; i < grid_width - 1; ++i) {
        for (int j = 1; j < grid_height - 1; ++j) {
            int idx = i + j*grid_width;
            dye[idx].red = dye[idx].red_next;
            dye[idx].blue = dye[idx].blue_next;
        }
    }
}

#endif