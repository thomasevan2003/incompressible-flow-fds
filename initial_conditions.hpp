#ifndef INITIAL_CONDITIONS_HPP
#define INITIAL_CONDITIONS_HPP

#include <vector>
#include "simulation.hpp"
#include "settings.hpp"

void initialize_grid(std::vector<Cell>& grid) {

    // wall positions
    for (int idx = 0; idx < grid_width*grid_height; ++idx) {
        grid[idx].wall = false;
        float x = (float)(idx % grid_width)/(float)(grid_width-1);
        float y = (float)(idx / grid_width)/(float)(grid_height-1);
        // top and bottom borders
        if (y < 0.01f || y > 0.99f) {
            grid[idx].wall = true;
        }
        if (y < 0.8f && x < 0.05f) {
            grid[idx].wall = true;
        }
        if (y > 0.2f && x > 0.95f) {
            grid[idx].wall = true;
        }
        if (x > 0.28f && x < 0.38f && y > 0.3f) {
            grid[idx].wall = true;
        }
        if (x > 0.62f && x < 0.72f && y < 0.7f) {
            grid[idx].wall = true;
        }
    }

}

#endif