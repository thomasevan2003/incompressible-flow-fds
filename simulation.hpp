#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "settings.hpp"
#include <vector>
#include <iostream>
#include <omp.h>
#include <random>

struct Cell {
    real p;
    real u;
    real v;
    real u_star;
    real v_star;
    real poisson;
    real p_next;
    bool wall;
};

void apply_inlet_velocity_boundary_condition(std::vector<Cell>& grid) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<real> dist(0.0, 1.0);
    #pragma omp parallel for
    for (int j = 0; j < grid_height; ++j) {
        int idx = j*grid_width;
        grid[idx].u = V_inlet*(1.0 + 0.05*dist(gen));
        grid[idx].v = 0.0;
    }
}

void apply_outlet_velocity_boundary_condition(std::vector<Cell>& grid) {
    #pragma omp parallel for
    for (int j = 0; j < grid_height; ++j) {
        int idx = j*grid_width + grid_width-1;
        grid[idx].u = V_inlet;
        grid[idx].v = 0.0;
    }
}

void compute_intermediate_velocity(std::vector<Cell>& grid) {
    #pragma omp parallel for
    for (int i = 1; i < grid_width-1; ++i) {
        for (int j = 1; j < grid_height-1; ++j) {
            int ix = i + j*grid_width;

            // indices of adjacent cells
            int wx = ix - 1;
            int ex = ix + 1;
            int nx = ix + grid_width;
            int sx = ix - grid_width;
            int nw = nx - 1;
            int ne = nx + 1;
            int sw = sx - 1;
            int se = sx + 1;

            // get intermediate velocity using viscosity and advection terms
            real laplace_u = (grid[wx].u + grid[ex].u + grid[nx].u + grid[sx].u - 4*grid[ix].u)/(dx*dx);
            real laplace_v = (grid[wx].v + grid[ex].v + grid[nx].v + grid[sx].v - 4*grid[ix].v)/(dx*dx);
            // get velocity gradients using upwind scheme
            real du_dx;
            if (grid[ix].u > 0.0) {
                du_dx = (grid[ix].u - grid[wx].u)/dx;
            } else {
                du_dx = (grid[ex].u - grid[ix].u)/dx;
            }
            real v_interp = ((real)0.25)*(grid[wx].v + grid[ix].v + grid[sx].v + grid[sw].v);
            real du_dy;
            if (v_interp > 0.0) {
                du_dy = (grid[ix].u - grid[sx].u)/dx;
            } else {
                du_dy = (grid[nx].u - grid[ix].u)/dx;
            }
            real u_interp = ((real)0.25)*(grid[sx].u + grid[se].u + grid[ix].u + grid[ex].u);
            real dv_dx;
            if (u_interp > 0.0) {
                dv_dx = (grid[ix].v - grid[wx].v)/dx;
            } else {
                dv_dx = (grid[ex].v - grid[ix].v)/dx;
            }
            real dv_dy;
            if (grid[ix].v > 0.0) {
                dv_dy = (grid[ix].v - grid[sx].v)/dx;
            } else {
                dv_dy = (grid[nx].v - grid[ix].v)/dx;
            }
            real du_dt = nu*laplace_u - (grid[ix].u*du_dx + v_interp*du_dy);
            real dv_dt = nu*laplace_v - (u_interp*dv_dx + grid[ix].v*dv_dy);
            grid[ix].u_star = grid[ix].u + dt*du_dt;
            grid[ix].v_star = grid[ix].v + dt*dv_dt;

            // no slip condition
            if (grid[ix].wall || grid[wx].wall || grid[sw].wall || 
                grid[sx].wall || grid[nw].wall || grid[nx].wall) {
                grid[ix].u_star = 0.0;
            }
            if (grid[sw].wall || grid[sx].wall || grid[se].wall ||
                grid[wx].wall || grid[ix].wall || grid[ex].wall) {
                grid[ix].v_star = 0.0;
            }

        }
    }
}

void compute_poisson_rhs(std::vector<Cell>& grid) {
    #pragma omp parallel for
    for (int i = 1; i < grid_width-1; ++i) {
        for (int j = 1; j < grid_height-1; ++j) {
            int ix = i + j*grid_width;
            grid[ix].poisson = rho/dt*(
                grid[ix+1].u_star - grid[ix].u_star + grid[ix+grid_width].v_star - grid[ix].v_star
            ) / dx;
        }
    }
}

void iterate_pressure_poisson(std::vector<Cell>& grid) {
    for (int offset = 0; offset < 2; ++offset) {
        #pragma omp parallel for
        for (int ix = 0; ix < grid_width*grid_height; ++ix) {

            int i = ix % grid_width;
            int j = ix / grid_width;

            // checkerboard iteration
            if ((i+j) % 2 == offset) {
                continue;
            }

            // skip last row and column to account for staggered grid
            if (i == grid_width-1 || j == grid_height-1) {
                continue;
            }
 
            // skip top and bottom
            if (j == 0 || j == grid_height - 1) {
                continue;
            }

            // inlet condition
            if (i == 0) {
                grid[ix].p_next = grid[ix+1].p;
                continue;
            }

            // outlet condition (dp/dn = 0)
            if (i == grid_width-2) {
                grid[ix].p_next = grid[ix-1].p;
                continue;
            }

            // wall boundary conditions (dp/dn = 0)
            if (grid[ix].wall) {
                if (!grid[ix-1].wall) {
                    grid[ix].p_next = grid[ix-1].p;
                } 
                if (!grid[ix+1].wall) {
                    grid[ix].p_next = grid[ix+1].p;
                } 
                if (!grid[ix-grid_width].wall) {
                    grid[ix].p_next = grid[ix-grid_width].p;
                } 
                if (!grid[ix+grid_width].wall) {
                    grid[ix].p_next = grid[ix+grid_width].p;
                }
                continue;
            }

            // poisson relaxation method
            grid[ix].p_next = ((real)0.25)*(
                grid[ix-1].p + grid[ix+1].p + grid[ix-grid_width].p + grid[ix+grid_width].p
                - grid[ix].poisson*(dx*dx)
            );
        }
        #pragma omp parallel for
        for (int ix = 0; ix < grid_width*grid_height; ++ix) {
            int i = ix % grid_width;
            int j = ix / grid_width;
            if ((i+j) % 2 == offset) {
                continue;
            }
            grid[ix].p = grid[ix].p_next;
        }
    }
}

void compute_final_velocity(std::vector<Cell>& grid) {
    #pragma omp parallel for
    for (int i = 1; i < grid_width-1; ++i) {
        for (int j = 1; j < grid_height-1; ++j) {
            int ix = i + j*grid_width;
            grid[ix].u = grid[ix].u_star - (dt/rho)*(grid[ix].p - grid[ix-1].p)/dx;
            grid[ix].v = grid[ix].v_star - (dt/rho)*(grid[ix].p - grid[ix-grid_width].p)/dx;
            // no slip condition
            int wx = ix - 1;
            int ex = ix + 1;
            int nx = ix + grid_width;
            int sx = ix - grid_width;
            int nw = nx - 1;
            int ne = nx + 1;
            int sw = sx - 1;
            int se = sx + 1;
            if (grid[ix].wall || grid[wx].wall || grid[sw].wall || 
                grid[sx].wall || grid[nw].wall || grid[nx].wall) {
                grid[ix].u = 0.0;
            }
            if (grid[sw].wall || grid[sx].wall || grid[se].wall ||
                grid[wx].wall || grid[ix].wall || grid[ex].wall) {
                grid[ix].v = 0.0;
            }
        }
    }
}

void step_fluid(std::vector<Cell>& grid) {

    apply_inlet_velocity_boundary_condition(grid);
    apply_outlet_velocity_boundary_condition(grid);
    compute_intermediate_velocity(grid);
    compute_poisson_rhs(grid);
    for (int k = 0; k < num_pressure_iterations; ++k) {
        iterate_pressure_poisson(grid);
    }
    compute_final_velocity(grid);

}

#endif