#include "settings.hpp"
#include "graphics.hpp"
#include "simulation.hpp"
#include <iostream>
#include "initial_conditions.hpp"
#include <filesystem>
#include <cstdio>
#include "Particles.hpp"
#include "dye.hpp"

int main() {{

    printf("DIFFUSION NUMBER : %e\n", diffusion_number);

    if (std::filesystem::exists("output")) {
        for (const auto& entry : std::filesystem::recursive_directory_iterator("output")) {
            if (std::filesystem::is_regular_file(entry.path())) {
                std::filesystem::remove(entry.path());
            }
        }
    }
    std::filesystem::create_directory("output");
    std::filesystem::create_directory("output/velocity");
    std::filesystem::create_directory("output/dye");

    std::vector<unsigned char> velocity_image(3*image_width*image_height);
    std::vector<unsigned char> dye_image(3*image_width*image_height);
    std::vector<Cell> grid(grid_width*grid_height);
    std::vector<Particle> particles(num_particles);
    std::vector<Dye> dye(grid_width*grid_height);
    initialize_grid(grid);
    initialize_particles(particles);
    for (int i = 0; i < total_frames; ++i) {
        printf("frame %d / %d (%f %)\n", i, total_frames, (float)i/(float)total_frames*100.0f);
        clear_image(dye_image);
        clear_image(velocity_image);
        draw_dye(dye_image, dye);
        draw_velocity(velocity_image, grid);
        for (int j = 0; j < steps_per_frame; ++j) {
            step_particles(particles, grid);
            draw_particles(velocity_image, particles);
            step_fluid(grid);
            step_dye(dye, grid);
        }
        draw_walls(velocity_image, grid);
        draw_walls(dye_image, grid);
        save_image(velocity_image, "velocity", i);
        save_image(dye_image, "dye", i);
    }

} std::cout << "END OF MAIN" << std::endl; }