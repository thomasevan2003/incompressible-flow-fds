# incompressible-flow-fds

This program generates a visualization of incompressible flow through a series of channels. The flow solution is computed using an explicit scheme on a staggered finite difference grid with pressure resolved using Chorin's projection method. Pressure is computed using iterative relaxation of the Poisson equation. The flow velocity is visualized using tracer particles superimposed on a velocity field plot as well as dye tracking from the inlet. The flow of dye across the domain is computed at each step with an upwind advection scheme using the current velocity solution. Computation on the grid is parallelized using OpenMP.

To use this program, compile and run main.cpp. The resulting images will be placed in the output folder. The program will take several hours to run. The video files can then be generated using tools such as ffmpeg.

https://github.com/user-attachments/assets/422d6a9c-0a11-45e5-bfe8-64bfde26468e

https://github.com/user-attachments/assets/26eaa9bb-7338-43d2-bde4-0743e7fc76cb
