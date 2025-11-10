# FPGA AESA Radar Control, Simulation, and Visualization

Welcome to my github repository for my masters project at SSU. Currently the repository houses Simulink files which I will be using to program my Artix-7 FPGA, and the Matlab scripts for initializing the variables that the Simulink blocks use.

**Tools and Software used in this project:**
- Simulink/Matlab 2024a
- Vitis Model Composer 2025.1
- Vivado 2025.1
- Artix-7 AC701 Development Board

Inside Main are the Matlab scripts, and the Outdated folder is simply older scripts and Simulink files in case I need to reference them later.
InitializationOfVariablesRev00x.m is a script that needs to be ran before opening the Simulink file, as it initializes all the variables contained in the blocks.
RadarSimulinkProjectRev00x.slx is the main Simulink file containing the block model
The .slxc files can be ignored.

Current Progress:
- Generating simulated signal (X)
- Range Processing (X)
- Velocity Calculation (X)
- Angle Estimation ( )
- Display in Matlab ( )
- Testing ( )

What is left to do:
- Review delays and latencies, rename blocks
- Write Matlab interface code
- Test

More will be added to this readme later.
