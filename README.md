# FPGA AESA Radar Control, Simulation, and Visualization
(Updated 2/9/2026)

Welcome to my github repository for my masters project at SSU. Currently the repository houses Simulink files which I will be using to program my Artix-7 FPGA, and the Matlab scripts for initializing the variables that the Simulink blocks use.

**Tools and Software used in this project:**
- Simulink/Matlab 2024a
- Vitis Model Composer 2025.1
- Vivado 2025.1
- ~~Artix-7 AC701 Development Board~~ 
(Artix board no longer required, pure simulation. Program still designed around board for future implementation. Keep in mind for the future, the license for compiling VMC into code for the board costs $2,000USD).

Inside Main are the Matlab scripts, and the Outdated folder is simply older scripts and Simulink files in case I need to reference them later.
InitializationOfVariablesRev00x.m is a script that needs to be ran before opening the Simulink file, as it initializes all the variables contained in the blocks.
RadarSimulinkProjectRev00x.slx is the main Simulink file containing the block model
The .slxc files can be ignored.
The last file, PostProcessingAndVisualizationRev00x.m is run at the completion of the simulation, and will handle the graphical representation of the radar sweep. It is not streaming, and would need to be adjusted if it were to be implemented to take serial data from the Artix-7 board in the future.

Current Progress:
- Generating simulated signal (X)
- Range Processing (X)
- Velocity Calculation () -> (Debugging)
- Angle Estimation ( )
- Display in Matlab (X)
- Testing ( )

What is left to do:
- Ensure velocity results match expected values
- Confirm that all data is being sent out to the workspace correctly. Identify and fix any other issues.
- Test

