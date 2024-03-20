===========================================
Instructions - rebkit-based solvers 
===========================================

Pre-requisites
--------------

For the basic installation, you need:<br>
-> a reasonably recent version of Matlab (anything above R2017a should be fine)<br>
-> a supported C compiler, see e.g. http://it.mathworks.com/support/compilers/R2016a/<br>
   (gcc on Unix and gcc/clang on Mac OS X work smoothly)


Basic Installation
------------------

Open Matlab and navigate to the root folder of the downloaded solver that you will be using (e.g. **3D_Part_FSI_NH**). Then type 

>> make

to compile the C-assembly routines and "mexify" some other files.


Running and Post-Processing the Example Provided 
------------------

With Matlab opened; in the same folder as the *make.m* file ran previously, type 

>> main_FSI

doing this begins the example simulation and creates a "Figures" folder where the simulation configuration is output periodically as a series of *.vtk* files. Additionaly, a "Results" folder is created where the code prints physical and coupling results at the end of each time step to a handful of output files. Namely; the resultant aerodynamic forces integrated over the defined FSI surface, the displacement and velocity results at the solid (Wid) and fluid (WidF) watchpoint nodes (user defined in *init_couple.m*), the required number of coupling iterations needed per time step, and the resdiual errors of the interfacial variables after each fixed-point iteration. 

Some basic post processing tools are also provided in the `Post Processing` folder which plots the evolution of the collected variables (force, displacement, velocity, residual errors, fixed-point iterations) over the course of the performed simulation. 

Additional Example Simulations Provided 
------------------
For each of the solvers there are additional example simulations provided that were used to validate the code against well established FSI benchmarks found in the literature. To run these files the user just needs to replace the ***CFD_data.m, CSM_data.m, and .mat*** files currently in the root folder with the corresponding files found in the respective `Examples/*` folder.

Next, the user simply updates the `init_data.m` file as described below to instruct the solver which .mat files to import and what names to give to the output files. These new names given to the output files also naturally requires eventually editting the post processing files when the time comes to import the new results. 

Running your own Simulation
------------------
<ins>***Mesh and .MAT file generation***</ins><br>
See this brief [tutorial](MESHING.md) on how to set up the **.mat** files needed to perform your own simulations. These .mat files define the geometric/spatial properties of the problem you intend to investigate. Pay special attention to how you define and identify the boundary and interface surfaces for each of the respective domains. 

<ins>***Editting init_data.m***</ins><br>
The only adjustments that are needed in this file to run your own simulations are 

*line 6 :* problemString = 'Tube';

*line 12 :* vtk_filename = 'Figures/Tube_';

*line 21 :* load('Tube_S.mat')

*line 28 :* load('Tube_F.mat')

These lines are where you indicate the names you chose for the domain specific .mat files that you will be importing into the solver as well as the prefix string that you assign to your output data and files.  

<ins>***Editting init_coupling.m***</ins><br>
This file is where the user defines all the coupling parameters they'd like to use for their specific simulation. Some of the setting that can be manipulated include;

   1. Watchpoint node IDs
   2. Variable output frequency 
   3. Type of coupling strategy
   4. Type of filtering strategy
   5. Type of extrapolator
   6. Maximum number of coupling iterations per time-step
   7. Relative and Absolute convergence threshold 
   8. Relaxation weighting
   9. Number of time-steps to reuse (Anderson Acceleration)
   10. Linear dependence filtering threshold (mainly for Anderson Acceleration)

given that the optimal setting of these parameters (from a performance, complexity, and cost standpoint) can be highly problem dependent when tackling particularly challenging FSI problems I plan to post a short video soon instructing the user on how best to do so. Until then it is advised that the user keep the settings as is (apart from defining desired watchpoint IDs of course) especially so if the user is uninterested in the details of coupling and would just like to study the physics of an FSI problem itself. 

<ins>***Editting NS_data.m and CSM_data.m***</ins><br>
These two files are the primary locations where you define all the domain specific parameters needed for your unqiue problem. Parameters that you are able to set include: 

1. material properties 
2. boundary conditions
3. initial conditions
4. structure's constitutive model to use
5. temporal and spatially varying properties
6. time step/integration properties
7. aerodynamic forces output variables

===========================================
Instructions - 2D rigid-body FSI solver
===========================================

Pre-requisites
--------------

For the basic installation, you just need a reasonably recent version of Matlab (anything above R2015a should be fine). There are no additional compiler or system requirements. 


Running the Example Provided 
------------------

Open Matlab and navigate to the root folder of the downloaded solver then type 

>> main_FSI

doing this begins the example simulation provided. 

The simulation configuration is output periodically as a series of *.plt* files, and the displacement and interface aerodynamic force results experienced by the rigid body are printed to a local *.oisd* file in the folder. 