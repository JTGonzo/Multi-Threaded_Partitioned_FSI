===========================================
redbKIT solver INSTALLATION
===========================================

Pre-requisites
--------------

For the basic installation, you need:
-> a reasonably recent version of Matlab (anything above R2017a should ve fine)
-> a supported C compiler, see e.g. http://it.mathworks.com/support/compilers/R2016a/
   (gcc on Unix and gcc/clang on Mac OS X work smoothly)


Basic Installation
------------------

Open Matlab and navigate to the root folder of the solver that you will be using (e.g. **2D_FSI_implicit_nonlinear**). Then type 

>> make

to compile the C-assembly routines and "mexify" some other files.


Running the Example Provided 
------------------

With Matlab opened; in the same folder as the *make.m* file, type 

>> main_FSI


Running your own Simulation
------------------