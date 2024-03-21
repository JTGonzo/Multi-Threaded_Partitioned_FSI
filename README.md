![alt text](mok_white.gif)
# DISCLAIMER!!
These codes were very recently uploaded (March 20th, 2024) from my system and have yet to be vetted for general distribution/use. I am currently working through commenting, cleaning and re-testing these solver which were put together approximately 1 year ago. Please allow me the next few weeks to clean up the codes so as to make them more instructive, easy to use, and bug-free. 

If the reader is still compelled to experiment I encourgae them to do so, however, `while this disclaimer is being displayed` I suggest that one exercises ***reasonable*** caution when using the partitoned solvers provided.

## Download, Installation and Running
-------

For detailed configuartion, installation and operation instructions, in particular for the redbKIT based solvers, please follow the steps contained in the brief [INSTALL.md](INSTALL.md) file provided.

## Strongly Coupled Partitioned FSI Solver
<ins>***2D & 3D redbKIT-based solvers***</ins><br>
This repository contains **MULTI-THREADED** 2D and 3D high-fidelity computational fluid-structure interaction solvers built (mainly) in MATLAB. They are primarily focused on solving two-field problems in a partitioned fashion, where the incompressible Navier-Stokes equations - written in the Arbitrary Lagrangian Eulerian (ALE) form - are coupled with the nonlinear elastodynamics structural equations. For 2D FSI problems the user has the option to choose either a linear elastic constitutive material model for the structural mechanics solver (**2D_Part_FSI_LE**) or a nonlinear St. Venant Kirchhoff material model (**2D_Part_FSI_SVK**). For 3D FSI problems, the user can additionally choose a nearly-incompressible NeoHookean constitutive material model for the structural mechanics solver (**3D_Part_FSI_NH**). 

These solvers all use a Galerkin least-squares finite element strategy for each set of domain specific equations. Here we exclusively use higher-order P<sub>2</sub> finite elements for the structural equations and stable Taylor-Hood P<sub>2</sub>/P<sub>1</sub> mixed finite elements for the fluid equations. A second-order backward difference scheme is adopted for the temporal discretization of the fluid domian whilst the generalized-&alpha; method is used for the structural domian. 

In all cases each code was derived from a source [redbKIT monolithic FSI](https://github.com/JTGonzo/Monolithic_FSI) solver. These solvers were decomposed into their baseline elements and were then strategically reassembled and further built upon such that we finally ended up with a **partitioned** numerical framework. These partitoned solvers (given above) served as a robust testbed for my subsequent investigation of ["*coupling  instability in low mass-ratio partitioned FSI simulations*"](https://jtgonzo.github.io/). From this effort, algorithmic solutions to stabilize and accelerate the fixed point coupling procedure comensurate with the partitioned approach [were recommended](https://github.com/JTGonzo/Multi-Threaded_Partitioned_FSI).  

Features unique to the partitioned strategy that have been included in each of the solvers provided include:<br>
* Coupling Methods:
  - Anderson Acceleration
  - Generalized Broyden
  - Broyden's Second Method
  - Multi-Vector Least Square 
  - Dynamic Aitken Relaxation 
  - Constant Relaxation 
  - Pure Gauss-Seidel Iterations

* Filtering Strategies:
  - QR1 filtering
  - QR2 (Classic Gram-Schmidt) filtering
  - QR2 (Modified Gram-Schmidt) filtering
  - QR2 (Householder Reflections) filtering
  - POD filtering

* Interface Mapping Techniques:
  - Nearest Neighbour Interpolation
  - Projection and Linear Interpolation
  - Radial Basis Function Interpolation

* Extrapolation Methods:
  - Constant
  - Linear 
  - Quadratic
  - Cubic

* [Work in Progress Tools](https://github.com/JTGonzo/Multi-Threaded_Partitioned_FSI):
  - Interface tracking 
  - Adaptive Regularization Least Squares
  - Selective Secant Inclusion
  - Truncated SVD filtering

## Further Clarifying Details
For absolute clarity, the [**redbKIT**](https://github.com/redbKIT/redbKIT) toolkit by itself is an ***AMAZING*** resource for (young) academics seeking an education in the finite element method and reduced-order modelling. As such, I **strongly** recommend the reader fork that repository and work through redbKIT's numerous well crafted example problems. `*I personally give this resource the lion's share of credit in helping me complete my graduate studies.*`

These solvers, however, are of a monolthic construction and thus require the linearization of one set of domain equations with respect to the other for the system of equations to be solvable. If the interested reader is more interested in the details of **partitioned** coupling, however, or would like to develop one of the respective domain solvers independently of the other domain they may find the tools found in this repository very useful.

Lastly, for those readers seeking a more digestible version of the rebKIT-based monolthic FSI solvers or the individual (rebKIT-based) domain-specific solvers I have additionally compiled what I believe to be useful resources in the following repositories: 
  1. [Monolithic Fluid Structure Interaction solvers](https://github.com/JTGonzo/Monolithic_FSI)
  2. [Single Domain (CFD & CSM) solvers](https://github.com/JTGonzo/Single_Physics_Solvers)

## redbKIT : a MATLAB(R) library for reduced-order modeling of parametrized PDEs

redbKIT is a MATLAB library for finite element simulation and reduced-order modeling of Partial Differential Equations developed at [EPFL](https://www.epfl.ch/) - [Chair of Modeling and Scientific Computing](http://cmcs.epfl.ch/)). 

In particular, it includes straightforward implementations of many of the algorithms presented in the companion book:

>[**[QMN16] A. Quarteroni, A. Manzoni, F. Negri. Reduced Basis Methods for Partial Differential Equations. An Introduction**, Springer, 2016.](http://www.springer.com/us/book/9783319154305#aboutBook)

For further details please visit [redbKIT website](http://redbkit.github.io/redbKIT/).

### License
-------

**redbKIT is distributed under BSD 2-clause license**

Copyright (c) 2015-2017, Ecole Polytechnique Fédérale de Lausanne (EPFL)
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


### **redbkit** Development
-------

redbKIT was developed and maintained by [`Federico Negri`](https://www.linkedin.com/in/negrifederico/). Paola Gervasio (Università degli Studi di Brescia) is gratefully acknowledged for granting the use of parts of the finite element code MLife.


## Contact
-------
If you have any questions regarding the contents of this repository or the use of these solvers, please feel to email me at <josetg@vt.edu>.
