# DISCLAIMER!!
These codes were very recently uploaded (March 20th, 2024) from my system and have yet to be vetted for general distribution/use. I am currently working through commenting, cleaning and re-testing these solver which were put together approximately 2 years ago. Please allow me the next few weeks to clean up the codes so as to make them more instructive, easy to use, and bug-free. 

`While this disclaimer is being displayed` **I suggest the user refrain from using the CFEI solvers** and exercise *some* caution when using the deconstructed **redbkit** solvers provided.

# Monolithic_FSI
This repository contains three 2D and 3D high-fidelity computational fluid-structure interaction solvers built (mainly) in MATLAB. They are primarily focused on solving two-fields problems monolithically, where the incompressible Navier-Stokes equations - written in the Arbitrary Lagrangian Eulerian (ALE) form - are coupled with the nonlinear elastodynamics structural equations. The numerical strategies employed in the above tools fall into the following three catagoriess: 

1. Fully Implicit Nonlinear FSI 
2. Semi-Implicit Geometry–Convective Explicit (GCE) 
3. Semi-Implicit Combined Field with Explicit Interface (CFEI)

they are each of the Galerkin least-squares finite element (FE) variety and are based on the Arbitrary Lagrangian–Eulerian (ALE) methodology for solving FSI problems. 

## Download, Installation and Running
-------

For detailed configuartion, installion and operation instructions, in particular for the redbKIT based solvers, please follow the steps contained in the brief [INSTALL.md](INSTALL.md) file provided.

## Solver Details 
<ins>***CFEI Solver***</ins><br>
The "Combined Field and Explicit Interface" scheme was developed at the [National University of Singapore](https://cde.nus.edu.sg/me/) by the graduate students of my superviosor [Dr. Rajeev Jaiman](https://scholar.google.com/citations?user=iofAU68AAAAJ&hl=en&oi=ao). This technique uses a finite-element spatial discretization for both the fluid and structural equations. It uses stabilized streamline-upwind Petrov-Galerkin (SUPG) P<sub>m</sub>/P<sub>m-1</sub>/P<sub>m</sub> iso-parametric finite elements to satisfy the inf–sup condition. Temporal discretization is achieved via the second-order bacward difference formula in each domain.  

Further reading on this approach can be found in 
>[**[RKJ16] J. Liu, R. Jaiman, P. Gurugubelli. A stable second-order scheme for fluid–structure interaction with strong added-mass effects **, Journal of Computational Physics, 2014.](https://doi.org/10.1016/j.jcp.2014.04.020)

<ins>***redbKIT Solver***</ins><br>
The redbKIT-based solvers presented above are bare-bones deconstructions (/simplification) of the original [toolkit](https://github.com/redbKIT/redbKIT). At face value they are exclusively the intellectual property of the original repository with a few original lines of functional code included for initialization, data output, and operational purposes. 

These solvers use a finite-element spatial discretization for both the fluid and structural equations. In the above  solvers specifically (as opposed to the original toolkit), we exclusively use higher-order P<sub>2</sub> finite elements for the structural equations and stable Taylor-Hood P<sub>2</sub>/P<sub>1</sub> mixed finite elements for the fluid equations. A second-order backward difference scheme is adopted for the temporal discretization of the fluid domian whilst the generalized-&alpha; method is used for the structural domian. 

As described in detail below; while the [**redbKIT**](https://github.com/redbKIT/redbKIT) toolkit should be the primary reference point for the interested reader, it is my belief that the condensed monolithic solvers (as presented here in their simplified intermediary forms) hold their own unique value to the wider student/researcher community. In particular, for the purposes of streamlined learning, efficient tool adoption, and additional/alternative technique development. 

## Further Clarifying Details
For absolute clarity, the [**redbKIT**](https://github.com/redbKIT/redbKIT) toolkit by itself is an ***AMAZING*** resource for (young) academics seeking an education on the finite element method and reduced-order modelling. As such I **strongly** recommend the reader fork that repository and work through redbKIT's numerous well crafted example problems. *I personally give this resource the lion's share of credit in helping me complete my graduate studies.*

Having said that, however, as is true for (almost) all **well developed / robust** software packages used to tackle complex problems in computational physics; its sophistication tends to be its own undoing for those novice students beginning their journey in the simulation sciences or those young numerical methods researchers who are focusing on a niche element of the problem/framework/technique and just want the machinery around the issue to be (first) working and (secondly) easily digestible for peace of mind, debugging, development etc. purposes. As such, the codes provided above seek to shallow the learning curve and expedite the toolkit's adoption/integration. 

In the case of my graduate studies, for example, having simplified the redbKIT solvers into just those most essential numerical elements needed for running a reasonably complex fluid-structure interaction simulation, I was **then** able to more easily/quickly deconstruct the monolthic system of equations into a partitioned framework for solving FSI problem. From this point I built out all the supporting tools needed to perform [partitioned simulations](https://github.com/JTGonzo/Partitioned_FSI) (interface mapping, iterative acceleration etc.) and only then was I was able to **BEGIN** my research into ["*stabilization strategies for low mass-ratio partitioned FSI simulations*"](https://jtgonzo.github.io/).  


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
