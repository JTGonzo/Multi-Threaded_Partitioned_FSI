%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri@epfl.ch>

%DATAFILE Fluid
P_bar = 1333.2;

% Source term
data.force{1} = @(x, y, z, t, param)(0.*x.*y.*z);
data.force{2} = @(x, y, z, t, param)(0.*x.*y.*z);
data.force{3} = @(x, y, z, t, param)(0.*x.*y.*z);

% Dirichlet
data.bcDir_t  = @(t)( 1.*(t<=0.003) + 0.*(t>0.003));

data.bcDir{1} = @(x, y, z, t, param)(0.*x.*y.*z); 
data.bcDir{2} = @(x, y, z, t, param)(0.*x.*y.*z); 
data.bcDir{3} = @(x, y, z, t, param)(0.*x.*y.*z); 

% Neumann
data.bcNeu{1} = @(x, y, z, t, param)(0.*x.*y.*z);
data.bcNeu{2} = @(x, y, z, t, param)(0.*x.*y.*z);
data.bcNeu{3} = @(x, y, z, t, param)(0.*x.*y.*z);

% Normal Pressure
data.bcPrex   = @(x, y, z, t, param)(data.bcDir_t(t)*P_bar.*(x==0) + 0.*x.*y.*z);
% data.bcPrex{2}   = @(x, y, z, t, param)(0.*x.*y.*z);
% data.bcPrex{3}   = @(x, y, z, t, param)(0.*x.*y.*z);

% initial condition
data.u0{1}    = @(x, y, z, t, param)(0.*x.*y.*z);
data.u0{2}    = @(x, y, z, t, param)(0.*x.*y.*z);
data.u0{3}    = @(x, y, z, t, param)(0.*x.*y.*z);

% BC flag
data.flag_dirichlet{1}    =  []; %1 5
data.flag_neumann{1}      =  [];
data.flag_FSinterface{1}  =  [3]; %6
data.flag_ring{1}         =  [1 2]; 
data.flag_ALE_fixed{1}    =  [1 2]; %1 5 2
data.flag_pressure{1}     =  [1]; 

data.flag_dirichlet{2}    =  [];
data.flag_neumann{2}      =  [];
data.flag_FSinterface{2}  =  [3]; %6
data.flag_ring{2}         =  [1 2]; 
data.flag_ALE_fixed{2}    =  [1 2];%1 5 4
data.flag_pressure{2}     =  []; 

data.flag_dirichlet{3}    =  [];
data.flag_neumann{3}      =  [];
data.flag_FSinterface{3}  =  [3]; %6
data.flag_ring{3}         =  [1 2]; 
data.flag_ALE_fixed{3}    =  [1 2]; %1 5 3
data.flag_pressure{3}     =  []; 

% Model parameters
data.dynamic_viscosity  =   0.003;
data.density  =   1000;
%data.Stabilization        =   'SUPG';
data.flag_resistance = [];
data.flag_absorbing = [];

% Nonlinear solver
data.NonLinearSolver.tol         = 1e-6; 
data.NonLinearSolver.maxit       = 30;
 
% Linear solver
data.LinearSolver.type              = 'backslash'; % MUMPS, backslash, gmres
data.LinearSolver.mumps_reordering  = 7;

% Preconditioner
data.Preconditioner.type         = 'None'; % AdditiveSchwarz, None, ILU

%% Time Setting
data.time.BDF_order  = 2;
data.time.t0         = 0;
data.time.dt         = 0.0001; 
data.time.tf         = 0.02;
data.time.nonlinearity  = 'semi-implicit';

%% Output options
data.Output.DragLift.computeDragLift = 1;
data.Output.DragLift.factor          = 2/(1.18*1e-3*51.3^2*1);
data.Output.DragLift.flag            = [5 6];
data.Output.DragLift.filename        = 'Results/AerodynamicForcesFSI_Final.txt';