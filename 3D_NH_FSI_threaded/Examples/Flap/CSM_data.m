%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% DATAFILE for Solid Domain %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Applied boundary conditions 
% Source term
data.force{1} = @(x, y, z, t, param)(0.*x.*y);
data.force{2} = @(x, y, z, t, param)(0.*x.*y);
data.force{3} = @(x, y, z, t, param)(0.*x.*y);

% 3D spatial Dirchlet boundary condition
data.bcDir{1} = @(x, y, z, t, param)(0.*x.*y);
data.bcDir{2} = @(x, y, z, t, param)(0.*x.*y);
data.bcDir{3} = @(x, y, z, t, param)(0.*x.*y);

% 3D Neumann boundary condition
data.bcNeu{1} = @(x, y, z, t, param)(0.*x.*y);
data.bcNeu{2} = @(x, y, z, t, param)(0.*x.*y);
data.bcNeu{3} = @(x, y, z, t, param)(0.*x.*y);

% Normal Pressure
data.bcPrex   = @(x, y, z, t, param)(0.*x.*y);

% 3D initial displacement & velocity condition
data.u0{1} = @(x, y, z, t, param)(0.*x.*y);
data.u0{2} = @(x, y, z, t, param)(0.*x.*y);
data.u0{3} = @(x, y, z, t, param)(0.*x.*y);

data.du0{1} = @(x, y, z, t, param)(0.*x.*y);
data.du0{2} = @(x, y, z, t, param)(0.*x.*y);
data.du0{3} = @(x, y, z, t, param)(0.*x.*y);

%% BC flags for physical gemomtry boundaries
% x-direction conditions
data.flag_dirichlet{1}    =  [6];
data.flag_neumann{1}      =  [7];
data.flag_FSinterface{1}  =  [5];
data.flag_pressure{1}     =  [];
data.flag_robin{1}     = [];
data.flag_clamp_points{1} = [];

% y-direction conditions
data.flag_dirichlet{2}    =  [6];
data.flag_neumann{2}      =  [7];
data.flag_FSinterface{2}  =  [5];
data.flag_pressure{2}     =  [];
data.flag_robin{2}     = [];
data.flag_clamp_points{2} = [];

% z-direction conditions
data.flag_dirichlet{3}    =  [6 7];
data.flag_neumann{3}      =  [];
data.flag_FSinterface{3}  =  [5];
data.flag_pressure{3}     =  [];
data.flag_robin{3}     = [];
data.flag_clamp_points{3} = [];

data.flag_dirichletNormal = [];
data.Output.ComputeVonMisesStress = false;

%% Domain material parameters/properties
% material parameters 
data.Young   = 4*10^6;
data.Poisson = 0.3;
data.Density = 3*10^3; %30; %
data.Material_Model   =  'NeoHookean'; 
data.model   = 'CSM';

data.ElasticCoefRobin = 1e+6;

%% Numerical scheme and solver properties/settings
% Compute structural stress
data.Output.ComputeVonMisesStress = false;

% NonLinear Solver
data.NonLinearSolver.tol               = 1e-6;
data.NonLinearSolver.maxit             = 25;

% Linear solver 
data.options.LinSolver.solver            = 'backslash';
data.LinearSolver.type                   =  'backslash';
data.LinearSolver.mumps_reordering  = 7;

% Preconditioner
data.Preconditioner.type         = 'None'; % AdditiveSchwarz, None, ILU

%% Simulation time settings
data.time.t0         = 0;
data.time.dt         = 0.001; 
data.time.tf         = 5;
data.time.gamma      = 1/2;
data.time.beta       = 1/4;
data.time.alpha_m    = 0;
data.time.alpha_f    = 0;