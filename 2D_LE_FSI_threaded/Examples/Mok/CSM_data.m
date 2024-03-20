%DATAFILE for FS1 - Solid

% Source term
data.force{1} = @(x, y, t, param)(0.*x.*y);
data.force{2} = @(x, y, t, param)(0.*x.*y);

% Dirichlet
data.bcDir{1} = @(x, y, t, param)(0.*x.*y); 
data.bcDir{2} = @(x, y, t, param)(0.*x.*y); 

% Neumann
data.bcNeu{1} = @(x, y, t, param)(0.*x.*y);
data.bcNeu{2} = @(x, y, t, param)(0.*x.*y);

% Normal Pressure
data.bcPrex   = @(x, y, t, param)(0.*x.*y);

% BC flag
data.flag_dirichlet{1}    =  [7];
data.flag_neumann{1}      =  [];
data.flag_FSinterface{1}  =  [6];
data.flag_pressure{1}     =  [];
data.flag_robin{1}     = [];
data.flag_clamp_points{1} = [];

data.flag_dirichlet{2}    =  [7];
data.flag_neumann{2}      =  [];
data.flag_FSinterface{2}  =  [6];
data.flag_pressure{2}     =  [];
data.flag_robin{2}     = [];
data.flag_clamp_points{2} = [];

data.u0{1} = @(x, y, t, param)(0.*x.*y);
data.u0{2} = @(x, y, t, param)(0.*x.*y);
data.du0{1} = @(x, y, t, param)(0.*x.*y);
data.du0{2} = @(x, y, t, param)(0.*x.*y);

% material parameters 
data.Young   = 2.3*10^6;
data.Poisson = 0.45;
data.Density = 1500;
data.Material_Model   =  'Linear'; 
data.model   = 'CSM';

data.ElasticCoefRobin = 1e+6;
data.flag_dirichletNormal = [];
data.Output.ComputeVonMisesStress = false;

% NonLinear Solver
data.NonLinearSolver.tol               = 1e-6;
data.NonLinearSolver.maxit             = 35;

% Linear solver 
data.options.LinSolver.solver            = 'backslash';
data.LinearSolver.type                   =  'backslash';
data.LinearSolver.mumps_reordering  = 7;

% Preconditioner
data.Preconditioner.type         = 'None'; % AdditiveSchwarz, None, ILU

% Time options
data.time.t0         = 0;
data.time.dt         = 0.01; 
data.time.tf         = 15;
data.time.gamma      = 1/2;
data.time.beta       = 1/4;
data.time.alpha_m    = 0;
data.time.alpha_f    = 0;
