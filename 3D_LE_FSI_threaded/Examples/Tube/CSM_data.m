% Source term
data.force{1} = @(x, y, z, t, param)(0.*x.*y.*z);
data.force{2} = @(x, y, z, t, param)(0.*x.*y.*z);
data.force{3} = @(x, y, z, t, param)(0.*x.*y.*z);

% Dirichlet
data.bcDir{1} = @(x, y, z, t, param)(0.*x.*y.*z);
data.bcDir{2} = @(x, y, z, t, param)(0.*x.*y.*z);
data.bcDir{3} = @(x, y, z, t, param)(0.*x.*y.*z);

% Neumann
data.bcNeu{1} = @(x, y, z, t, param)(0.*x.*y.*z);
data.bcNeu{2} = @(x, y, z, t, param)(0.*x.*y.*z);
data.bcNeu{3} = @(x, y, z, t, param)(0.*x.*y.*z);

% Normal Pressure
data.bcPrex   = @(x, y, z, t, param)(0.*x.*y.*z);

% BC flag
data.flag_dirichlet{1}    =  [4 5];
data.flag_neumann{1}      =  [];
data.flag_FSinterface{1}  =  [3];
data.flag_pressure{1}     =  [];
data.flag_robin{1}     = [];
data.flag_clamp_points{1} = [];

data.flag_dirichlet{2}    =  [4 5];
data.flag_neumann{2}      =  [];
data.flag_FSinterface{2}  =  [3];
data.flag_pressure{2}     =  [];
data.flag_robin{2}     = [];
data.flag_clamp_points{2} = [];

data.flag_dirichlet{3}    =  [4 5];
data.flag_neumann{3}      =  [];
data.flag_FSinterface{3}  =  [3];
data.flag_pressure{3}     =  [];
data.flag_robin{3}     = [];
data.flag_clamp_points{3} = [];

data.u0{1} = @(x, y, z, t, param)(0.*x.*y.*z);
data.u0{2} = @(x, y, z, t, param)(0.*x.*y.*z);
data.u0{3} = @(x, y, z, t, param)(0.*x.*y.*z);

data.du0{1} = @(x, y, z, t, param)(0.*x.*y.*z);
data.du0{2} = @(x, y, z, t, param)(0.*x.*y.*z);
data.du0{3} = @(x, y, z, t, param)(0.*x.*y.*z);

% material parameters 
data.Young   = 3*10^5;
data.Poisson = 0.3;
data.Density = 1200;
data.Material_Model   = 'Linear';
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
data.time.dt         = 0.0001; 
data.time.tf         = 0.02;
data.time.gamma      = 1/2;
data.time.beta       = 1/4;
data.time.alpha_m    = 0;
data.time.alpha_f    = 0;