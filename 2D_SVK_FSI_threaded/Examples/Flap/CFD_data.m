%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% DATAFILE for Fluid Domain %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U_bar = 10; % mean velocity 

%% Applied boundary conditions 
% Source term
data.force{1} = @(x, y, t, param)(0.*x.*y);
data.force{2} = @(x, y, t, param)(0.*x.*y);

% temporal Dirichlet function
data.bcDir_t  = @(t)( 0.5*(1-cos(pi/0.2*t)).*(t<0.0) + 1 * (t>=0.0));
%data.bcDir_t  = @(t)( 0.5*(1-cos(pi/1*t)).*(t<1) + 1 * (t>=1));

% 2D spatial Dirchlet boundary condition
%data.bcDir{1} = @(x, y, t, param)(data.bcDir_t(t) * 1.5*U_bar*y.*(H-y)/(H/2)^2.*(x==0) + 0.*x.*y); 
data.bcDir{1} = @(x, y, t, param)(data.bcDir_t(t) * U_bar.*(x==0) + 0.*x.*y); 
data.bcDir{2} = @(x, y, t, param)(0.*x.*y); 

% 2D Neumann boundary condition
data.bcNeu{1} = @(x, y, t, param)(0.*x.*y);
data.bcNeu{2} = @(x, y, t, param)(0.*x.*y);

% 2D initial velocity condition
data.u0{1}    = @(x, y, t, param)(0.*x.*y);
data.u0{2}    = @(x, y, t, param)(0.*x.*y);

%% BC flags for physical gemomtry boundaries
% x-direction conditions
data.flag_dirichlet{1}   =  [2 4]; % numbers indicate physical surfaces of mesh
data.flag_neumann{1}     =  [3];
data.flag_FSinterface{1} =  [5];
data.flag_ring{1}        =  [1]; 
data.flag_ALE_fixed{1}   =  [2 4 3];
data.flag_pressure{1}    =  [];
%data.flag_ring{1}       =  [];

% y-direction conditions
data.flag_dirichlet{2}   =  [2 4];
data.flag_neumann{2}     =  [3];
data.flag_FSinterface{2} =  [5];
data.flag_ring{2}        =  [1]; 
data.flag_ALE_fixed{2}   =  [2 4 3];
data.flag_pressure{2}    =  [];
%data.flag_ring{2}       =  [];

data.flag_resistance  = [];
data.flag_absorbing   = [];

%% Domain material parameters/properties
data.dynamic_viscosity    =   1;
data.density              =   1;

%% Numerical scheme and solver properties/settings
% Fluid Pressure Stabilization
% data.Stabilization        =   'SUPG';

% Nonlinear solver
data.NonLinearSolver.tol         = 1e-6; 
data.NonLinearSolver.maxit       = 30;
 
% Linear Solver
data.LinearSolver.type              = 'backslash'; % MUMPS, backslash, gmres
data.LinearSolver.mumps_reordering  = 7;

% Preconditioner
data.Preconditioner.type         = 'None'; % AdditiveSchwarz, None, ILU

%% Simulation time settings
data.time.BDF_order  = 2;
data.time.t0         = 0;
data.time.dt         = 0.001; %0.00165; 
data.time.tf         = 5;
data.time.nonlinearity  = 'semi-implicit';

%% Lift/Drag output options
data.Output.DragLift.computeDragLift = 1;
data.Output.DragLift.factor          = 2/(1*10^2*0.1);  % (2 / (r A v^2))
data.Output.DragLift.flag            = [5];             % physical surface to compute lift/drag 
data.Output.DragLift.filename        = 'Results/AeroForces.txt';
