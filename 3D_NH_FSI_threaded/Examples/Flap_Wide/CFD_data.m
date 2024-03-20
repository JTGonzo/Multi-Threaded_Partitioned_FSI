%DATAFILE Fluid
U_bar = 10;

% Source term
data.force{1} = @(x, y, z, t, param)(0.*x.*y);
data.force{2} = @(x, y, z, t, param)(0.*x.*y);
data.force{3} = @(x, y, z, t, param)(0.*x.*y);

% Dirichlet
data.bcDir_t  = @(t)( 0.5*(1-cos(pi/0.2*t)).*(t<0.0) + 1 * (t>=0.0));

data.bcDir{1} = @(x, y, z, t, param)(data.bcDir_t(t) * U_bar.*(x==0) + 0.*x.*y); 
data.bcDir{2} = @(x, y, z, t, param)(0.*x.*y); 
data.bcDir{3} = @(x, y, z, t, param)(0.*x.*y); 

% Neumann
data.bcNeu{1} = @(x, y, z, t, param)(0.*x.*y);
data.bcNeu{2} = @(x, y, z, t, param)(0.*x.*y);
data.bcNeu{3} = @(x, y, z, t, param)(0.*x.*y);

% initial condition
data.u0{1}    = @(x, y, z, t, param)(0.*x.*y);
data.u0{2}    = @(x, y, z, t, param)(0.*x.*y);
data.u0{3}    = @(x, y, z, t, param)(0.*x.*y);

% BC flag
data.flag_dirichlet{1}    =  [1];
data.flag_neumann{1}      =  [2 3 4 5];
data.flag_FSinterface{1}  =  [6];
data.flag_ring{1}         =  [1]; 
data.flag_ALE_fixed{1}    =  [1 2];

data.flag_dirichlet{2}    =  [1 4 5];
data.flag_neumann{2}      =  [2 3];
data.flag_FSinterface{2}  =  [6];
data.flag_ring{2}         =  [1]; 
data.flag_ALE_fixed{2}    =  [1 4 5];

data.flag_dirichlet{3}    =  [1 3];
data.flag_neumann{3}      =  [2 4 5];
data.flag_FSinterface{3}  =  [6];
data.flag_ring{3}         =  [1]; 
data.flag_ALE_fixed{3}    =  [1 3];

% Model parameters
data.dynamic_viscosity    =   1;
data.density              =   1;
data.Stabilization        =   'SUPG';

% Nonlinear solver
data.NonLinearSolver.tol         = 1e-6; 
data.NonLinearSolver.maxit       = 25;
 
% Linear Solver
data.LinearSolver.type              = 'backslash'; % MUMPS, backslash, gmres
data.LinearSolver.mumps_reordering  = 7;

% Preconditioner
data.Preconditioner.type         = 'None'; % AdditiveSchwarz, None, ILU

%% Time Setting
data.time.BDF_order  = 2;
data.time.t0         = 0;
data.time.dt         = 0.005; 
data.time.tf         = 2;
data.time.nonlinearity  = 'semi-implicit';

%% Output options
data.Output.DragLift.computeDragLift = 1;
data.Output.DragLift.factor          = 2/(1*10^2*0.1);
data.Output.DragLift.flag            = [5];
data.Output.DragLift.filename        = 'Results/AerodynamicForcesFSI_Final.txt';