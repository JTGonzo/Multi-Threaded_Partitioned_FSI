%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  Initialize Fluid solver variable spaces, export files %%%%%%%
%%%%%% integration parameters, reference mesh and soln' tractions %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize time advance parameters
BDF_orderF = DATA.Fluid.time.BDF_order;
t0        = DATA.Fluid.time.t0;
dt        = DATA.Fluid.time.dt;
tf        = DATA.Fluid.time.tf;
t         = DATA.Fluid.time.t0;
k_t       = 0;

%% Initialize variable space
% read initial velocity condition
v0  = [];
for k = 1 : FE_SPACE_v.numComponents
    v0  = [v0; DATA.Fluid.u0{k}(  MESH.Fluid.nodes(1,:), MESH.Fluid.nodes(2,:), t0, param )'];
end

% pressure variables
p0 = zeros(FE_SPACE_p.numDof,1);

% full fluid domain solution space
uF = [v0; p0];

%% Initialize export data
% export initial condition in vtk format
if ~isempty(vtk_filename)
    CFD_export_solution(MESH.dim, uF(1:FE_SPACE_v.numDof), uF(1+FE_SPACE_v.numDof:end), ...
        MESH.Fluid.vertices, MESH.Fluid.elements, MESH.Fluid.numNodes, [vtk_filename,'Fluid'], 0);
end

%% Initialize time integration scheme
TimeAdvanceF = BDF_TimeAdvance( BDF_orderF );  % set BDF order

% provide intial domain velocity solution
TimeAdvanceF.Initialize( v0 );
for bd = 2 : BDF_orderF
    TimeAdvanceF.Append( v0 );
end

%% Initialize aerodynamic forces output
Fx(k_t+1)  = 0;
Fy(k_t+1)  = 0;
Fz(k_t+1)  = 0;
dofs_drag    = [];
p_dofs    = [];

% Define the interface/boundary nodes to export lift/drag data
    for j = 1 : length(DATA.Fluid.Output.DragLift.flag)
        Dirichlet_side         = find(MESH.Fluid.boundaries(MESH.Fluid.bc_flag_row,:) == DATA.Fluid.Output.DragLift.flag(j));
        Dirichlet_side         = unique(Dirichlet_side);
        Dirichlet_dof          = MESH.Fluid.boundaries(1:MESH.Fluid.numBoundaryDof,Dirichlet_side);
        dofs_drag              = [dofs_drag; Dirichlet_dof(:)];
    end
dofs_drag = unique(dofs_drag);
dofs_drag = setdiff(dofs_drag, intersect(dofs_drag,MESH.Fluid.Dirichlet_dof_c{1}));

for i = 1:length(dofs_drag)
    if dofs_drag(i) < MESH.Fluid.numVertices
        p_dofs =  [p_dofs; dofs_drag(i)];       
    end
end

% lift and drag output file
filename5 = sprintf('Results/%s_forces.txt',problemString);   % iterations residual info data file
fileDragLift = fopen(filename5, 'w+');
fprintf(fileDragLift, '%1.4e  %1.4e  %1.4e  %1.4e', t, Fx(k_t+1), Fy(k_t+1), Fz(k_t+1));

%% initialize algebraic Navier-Stokes LHS & RHS matricices/vector spaces
C_NS = 0*sparse(FE_SPACE_v.numDof + FE_SPACE_p.numDof, FE_SPACE_v.numDof + FE_SPACE_p.numDof);
F_NS = 0*sparse(FE_SPACE_v.numDof + FE_SPACE_p.numDof, 1);
trac = zeros(size(F_NS));

%% Initialize the ALE mes solution space
ALE_velocity = zeros(MESH.Fluid.numNodes*dim, 1);
d_Fn         = zeros(MESH.Fluid.numNodes*dim, 1);
