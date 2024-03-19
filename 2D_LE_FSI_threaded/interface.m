function [MESH] = interface(DATA, MESH, FE_SPACE_v, FE_SPACE_p, FE_SPACE_s, mapper, dim)
%% Generates mapping between solid and fluid interface 
%[MESH] = imap(DATA, MESH);     % assumes co-location of fluid/solid interface nodes
[MESH] = imap2(DATA, MESH, mapper);

%% Define fluid and solid interface arrays for solution space
for k = 1 : dim
    MESH.ndof_interface{k} = length(MESH.Interface_FSmap{k});
end

F_interface_dofs = [];
S_interface_dofs = [];
for k = 1 : MESH.dim
    F_interface_dofs  = [F_interface_dofs; MESH.Fluid.numNodes*(k-1)+MESH.Fluid.dof_interface{k}];
    S_interface_dofs  = [S_interface_dofs; MESH.Solid.numNodes*(k-1)+MESH.Solid.dof_interface{k}];
end

%% Define fluid interface indices in global and internal numbering
tmp = zeros(FE_SPACE_v.numDof+FE_SPACE_p.numDof,1);
for i = 1 : dim
    tmp([FE_SPACE_v.numDofScalar*(i-1)+MESH.Fluid.dof_interface{i}]) = 1; %which of the domain nodes are interfcae nodes
end

MESH.Fluid.Gamma     = find(tmp(MESH.Fluid.internal_dof));
MESH.Fluid.II        = setdiff(1:length(MESH.Fluid.internal_dof),MESH.Fluid.Gamma);

MESH.Fluid.II_global = setdiff(MESH.Fluid.internal_dof,MESH.Fluid.internal_dof(MESH.Fluid.Gamma));
MESH.Fluid.Gamma_global = F_interface_dofs;

%% Define solid interface indices in global and internal numbering
tmp = zeros(FE_SPACE_s.numDof,1);
for i = 1 : dim
    tmp([FE_SPACE_s.numDofScalar*(i-1)+MESH.Solid.dof_interface{i}]) = 1; %which of the domain nodes are interfcae nodes
end

MESH.Solid.Gamma     = find(tmp(MESH.Solid.internal_dof));
MESH.Solid.II        = setdiff(1:length(MESH.Solid.internal_dof),MESH.Solid.Gamma);

MESH.Solid.II_global = setdiff(MESH.Solid.internal_dof,MESH.Solid.internal_dof(MESH.Solid.Gamma));
MESH.Solid.Gamma_global = S_interface_dofs;

% All free nodes for both domains
MESH.internal_dof    = [MESH.Fluid.II MESH.Fluid.Gamma' ...
            length(MESH.Fluid.internal_dof)+1:(length(MESH.Fluid.internal_dof)+length(MESH.Solid.II))];      
end