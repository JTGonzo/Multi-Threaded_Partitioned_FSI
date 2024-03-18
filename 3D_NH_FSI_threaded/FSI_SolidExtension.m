function [d_F] = FSI_SolidExtension(MESH, Displacement, Solid_Extension)

% initialize variables and size of solution space 
d_F = zeros(MESH.Fluid.numNodes * MESH.dim, 1);

internal_dofs = Solid_Extension.internal_dofs;

interface_dofs  = [];
d_S             = [];

% transfer the displacement of the solid interface nodes to the
% corresponding fluid mesh nodes at interface
for k = 1 : MESH.dim
    d_inc = zeros(MESH.ndof_interface{k},1);
    interface_dofs    = [interface_dofs; MESH.Fluid.numNodes*(k-1)+MESH.Fluid.dof_interface{k}];
    tmp  = Displacement(MESH.Solid.numNodes*(k-1)+MESH.Solid.dof_interface{k});
    if MESH.mapper == 1
        d_inc =  tmp(MESH.Interface_SFmap{k});
    else
        for i = 1:length(MESH.coeff_SFmap{k}(1,:))    
             d_add = MESH.coeff_SFmap{k}(:,i).*tmp(MESH.Interface_SFmap{k}(:,i));  
             d_inc = d_inc + d_add;
        end
    end           
    d_S  = [d_S; d_inc];
end

% per the interface disp. compute the global mesh displacement at free nodes
F   = - Solid_Extension.matrix(internal_dofs,interface_dofs) * d_S;

x   = Solid_Extension.U \ (Solid_Extension.L \ F(Solid_Extension.perm));

d_F(internal_dofs) = x(Solid_Extension.invp);

% enforce the mesh displacement at the FSI boundary 
for k = 1 : MESH.dim   
    d_inc = zeros(MESH.ndof_interface{k},1);
    tmp     = Displacement(MESH.Solid.numNodes*(k-1)+MESH.Solid.dof_interface{k});
    if MESH.mapper == 1
        d_F(MESH.Fluid.numNodes*(k-1)+MESH.Fluid.dof_interface{k}) = tmp(MESH.Interface_SFmap{k}); 
    else
        for i = 1:length(MESH.coeff_SFmap{k}(1,:))    
             d_add = MESH.coeff_SFmap{k}(:,i).*tmp(MESH.Interface_SFmap{k}(:,i));  
             d_inc = d_inc + d_add;
        end
        d_F(MESH.Fluid.numNodes*(k-1)+MESH.Fluid.dof_interface{k}) = d_inc;
    end
end

d_F = reshape(d_F, MESH.Fluid.numNodes,  MESH.dim)';
end