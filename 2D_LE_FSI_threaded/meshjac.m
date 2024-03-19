function [MESH, DATA, Solid_Extension] = meshjac(DATA, MESH, FE_SPACE_v, FE_SPACE_g, dim)
% Define the elastic mesh material properties
DATA.Geometry                  = DATA.Solid;
DATA.Geometry.Material_Model   = 'SEMMT';
DATA.Geometry.Stiffening_power = 0.8;

% assemble the geometry problem Jacobian matrix
MeshMotionAssembler = CSM_Assembler( MESH.Fluid, DATA.Geometry, FE_SPACE_g );
Solid_Extension.matrix  = MeshMotionAssembler.compute_jacobian( zeros(FE_SPACE_g.numDof, 1) );

% isolate solid-extension internal DoFs (excluding FSI interface)
internal_dofs_HE = [];
for k = 1 : dim
    internal_dofs_HE       = [ internal_dofs_HE (k-1)*FE_SPACE_v.numDofScalar+setdiff(1:FE_SPACE_v.numDofScalar, ...
        [MESH.Fluid.dof_interface{k}; MESH.ALE_dirichlet{k}])];
end
Solid_Extension.internal_dofs = internal_dofs_HE;

% compute LU factorization of mesh Jacobian (for internal dofs) and store it
[Solid_Extension.L , Solid_Extension.U,...
    Solid_Extension.perm , q ]   = lu(Solid_Extension.matrix(internal_dofs_HE,internal_dofs_HE), 'vector');
Solid_Extension.invp             = 0*q ;
Solid_Extension.invp(q)          = 1:length(q);

end
