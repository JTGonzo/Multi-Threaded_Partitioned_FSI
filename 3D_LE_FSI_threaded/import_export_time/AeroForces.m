function [Z] = AeroForces(MESH, DATA, k_t, t, u, C_NS, F_NS, FE_SPACE_v, FE_SPACE_p, dofs_drag, fileDragLift)
%% Compute Aerodynamic Forces on boundary faces
% identify the solution nodes used to compute forces
Z              = zeros(FE_SPACE_v.numDofScalar,1);
Z(dofs_drag)   = 1;

% compute x-force (drag)
W               = zeros(FE_SPACE_v.numDof+FE_SPACE_p.numDof,1);
W(1:FE_SPACE_v.numDofScalar)        = Z;
AeroF_x(k_t+1) = DATA.Fluid.Output.DragLift.factor*(W'*(-C_NS*u + F_NS));

% compute y-force (lift)
W               = zeros(FE_SPACE_v.numDof+FE_SPACE_p.numDof,1);
W(FE_SPACE_v.numDofScalar+(1:FE_SPACE_v.numDofScalar))  = Z;
AeroF_y(k_t+1)  = DATA.Fluid.Output.DragLift.factor*(W'*(-C_NS*u  + F_NS));

% compute z-force (twist)
if MESH.dim == 3
	W               = zeros(FE_SPACE_v.numDof+FE_SPACE_p.numDof,1);
	W(2*FE_SPACE_v.numDofScalar+(1:FE_SPACE_v.numDofScalar))  = Z;
	AeroF_z(k_t+1)  = DATA.Fluid.Output.DragLift.factor*(W'*(-C_NS*u  + F_NS));
else
	AeroF_z(k_t+1) = 0.0;
end

% print and output coefficients
%fprintf('\n *** F_x = %e, F_y = %e, F_z = %e *** \n',  AeroF_x(k_t+1), AeroF_y(k_t+1), AeroF_z(k_t+1));
fprintf(fileDragLift, '\n%1.4e  %1.4e  %1.4e  %1.4e', t, AeroF_x(k_t+1), AeroF_y(k_t+1), AeroF_z(k_t+1));
