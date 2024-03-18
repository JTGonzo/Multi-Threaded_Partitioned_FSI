function [A_in, F_in, u_D, u_gamma] =  CFD_ApplyBC(A, F, FE_SPACE, FE_SPACE_p, MESH, DATA, t, TimeAdvanceS, uS, Couple, nLiter, zero_Dirichlet)
%% Initialize some function parameters
if nargin < 7
    t = [];
end

if isempty(A)
    A = sparse(FE_SPACE.numDof + FE_SPACE_p.numDof, FE_SPACE.numDof + FE_SPACE_p.numDof);
end

if isempty(F)
    F = sparse(FE_SPACE.numDof + FE_SPACE_p.numDof, 1);
end

if nargin < 12
    zero_Dirichlet = 0;
end

param = DATA.param;
u_D = [];
u_gamma = [];

%% Pressure condition
for k = 1 : MESH.Fluid.dim
	if ~isempty(MESH.Fluid.Pressure_side{k})
		% evaluate the surface quadrature nodes and basis function
        [quad_points, wi, phi] = surf_quad( );           
        
        csi = quad_points(1,:);
        eta = quad_points(2,:);          
		eta1           =  1-csi-eta;
		nqn            =  length(wi);
		
		nbn         = MESH.Fluid.numBoundaryDof;
		
		for flag = 1 : length(DATA.flag_pressure{k})
			
			nof         = length(MESH.Fluid.Pressure_side_CompFlag{k,flag});

			Rrows       = zeros(nbn*nof,1);
			Rcoef       = Rrows;
			
            % compute global coordinates of quadrature nodes 
			xlt = zeros(nof,nqn); ylt = xlt; zlt = xlt;
			coord_ref = [eta1; csi; eta];
			for j = 1 : 2
				dof = MESH.Fluid.boundaries(j,MESH.Fluid.Pressure_side_CompFlag{k,flag});
				vtemp = MESH.Fluid.vertices(1,dof);
				xlt = xlt + vtemp'*coord_ref(j,:);
				vtemp = MESH.Fluid.vertices(2,dof);
				ylt = ylt + vtemp'*coord_ref(j,:);
				vtemp = MESH.Fluid.vertices(3,dof);
				zlt = zlt + vtemp'*coord_ref(j,:);
            end
			
            % Evaluate/import applied pressures at the quadrature nodes
			if length(DATA.bcPrex) == 1
				pressure = DATA.bcPrex(xlt,ylt,zlt,t,param);
			else
				pressure = DATA.bcPrex{DATA.flag_pressure{k}(flag)}(xlt,ylt,zlt,t,param);
            end			
			one      = ones(nof,nqn);
			pressure = pressure.*one;
			
            % compute the boundary elemet length  
			x    =  MESH.Fluid.vertices(1,MESH.Fluid.boundaries(1:3, MESH.Fluid.Pressure_side_CompFlag{k,flag}));
			y    =  MESH.Fluid.vertices(2,MESH.Fluid.boundaries(1:3, MESH.Fluid.Pressure_side_CompFlag{k,flag}));
			z    =  MESH.Fluid.vertices(3,MESH.Fluid.boundaries(1:3, MESH.Fluid.Pressure_side_CompFlag{k,flag}));
			
			areav = cross(  [x(2:3:end)-x(1:3:end);  y(2:3:end)-y(1:3:end);  z(2:3:end)-z(1:3:end)], ...
				[x(3:3:end)-x(1:3:end);  y(3:3:end)-y(1:3:end);  z(3:3:end)-z(1:3:end)]);
			
            % compute the boundary element pressure force acting on the face normal
			for l = 1 : nof				
				area   = 0.5*norm(areav(:,l));
				detjac = 2*area;
				
				face = MESH.Fluid.Pressure_side_CompFlag{k,flag}(l);
				
				pressure_loc  = pressure(l,:).*wi;
				pressure_loc  = pressure_loc(1,:)';
				
				Rrows(1+(l-1)*nbn:l*nbn)    = MESH.Fluid.boundaries(1:nbn,face);
				Rcoef(1+(l-1)*nbn:l*nbn)    = MESH.Fluid.Normal_Faces(k,face)*detjac*phi*pressure_loc;
            end
            %  sparse vector of boundary pressure forces
			F = F + fsparse(Rrows+(k-1)*MESH.Fluid.numNodes,1,Rcoef,[FE_SPACE.numDof + FE_SPACE_p.numDof,1]);		
			%F = F + sparse(Rrows+(k-1)*MESH.Fluid.numNodes,1,Rcoef,FE_SPACE.numDof + FE_SPACE_p.numDof,1);			
		end
	end
end

%% Dirichlet condition
% extract coordinates of dirichlet nodes and evaluate the BC
for k = 1 : MESH.Fluid.dim
      x  = MESH.Fluid.nodes(1,MESH.Fluid.Dirichlet_dof_c{k});
      y  = MESH.Fluid.nodes(2,MESH.Fluid.Dirichlet_dof_c{k});
      z  = MESH.Fluid.nodes(3,MESH.Fluid.Dirichlet_dof_c{k});
      u_Dirichlet{k} = DATA.bcDir{k}(x,y,z,t,param);
      u_D = [u_D; u_Dirichlet{k}'];
end

% either a zero or non-zero enforced condition
u_D  = u_D * (1 - zero_Dirichlet);

% extract the solid interface velocity
if (nLiter == 1) && (Couple.extra == 1) 
    duS = TimeAdvanceS.velocityNL1( );
else
    duS = TimeAdvanceS.velocity(uS);
end
for k = 1 : MESH.Fluid.dim
    u_inc = zeros(MESH.ndof_interface{k},1);
    tmp  = duS(MESH.Solid.numNodes*(k-1)+MESH.Solid.dof_interface{k});
    if MESH.mapper == 1
            u_inc =  tmp(MESH.Interface_SFmap{k});
    else
        for i = 1:length(MESH.coeff_SFmap{k}(1,:))
            u_add = MESH.coeff_SFmap{k}(:,i).*tmp(MESH.Interface_SFmap{k}(:,i));
            u_inc = u_inc + u_add;           
        end
    end
    u_gamma = [u_gamma; u_inc];
end

% boundary and interface Dirichlet conditions
u_BC = [u_D; u_gamma];

%% Apply the BCs to the LHS/RHS matrices
% output the reduced LHS and RHS NS-matrices for the free internal nodes 
F_in = F(MESH.Fluid.II_global) - A(MESH.Fluid.II_global,[MESH.Fluid.Dirichlet_dof; MESH.Fluid.Gamma_global])*u_BC;
A_in = A(MESH.Fluid.II_global,MESH.Fluid.II_global);

end
