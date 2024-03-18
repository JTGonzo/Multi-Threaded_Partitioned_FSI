function [A_in, F_in, u_D] =  CSM_ApplyBC(A, F, FE_SPACE, MESH, DATA, t, zero_Dirichlet)
%% Initialize some function parameters
if nargin < 6
    t = [];
end

if isempty(A)
    A = sparse(MESH.numNodes*MESH.dim, MESH.numNodes*MESH.dim);
end

if isempty(F)
    F = sparse(MESH.numNodes*MESH.dim, 1);
end

if nargin < 7
    zero_Dirichlet = 0;
end

param = DATA.param;

u_D = [];
        
%% Neumann condition
   for k = 1 : MESH.dim
       if ~isempty(MESH.Neumann_side{k})
          % evaluate the surface quadrature nodes and basis function      
          [quad_points, wi, phi] = surf_quad( );
           
           csi = quad_points(1,:);
           eta = quad_points(2,:);
           eta1           =  1-csi-eta;
           nqn            =  length(wi);
                
           nof         = length(MESH.Neumann_side{k});
           nbn         = MESH.numBoundaryDof;
                
           Rrows       = zeros(nbn*nof,1);
           Rcoef       = Rrows;
           
           % compute global coordinates of quadrature nodes 
           xlt = zeros(nof,nqn); ylt = xlt; zlt = xlt;
           coord_ref = [eta1; csi; eta];          
           for j = 1 : 2
               dof = MESH.boundaries(j,MESH.Neumann_side{k});
               vtemp = MESH.vertices(1,dof);
               xlt = xlt + vtemp'*coord_ref(j,:);
               vtemp = MESH.vertices(2,dof);
               ylt = ylt + vtemp'*coord_ref(j,:);
               vtemp = MESH.vertices(3,dof);
               zlt = zlt + vtemp'*coord_ref(j,:);
           end
           
           % Evaluate/import applied neumann at the quadrature nodes
            u_Neumann = DATA.bcNeu{k}(xlt,ylt,zlt,t,param);
            one       = ones(nof,nqn);
            u_Neumann = u_Neumann.*one;
            
            % compute the boundary element surface area 
            x =  MESH.vertices(1,MESH.boundaries(1:3, MESH.Neumann_side{k}));
            y =  MESH.vertices(2,MESH.boundaries(1:3, MESH.Neumann_side{k}));
            z =  MESH.vertices(3,MESH.boundaries(1:3, MESH.Neumann_side{k}));
                
            areav = cross(  [x(2:3:end)-x(1:3:end);  y(2:3:end)-y(1:3:end);  z(2:3:end)-z(1:3:end)], ...
                  [x(3:3:end)-x(1:3:end);  y(3:3:end)-y(1:3:end);  z(3:3:end)-z(1:3:end)]);
            
            % compute the boundary element pressure force acting on the face normal   
            for l = 1 : nof  
                area   = 0.5*norm(areav(:,l));
                detjac = 2*area;
                    
                face = MESH.Neumann_side{k}(l);
                    
                u_Neumann_loc  = u_Neumann(l,:).*wi;
                u_Neumann_loc  = u_Neumann_loc(1,:)';
                    
                Rrows(1+(l-1)*nbn:l*nbn) = MESH.boundaries(1:nbn,face);
                Rcoef(1+(l-1)*nbn:l*nbn) = detjac*phi*u_Neumann_loc;
            end
             %  sparse vector of boundary nuemann condition
             F = F + sparse(Rrows+(k-1)*MESH.numNodes,1,Rcoef,MESH.dim*MESH.numNodes,1);                
        end
   end        
                
%% Dirichlet condition
% extract coordinates of dirichlet nodes and evaluate the BC
   for k = 1 : 3           
       x  = MESH.nodes(1,MESH.Dirichlet_dof_c{k});
       y  = MESH.nodes(2,MESH.Dirichlet_dof_c{k});
       z  = MESH.nodes(3,MESH.Dirichlet_dof_c{k});
       u_Dirichlet{k} = DATA.bcDir{k}(x,y,z,t,param);                
       u_D = [u_D; u_Dirichlet{k}'];
   end

% either a zero or non-zero enforced condition   
u_D  = u_D * (1 - zero_Dirichlet);

%% Apply the BCs to the LHS/RHS matrices
% output the reduced LHS and RHS matrices for the free internal nodes 
F_in = F(MESH.internal_dof) - A(MESH.internal_dof,MESH.Dirichlet_dof)*u_D;    
A_in = A(MESH.internal_dof,MESH.internal_dof);
    
end
