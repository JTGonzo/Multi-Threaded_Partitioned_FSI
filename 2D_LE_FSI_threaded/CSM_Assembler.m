classdef CSM_Assembler < handle
    
    properties (GetAccess = public, SetAccess = protected)
        %% Set global handle variables
        M_MESH;
        M_DATA;
        M_FE_SPACE;
        M_subdomain;
        M_MaterialModel;
        M_MaterialParam;
    end
   
    methods
        %==========================================================================
        %% Initialize solid solver variables/space
        function obj = CSM_Assembler( MESH, DATA, FE_SPACE )
            obj.M_MESH      = MESH;
            obj.M_DATA      = DATA;
            obj.M_FE_SPACE  = FE_SPACE;
            obj.M_MaterialModel = DATA.Material_Model;
            obj = SetMaterialParameters(obj);      
        end
        
        %==========================================================================
        %% Set Constitutive Law Material Parameters
        function obj = SetMaterialParameters( obj )
            switch obj.M_MaterialModel
                case {'Linear', 'StVenantKirchhoff', 'NeoHookean'}
                    obj.M_MaterialParam = [obj.M_DATA.Young obj.M_DATA.Poisson];
                    
                case {'RaghavanVorp'}
                    obj.M_MaterialParam = [obj.M_DATA.Alpha obj.M_DATA.Beta obj.M_DATA.Bulk];
                     
                case 'SEMMT'
                    obj.M_MaterialParam = [obj.M_DATA.Young obj.M_DATA.Poisson obj.M_DATA.Stiffening_power];
            end    
        end
        
        %==========================================================================
        %% Compute Volumetric Forces
        function F_ext = compute_volumetric_forces( obj, t )
            
            if nargin < 2 || isempty(t)
                t = [];
            end
            
            % Computations of all quadrature nodes in the elements
            coord_ref = obj.M_MESH.chi;
            
			% compute global coordinates of quadrature nodes 
			x = zeros(obj.M_MESH.numElem,obj.M_FE_SPACE.numQuadNodes); y = x;
			for j = 1 : 3
				i = obj.M_MESH.elements(j,:);
				vtemp = obj.M_MESH.vertices(1,i);
				x = x + vtemp'*coord_ref(j,:);
				vtemp = obj.M_MESH.vertices(2,i);
				y = y + vtemp'*coord_ref(j,:);
			end
                     
			% Evaluate/import external forces at the quadrature nodes
			for k = 1 : obj.M_MESH.dim
				f{k}  = obj.M_DATA.force{k}(x,y,t,obj.M_DATA.param);
            end                    

            
            F_ext = [];
            for k = 1 : obj.M_MESH.dim
                
                [rowF, coefF] = CSM_assembler_ExtForces(f{k}, obj.M_MESH.elements, obj.M_FE_SPACE.numElemDof, ...
                    obj.M_FE_SPACE.quad_weights, obj.M_MESH.jac, obj.M_FE_SPACE.phi);
                
                %  sparse vector of external forces
                F_ext    = [F_ext; fsparse(rowF, 1, coefF, [obj.M_MESH.numNodes, 1])];				
                %F_ext    = [F_ext; sparse(rowF, 1, coefF, obj.M_MESH.numNodes, 1)]; 
            end 
        end
        
        %==========================================================================
        %% Compute Surface Forces
        function F = compute_surface_forces( obj, t )
            
            if nargin < 2 || isempty(t)
                t = [];
            end
            
            F = sparse(obj.M_MESH.numNodes*obj.M_MESH.dim, 1);

            % Pressure condition
            for k = 1 : obj.M_MESH.dim
                if ~isempty(obj.M_MESH.Pressure_side{k})
                    % compute the quadrature nodes of the surface element    
                    [csi,wi,phi]  =  bound_quad(obj.M_FE_SPACE.quad_order, 0, 1);
                    eta            =  1 - csi;
                    nqn            =  length(csi);

                    nof         = length(obj.M_MESH.Pressure_side{k});
                    nbn         = obj.M_MESH.numBoundaryDof;

                    Rrows       = zeros(nbn*nof,1);
                    Rcoef       = Rrows;
                    
                    % compute global coordinates of quadrature nodes 
                    xlt = zeros(nof,nqn); ylt = xlt;
                    coord_ref = [eta; csi];
                    for j = 1 : 2
                        dof = obj.M_MESH.boundaries(j,obj.M_MESH.Pressure_side{k});
                        vtemp = obj.M_MESH.vertices(1,dof);
                        xlt = xlt + vtemp'*coord_ref(j,:);
                        vtemp = obj.M_MESH.vertices(2,dof);
                        ylt = ylt + vtemp'*coord_ref(j,:);
                    end
                    
                    % Evaluate/import applied pressures at the quadrature nodes
                    pressure = obj.M_DATA.bcPrex(xlt,ylt,t,obj.M_DATA.param);
                    one       = ones(nof,nqn);
                    pressure = pressure.*one;
                    
                    % compute the boundary elemet length
                    x    =  obj.M_MESH.vertices(1,obj.M_MESH.boundaries(1:2, obj.M_MESH.Pressure_side{k}));
                    y    =  obj.M_MESH.vertices(2,obj.M_MESH.boundaries(1:2, obj.M_MESH.Pressure_side{k}));
                    side_length = sqrt((x(2:2:end)-x(1:2:end-1)).^2+(y(2:2:end)-y(1:2:end-1)).^2);
                    
                    % compute the boundary element pressure force acting on
                    % the face normal
                    for l = 1 : nof
                        face = obj.M_MESH.Pressure_side{k}(l);

                        pressure_loc  = pressure(l,:).*wi;
                        pressure_loc  = pressure_loc(1,:)';

                        Rrows(1+(l-1)*nbn:l*nbn)    = obj.M_MESH.boundaries(1:nbn,face);
                        Rcoef(1+(l-1)*nbn:l*nbn)    = obj.M_MESH.Normal_Faces(k,face)*side_length(l)*phi*pressure_loc;
                    end
                    %  sparse vector of surface pressure forces
                    F = F + sparse(Rrows+(k-1)*obj.M_MESH.numNodes,1,Rcoef,obj.M_MESH.dim*obj.M_MESH.numNodes,1);
                end
            end           
        end
        
        %==========================================================================
        %% Compute External Forces
        function F_ext = compute_external_forces( obj, t )
            % summation of volumetric and surface forces
            if nargin < 2 || isempty(t)
                t = [];
            end
            
            F_ext = compute_volumetric_forces( obj, t ) + compute_surface_forces( obj, t );
        end
        
        %==========================================================================
        %% Compute mass matrix
        function [M] = compute_mass( obj )
            
            [rowM, colM, coefM] = Mass_assembler_C_omp(obj.M_MESH.dim, obj.M_MESH.elements, obj.M_FE_SPACE.numElemDof, ...
                obj.M_FE_SPACE.quad_weights, obj.M_MESH.jac, obj.M_FE_SPACE.phi);
            
            % Build sparse matrix
            M_scalar   = fsparse(rowM, colM, coefM, [obj.M_MESH.numNodes, obj.M_MESH.numNodes]);			
            %M_scalar   = sparse(rowM, colM, coefM, obj.M_MESH.numNodes, obj.M_MESH.numNodes);
            M          = [];
            for k = 1 : obj.M_FE_SPACE.numComponents
                M = blkdiag(M, M_scalar);
            end
        end

        %==========================================================================
        %% Compute internal forces
        function [F_in] = compute_internal_forces(obj, U_h)
            
            [rowG, coefG] = ...
                CSM_assembler_C_omp(obj.M_MESH.dim, [obj.M_MaterialModel,'_forces'], obj.M_MaterialParam, full( U_h ), ...
                obj.M_MESH.elements, obj.M_FE_SPACE.numElemDof, ...
                obj.M_FE_SPACE.quad_weights, obj.M_MESH.invjac, obj.M_MESH.jac, obj.M_FE_SPACE.phi, obj.M_FE_SPACE.dphi_ref);
            
            % Build sparse matrix and vector
            F_in    = fsparse(rowG, 1, coefG, [obj.M_MESH.numNodes*obj.M_MESH.dim, 1]);			
            %F_in    = sparse(rowG, 1, coefG, obj.M_MESH.numNodes*obj.M_MESH.dim, 1);          
        end

        %==========================================================================
        %% Compute internal forces Jacobian
        function [dF_in] = compute_jacobian(obj, U_h)
            
            if nargin < 2 || isempty(U_h)
                U_h = zeros(obj.M_MESH.dim*obj.M_MESH.numNodes,1);
            end

            [rowdG, coldG, coefdG] = ...
                CSM_assembler_C_omp(obj.M_MESH.dim, [obj.M_MaterialModel,'_jacobian'], obj.M_MaterialParam, full( U_h ), ...
                obj.M_MESH.elements, obj.M_FE_SPACE.numElemDof, ...
                obj.M_FE_SPACE.quad_weights, obj.M_MESH.invjac, obj.M_MESH.jac, obj.M_FE_SPACE.phi, obj.M_FE_SPACE.dphi_ref);
            
            % Build sparse matrix and vector
            dF_in   = fsparse(rowdG, coldG, coefdG, [obj.M_MESH.numNodes*obj.M_MESH.dim, obj.M_MESH.numNodes*obj.M_MESH.dim]);  			
            %dF_in   = sparse(rowdG, coldG, coefdG, obj.M_MESH.numNodes*obj.M_MESH.dim, obj.M_MESH.numNodes*obj.M_MESH.dim);
        end
                
        %==========================================================================
        %% Assemble Robin Condition: Pn + K d = 0 on \Gamma_Robin, with K = ElasticCoefRobin
        function [A] = assemble_ElasticRobinBC(obj)
            
            A = sparse(obj.M_MESH.numNodes*obj.M_MESH.dim, obj.M_MESH.numNodes*obj.M_MESH.dim);
                       
      end
    end
end
