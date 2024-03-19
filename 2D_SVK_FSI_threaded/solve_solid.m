function [uS_t] = solve_solid(FE_SPACE_s, FE_SPACE_v, MESH, DATA, SolidModel, TimeAdvanceS, uS_n, t, Coef_Mass, M_s, A_robin, trac)
%% Newton Method/Iter Parameters
    tol        = DATA.Solid.NonLinearSolver.tol;
    resRelNorm = tol + 1;
    incrNorm   = tol + 1;
    maxIter    = DATA.Solid.NonLinearSolver.maxit;

%% Initialize iterated solid variable and increment     
    [~, ~, u_D]   =  CSM_ApplyBC([], [], FE_SPACE_s, MESH.Solid, DATA.Solid, t);
    dU            = zeros(MESH.Solid.numNodes*MESH.dim,1);  %| solution increment
    Us_k          = uS_n(:,end);                            %| iterated solid solution
    Us_k(MESH.Solid.Dirichlet_dof) = u_D;                   %| apply Dirichlet conditions
    
    Csi = TimeAdvanceS.RhsContribute( );

%% Compute/apply internal and external forces   
    % Assemble matrix and right-hand side
    %F_ext = SolidModel.compute_volumetric_forces( (1 - TimeAdvanceS.M_alpha_f) * t + TimeAdvanceS.M_alpha_f * (t-dt) );
    %F_ext = SolidModel.compute_external_forces( (1 - TimeAdvanceS.M_alpha_f) * t + TimeAdvanceS.M_alpha_f * (t-dt) );
    
    % directly transfer te boundary node tractions to the solid interface     
    F_ext = zeros(FE_SPACE_s.numDof,1);
    for k = 1 : MESH.Fluid.dim
        f_inc = zeros(MESH.ndof_interface{k},1);
        tmp = trac(FE_SPACE_v.numDofScalar*(k-1)+MESH.Fluid.dof_interface{k});  
        if MESH.mapper == 1
            F_ext(FE_SPACE_s.numDofScalar*(k-1)+MESH.Solid.dof_interface{k}) = tmp(MESH.Interface_FSmap{k});
        else
            for i = 1:length(MESH.coeff_FSmap{k}(1,:))
                f_add = MESH.coeff_FSmap{k}(:,i).*tmp(MESH.Interface_FSmap{k}(:,i)); 
                f_inc = f_inc + f_add;
            end
            F_ext(FE_SPACE_s.numDofScalar*(k-1)+MESH.Solid.dof_interface{k}) = f_inc;
        end
    end
    
    k = 1;
     
    % compute internal solid forces
    F_in  = SolidModel.compute_internal_forces( (1 - TimeAdvanceS.M_alpha_f) * Us_k + TimeAdvanceS.M_alpha_f * uS_n );

%%  construct the elasticity algebraic formulation and apply boundary conditions
    % RHS residual array (force balance)
    Residual  = Coef_Mass * M_s * Us_k + F_in - F_ext - M_s * Csi ...
                + A_robin * ((1 - TimeAdvanceS.M_alpha_f) * Us_k + TimeAdvanceS.M_alpha_f * uS_n);
    
    % solid model Jacobian
    dF_in     = SolidModel.compute_jacobian( (1 - TimeAdvanceS.M_alpha_f) * Us_k + TimeAdvanceS.M_alpha_f * uS_n );
    
    % LHS Matrix
    Jacobian  = Coef_Mass * M_s + (1 - TimeAdvanceS.M_alpha_f) * dF_in + A_robin * (1 - TimeAdvanceS.M_alpha_f);
    
    % Apply boundary conditions
    [C_STR, F_STR] =  CSM_ApplyBC(Jacobian, -Residual, FE_SPACE_s, MESH.Solid, DATA.Solid, t, 1);
    
    % initial residual norm
    res0Norm = norm(F_STR);
 
%% Solve the linear system via Newton Iterations
    % fprintf('\n============ Start Newton Iterations ============\n\n');
    % iterate until the relative residual norm convergence condition is met
    while (k <= maxIter && incrNorm > tol && resRelNorm > tol)
        
        % Solve the linearized system of equation for the deformation increment 
        dU(MESH.Solid.internal_dof) = C_STR \ F_STR;
        
        % update the iterated displacement solution
        Us_k     = Us_k + dU;
        
        % compute the current iteration's relative increment norm
        incrNorm = norm(dU)/norm(Us_k);
        
        % re-assemble LHS matrix and right-hand side residual array
        F_in  = SolidModel.compute_internal_forces( (1 - TimeAdvanceS.M_alpha_f) * Us_k + TimeAdvanceS.M_alpha_f * uS_n );
        dF_in = SolidModel.compute_jacobian( (1 - TimeAdvanceS.M_alpha_f) * Us_k + TimeAdvanceS.M_alpha_f * uS_n );
        
        Residual  = Coef_Mass * M_s * Us_k + F_in - F_ext - M_s * Csi ...
                    + A_robin * ((1 - TimeAdvanceS.M_alpha_f) * Us_k + TimeAdvanceS.M_alpha_f * uS_n);
            
        Jacobian  = Coef_Mass * M_s + (1 - TimeAdvanceS.M_alpha_f) * dF_in + A_robin * (1 - TimeAdvanceS.M_alpha_f);
        
        % apply boundary conditions
        [C_STR, F_STR]   =  CSM_ApplyBC(Jacobian, -Residual, FE_SPACE_s, MESH.Solid, DATA.Solid, t, 1);
        
        % evaluate the relative residual norm
        resRelNorm = norm(F_STR) / res0Norm;
        
        %fprintf('\n **** Iteration  k = %d:  norm(dU)/norm(Uk) = %1.2e, Residual Rel Norm = %1.2e \n\n',k,full(incrNorm), full(norm(resRelNorm)));
        k = k + 1;   % update Newton increment count    
    end
    
    % update/output the domain solution
    uS_t = Us_k;
end