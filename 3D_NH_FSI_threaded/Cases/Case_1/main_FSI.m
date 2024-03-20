% This solver was largely built starting from the redbkit monolithic FSI code which was
% deconstructed to its baseline elements, reassembled and added to (as needed), and used 
% to solve specific FSI benchmarks in a Partitioned fashion for subsequent coupling analysis

clear all; close all; clc;

%% Initialize problem and solver data 
init_data;

%% Build finite element meshes and initialize boundary conditions

[MESH.Fluid, FE_SPACE_v, FE_SPACE_p] = build_F( dim, meshFluid.elements, meshFluid.vertices, meshFluid.boundaries, fem_F{1}, quad_order, DATA.Fluid, 'CFD', meshFluid.rings );
[MESH.Solid, FE_SPACE_s ] = build_S( dim, meshSolid.elements, meshSolid.vertices, meshSolid.boundaries, fem_S, quad_order, DATA.Solid, 'CSM', meshSolid.rings );
[FE_SPACE_g] = build_G( MESH.Fluid, fem_F{1}, dim, quad_order );%

MESH.Fluid.internal_dof_c{MESH.dim+1} = 1:FE_SPACE_p.numDof;

%% Definte interface mapping and segregate domains (internal/interface)
[MESH] = interface(DATA, MESH, FE_SPACE_v, FE_SPACE_p, FE_SPACE_s, MESH.mapper, dim);

% undisturbed/orginal node coordinates for mesh updating
Fluid_ReferenceNodes = MESH.Fluid.nodes(1:dim,:);

%% Initialize fluid and solid domain variables
init_fluid;
init_solid;

%% Initialize cpupling and acceleration scheme variables
init_couple;

%% Initalize Geometry Jacobian
[MESH, DATA, Solid_Extension] = meshjac(DATA, MESH, FE_SPACE_v, FE_SPACE_g, dim);

tol        = DATA.Solid.NonLinearSolver.tol;
maxIter    = DATA.Solid.NonLinearSolver.maxit;

uS_n = u0;
t_iters = 0;

%% Time Loop
while ( t < tf )
    
    iter_time = tic;
    cellfun(@increase_time,objects); % update current step info for all sub-functions 
        
    t   = t   + dt;
    k_t = k_t + 1;
    
    fprintf('\n=========================================================================')
    fprintf('\n==========  t0 = %2.4f  t = %2.4f  tf = %2.4f\n', t0, t, tf );
    
    % Enter the non-linear FP iteration loop    
    for nLiter = 1:Couple.imax
        
        t_iters = t_iters + 1;
        
        % accelerate the FP convergence of the approximate solid displacement 
        [uS, omega] = accelerator( MESH, uS, r, nLiter, extrapolator, model, Couple, omega, Acc_flag); 
         
        % apply Dirichlet BC to solid nodes
        [~, ~, u_D]   =  CSM_ApplyBC([], [], FE_SPACE_s, MESH.Solid, DATA.Solid, t, 1);
        uS(MESH.Solid.Dirichlet_dof) = u_D;
        
        %clear MESH.Fluid.vertices MESH.Fluid.nodes MESH.Fluid.jac MESH.Fluid.invjac MESH.Fluid.h
        % move the ALE mesh per the solid displacement and compute the ALE velocity
        if (nLiter == 1) && (Couple.extra == 1)     
             [~, ~, ALE_velocity] = movemesh(MESH, Fluid_ReferenceNodes, Solid_Extension, TimeAdvanceS, uS, d_Fn, dim, dt, nLiter, ALE_velocity, Couple);
        else
             [MESH, d_Fn_t, ALE_velocity] = movemesh(MESH, Fluid_ReferenceNodes, Solid_Extension, TimeAdvanceS, uS, d_Fn, dim, dt, nLiter, ALE_velocity, Couple);
        end

		uF_l = uF;        
        trac_l = trac;
        
        % solve the Navier-Stokes equations (whatever formulation)
        [uF, C_NS, F_NS, trac] =  solve_fluid(MESH, DATA, FE_SPACE_v, FE_SPACE_p, TimeAdvanceF, TimeAdvanceS, t, dt, ALE_velocity, uS, Couple, nLiter);
        
        % solve the solid elasticity equations (whatever formulation/constitutive model)
        [uS_t] = solve_solid(FE_SPACE_s, FE_SPACE_v, MESH, DATA, SolidModel, TimeAdvanceS, uS_n, t, Coef_Mass, M_s, A_robin, trac);   
        
        % compute the FP displacement solution residual
        r =  uS_t - uS;
		rF = uF_l - uF; 		
        rT = trac_l - trac;
        
        % Retain iteration residual data
        if ((Acc_flag >= 3) && (Acc_flag ~= 4))
            if Acc_flag == 3
                model.add(r(MESH.Solid.Gamma_global))
            else
                model.add(r(MESH.Solid.Gamma_global), uS_t(MESH.Solid.Gamma_global))  
            end
        end
        
       % Export to solid displacement (VTK)
        if ~isempty(vtk_filename2)
            CSM_export_solution(MESH.dim, r, MESH.Solid.vertices, ...
                MESH.Solid.elements, MESH.Solid.numNodes, [vtk_filename2, 'Residual'], t_iters, 'DisplacementResidual'); 
        end
        
        % Check convergence criteria
        convergence.add(r,r(MESH.Solid.Gamma_global),rF(MESH.Fluid.Gamma_global),rF(dim*MESH.Fluid.numNodes+p_dofs),rT(MESH.Fluid.Gamma_global));
        if (convergence.is_satisfied())
            break;
        end                
    end

%% Store current time-step data and update time-integration variables
    if Acc_flag == 0
        uS = uS_t;
    end
    
    uS_n = uS;
    d_Fn = d_Fn_t;
    
    % update time advance variables for integration 
    TimeAdvanceS.Update( uS );
    TimeAdvanceF.Append( uF(1:FE_SPACE_v.numDof) );
    
    % Extrapolate initial estimate for next time step
    extrapolator.add(uS(MESH.Solid.Gamma_global)); 
    
%% Export Domain and Iteration Data 
   
    if (mod(k_t,Couple.probeOutFreq)==0)       
        % Export FP iterations data
        fprintf(fileId1,'%d %d\n',k_t,nLiter);  
        
        % Export watchpoint physical data
        fprintf(fileId2,'%d %d %d %d\n',t,uS(Wid,1),uS(MESH.Solid.numNodes+Wid,1),uS(2*MESH.Solid.numNodes+Wid,1));
        
        % Export watchpoint physical data
        fprintf(fileId3,'%d %d %d %d\n',t,uF(WidF,1),uF(MESH.Fluid.numNodes+WidF,1),uF(2*MESH.Fluid.numNodes+WidF,1));
        
        % Export lift and drag coefficients
        [~] = AeroForces(MESH, DATA, k_t, t, uF, C_NS, F_NS, FE_SPACE_v, FE_SPACE_p, dofs_drag, fileDragLift);
        
    end
    if (mod(k_t,Couple.probeOutFreq2)==0) 
        % Export to solid displacement (VTK)
        if ~isempty(vtk_filename)
            CSM_export_solution(MESH.dim, uS, MESH.Solid.vertices, ...
                MESH.Solid.elements, MESH.Solid.numNodes, [vtk_filename, 'Solid'], k_t); %Us_k
        end
        
        % Fluid velocity and pressure on deformed mesh    
        if ~isempty(vtk_filename)
            CFD_export_solution(dim, uF(1:FE_SPACE_v.numDof), uF(1+FE_SPACE_v.numDof:end), ...
                MESH.Fluid.vertices, MESH.Fluid.elements, MESH.Fluid.numNodes, [vtk_filename,'Fluid'], k_t);
        end
    end
    
    % view time step duration for monitoring 
    iter_time = toc(iter_time);
    fprintf('\n-------------- iter time: %3.2f s -----------------',iter_time);
end

% close all output files
fclose('all'); 
