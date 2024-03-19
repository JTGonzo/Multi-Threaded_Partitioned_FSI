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

%% Initialize Linear/Non-Linear Solver Settings
tol        = DATA.Solid.NonLinearSolver.tol;
maxIter    = DATA.Solid.NonLinearSolver.maxit;

uS_n = u0; 
uFprev = uF;

%% Time Loop
while ( t < tf )
    
    iter_time = tic;
    
    t   = t   + dt;
    k_t = k_t + 1;
    
    fprintf('\n=========================================================================')
    fprintf('\n==========  t0 = %2.4f  t = %2.4f  tf = %2.4f\n', t0, t, tf );
    
    % Enter the non-linear FP iteration loop    
    for nLiter = 1:Couple.imax

        if nLiter > 1     
            [MESH, d_Fn_t, ALE_velocity] = movemesh(MESH, Fluid_ReferenceNodes, Solid_Extension, TimeAdvanceS, uS, d_Fn, dim, dt, nLiter, ALE_velocity);
        else
            [~, ~, ALE_velocity] = movemesh(MESH, Fluid_ReferenceNodes, Solid_Extension, TimeAdvanceS, uS, d_Fn, dim, dt, nLiter, ALE_velocity);
        end
        
        % solve the Navier-Stokes equations (whatever formulation)
        [uF, C_NS, F_NS, trac] =  solve_fluid(MESH, DATA, FE_SPACE_v, FE_SPACE_p, TimeAdvanceF, TimeAdvanceS, t, dt, ALE_velocity, uS);
               
        % Check convergence criteria
%         if k_t > 2
%             convg = norm(uF - uFprev)/norm(uFprev);
%             fprintf('NS: %e, ',convg);
%             if convg < Couple.rtol
%                 break
%             end
%         else 
            if nLiter > 4
                break
            end
%         end                
    end

%% Store current time-step data and update time-integration variables
    d_Fn = d_Fn_t;
    uFprev = uF;
    
    % update time advance variables for integration 
    TimeAdvanceS.Update( uS );
    TimeAdvanceF.Append( uF(1:FE_SPACE_v.numDof) );
         
%% Export Domain and Iteration Data 
   
    if (mod(k_t,Couple.probeOutFreq)==0)       
        
        % Export watchpoint physical data
        fprintf(fileId3,'%d %d %d\n',t,uF(WidF,1),uF(MESH.Fluid.numNodes+WidF,1));
        
        % Export lift and drag coefficients
        [~] = AeroForces(MESH, DATA, k_t, t, uF, C_NS, F_NS, FE_SPACE_v, FE_SPACE_p, dofs_drag, fileDragLift);
               
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
fclose(fileDragLift);
fclose(fileId1);
fclose(fileId2);
fclose(fileId3);