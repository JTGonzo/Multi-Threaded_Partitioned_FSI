function [uF, C_NS, F_NS, trac] =  solve_fluid(MESH, DATA, FE_SPACE_v, FE_SPACE_p, TimeAdvanceF, TimeAdvanceS, t, dt, ALE_velocity, uS, Couple, nLiter)
%% BDF scheme to approximate fluid velocity and extrapolate convective velocity
    % BDF approx.
    v_BDF = TimeAdvanceF.RhsContribute( );
    u_BDF = [v_BDF; zeros(FE_SPACE_p.numDof,1)];
    alphaF = TimeAdvanceF.GetCoefficientDerivative();
    
    % extrapolated convective velocity
    v_extrapolated = TimeAdvanceF.Extrapolate();

%% initialize solution space and NS-solver parameters
    Uf_k  = zeros(FE_SPACE_v.numDof + FE_SPACE_p.numDof,1);
    
    % set material properties and FE info
    FluidModel = CFD_Assembler( MESH.Fluid, DATA.Fluid, FE_SPACE_v, FE_SPACE_p );

%%  compile the solver matrices  
    % external or body forces
    %F_gravity = FluidModel.compute_external_forces();
    
    % divergence and stiffness matrices
    [A_Stokes] = FluidModel.compute_Stokes_matrix();
    
    % velocity and pressure block mass matrices 
    Mv = FluidModel.compute_mass_velocity();
    Mp = FluidModel.compute_mass_pressure();
    M_f  = blkdiag(DATA.Fluid.density * Mv, 0*Mp);
    
    % fluid convective matrix
    [C1] = FluidModel.compute_convective_Oseen_matrix( v_extrapolated - ALE_velocity);

%%  construct the NS algebraic formulation and apply boundary conditions
    F_NS = 1/dt * M_f * u_BDF ;  % RHS vector
    %F_NS = 1/dt * M_f * u_BDF - F_gravity;  %| RHS vector
    C_NS = alphaF/dt * M_f + A_Stokes + C1;  %| LHS matrix
                      
    [C_NS_in, F_NS_in, v_D, v_gamma] =  CFD_ApplyBC(C_NS, F_NS, FE_SPACE_v, FE_SPACE_p, MESH, DATA.Fluid, t, TimeAdvanceS, uS, Couple, nLiter);
    
%%  solve the linear system of equation for free nodes and enforce boundary solns'    
    Uf_k(MESH.Fluid.II_global) = C_NS_in \ F_NS_in;  %| free internal nodes
    Uf_k(MESH.Fluid.Dirichlet_dof) = v_D;            %| boundary nodes   
    Uf_k(MESH.Fluid.Gamma_global) = v_gamma;         %| interface nodes

%% update the fluid variables and compute the discrete nodal tractions    
    uF = Uf_k;  
    
    % discrete nodal forces
    trac = -C_NS*uF + F_NS;
    
%     Tx_f = trac(MESH.Fluid.dof_interface{1});
%     Ty_f = trac(FE_SPACE_v.numDofScalar+MESH.Fluid.dof_interface{2});  

end