%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  Initialize Solid solver variable spaces, export files %%%%%%%
%%%%%% integration parameters, reference mesh and soln' tractions %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize variable space
% read displacement and velocity initial condition
u0  = [];
du0 = [];
for k = 1 : FE_SPACE_s.numComponents
      u0  = [u0; DATA.Solid.u0{k}(  MESH.Solid.nodes(1,:), MESH.Solid.nodes(2,:), MESH.Solid.nodes(3,:), t0, param )'];
      du0 = [du0; DATA.Solid.du0{k}( MESH.Solid.nodes(1,:), MESH.Solid.nodes(2,:), MESH.Solid.nodes(3,:), t0, param )'];
end

uS = u0;
uS_t = u0; % define a temporary solid variable for FP iterations

%% Initialize export data
% export initial condition (if it's the case)
if ~isempty(vtk_filename)
    CSM_export_solution(MESH.dim, u0, MESH.Solid.vertices, MESH.Solid.elements, MESH.Solid.numNodes, [vtk_filename,'Solid'], 0);
end

%% Create Solid Assembler Object
SolidModel = CSM_Assembler( MESH.Solid, DATA.Solid, FE_SPACE_s );

% Assemble mass matrix
M_s  =  SolidModel.compute_mass() * DATA.Solid.Density;

% Assemble Robin BC (if applicable)
A_robin = SolidModel.assemble_ElasticRobinBC();

%% Cpompute initial domain acceleration
F_ext_0  = SolidModel.compute_volumetric_forces(t0);   % external forces

F_in_0   =  SolidModel.compute_internal_forces( u0 ) + A_robin * u0; % internal forces

d2u0 = M_s \ (F_ext_0 - F_in_0); % solid acceleration

%% Initialize time advance scheme
TimeAdvanceS = GeneralizedAlpha_TimeAdvance( DATA.Solid.time.beta, DATA.Solid.time.gamma, DATA.Solid.time.alpha_m, DATA.Solid.time.alpha_f, dt );

Coef_Mass = TimeAdvanceS.MassCoefficient( );

TimeAdvanceS.Initialize( u0, du0, d2u0 );
