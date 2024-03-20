%% Definte paths, folders and outputs locations
[~,~,~] = mkdir('Figures');
[~,~,~] = mkdir('Results');

wrkDir = './' ;
problemString = 'TB_3DL1T1Z1'; 
addpath('../../')
addpath('../../C_Files/')
addpath('../../Mappers/')
addpath('../../Accelerators/')
addpath('../../Accelerators/Filters/')
addpath('../../import_export_time/')
vtk_filename = 'Figures/TB_3DL1T1Z1'; % []; % 
vtk_filename2 = []; 

t = [];
param = [];
dim  =  3;
MESH.dim  = dim; % problem dimensionality

%% Load problem mesh and element properties
% Solid mesh
%load('Tube_S_fine.mat')
load('../../../Geometry_Files/Tube_MATS/TB_3DL1T1Z1_S.mat')
meshSolid.boundaries = boundaries;
meshSolid.elements = elements;
meshSolid.vertices = vertices;
meshSolid.rings = rings;

% Fluid mesh
%load('Tube_F_fine.mat')
load('../../../Geometry_Files/Tube_MATS/TB_3DL1T1Z1_F.mat')
meshFluid.boundaries = boundaries;
meshFluid.elements = elements;
meshFluid.vertices = vertices;
meshFluid.rings = rings;

use_SUPG = false;
quad_order = 5;     % for 2D problem quad_order = 4
fem_F = {'P2', 'P1'};
fem_S = 'P2';
MESH.mapper = 1;                  %| NN = 1; Linear = 2; RBF =3

%% Load problem, solver and domain variable details
data_file_F = 'CFD_data';
data_file_S = 'CSM_data';

% read fluid domain parameters/settings
eval(data_file_F);
data_fields_name = fieldnames(data);

for i = 1 : length(data_fields_name)    
     eval(['DATA.Fluid.',data_fields_name{i}, '=', 'data.',data_fields_name{i},';']);
end
DATA.Fluid.param = param;

clear data_fields_name data

% read solid domain parameters/settings
eval(data_file_S);
data_fields_name = fieldnames(data);

for i = 1 : length(data_fields_name)   
    eval(['DATA.Solid.',data_fields_name{i}, '=', 'data.',data_fields_name{i},';']);
end
DATA.Solid.param = param;
