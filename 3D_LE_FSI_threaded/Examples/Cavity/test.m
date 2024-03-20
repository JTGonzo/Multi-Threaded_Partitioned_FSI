close all; clear all; clc;

load('Cavity_S.mat') %'flap_S.mat' 'Mok_S.mat' 'Cavity_S.mat'
meshSolid.boundaries = boundaries;
meshSolid.elements = elements;
meshSolid.vertices = vertices;
meshSolid.rings = rings;

% Fluid mesh
load('Cavity_F.mat') ; %'flap_F.mat' 'Mok_S.mat' 'Cavity_S.mat'
meshFluid.boundaries = boundaries;
meshFluid.elements = elements;
meshFluid.vertices = vertices;
meshFluid.rings = rings;

aa = 1;