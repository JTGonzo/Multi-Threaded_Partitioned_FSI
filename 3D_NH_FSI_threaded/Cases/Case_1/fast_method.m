loca = pwd;

cd ../../stenglib/

s = pwd;
addpath([s '/Tensor']);
addpath([s '/Fast']);
addpath([s '/Utils']);
addpath([s '/Scicomp']);
addpath([s '/Misc']);

cd Fast/
startup
make('openmp',true,'fsparseonly',1)

cd(loca)

main_FSI
