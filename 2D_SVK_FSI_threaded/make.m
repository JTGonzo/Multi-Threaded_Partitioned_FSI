function makenew(hasOpenMP, sources)

if nargin < 1 || isempty(hasOpenMP)
    hasOpenMP = 0;
end

%% Compile C code

if hasOpenMP
    fprintf('\nCompiling with openmp enabled\n');
    Flags = 'CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"';

else
    fprintf('\nCompiling without openmp\n');
    Flags = '';
end

source_files{1} = {'C_Files/','Mass_assembler_C_omp.c'};
dependencies{1} = {};
source_files{2} = {'C_Files/','CSM_assembler_ExtForces.c'};
dependencies{2} = {};
source_files{3} = {'C_Files/','CSM_assembler_C_omp.c'};
dependencies{3} = {'Tools.c', 'StVenantKirchhoffMaterial.c', 'SEMMTMaterial.c'};
source_files{4} = {'C_Files/','CFD_assembler_C_omp.c'};
dependencies{4} = {'Tools.c'};
source_files{5} = {'C_Files/','CFD_assembler_ExtForces.c'};
dependencies{5} = {};

sources = 1:length(source_files);    
n_sources = length(sources); 
k = 0;

for i = sources
    
    k = k + 1;
    
    fprintf('\n ** Compiling source Nr. %d; %d of %d \n', i, k, n_sources)
    file_path = [pwd, '/', source_files{i}{1}];
    file_name = source_files{i}{2};
    
    all_dep = '';
    if ~isempty( dependencies{i} )
        for j = 1 : length( dependencies{i})
            
            this_dep =  sprintf('%s%s', file_path, dependencies{i}{j});
            all_dep = [all_dep, ' ', this_dep]; 
        end
        
    end
    
    mex_command = sprintf( 'mex %s%s %s %s -outdir %s', file_path, file_name, all_dep, Flags, file_path);
    eval( mex_command );
end

end