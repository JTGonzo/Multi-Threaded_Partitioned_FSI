function CSM_export_solution(dim, u, vertices, elements, nov, outputFileName, iter, variableName)
%% Define output variables, labelling and formatting 
titleData = 'CSM_solution';
if nargin < 7
      iter = -1;
end

if nargin < 8
      variableName = 'StructureDisplacement';
end

novP1 = size(vertices,2);
    
      data = struct('iteration', {iter},...
            'vertices', {vertices'},...
            'elements', {elements(1:3,:)'},...
            'outputFile', {outputFileName},...
            'title', {titleData},...
            'variableName',{{variableName}},...
            'variableType',{{'VECTORS'}},...
            'variableData',{{[full(u(1:novP1)) full(u(nov+1:novP1+nov))  0*full(u(1:novP1))]'}});
      
timewrite = tic;

%% Define output space size, file formatting and naming
%  Get the size of your mesh
[ nodeNum, dim ] = size ( data.vertices );
[ elementNum,  nodePerElement ] = size ( data.elements );

%  Open the output file.
if data.iteration == -1
    output_filename = sprintf('%s.vtk', data.outputFile);
else
    output_filename = sprintf('%s%04d.vtk', data.outputFile,data.iteration);
end

if ( isempty ( output_filename ) )
    output_filename = 'solution.vtk';
end

outputSolution = fopen ( output_filename, 'w' );

%%  Transpose or otherwise modify the data.
coorz = zeros ( 3, nodeNum );
if size(data.vertices,2) == 3
    coorz(1:3,1:nodeNum) = (data.vertices(1:nodeNum,1:3))';
else
    coorz(1:2,1:nodeNum) = (data.vertices(1:nodeNum,1:2))';
end

elementNode = zeros ( nodePerElement, elementNum );
elementNode(1:nodePerElement,1:elementNum) = ...
    data.elements(1:elementNum,1:nodePerElement)' - 1;

%%  Write the data.
vtk_headerfile ( outputSolution, data.title, nodeNum, elementNum, ...
    nodePerElement, coorz, elementNode, output_filename);
for i=1:size(data.variableName,2)
    vtk_write ( outputSolution, nodeNum, data.variableName{i}, data.variableType{i}, data.variableData{i}, output_filename);
end
fclose ( outputSolution );

timewrite = toc(timewrite);
fprintf ( 1, '\n' );
fprintf ( 1, '  The data was written to "%s" in %1.3f seconds\n', output_filename, timewrite);

return
end

function vtk_headerfile ( outputSolution, title, nodeNum, elementNum, ...
    nodePerElement, coorz, elementNode, filename)
%*****************************************************************************
%
%% VTK_HEADERFILE write the header of the VTK file.
fprintf ( outputSolution, '# vtk DataFile Version 3.0\n' );
fprintf ( outputSolution, '%s\n', title );
fprintf ( outputSolution, 'BINARY\n' );
fprintf ( outputSolution, '\n' );
fprintf ( outputSolution, 'DATASET UNSTRUCTURED_GRID\n' );

fprintf ( outputSolution, 'POINTS %d float\n', nodeNum );
fwrite(outputSolution, coorz(1:3,:),'float','b');

cell_size = elementNum * ( nodePerElement + 1 );

fprintf ( outputSolution, '\n' );
fprintf ( outputSolution, 'CELLS  %d  %d\n', elementNum, cell_size );

fwrite(outputSolution, [nodePerElement*ones(1,elementNum); elementNode],'int','b');

fprintf ( outputSolution, '\n' );
fprintf ( outputSolution, 'CELL_TYPES %d\n', elementNum );

fwrite(outputSolution, [5*ones(elementNum,1)]','int','b');

fprintf ( outputSolution, '\n' );
fprintf ( outputSolution, 'POINT_DATA %d\n', nodeNum );
return
end

function vtk_write ( outputSolution, nodeNum, variableName, variableType, variableData, filename)
%*****************************************************************************80
%
%% VTK_WRITE vector or scalar data to a VTK file.
if(variableType=='SCALARS')
    fprintf ( outputSolution, 'SCALARS %s float\n', variableName );
    fprintf ( outputSolution, 'LOOKUP_TABLE default\n' );    
    fwrite(outputSolution, variableData(1,:),'float','b');
    fprintf ( outputSolution, '\n' );

else
    fprintf ( outputSolution, 'VECTORS %s float\n',variableName );
    fwrite(outputSolution, variableData(1:3,:),'float','b');
    fprintf ( outputSolution, '\n' );
    return
end
end
