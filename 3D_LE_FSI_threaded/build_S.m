function [ MESH, FE_SPACE] = build_S( dim, elements, vertices, boundaries, fem, quad_order, DATA, ...
    model, rings)

%% Load domain specific MESH and element data 
MESH.dim         = dim;
MESH.fem         = fem;
MESH.vertices    = vertices;
MESH.boundaries  = boundaries;
MESH.elements    = elements;
MESH.numVertices = size(vertices,2);
MESH.quad_order  = quad_order;

if nargin > 8
    MESH.rings    = rings;
end

%% Build higher order mesh  
%begin P1 to P2 convert
    [~,nov]     =  size(vertices);
    [~,noe]     =  size(elements);
    nside       =  nov;

% create space for new nodes in the element connectivity array      
    elements = [elements(1:4,:); zeros(6,noe); elements(5,:)];

% compute adjacency of verticies
    noe = size(elements,2);
    nov = size(vertices,2);

    nln = 4;        % original number of vertices per triangular element

    nln2 = nln^2;

    [X,Y] = meshgrid(1:nln,1:nln);

    rr = X(:)';
    tt = Y(:)';
    cc = ones(1,nln2);

% Global array for all element node ids with extra spaces to be filled      
    row  = zeros(nln2*noe,1);
    col  = row;
    coef = row;

    iii  = 1:nln2;

% initializing element array with current element node ids      
    for i = 1 : noe

        p = elements(1:nln,i);

        row(iii)  = p(rr);
        col(iii)  = p(tt);
        coef(iii) = cc;

        iii = iii + nln2;

    end

% assemble sparse global eleemnt/node id matrix  
    a = sparse(row,col,coef,nov,nov); 
    [ii,jj,vv] = find(a);
    a = sparse(ii, jj, vv*0 - 1, nov, nov);

%identify the adjacent node pairs within each element     
    for ie = 1:noe
          i = elements(1,ie);%0
          j = elements(2,ie);%1
          k = elements(3,ie);%2
          h = elements(4,ie);%3

          % 4, mid point of 0-1
          l1 = a(i,j);
          if l1 == -1
              % give new node ID for added element nodes 
                nside = nside + 1;
                a(i,j) = nside;
                a(j,i) = nside;
                
                % assign the new node ID to the elements connectivity array
                elements(5,ie) = nside;
                
                % compute the spatial coordinate of the new added node  
                vertices(1:3,nside) = (vertices(1:3,i)+vertices(1:3,j))*0.5;
          else
                elements(5,ie) = l1;
          end

% repeat for other prismatic edges 
          % 5, mid point of 1-2
          l2 = a(j,k);
          if l2 == -1
                nside = nside + 1;
                a(j,k) = nside;
                a(k,j) = nside;
                elements(6,ie) = nside;
                vertices(1:3,nside) = (vertices(1:3,j)+vertices(1:3,k))*0.5;
          else
                elements(6,ie) = l2;
          end

          % 6, mid point of 2-0
          l3 = a(k,i);
          if l3 == -1
                nside = nside + 1;
                a(k,i) = nside;
                a(i,k) = nside;
                elements(7,ie) = nside;
                vertices(1:3,nside) = (vertices(1:3,k)+vertices(1:3,i))*0.5;
          else
                elements(7,ie) = l3;
          end

          % 7, mid point of 0-3
          l4 = a(i,h);
          if l4 == -1
                nside = nside + 1;
                a(h,i) = nside;
                a(i,h) = nside;
                elements(8,ie) = nside;
                vertices(1:3,nside) = (vertices(1:3,h)+vertices(1:3,i))*0.5;
          else
                elements(8,ie) = l4;
          end

          % 8, mid point of 1-3
          l5 = a(j,h);
          if l5 == -1
                nside = nside + 1;
                a(h,j) = nside;
                a(j,h) = nside;
                elements(9,ie) = nside;
                vertices(1:3,nside) = (vertices(1:3,h)+vertices(1:3,j))*0.5;
          else
                elements(9,ie) = l5;
          end

          % 9, mid point of 2-3
          l6 = a(k,h);
          if l6 == -1
                nside = nside + 1;
                a(h,k) = nside;
                a(k,h) = nside;
                elements(10,ie) = nside;
                vertices(1:3,nside) = (vertices(1:3,k)+vertices(1:3,h))*0.5;
          else
                elements(10,ie) = l6;
          end


    end

% add the new node IDs to the boundary elements connectivity arrays    
    [n,nside]=size(boundaries);
    for i = 1 : nside
          boundaries(4,i) = a(boundaries(1,i),boundaries(2,i));
          boundaries(5,i) = a(boundaries(2,i),boundaries(3,i));
          boundaries(6,i) = a(boundaries(3,i),boundaries(1,i));
    end

% add 3rd node index for ring conditions 
    [~,nrings]=size(rings);
    for i = 1 : nrings
        rings(3,i) = a(rings(1,i),rings(2,i));
    end

% end P1 to P2 element conversion and redeifine domain mesh properties
    MESH.elements = elements;
    MESH.nodes = vertices;
    MESH.boundaries = boundaries;
    MESH.rings = rings;
        
%% Update Mesh data with geometrical maps
% define element properties for 3D P2 elements
MESH.numElemDof = 10;       % P1 = 4
MESH.numBoundaryDof = 6;    % P1 = 3
MESH.numRingsDof = 3;       % P1 = 2
MESH.numNodes  = size(MESH.nodes,2);
MESH.numElem   = size(MESH.elements,2);

% Compute geometrical map (ref to physical elements) information  
noe    = size(MESH.elements,2);
MESH.jac = zeros(1,noe);
MESH.invjac = zeros(noe, dim, dim);
MESH.h = zeros(1,noe);

% Corner point indices
a1 = MESH.elements(1,:); a2 = MESH.elements(2,:); a3 = MESH.elements(3,:); a4 = MESH.elements(4,:);
       
% Lengths of prism sides
s21x = MESH.vertices(1,a2)-MESH.vertices(1,a1);
s31x = MESH.vertices(1,a3)-MESH.vertices(1,a1);
s41x = MESH.vertices(1,a4)-MESH.vertices(1,a1);
s21y = MESH.vertices(2,a2)-MESH.vertices(2,a1);
s31y = MESH.vertices(2,a3)-MESH.vertices(2,a1);
s41y = MESH.vertices(2,a4)-MESH.vertices(2,a1);
s21z = MESH.vertices(3,a2)-MESH.vertices(3,a1);
s31z = MESH.vertices(3,a3)-MESH.vertices(3,a1);
s41z = MESH.vertices(3,a4)-MESH.vertices(3,a1);
pzkI   = s31y.*s41z-s41y.*s31z;
pzkII  = s31x.*s41z-s41x.*s31z;
pzkIII = s31x.*s41y-s41x.*s31y;

% Volumes (multiplied for 6)
volumes = s21x.*pzkI-s21y.*pzkII+s21z.*pzkIII;
        
uno_su_volumes = 1./volumes;

% Jacobian elements
dcdx = pzkI.*uno_su_volumes;
dcdy = -pzkII.*uno_su_volumes;
dcdz = pzkIII.*uno_su_volumes;
        
dedx = (s41y.*s21z-s21y.*s41z).*uno_su_volumes;
dedy = (s21x.*s41z-s41x.*s21z).*uno_su_volumes;
dedz = (s21y.*s41x-s21x.*s41y).*uno_su_volumes;
dtdx = (s21y.*s31z-s31y.*s21z).*uno_su_volumes;
dtdy = (s21z.*s31x-s21x.*s31z).*uno_su_volumes;
dtdz = (s21x.*s31y-s31x.*s21y).*uno_su_volumes;
        
s23x = MESH.vertices(1,a2)-MESH.vertices(1,a3);
s43x = MESH.vertices(1,a4)-MESH.vertices(1,a3);
s24x = MESH.vertices(1,a2)-MESH.vertices(1,a4);
s23y = MESH.vertices(2,a2)-MESH.vertices(2,a3);
s43y = MESH.vertices(2,a4)-MESH.vertices(2,a3);
s24y = MESH.vertices(2,a2)-MESH.vertices(2,a4);
s23z = MESH.vertices(3,a2)-MESH.vertices(3,a3);
s43z = MESH.vertices(3,a4)-MESH.vertices(3,a3);
s24z = MESH.vertices(3,a2)-MESH.vertices(3,a4);

MESH.h = max([s31x.^2+s31y.^2+s31z.^2;s21x.^2+s21y.^2+s21z.^2;...
            s41x.^2+s41y.^2+s41z.^2;s23x.^2+s23y.^2+s23z.^2;...
            s43x.^2+s43y.^2+s43z.^2;s24x.^2+s24y.^2+s24z.^2]);

% Determinant of the Jacobian inverse of the linear transformation      
MESH.invjac(:,1,1) = dcdx;
MESH.invjac(:,1,2) = dedx;
MESH.invjac(:,1,3) = dtdx;
        
MESH.invjac(:,2,1) = dcdy;
MESH.invjac(:,2,2) = dedy;
MESH.invjac(:,2,3) = dtdy;
        
MESH.invjac(:,3,1) = dcdz;
MESH.invjac(:,3,2) = dedz;
MESH.invjac(:,3,3) = dtdz;

% Determinant of the Jacobian of the linear transformation
MESH.jac = abs(volumes);

%% Compute quadrature nodes and weights on the reference element
% quadrature for 3D elements, quad_order = 5
ax = [0 1 0 0];
ay = [0 0 1 0];
az = [0 0 0 1];
volume = 1/6;

w  = [16/135 (2665+14*sqrt(15))/37800*ones(1,4) (2665-14*sqrt(15))/37800*ones(1,4) 20/378*ones(1,6)];
s1 = (7-sqrt(15))/34;
t1 = (13+3*sqrt(15))/34;
s2 = (7+sqrt(15))/34;
t2 = (13-3*sqrt(15))/34;
u  = (10-2*sqrt(15))/40;
v  = (10+2*sqrt(15))/40;

x  = [1/4 s1 s1 s1 t1 s2 s2 s2 t2 u u v v u v];% Nodi di Gauss per il tetraedro di riferimento
y  = [1/4 s1 s1 t1 s1 s2 s2 t2 s2 u v v u v u];
z  = [1/4 s1 t1 s1 s1 s2 t2 s2 s2 v v u u u v];

w = w*volume;
quad_nodes = [x;y;z];
%% Evaluate P1 geometrical mapping basis functions in the quad points
% define finite element basis functions for 3D P1 elements
phi   = [];

x = quad_nodes(1,:);
y = quad_nodes(2,:);
z = quad_nodes(3,:);

phi(1,:) = 1-x-y-z;
phi(2,:) = x;
phi(3,:) = y;
phi(4,:) = z;
        
MESH.chi =  phi;

%% Generate mesh normals
[MESH.Normal_Faces] = ComputeSurfaceNormals3D(MESH.boundaries(1:3,:),MESH.vertices(1:3,:), elements(1:4,:));

%% Update MESH data with BC information
if nargin >= 7 && ~isempty(DATA)   
    MESH_parser = inputParser;
    MESH_parser.KeepUnmatched = true;

    addParameter(MESH_parser,'nodes',0);
    addParameter(MESH_parser,'boundaries',0);
    addParameter(MESH_parser,'elements',0);
    addParameter(MESH_parser,'numNodes',0);
    addParameter(MESH_parser,'dim',0);
    
    parse(MESH_parser,MESH);
    bc_flag_row = 12; %connectivity row identifies what physcial face the node belongs to (for 2D = 5)
    
    MESH.Dirichlet_dof = [];
    MESH.internal_dof  = [];
        
        for d = 1 : MESH.dim
            % initialize boundary flags per type, face and dimension
            type_Dirichlet = DATA.flag_dirichlet{d};
            type_Neumann   = DATA.flag_neumann{d};
            type_Pressure  = DATA.flag_pressure{d};
            type_Robin     = DATA.flag_robin{d};
            type_clamp_points = DATA.flag_clamp_points{d};
            
            if isempty(type_Dirichlet) && isempty(type_Neumann) && isempty(type_Pressure) && isempty(type_Robin) && isempty(type_clamp_points) 
                error(['No boundary conditions are imposed on component ', num2str(d)]);
            end
            
            %% Find Dirichlet dofs (if any)
            % Identifies Dirichlet nodes/indicies of the domain
            nDir                    = length(type_Dirichlet);
            Dirichlet_side          = [];

            for kk = 1 : nDir
                this_Dirichlet_side     = find(MESH.boundaries(bc_flag_row,:) == type_Dirichlet(kk));
                this_Dirichlet_dof      = MESH.boundaries(1:MESH.numBoundaryDof, unique( this_Dirichlet_side ) );
                MESH.DiriDof_CompFlag{d,kk}  = unique(this_Dirichlet_dof(:));
                Dirichlet_side          = [Dirichlet_side, this_Dirichlet_side];
            end
            Dirichlet_side             = unique(Dirichlet_side);
            Dirichlet_dof              = MESH.boundaries(1:MESH.numBoundaryDof,Dirichlet_side);
            Dirichlet_dof              = unique( Dirichlet_dof(:) );
            
            % locations where 2 dirichlet boundary conditions intersect
            dir_ringDofs = type_clamp_points;
            MESH.clamp_points{d} = dir_ringDofs;
                              
            MESH.Dirichlet_dof_c{d}    = unique([Dirichlet_dof; dir_ringDofs]);
            MESH.internal_dof_c{d}     = setdiff([1:MESH.numNodes]',MESH.Dirichlet_dof_c{d});
            
            % Dirichlet and internal indicies in terms of global solution space
            MESH.Dirichlet_dof = [MESH.Dirichlet_dof;  (d-1)*MESH.numNodes+MESH.Dirichlet_dof_c{d}];
            MESH.internal_dof  = [MESH.internal_dof; (d-1)*MESH.numNodes+MESH.internal_dof_c{d}];            
                        
            %% Find and apply Neumann boundaries (if any)
            if ~isempty(type_Neumann)
                nNeu         = length(type_Neumann);
                Neumann_side = [];
                for k = 1 : nNeu
                    Neumann_side = [Neumann_side,find(MESH.boundaries(bc_flag_row,:) == type_Neumann(k))];
                end
                MESH.Neumann_side{d} = unique(Neumann_side);
            else
                MESH.Neumann_side{d} = [];
            end
            
            %% Find Robin boundaries (if any)
            MESH.Robin_side{d} = [];
            
            %% Find and apply Pressure boundaries (if any)
            if ~isempty(type_Pressure)
                nPrex       = length(type_Pressure);
                Pressure_side = [];
                for kk = 1 : nPrex
                    this_Pressure_side = find(MESH.boundaries(bc_flag_row,:) == type_Pressure(kk));
                    MESH.Pressure_side_CompFlag{d,kk} = unique(this_Pressure_side);
                    Pressure_side = [Pressure_side, this_Pressure_side];
                end
                MESH.Pressure_side{d} = unique(Pressure_side);
            else
                MESH.Pressure_side{d} = [];
            end
        end
        
            MESH.DirichletNormal_dof = [];
            MESH.DirichletNormal_N   = [];
            MESH.DirichletNormal_R   = speye( MESH.dim*MESH.numNodes, MESH.dim*MESH.numNodes );

    MESH.bc_flag_row = bc_flag_row;
end

%% Define the finite element properties of the fluid variable spaces
% for the pressure solution elements (P1 elements)
    numElemDof = 10;
    numBoundaryDof  = 6;
    numComponents = dim;
    quad_weights = w;

% quadrature properties and weights 
    FE_SPACE.dim              = MESH.dim;
    FE_SPACE.fem              = fem;
    FE_SPACE.numComponents    = numComponents;
    FE_SPACE.numElemDof       = numElemDof;
    FE_SPACE.numBoundaryDof   = numBoundaryDof;

% variable solution space size   
    FE_SPACE.numDof           = numComponents * MESH.numNodes;
    FE_SPACE.numDofScalar     = MESH.numNodes;

% store quadrature nodes and weights on the reference element
    FE_SPACE.quad_order    = quad_order;
    FE_SPACE.quad_nodes    = quad_nodes;
    FE_SPACE.quad_weights  = quad_weights;
    FE_SPACE.numQuadNodes  = length(FE_SPACE.quad_nodes);

% finite element basis functions for P2 displacement elements     
    phi   = [];
    dphi   = [];
    dphix = [];
    dphiy = [];
    dphiz = [];

    phi(1,:) = (1-x-y-z).*(1-2*x-2*y-2*z);
    phi(2,:) = x.*(2*x-1);
    phi(3,:) = y.*(2*y-1);
    phi(4,:) = z.*(2*z-1);
    phi(5,:) = 4*x.*(1-x-y-z);
    phi(6,:) = 4*x.*y;
    phi(7,:) = 4*y.*(1-x-y-z);
    phi(8,:) = 4*z.*(1-x-y-z);
    phi(9,:) = 4*x.*z;
    phi(10,:) = 4*y.*z;
    
    dphix(1,:) = -3+4*x+4*y+4*z;
    dphix(2,:) =  4*x-1;
    dphix(3,:) =  0+0.*x;
    dphix(4,:) =  0+0.*x;
    dphix(5,:) = 4-8*x-4*y-4*z;
    dphix(6,:) =  4*y;
    dphix(7,:) =  -4*y;
    dphix(8,:) =  -4*z;
    dphix(9,:) =  4*z;
    dphix(10,:) =  0+0.*x;
    
    dphiy(1,:) = -3+4*x+4*y+4*z;
    dphiy(2,:) =  0+0.*x;
    dphiy(3,:) =  -1+4*y;
    dphiy(4,:) =  0+0.*x;
    dphiy(5,:) = -4*x;
    dphiy(6,:) =  4*x;
    dphiy(7,:) =  4-4*x-8*y-4*z;
    dphiy(8,:) =  -4.*z;
    dphiy(9,:) =  0+0.*x;
    dphiy(10,:) =  4*z;
    
    dphiz(1,:) = -3+4*x+4*y+4*z;   
    dphiz(2,:) =  0+0.*x;
    dphiz(3,:) =  0+0.*x;
    dphiz(4,:) =  -1+4*z;
    dphiz(5,:) = -4*x;
    dphiz(6,:) =  0+0.*x;
    dphiz(7,:) =  -4*y;
    dphiz(8,:) =  4-4*x-4*y-8*z;
    dphiz(9,:) =  4*x;
    dphiz(10,:) =  4*y;
    
    dphi(:,:,1) = dphix;
    dphi(:,:,2) = dphiy;
    dphi(:,:,3) = dphiz;

% add basis functions to element variable space     
    FE_SPACE.phi = phi;
    FE_SPACE.dphi_ref = dphi;  
end
