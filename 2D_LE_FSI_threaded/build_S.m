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
    elements = [elements(1:3,:); zeros(3,noe); elements(4,:)];

% compute adjacency of verticies
    noe = size(elements,2);
    nov = size(vertices,2);

    nln = 3; % original number of vertices per triangular element
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
        i = elements(1,ie);
        j = elements(2,ie);
        k = elements(3,ie);

        l1 = a(i,j);
        if l1 == -1
           % give new node ID for added element nodes 
           nside = nside + 1;
           a(i,j) = nside;
           a(j,i) = nside;
           
           % assign the new node ID to the elements connectivity array
           elements(4,ie) = nside;
           
           % compute the spatial coordinate of the new added node  
           vertices(1,nside) = (vertices(1,i)+vertices(1,j))*0.5;
           vertices(2,nside) = (vertices(2,i)+vertices(2,j))*0.5;
       else 
           elements(4,ie) = l1;
       end

% repeat for other 2 triangle edges         
       l2 = a(j,k);
       if l2 == -1 
           nside = nside + 1;
           a(j,k) = nside;
           a(k,j) = nside;
           elements(5,ie) = nside;
           vertices(1,nside) = (vertices(1,j)+vertices(1,k))*0.5;
           vertices(2,nside) = (vertices(2,j)+vertices(2,k))*0.5;
       else 
           elements(5,ie) = l2;
       end

       l3 = a(k,i);
       if l3 == -1 
           nside = nside + 1;
           a(k,i) = nside;
           a(i,k) = nside;
           elements(6,ie) = nside;
           vertices(1,nside) = (vertices(1,k)+vertices(1,i))*0.5;
           vertices(2,nside) = (vertices(2,k)+vertices(2,i))*0.5;
       else
           elements(6,ie) = l3;
       end
    end

% add the new node IDs to the boundary elements connectivity arrays    
    [n,nside] = size(boundaries);
    for i = 1 : nside
        boundaries(3,i) = a(boundaries(1,i),boundaries(2,i));
    end

% end P1 to P2 element conversion and redeifine domain mesh properties
    MESH.elements = elements;
    MESH.nodes = vertices;
    MESH.boundaries = boundaries;
            
%% Update Mesh data with geometrical maps
% define element properties for 2D P2 elements
MESH.numElemDof = 6;      % P1 = 3
MESH.numBoundaryDof = 3;  % P1 = 2
MESH.numRingsDof = 1;     % P1 = 1
MESH.numNodes = size(MESH.nodes,2);
MESH.numElem = size(MESH.elements,2);

% Compute geometrical map (ref to physical elements) information  
noe    = size(MESH.elements,2);
MESH.jac = zeros(1,noe);
MESH.invjac = zeros(noe, dim, dim);
MESH.h = zeros(1,noe);

% Corner point indices
a1 = MESH.elements(1,:); a2 = MESH.elements(2,:); a3 = MESH.elements(3,:);
        
% Triangle sides
s13x = MESH.vertices(1,a1)-MESH.vertices(1,a3);
s31y = MESH.vertices(2,a3)-MESH.vertices(2,a1);
s32x = MESH.vertices(1,a3)-MESH.vertices(1,a2);
s23y = MESH.vertices(2,a2)-MESH.vertices(2,a3);
s21x = MESH.vertices(1,a2)-MESH.vertices(1,a1);
s21y = MESH.vertices(2,a2)-MESH.vertices(2,a1);
        
% Determinant of the Jacobian matrix with sign (two times the area with sign of T)
MESH.jac = s13x.*s23y-s31y.*s32x;
        
uno_su_detjac = 1./MESH.jac;
        
% Jacobian elements
dcdx = s31y.*uno_su_detjac;
dcdy = s13x.*uno_su_detjac;
dedx = -s21y.*uno_su_detjac;
dedy = s21x.*uno_su_detjac;
MESH.h = max([s13x.^2+s31y.^2;s32x.^2+s23y.^2;s21x.^2+s21y.^2]);
MESH.h = sqrt(MESH.h);

% Determinant of the Jacobian inverse of the linear transformation
MESH.invjac = zeros(size(elements,2), 2, 2);
                
MESH.invjac(:,1,1) = dcdx;
MESH.invjac(:,1,2) = dedx;
MESH.invjac(:,2,1) = dcdy;
MESH.invjac(:,2,2) = dedy;
              
% Determinant of the Jacobian of the linear transformation
MESH.jac = abs(MESH.jac); 

%% Compute quadrature nodes and weights on the reference element
% quadrature for 2D elements, quad_order = 4
ax = [0 1 0];
ay = [0 0 1];
area = 0.5;

w = [0.223381589678011.*ones(1,3) 0.109951743655322.*ones(1,3)];
p1 = 0.108103018168070; p2 = 0.445948490915965;
x = [p1*ax(1)+p2*ax(2)+p2*ax(3) p2*ax(1)+p1*ax(2)+p2*ax(3) p2*ax(1)+p2*ax(2)+p1*ax(3)];
y = [p1*ay(1)+p2*ay(2)+p2*ay(3) p2*ay(1)+p1*ay(2)+p2*ay(3) p2*ay(1)+p2*ay(2)+p1*ay(3)];

p1 = 0.816847572980459; p2 = 0.091576213509771;
x = [x, p1*ax(1)+p2*ax(2)+p2*ax(3) p2*ax(1)+p1*ax(2)+p2*ax(3) p2*ax(1)+p2*ax(2)+p1*ax(3)];
y = [y, p1*ay(1)+p2*ay(2)+p2*ay(3) p2*ay(1)+p1*ay(2)+p2*ay(3) p2*ay(1)+p2*ay(2)+p1*ay(3)];

w = w*area;
quad_nodes = [x;y];

%% Evaluate P1 geometrical mapping basis functions in the quad points
% define finite element basis functions for 2D P1 elements
phi   = [];
dphix = [];
dphiy = [];
        
x = quad_nodes(1,:);
y = quad_nodes(2,:);

phi(1,:) = 1-x-y;
phi(2,:) = x;
phi(3,:) = y;

dphix(1,:) = -1+0.*x;
dphix(2,:) =  1+0.*x;
dphix(3,:) =  0+0.*x;
dphiy(1,:) = -1+0.*x;
dphiy(2,:) =  0+0.*x;
dphiy(3,:) =  1+0.*x;               
dphi(:,:,1) = dphix;
dphi(:,:,2) = dphiy;                    
                
% Evaluate P1 geometrical mapping basis functions in the quad points
MESH.chi =  phi;

%% Generate mesh normals
if strcmp( model, 'CSM') || strcmp( model, 'CFD')

% begin mesh normal computation    
    boundaries = MESH.boundaries(1:2,:);
    vertices = MESH.vertices(1:2,:);
    elements = MESH.elements(1:3,:);
    
    [~,nside]    = size(boundaries);  % number of boundary elements

    normalf      = zeros(2,nside);

    elements_row = [elements(1,:) elements(2,:) elements(3,:)];
    noe          = size(elements,2);

    for iside = 1:nside
         % identify the indicies of each boundary element
          i1 = boundaries(1,iside);
          i2 = boundaries(2,iside);
          
          % boundary length vector
          l  = sqrt((vertices(1,i1)-vertices(1,i2))^2+(vertices(2,i1)-vertices(2,i2))^2);

          % coordinates of the vertices i1 and i2
          x1 = vertices(1,i1);
          x2 = vertices(1,i2);
          y1 = vertices(2,i1);
          y2 = vertices(2,i2);

          % normal of the edge
          n1 =   (y2-y1)/l;
          n2 = - (x2-x1)/l;

          % find the element to which the edge (i1,i2) belongs
          ie = find(elements_row == i1);
          
          for kk = 1:length(ie)  
                ie_tmp = ie(kk);
                if mod(ie_tmp,noe) == 0
                      ie_tmp = noe;
                else
                      ie_tmp = mod(ie_tmp,noe);
                end

                i3 = setdiff(elements(1:3,ie_tmp),[i1 i2]);
                if length(i3) == 1
                      break;
                end
          end

          % a rough check fo search logic
          if length(i3) ~=1
                disp('error')
          end

          % coordinates of the third nodes of the elements to which 
          % the edge (i1,i2) belongs 
          x3  = vertices(1,i3);
          y3  = vertices(2,i3);

          % check if the normal is outward or not
          if (n1*(x3 - x2) + n2*(y3-y2)) > 0
                n1 = - n1;
                n2 = - n2;
          end
          
          % define and assign normals of boundary faces
          normalf(:,iside) = [n1; n2];
          MESH.Normal_Faces = normalf;
    end
end

%% Update Mesh data with BC information
if nargin >= 7 && ~isempty(DATA)
    MESH_parser = inputParser;
    MESH_parser.KeepUnmatched = true;
    
    addParameter(MESH_parser,'nodes',0);
    addParameter(MESH_parser,'boundaries',0);
    addParameter(MESH_parser,'elements',0);
    addParameter(MESH_parser,'numNodes',0);
    addParameter(MESH_parser,'dim',0);
    
    parse(MESH_parser,MESH);
    bc_flag_row = 5; %connectivity row identifies what physcial face the node belongs to (for 3D = 12)
    
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
% for the displacement solution elements (P2 elements)
    numElemDof = 6;
    numBoundaryDof  = 3;
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

    phi(1,:) = (1-x-y).*(1-2*x-2*y);
    phi(2,:) = x.*(-1+2*x);
    phi(3,:) = y.*(-1+2*y);
    phi(4,:) = 4*x.*(1-x-y);
    phi(5,:) = 4*x.*y;
    phi(6,:) = 4*y.*(1-x-y);
    
    dphix(1,:) = -3+4*x+4*y;
    dphix(2,:) = -1+4*x;
    dphix(3,:) = 0 + 0.*x;
    dphix(4,:) = 4-8*x-4*y;
    dphix(5,:) = 4*y;
    dphix(6,:) = -4*y;
    
    dphiy(1,:) = -3+4*x+4*y;
    dphiy(2,:) = 0+0.*x;
    dphiy(3,:) = -1+4*y;
    dphiy(4,:) = -4*x;
    dphiy(5,:) = 4*x;
    dphiy(6,:) = 4-4*x-8*y;

    dphi(:,:,1) = dphix;
    dphi(:,:,2) = dphiy;    

% add basis functions to variable space    
    FE_SPACE.phi = phi;
    FE_SPACE.dphi_ref = dphi; 
end

