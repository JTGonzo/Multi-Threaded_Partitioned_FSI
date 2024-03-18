function [ FE_SPACE ] = build_G( MESH, fem, numComponents, quad_order )
%% Compute quadrature nodes and weights on the reference element   
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

%% Define the finite element properties of the geometry variable spaces
% for the geometry displacement solution elements (P2 elements)
    numElemDof = 10;
    numBoundaryDof  = 6;
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
