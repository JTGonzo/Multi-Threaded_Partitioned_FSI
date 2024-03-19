function [ FE_SPACE ] = build_G( MESH, fem, numComponents, quad_order )
%% Compute quadrature nodes and weights on the reference element    
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

%% Define the finite element properties of the geometry variable spaces
% for the geometry displacement solution elements (P2 elements)
    numElemDof = 6;
    numBoundaryDof  = 3;
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

% add basis functions to element variable space       
    FE_SPACE.phi = phi;
    FE_SPACE.dphi_ref = dphi;   
end
