% NP weighs and the nodes of the corresponding Gauss-Legendre 
% quadrature formula in the reference interval (a,b) and corresponding
% boundary basis functions
function [x,w,phi] = bound_quad(np,a,b)
%% Initialize the variable space
n = np - 1;
phi = [];

x   = zeros(np,1);
w   = zeros(np,1);
jac = zeros(np);

%% Compute the GL quadrature nodes and weights
% done on the reference interval (-1,1)
k   = [1:n];
v   = (k)./(sqrt(4*(k.^2)-1));
jac = jac + diag(v,1) + diag(v,-1);

[w1,x] = eig(jac);
norm2  = sqrt(diag(w1'*w1));    
w1     = (2*w1(1,:)'.^2)./norm2;  
x      = diag(x);		    

[x,ip] = sort(x);

for i = 1:np
    w(i) = w1(ip(i));
end

%% Transform the GL nodes and weights to general interval (a,b)
if nargin == 3
    bma = (b-a)*.5;
    bpa = (b+a)*.5;
    x   = bma*x + bpa;
    w   = w*bma;
end

% GL quadrature and weights
x = x';
w = w';

%% Compute the basis functions for the surface quadrature nodes
% Only keep the x direction quadrature 
X = [x; 0*x];

x = X(1,:);
y = X(2,:);

% basis functions
phi(1,:) = (1-x-y).*(1-2*x-2*y);
phi(2,:) = x.*(-1+2*x);
phi(3,:) = y.*(-1+2*y);
phi(4,:) = 4*x.*(1-x-y);
phi(5,:) = 4*x.*y;
phi(6,:) = 4*y.*(1-x-y);

% keep the planar functions only
BDOF = [1 2 4];

phi   = phi(BDOF,:);

end