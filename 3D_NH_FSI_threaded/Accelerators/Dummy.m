% Gauss Seidel Iterations
classdef Dummy < handle
%% Set in-function fixed (internal) and dynamic (interface) variable values 
    properties (SetAccess=private)
        xr = []; 
        V = [];
        T  = [];
        n  = 0;
    end

 %% Secant Jacobian computation
    methods
        function M = Dummy()
        end
        
        function add(M,xr)
            M.xr = xr;             %| initial difference input
            M.V = [M.xr M.V];      %| Update residual matrix
            M.T = [M.n M.T];       %| Update time catalog
        end
 
%% No Jacobian 
       function [dyn]= predict(M,r) 
                dyn = ones(size(r),1);                             
       end

%% Function maintenance and monitoring
        function clear(M)
            M.xr = [];                            % clear coupling iteration matrices
            M.V = [];          
            M.T = [];
        end
        
        function increase_time(M)       % increase time step
            M.n = M.n + 1;
            M.clear();
        end

       function closefile(M)
       end
        
       function LS = ready(M)             
            LS = true;
        end
    end
end