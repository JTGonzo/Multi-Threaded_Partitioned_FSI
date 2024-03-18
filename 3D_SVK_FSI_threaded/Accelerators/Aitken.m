% Aitken Relaxation
classdef Aitken < handle
%% Set in-function fixed (internal) and dynamic (interface) variable values 
    properties (SetAccess=private)
        xr = []; 
        V = [];
        T  = [];
        n  = 0;
    end

%% Function initialization and interface-data retention
    methods
        function M = Aitken()
        end
        
        function add(M,xr)
            M.xr = xr;             %| initial difference input
            M.V = [M.xr M.V];      %| Update residual matrix
            M.T = [M.n M.T];       %| Update time catalog
        end
 
%% Basic Aitken computation 
       function [dyn]= predict(M,r) 
                dyn = dot(M.V(:,2)',(r-M.V(:,2)))/dot((r-M.V(:,2))',(r-M.V(:,2))); 
                % dynamic relaxation term                           
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
       
       function LS = ready(M)             % check to see if information is available
            LS = true;
        end
    end
end