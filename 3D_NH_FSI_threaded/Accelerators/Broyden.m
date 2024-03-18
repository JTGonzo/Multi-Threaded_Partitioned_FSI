% Rank-2 Broyden
classdef Broyden < handle
%% Set in-function fixed (internal) and Jhatamic (interface) variable values 
    properties (SetAccess=private)
        r = []; 
        xt = [];
        V  = [];
        W  = [];
        T  = [];
        Jhat = [];
        n  = 0;
    end

    properties (SetAccess = immutable)
        reuse;
    end    
%% Function initialization and interface-data retention
    methods
        function M = Broyden(reuse)
            M.reuse = reuse;              %| number of time steps to reuse
        end
        
        function add(M,r,xt)
            if (~isempty(M.r))
                M.V = [r - M.r M.V];     %| Update residual difference matrix
                M.W = [xt - M.xt M.W];   %| Update solution difference matrix
                M.T = [M.n M.T];         %| Update info time identifier 
                
                if M.n > 1 
                    num = M.W(:,1) - M.Jhat*M.V(:,1);   
                    den = M.V(:,1)'*M.V(:,1);
                    Jup = (num/den)*M.V(:,1)';          %| least-change secant update
                    M.Jhat = M.Jhat + Jup;              %| Update approx. interface Jacobian           
                end
            end
            M.r = r;                     %| past residual update
            M.xt = xt;                   %| past intermediary  solution                
        end

%% Secant Jacobian computation
        function dx = predict(M,r)
            
            if ((~isempty(M.V))|| (~isempty(M.Jhat)))
                if M.n == 1
                    [Q,R] = qr(M.V,0);
                    c = R\(Q'*r);
                    dx = M.W*c;
                else    
                    dx = M.Jhat*r;        %| Quasi-newton increment
                end
            else
                if (isempty(M.xt))                  %| Min. data from 2 preceeding iterations                                            
                    error('Model cannot predict without any data');
                else
                    dx = zeros(size(M.xt));
                end
            end
        end

 %% Function maintenance 
       function clear(M)
           if M.n == 2
                temp = M.V'*M.V; 
                Z = LUinverse(temp)*M.V';            
                M.Jhat = M.W*Z;           
           end
            
            M.r = [];                 %| clear coupling iteration matrices
            M.xt = [];
            M.W  = [];
            M.V  = [];          
            M.T  = [];           
        end
        
        function increase_time(M)       
            M.n = M.n + 1;            %| increase time step
            M.clear();                %| clear retained info
        end
        
        function closefile(M)
        end 
       
        function LS = ready(M)
            LS =~isempty(M.Jhat);   %| data must be available for Jacobian approx.
        end
    end
end