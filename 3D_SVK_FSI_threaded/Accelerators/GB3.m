 % Rank-k Generalized Broyden Update
classdef GB3 < handle
    properties (SetAccess=private)
        r = []; 
        xt = [];
        Vprev = [];
        Wprev = []; 
        Tprev = []; 
        Jhat_prev = []; 
        V = [];
        W = [];        
        T  = [];
        B = [];        
        Z  = [];
        n  = 0;
        file;
        count = zeros(1,3);
    end
    
    properties (SetAccess=immutable)
        small;
        reuse;
        problemString;
    end
    
%% Function initialization and interface-data retention
    methods
        function M = GB3(small,reuse,count,problemString)
            M.small = small;      % Try range of regularization schemes
            M.reuse = reuse;      % 0 for gen. broyden
            
            filename6 = sprintf('Results/%s_filter.txt',problemString);
            M.file = fopen(filename6,'w');
            M.count = count;
        end
        
        function add(M,r,xt)
            if (~isempty(M.r))
                M.V = [r - M.r M.V];     %| Update residual difference matrix
                M.W = [xt - M.xt M.W];   %| Update solution difference matrix
                M.T = [M.n M.T];         %| Update info time identifier 
                
                % basic info filtering
                [M.V, M.W, M.T, ~, ~, ~, ~, ~,M.count] = filt_QR1(M.V, M.W, M.T, M.Vprev, M.Wprev, M.Tprev, M.small,M.count);                 

%                 [Q,R] = qr(M.V,0);
%                 for i = 1:length(Q(1,:))
%                     Z(i) = R\Q(:,1)';
%                 end
                if M.n > 1
                    temp = M.V'*M.V; 
                    Z = LUinverse(temp)*M.V';  

                    M.Z = Z;
                    M.B = (M.W - M.Jhat_prev*M.V);
                end
                               
            end            
            M.r = r;                     %| past residual update
            M.xt = xt;                   %| past intermediary  solutio            
        end
            
%% Secant Jacobian computation        
        function dx = predict(M,r) 
           
            if ((~isempty(M.V))||(~isempty(M.Jhat_prev)))
                if M.n == 1
                    [Q,R] = qr(M.V,0);
                    c = R\(Q'*r);
                    dx = M.W*c;
                else 
                    if (isempty(M.V))
                        dx = M.Jhat_prev*r;
                    else
                        Jup = M.Z*r;
                        dx = M.Jhat_prev*r + M.B*Jup;     %| Update approx. interface Jacobian             
                    end
                end
            else
                if (isempty(M.xt))
                    error('Model cannot predict without any data');
                else
                    dx = zeros(size(M.yr));
                end
            end

        end

%% Function maintenance and monitoring
         function clear(M)
          
          if M.n == 2
                temp = M.V'*M.V; 
                Z = LUinverse(temp)*M.V';            
                M.Jhat_prev = M.W*Z;           
           end
            
           if M.n > 2  
%                 [Q,R] = qr(M.V,0);
%                 for i = 1:length(Q(1,:))
%                     Z(i) = R\Q(:,1)';
%                  end   
                 temp = M.V'*M.V; 
                 Z = LUinverse(temp)*M.V';  
                 Jup = (M.W - M.Jhat_prev*M.V)*Z;
             
                 M.Jhat_prev = M.Jhat_prev + Jup;    %| update new approx. Jacobian
           end
                          
            % clear coupling iteration matrices
            M.r  = [];      
            M.xt = [];
            M.V  = [];
            M.W  = [];
            M.T  = [];
            M.count = zeros(1,3);
        end

        function increase_time(M)
            M.count(:,3) = size(M.V,2);
            format = '%4i %4i %4i %4i\n';           %| data format         
            data = [M.n M.count(:,1) M.count(:,2) M.count(:,3)];           %| time-step, coupling iteration, residual 
            fprintf(M.file,format,data);    %| write results to file
                
            M.n = M.n + 1;      %| increase time step
            M.clear();          %| update retained info
        end
        
       function closefile(M)
            fclose(M.file);
       end
                
        function LS = ready(M)
            LS =~isempty([M.V M.Vprev]);   %| data must be available for Jacobian approx.
        end
    end
end