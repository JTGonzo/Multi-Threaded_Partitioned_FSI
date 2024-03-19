 % Rank-k Generalized Broyden Update
classdef GB < handle
    properties (SetAccess=private)
        r = []; 
        xt = [];
        Jhat = [];
        Jhat_prev = []; 
        Vprev = [];
        Wprev = []; 
        Tprev = [];
        Vtemp = []; 
        Wtemp = []; 
        V = [];
        W = [];        
        T  = [];
        n  = 0;
        file;
        count = zeros(1,3);
    end
    
    properties (SetAccess=immutable)
        small;
        reuse;
        filter;
        problemString;
    end
    
%% Function initialization and interface-data retention
    methods
        function M = GB(small, reuse, filter, count,problemString)
            M.small = small;      % Try range of regularization schemes
            M.reuse = reuse;      % 0 for gen. broyden
            M.filter = filter;            %| what type of filter to use

            filename6 = sprintf('Results/%s_filter.txt',problemString);
            M.file = fopen(filename6,'w');
            M.count = count;
        end
        
        function add(M,r,xt)
            if (~isempty(M.r))
                M.V = [r - M.r M.V];     %| Update residual difference matrix
                M.W = [xt - M.xt M.W];   %| Update solution difference matrix
                M.T = [M.n M.T];         %| Update info time identifier 
                
                if M.n > 1
                    if M.filter == 3
                        den = M.Vtemp'*M.Vtemp;
                        Z = LUinverse(den)*M.Vtemp';  
                        Jup = (M.Wtemp - M.Jhat_prev*M.Vtemp)*Z;
                        M.Jhat = M.Jhat_prev + Jup;                         
                    else
                        den = M.V'*M.V;
                        Z = LUinverse(den)*M.V';  
                        Jup = (M.W - M.Jhat_prev*M.V)*Z;
                        M.Jhat = M.Jhat_prev + Jup;     %| Update approx. interface Jacobian
                    end
                end
            end            
            M.r = r;                     %| past residual update
            M.xt = xt;                   %| past intermediary  solutio            
        end
            
%% Secant Jacobian computation        
        function dx = predict(M,r) 
            
            %% Filter the retained information
           if M.filter ~= 0
               if ~isempty(M.V)
                    if M.filter == 1    %| QR1 -Filter
                        [M.V, M.W, M.T, ~, ~, ~, Q, R, M.count] = filt_QR1(M.V, M.W, M.T, M.Vprev, M.Wprev, M.Tprev, M.small, M.count);
                    elseif M.filter == 2 %| QR2 -Filter
                        vt = [M.V M.Vprev];
                        if length(vt(1,:))>= 2
                            [M.V, M.W, M.T, ~, ~, ~, Q, R, M.count] = filt_NM(M.V, M.W, M.T, M.Vprev, M.Wprev, M.Tprev, M.small, M.count);                   
                        else                     
                            vt = [M.V M.Vprev]; 
                            [Q,R] = qr(vt,0);
                        end
                    else     %| POD -Filter
                       [M.V, M.W, M.T, ~, ~, ~, M.Vtemp, M.Wtemp, Q, R, M.count, xmodes] = filt_POD(M.V, M.W, M.T, M.Vprev, M.Wprev, M.Tprev, M.small, M.count, M.n);
                       M.count(:,1) = xmodes;
                    end
               end
           else
                vt = [M.V M.Vprev]; 
               [Q,R] = qr(vt,0);
           end
            
            %% Solve the least-squares porblem and compute approx. Jacobian
            if ((~isempty(M.V))|| (~isempty(M.Jhat)))
                if M.n == 1
                    c = R\(Q'*r);
                        if M.filter == 3
                            dx = M.Wtemp*c;
                        else
                            dx = M.W*c;
                        end                
                else    
                    dx = M.Jhat*r;        %| Quasi-newton increment
                end
            else
                if (isempty(M.xt))
                    error('Model cannot predict without any data');
                else
                    dx = zeros(size(M.xt));
                end
            end
        end

%% Function maintenance and monitoring
         function clear(M)
            if M.n == 2
                if M.filter == 3
                    temp = M.Vtemp'*M.Vtemp; 
                    Z = LUinverse(temp)*M.Vtemp';            
                    M.Jhat = M.Wtemp*Z;
                else
                    temp = M.V'*M.V; 
                    Z = LUinverse(temp)*M.V';            
                    M.Jhat = M.W*Z;
                end
            end
            
            M.Jhat_prev = M.Jhat;        %| save past time-step approx. Jacobian
            
            if M.reuse == 0              %| no past information is kept
                 M.Vprev = [ ];
                 M.Wprev = [ ];
                 M.Tprev = [ ];
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
            LS =~isempty(M.Jhat_prev);   %| data must be available for Jacobian approx.
        end
    end
end