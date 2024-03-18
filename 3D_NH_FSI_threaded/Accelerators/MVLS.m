% Generalized Broyden with Least-Squares Update
classdef MVLS < handle
    properties (SetAccess=private)
        n = 0;
        r = []; 
        xt = [];
        V = [];
        W = [];
        T = [];
        Wprev = [];
        Qprev = [];
        Rprev = [];
        Vtemp = []; 
        Wtemp = []; 
        Dummy = [];
        file;
        count = zeros(1,3);
    end
    
    properties (SetAccess=immutable)
        small;
        limit;
        filter;
        problemString;
    end

%% Function initialization and interface-data retention    
    methods
        function M = MVLS(small, limit, filter, count, problemString)
            M.small = small;         %| smallest floating point number
            M.limit = limit;         %| number of time steps to reuse
            M.filter = filter;    %| what type of filter to use
            
            filename6 = sprintf('Results/%s_filter.txt',problemString);
            M.file = fopen(filename6,'w');
            M.count = count;
        end
        
        function add(M,r,xt)
            if (~isempty(M.r))
                M.V = [r - M.r M.V];     %| Update residual difference matrix
                M.W = [xt - M.xt M.W];   %| Update solution difference matrix
                M.T = [M.n M.T];         %| Update info time identifier    
            end            
            M.r = r;                     %| past residual update
            M.xt = xt;                   %| past intermediary  solution             
        end

%% Secant Jacobian computation           
        function dx = predict(M,r)                       
            
            if ((isempty(M.V))&&(isempty(M.Wprev)))                                        
                    error('Model cannot predict without any data');
            end
            
            %% Current time-steps secant contributions        
            if (~isempty(M.V)) 
               if M.filter ~= 0
                    if M.filter == 1    %| QR1 -Filter
                        [M.V, M.W, M.T, ~, ~, ~, Q, R, M.count] = filt_QR1(M.V, M.W, M.T, M.Dummy, M.Dummy, M.Dummy, M.small, M.count);
                    elseif M.filter == 2 %| QR2 -Filter
                        vt = [M.V M.Dummy];
                        if length(vt(1,:))>= 2
                            [M.V, M.W, M.T, ~, ~, ~, Q, R, M.count] = filt_NM(M.V, M.W, M.T, M.Dummy, M.Dummy, M.Dummy,  M.small, M.count);                   
                        else                     
                            vt = [M.V M.Dummy]; 
                            [Q,R] = qr(vt,0);
                        end
                    else     %| POD -Filter
                       [M.V, M.W, M.T, ~, ~, ~, M.Vtemp, M.Wtemp, Q, R, M.count, xmodes] = filt_POD(M.V, M.W, M.T, M.Dummy, M.Dummy, M.Dummy, M.small, M.count, M.n);
                       M.count(:,1) = xmodes;
                    end
               else
                    vt = [M.V M.Dummy]; 
                   [Q,R] = qr(vt,0);
               end
                
                b = Q'*r;
                c = R\b;               %| least-squares weights 
                
                if M.filter == 3
                    dx = M.Wtemp*c;            %| approx. inverse Jacobian
                else
                    dx = M.W*c;
                end
                
                r = r - Q*b;            
            else
                dx = zeros(length(r),1);
            end
            
            %% Incorporate recursive Jacobian contributions
            if length(M.Wprev)>= 1             
                i = 0;
                if M.n-1 <= M.limit
                    kn = M.n-1;
                else
                    kn = M.limit;
                end

                while ((norm(r)>M.small) && (i<min(M.limit,length(M.Wprev))))
                    qq = M.Qprev{kn};
                    b = qq'*r;
                    c = M.Rprev{kn}\b;        %| past step leas-squares weights 
                    dx = dx + M.Wprev{kn}*c;  %| past time step contributions                
                    r = r - qq*b;             %| remove accounted residual info
                    kn = kn - 1;
                    i = i + 1;
                end
            end
        end

%% Function maintenance and monitoring        
        function clear(M)
           if M.n > 1 
               format = '%4i %4i %4i %4i\n';           %| data format         
               data = [M.n-1 M.count(:,1) M.count(:,2) M.count(:,3)];           %| time-step, coupling iteration, residual 
               fprintf(M.file,format,data);    %| write results to file
           end
           
            M.r  = [];      %} clear current coupling iteration matrices
            M.xt = [];
            M.V  = [];
            M.W  = [];
            M.T  = [];
            M.count = zeros(1,3);
        end

        function increase_time(M)
            M.n = M.n+1;
            
            if M.n > 1
                [Q,R] = qr(M.V,0);                     
                %| append coupling iteration matrices 
                if M.n-1 <= M.limit
                    M.Qprev{M.n-1} = Q;
                    M.Rprev{M.n-1} = R;
                    M.Wprev{M.n-1} = M.W;
                    
                   for i=1:M.n-1
                        M.count(:,3) = M.count(:,3) + size(M.Wprev{i},2);
                    end
                else
                %| remove outdated difference information
                    M.Qprev{1} = [];   
                    M.Rprev{1} = [];
                    M.Wprev{1} = [];
                    
                    for i=2:M.limit
                        M.Qprev{i-1} = M.Qprev{i};
                        M.Rprev{i-1} = M.Rprev{i};
                        M.Wprev{i-1} = M.Wprev{i};
                    end
                    
                    M.Qprev{M.limit} = Q;
                    M.Rprev{M.limit} = R;
                    M.Wprev{M.limit} = M.W;
                    
                    for i=1:M.limit
                        M.count(:,3) = M.count(:,3) + size(M.Wprev{i},2);
                    end
                end
            end     
            M.clear();
        end
        
        function closefile(M)
            fclose(M.file);
        end
       
        function LS = ready(M)
            LS =~isempty([M.V M.Wprev]); %| data must be available for Jacobian approx. 
        end
    end
end