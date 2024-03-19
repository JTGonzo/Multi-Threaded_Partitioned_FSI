% Anderson Acceleration
classdef AA_secant_select < handle
%% Set in-function fixed (internal) and dynamic (interface) variable values 
    properties (SetAccess=private)
        r = []; 
        xt = [];
        
        Vprev = [];
        Wprev = []; 
        Tprev = []; 
        T2prev = [];
        Kprev = []; 
        
        nRprev = [];
        nXprev = [];
        nVprev = []; 
        nWprev = [];

        V = [];
        W = [];        
        T  = [];
        T2  = [];
        K  = [];
        
        nR = [];
        nX = [];
        nW = [];
        nV = [];  
        
        n  = 0;
        k = 0;
        count = zeros(1,3);
        finished_primal = false;
        file;
    end
    
    properties (SetAccess = immutable)
        small;
        reuse;
        filter;
        adaptv;
        maxSec;
        problemString;
    end
    
%% Function initialization and interface-data retention
    methods
        function M = AA_secant_select(small, reuse, maxS, filter, adaptv, count, problemString)
            M.small = small;              %| smallest floating point number
            M.reuse = reuse;              %| number of time steps to reuse
            M.filter = filter;            %| what type of filter to use
            M.adaptv = adaptv;            %| what type of filter to use
            M.maxSec = maxS;
            
            filename6 = sprintf('Results/%s_filter.txt',problemString);
            M.file = fopen(filename6,'w');
            M.count = count;
        end
        
        function add(M,r,xt)
            M.k = M.k + 1;
            if M.n > 2
                if size(M.V,2) < 20
                    if (~isempty(M.r))
                        M.V = [r - M.r M.V];            %| Update residual difference matrix
                        M.W = [xt - M.xt M.W];          %| Update solution difference matrix
                        M.T = [M.n M.T];                %| Update info time identifier    
                        M.nW = [norm(M.W(:,1)) M.nW];   %| iterate difference vectorial property
                        M.nV = [norm(M.V(:,1)) M.nV];   %| residual difference vectorial property
                    else 
                        M.nW = [0 M.nW];
                        M.nV = [0 M.nV];                    
                    end  
                else
                     M.V = [r - M.r M.V]; 
                     M.W = [xt - M.xt M.W];   
                     M.T = [M.n M.T];     
                                      
                     M.V(:,end-6:end-2) = [];
                     M.W(:,end-6:end-2) = [];  
                     M.T(:,end-6:end-2) = [];
                     
                     M.T2(:,end-7:end-3) = [];
                     M.K(:,end-7:end-3)  = [];
                     M.nR(:,end-7:end-3) = [];
                     M.nX(:,end-7:end-3) = []; 
                     M.nW(:,end-7:end-3) = []; 
                     M.nV(:,end-7:end-3) = []; 
                end
            else   
                if (~isempty(M.r))
                    M.V = [r - M.r M.V];            
                    M.W = [xt - M.xt M.W];          
                    M.T = [M.n M.T];                 
                    M.nW = [norm(M.W(:,1)) M.nW];  
                    M.nV = [norm(M.V(:,1)) M.nV];   
                else 
                    M.nW = [0 M.nW];
                    M.nV = [0 M.nV];                    
                end  
            end    
            M.r = r;                     %| past residual update
            M.xt = xt;                   %| past intermediary  solution  
            
            M.T2 = [M.n M.T2];
            M.K = [M.k M.K];             %| Update info iteration identifier
            M.nR = [norm(r) M.nR];       %| residual norm storage
            M.nX = [norm(xt) M.nX];
        end

%% Secant Jacobian computation
        function [J] = predict(M, r)
            
            if M.n > 6
                [M.Vprev, M.Wprev, M.Tprev, M.T2prev, M.Kprev, M.nRprev, M.nXprev, M.nVprev, M.nWprev, M.V, M.W, ...
                    M.T, M.T2, M.K, M.nR, M.nX, M.nV, M.nW, Q, R, M.count, idx] = filt_QRSC(M, M.adaptv, M.small, M.count, M.k, r);              
            else               
                vt = [M.V M.Vprev];
                [Q,R] = qr(vt,0);
            end        
             
            %% Solve the least-squares porblem and compute approx. Jacobian
            if ((~isempty(M.V))||(~isempty(M.Vprev)))          

                c = R\(Q'*r);           % least-squares weights
                 
                wt = [M.W M.Wprev];                
                if M.n > 6
                    wt = wt(:,idx);
                end
                
                J = wt*c;               
            else
                if (isempty(M.xt))  %| Min. data from 2 preceeding iterations needed                                          
                    error('Model cannot predict without any data');                  
                else
                    J = zeros(size(M.xt));          
                end
            end       
        end
        
%% Function maintenance and monitoring
             function clear(M)
                 
                if M.reuse == 0                 %| no past information is kept
                     M.Vprev  = [ ];
                     M.Wprev  = [ ];
                     M.Tprev  = [ ];
                     M.T2prev = [ ];
                     M.Kprev  = [ ];
                     M.nRprev = [ ]; 
                     M.nXprev = [ ];
                     M.nVprev = [ ];                 
                     M.nWprev = [ ];  

                     M.count(:,3) = size(M.V,2);
                else 
                     M.Vprev = [M.V M.Vprev];   %| append coupling iteration matrices
                     M.Wprev = [M.W M.Wprev];
                     M.Tprev = [M.T M.Tprev];

                     M.T2prev = [M.T2 M.T2prev];
                     M.Kprev  = [M.K M.Kprev];
                     M.nRprev = [M.nR M.nRprev];          
                     M.nXprev = [M.nX M.nXprev];
                     M.nVprev = [M.nV M.nVprev];                 
                     M.nWprev = [M.nW M.nWprev];   

                     M.count(:,3) = size(M.Vprev,2);
                     
                     while size(M.Vprev,2) > M.maxSec
                        last1 = find(M.Tprev == min(M.Tprev));
                        last2 = find(M.T2prev == min(M.T2prev)); 

                        M.Vprev(:,last1) = [];      
                        M.Wprev(:,last1) = [];
                        M.Tprev(:,last1) = [];
                        
                        M.Kprev(:,last2) = [];                    
                        M.T2prev(:,last2) = [];                    
                        M.nRprev(:,last2) = [];
                        M.nXprev(:,last2) = [];
                        M.nVprev(:,last2) = [];
                        M.nWprev(:,last2) = [];
                        
                     end                     
                     
                     old = find(M.Tprev < M.n - M.reuse);   %| find outdated diff. information indicies
                     old2 = find(M.T2prev < M.n - M.reuse);   %| find outdated diff. information indicies

                     if (~isempty(old))                     %| remove outdated difference information   
                        M.Vprev(:,old) = [];      
                        M.Wprev(:,old) = [];
                        M.Tprev(:,old) = [];

                        M.Kprev(:,old2) = [];                    
                        M.T2prev(:,old2) = [];                    
                        M.nRprev(:,old2) = [];
                        M.nXprev(:,old2) = [];
                        M.nVprev(:,old2) = [];
                        M.nWprev(:,old2) = [];
                     end                 
                end
            
            if M.n > 1
                format = '%4i %4i %4i %4i\n';           %| data format         
                data = [M.n-1 M.count(:,1) M.count(:,2) M.count(:,3)];           %| time-step, coupling iteration, residual 
                fprintf(M.file,format,data);    %| write results to file
            end

            % clear coupling iteration matrices            
            M.r  = [];       
            M.xt = [];            
            M.V  = [];
            M.W  = [];
            M.T  = [];           
            M.K  = [];
            M.T2  = [];            
            M.nR  = [];
            M.nX  = [];
            M.nV  = [];
            M.nW  = [];
            
            M.count = zeros(1,3);
        end
        
        function increase_time(M)  
            M.k = 0;
            M.n = M.n + 1;                  %| increase time step
            M.clear();                      %| update retained info
        end
        
        function closefile(M)
            fclose(M.file);
        end
        
        function LS = ready(M)             
            LS =~ isempty([M.V M.Vprev]);   %| data must be available for Jacobian approx. 
        end
    end
end