% Anderson Acceleration
classdef AA_more_filters < handle
%% Set in-function fixed (internal) and dynamic (interface) variable values 
    properties (SetAccess=private)
        r = []; 
        xt = [];
        Vprev = [];
        Wprev = []; 
        Tprev = []; 
        V = [];
        W = [];        
        T  = [];
        n  = 0;
        k = 0;
        count = zeros(1,3);
        finished_primal = false;
        file;
        file2;
        file3;
    end
    
    properties (SetAccess = immutable)
        small;
        reuse;
        filter;
        problemString;
    end
    
%% Function initialization and interface-data retention
    methods
        function M = AA_more_filters(small, reuse, filter, count, problemString)
            M.small = small;              %| smallest floating point number
            M.reuse = reuse;              %| number of time steps to reuse
            M.filter = filter;            %| what type of filter to use
            
            filename6 = sprintf('Results/%s_filter.txt',problemString);
            M.file = fopen(filename6,'w');
            M.count = count;
            
            filename7 = sprintf('Results/%s_eign.txt',problemString);
            M.file2 = fopen(filename7,'w');
            
            filename8 = sprintf('Results/%s_step.txt',problemString);
            M.file3 = fopen(filename8,'w');
        end
        
        function add(M,r,xt)
            M.k = M.k + 1;
            if M.n > 2
                if size(M.V,2) < 20
                    if (~isempty(M.r))
                        M.V = [r - M.r M.V];     %| Update residual difference matrix
                        M.W = [xt - M.xt M.W];   %| Update solution difference matrix
                        M.T = [M.n M.T];         %| Update info time identifier    
                    end  
                else
                     M.V = [r - M.r M.V]; 
                     M.W = [xt - M.xt M.W];   
                     M.T = [M.n M.T];         

                     M.V(:,end) = [];
                     M.W(:,end) = [];  
                     M.T(:,end) = [];
                end
            else   
                if (~isempty(M.r))
                    M.V = [r - M.r M.V];     %| Update residual difference matrix
                    M.W = [xt - M.xt M.W];   %| Update solution difference matrix
                    M.T = [M.n M.T];         %| Update info time identifier    
                end  
            end   
            M.r = r;                     %| past residual update
            M.xt = xt;                   %| past intermediary  solution               
        end

%% Secant Jacobian computation
        function [J] = predict(M, r)
            
            
           %% Filter the retained information
           if M.filter ~= 0
                if M.filter == 1             %| QR1 -Filter
                    [M.V, M.W, M.T, M.Vprev, M.Wprev, M.Tprev, Q, R, M.count] = filt_QR1(M.V, M.W, M.T, M.Vprev, M.Wprev, M.Tprev, M.small,M.count);
                elseif M.filter == 2         %| QR2 -Filter
                    vt = [M.V M.Vprev];
                    if length(vt(1,:))>= 2
                        [M.V, M.W, M.T, M.Vprev, M.Wprev, M.Tprev, Q, R, M.count] = filt_NM(M.V, M.W, M.T, M.Vprev, M.Wprev, M.Tprev, M.small,M.count);                   
                    else                     
                        vt = [M.V M.Vprev]; 
                        [Q,R] = qr(vt,0);
                    end
                elseif M.filter == 3      %| POD -Filter
                   [M.V, M.W, M.T, M.Vprev, M.Wprev, M.Tprev, Vtemp, Wtemp, Q, R, M.count, xmodes] = filt_POD(M.V, M.W, M.T, M.Vprev, M.Wprev, M.Tprev, M.small, M.count, M.n);
                    M.count(:,1) = xmodes;
                elseif M.filter == 4      %| MGS -Filter
                    vt = [M.V M.Vprev];
                    if length(vt(1,:))>= 2
                        [M.V, M.W, M.T, M.Vprev, M.Wprev, M.Tprev, Q, R, M.count] = filt_NEW1(M.V, M.W, M.T, M.Vprev, M.Wprev, M.Tprev, M.small,M.count);                   
                    else                     
                        vt = [M.V M.Vprev]; 
                        [Q,R] = qr(vt,0);
                    end
                else                      %| Householder -Filter                              
                    vt = [M.V M.Vprev];
                    if length(vt(1,:))>= 2
                        [M.V, M.W, M.T, M.Vprev, M.Wprev, M.Tprev, Q, R, M.count] = filt_NEW2(M.V, M.W, M.T, M.Vprev, M.Wprev, M.Tprev, M.small,M.count);                   
                    else                     
                        vt = [M.V M.Vprev]; 
                        [Q,R] = qr(vt,0);
                    end
                end
           else
                vt = [M.V M.Vprev]; 
               [Q,R] = qr(vt,0);
            end
           
           %[~] = collector(M.V, M.Vprev, Q, R, M.n, M.k, M.file2); 
           
            %% Solve the least-squares porblem and compute approx. Jacobian
            if ((~isempty(M.V))||(~isempty(M.Vprev)))          
                c = R\(Q'*r);           % least-squares weights
                if M.filter == 3
                    J = Wtemp*c;
                else
                    wt = [M.W M.Wprev];
                    J = wt*c;
                end
            else
                if (isempty(M.xt))  %| Min. data from 2 preceeding iterations needed                                          
                    error('Model cannot predict without any data');                  
                else
                    J = zeros(size(M.xt));          
                end
            end
            
             if M.n >= 1  
                 format3 = '%3i %3i %11.4e %11.4e %11.4e %11.4e\n';                 
                 data3 = [M.n M.k norm(J) sum(c) min(c) max(c)]; 
                 fprintf(M.file3,format3,data3);    
             end
        end
        
%% Function maintenance and monitoring
         function clear(M)
             
            if M.reuse == 0                 %| no past information is kept
                 M.Vprev = [ ];
                 M.Wprev = [ ];
                 M.Tprev = [ ];
                 
                 M.count(:,3) = size(M.V,2);
            else 
                 M.Vprev = [M.V M.Vprev];   %| append coupling iteration matrices
                 M.Wprev = [M.W M.Wprev];
                 M.Tprev = [M.T M.Tprev];
                 
                  M.count(:,3) = size(M.Vprev,2);
                  
                  if M.filter == 3
                    limit = 150;
                  else
                    limit = 100; 
                  end
                  
                  while size(M.Vprev,2) > limit
                    last = find(M.Tprev == M.Tprev(end)); 
                    
                    M.Vprev(:,last) = [];      
                    M.Wprev(:,last) = [];
                    M.Tprev(:,last) = [];
                  end
                      
                 old = find(M.Tprev < M.n - M.reuse);   %| find outdated diff. information indicies
                 
                 if (~isempty(old))                     %| remove outdated difference information   
                    M.Vprev(:,old) = [];      
                    M.Wprev(:,old) = [];
                    M.Tprev(:,old) = [];
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
            M.count = zeros(1,3);
        end
        
        function increase_time(M)  
            M.k = 0;
            M.n = M.n + 1;                  %| increase time step
            M.clear();                      %| update retained info
        end
        
        function closefile(M)
            fclose(M.file);
            fclose(M.file2);
            fclose(M.file3);
        end
        
        function LS = ready(M)             
            LS =~ isempty([M.V M.Vprev]);   %| data must be available for Jacobian approx. 
        end
    end
end