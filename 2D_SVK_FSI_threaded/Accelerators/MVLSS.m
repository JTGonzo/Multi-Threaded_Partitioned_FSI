% Behr MVLS Techniques
classdef MVLSS < handle
    properties (SetAccess=private)
        n = 0;
        r = []; 
        xt = [];
        V = [];
        W = [];
        Z = [];
        b = [];
        B = [];
        T = [];
        Wprev = [];
        Vprev = [];
        Zprev = [];
        Dummy = [];
        file;
        count = zeros(1,3);
    end
    
    properties (SetAccess=immutable)
        rthresh;
        small;
        limit;
        filter;
        problemString;
    end

%% Function initialization and interface-data retention    
    methods
        function M = MVLSS(rthresh, small, limit, filter, count, problemString)
            M.rthresh = rthresh;
            M.small = small;         %| smallest floating point number
            M.limit = limit;         %| number of time steps to reuse
            M.filter = filter;       %| what type of filter to use
            
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
                 
            if ((isempty(M.Vprev))&&(~isempty(M.V))) 
                vt = [M.V M.Vprev];
                [Q,R] = qr(vt,0);
                c = R\(Q'*r);
                dx = M.W*c;
            else
                if (isempty(M.xt))
                    error('Model cannot predict without any data');
                else
                    dx = zeros(size(M.xt));
                end
            end
            
            if(~isempty(M.Vprev))
                if length(M.Vprev)>= 1
                    i = 0;
                    if M.n-1 <= M.limit
                        kn = M.n-1;
                    else
                        kn = M.limit;
                    end
                end
                
                if(isempty(M.V))
                    dx = zeros(length(r),1);
                    while i<min(M.limit,length(M.Vprev))
                        dx = dx + M.Wprev{kn}*M.Zprev{kn}*r;
                        r = r + M.Vprev{kn}*M.Zprev{kn}*r;  
                        kn = kn - 1;
                        i = i + 1; 
                    end
                    M.b = dx;
                else
                    if M.filter ~= 0
                        if M.filter == 1    %| QR1 -Filter
                            [M.V, M.W, M.T, ~, ~, ~, Q, R, M.count] = filt_QR1(M.V, M.W, M.T, M.Dummy, M.Dummy, M.Dummy, M.small, M.count);
                        else
                            vt = [M.V M.Dummy];
                            if length(vt(1,:))>= 2
                                [M.V, M.W, M.T, ~, ~, ~, Q, R, M.count] = filt_NM(M.V, M.W, M.T, M.Dummy, M.Dummy, M.Dummy,  M.small, M.count);                   
                            else                     
                                vt = [M.V M.Dummy]; 
                                [Q,R] = qr(vt,0);
                            end
                        end
                    else
                        vt = [M.V M.Dummy]; 
                       [Q,R] = qr(vt,0);
                    end

                    idx = zeros(length(r),1);
                    alpha = r;

                    while i<min(M.limit,length(M.Vprev))
                        idx = idx + M.Wprev{kn}*M.Zprev{kn}*alpha;
                        alpha = alpha + M.Vprev{kn}*M.Zprev{kn}*alpha;  
                        kn = kn - 1;
                        i = i + 1; 
                    end     

                    btemp = M.W(:,1) - M.b + idx;  
                    M.B = [btemp M.B];
                    M.b = idx;

                    c = R\(Q'*r);               
                    dx = -M.B*c;
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
           
            M.r  = [];      %| clear current coupling iteration matrices
            M.xt = [];
            M.V  = [];
            M.W  = [];
            M.Z  = [];
            M.T  = [];
            M.b  = [];
            M.B  = [];
            M.count = zeros(1,3);
        end

        function increase_time(M)
            M.n = M.n+1;

            if M.n > 1
                temp = M.V'*M.V; 
                M.Z = LUinverse(temp)*M.V';
                
                %| append coupling iteration matrices 
                if M.n-1 <= M.limit                    
                    M.Vprev{M.n-1} = M.V;
                    M.Zprev{M.n-1} = M.Z;
                    M.Wprev{M.n-1} = M.W;
                    
                    for i=1:M.n-1
                        M.count(:,3) = M.count(:,3) + size(M.Wprev{i},2);
                    end
                else
                %| remove outdated difference information
                    M.Vprev{1} = [];   
                    M.Zprev{1} = [];
                    M.Wprev{1} = [];
                    
                    for i=2:M.limit
                        M.Vprev{i-1} = M.Vprev{i};
                        M.Zprev{i-1} = M.Zprev{i};
                        M.Wprev{i-1} = M.Wprev{i};
                    end
                    
                    M.Vprev{M.limit} = M.V;
                    M.Zprev{M.limit} = M.Z;
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
            LS =~isempty([M.V M.Vprev]); %| data must be available for Jacobian approx. 
        end
    end
end