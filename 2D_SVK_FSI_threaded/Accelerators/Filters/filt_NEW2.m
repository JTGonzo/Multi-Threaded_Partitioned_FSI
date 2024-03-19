function [V, W, T, Vprev, Wprev, Tprev, Q, Re, count] = filt_NEW2(V, W, T, Vprev, Wprev, Tprev, small, count)
% filter interface data being collected   
    vt = [V Vprev];                   %| current and previous residual difference data
    
    singular = true;  
    
    H = @(z,x) x - z*(z'*x);
    sig = @(z) sign(z) + (z==0);
        
    % keep looping so long as singular condiiton is true         
    while (singular && (~isempty(vt)))
                
        [m,n] = size(vt);
        U = zeros(m,n);
        R = vt;
        
        for j = 1:min(m,n)
                           
            nu = norm(R(j:m,j));
            if nu ~= 0
                u = R(j:m,j)/nu;
                u(1) = u(1) + sig(u(1));
                u = u/sqrt(abs(u(1)));
            else
                u = R(j:m,j);
                u(1) = sqrt(2);
            end
                
            U(j:m,j) = u;
            
            R(j:m,j:n) = H(u,R(j:m,j:n));
            
            R(j+1:m,j) = 0;
            
            if j > 1
                v_bar = R(j,j);

                if  v_bar < small*norm(vt(:,j))
                    if isempty(V)
                       Vprev(:,j)=[];       %| eliminate singluar columns from past steps				
                       Wprev(:,j)=[];
                       Tprev(:,j)=[];
                       warning('MATLAB:DataRemoved','Model removed data');
                       count(:,1) = count(:,1)+1;
                    else
                        if j <= length(V(1,:))  
                            V(:,j)=[];           %| eliminate singluar columns from current data
                            W(:,j)=[];
                            T(:,j)=[];
                            warning('MATLAB:DataRemoved','Model removed data'); 
                            count(:,1) = count(:,1)+1;
                        else                   
                            cols = length(V(1,:));
                            m = j - cols;
                            Vprev(:,m)=[];      %| eliminate singluar columns from past steps
                            Wprev(:,m)=[];
                            Tprev(:,m)=[];
                            warning('MATLAB:DataRemoved','Model removed data');
                            count(:,1) = count(:,1)+1;
                        end
                    end
                    vt = [V Vprev]; 
                    break               
                end

                if j == length(vt(1,:))
                    singular=false;
                end
            end
        end
    end
        
%% Just want the "economy sized QR  
        Re = zeros(min(m,n));         
        for i = 1:min(m,n)            
            Re(i,:) = R(i,:);
        end
        
        Q = eye(size(U));
        for k = n:-1:1           
            Q = H(U(:,k),Q);
        end
        
%% Drop columns (info) if they exceed row (interface dofs) size
    if ~isempty(V)
        while length(vt(1,:)) > length(vt(:,1))
             if length(V(1,:)) > length(V(:,1))
                    V(:,end) = [];   
                    W(:,end) = [];   
                    T(:,end) = [];
                    count(:,2) = count(:,2)+1;
             else 
                    Vprev(:,end) = [];   
                    Wprev(:,end) = [];   
                    Tprev(:,end) = [];
                    count(:,2) = count(:,2)+1;
             end
             vt = [V Vprev];
        end
    end
end