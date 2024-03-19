function [V, W, T, Vprev, Wprev, Tprev, Q, R, count] = filt_NEW1(V, W, T, Vprev, Wprev, Tprev, small, count)
% filter interface data being collected   
    vt = [V Vprev];                   %| current and previous residual difference data
    
    singular = true;  
    
    % keep looping so long as singular condiiton is true         
    while (singular && (~isempty(vt)))
                
        [~,n] = size(vt);
        Q = vt;
        R = zeros(n,n);

        for j = 1:n
            R(j,j) = norm(Q(:,j));
            
            if (R(j,j) == 0)
                error('linearly dependent columns');
            end
            
            v_bar = Q(:,j);
            Q(:,j) = Q(:,j)/R(j,j);
            
            if norm(v_bar)< small*norm(vt(:,j))
                if isempty(V)
                   Vprev(:,j)=[];       %| eliminate singluar columns from past steps				
                   Wprev(:,j)=[];
                   Tprev(:,j)=[];
				   %warning('MATLAB:DataRemoved','Model removed data');
                   count(:,1) = count(:,1)+1;
                else
                    if j <= length(V(1,:))  
                        V(:,j)=[];           %| eliminate singluar columns from current data
                        W(:,j)=[];
                        T(:,j)=[];
                        %warning('MATLAB:DataRemoved','Model removed data'); 
                        count(:,1) = count(:,1)+1;
                    else                   
                        cols = length(V(1,:));
                        m = j - cols;
                        Vprev(:,m)=[];      %| eliminate singluar columns from past steps
                        Wprev(:,m)=[];
                        Tprev(:,m)=[];
                        %warning('MATLAB:DataRemoved','Model removed data');
                        count(:,1) = count(:,1)+1;
                    end
                end
                vt = [V Vprev]; 
                break                
            else
                % update the elements of the QR-decomposition
                R(j,j+1:n) = Q(:,j)'*Q(:,j+1:n);
                Q(:,j+1:n) = Q(:,j+1:n) - Q(:,j)*R(j,j+1:n);
            end
            
            if j == length(vt(1,:))
                singular=false;
            end
        end
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