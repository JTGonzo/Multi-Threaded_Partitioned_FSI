function [V, W, T, Vprev, Wprev, Tprev, Q, R, count] = filt_NM(V, W, T, Vprev, Wprev, Tprev, small, count)
% filter interface data being collected   
    vt = [V Vprev];                   %| current and previous residual difference data
    
    singular = true;    
    % keep looping so long as singular condiiton is true         
    while (singular && (~isempty(vt))) 
        [Q,R] = qr(vt,0);
        R(1,1) = norm(vt(:,1));
		Q(:,1) = vt(:,1)./R(1,1);         
        for i = 2:length(vt(1,:))
            v_bar = vt(:,i);
            
            for j = 1:i-1
				R(j,i) = dot(Q(:,j)',v_bar);    %| update QR-elements 
				v_bar = v_bar - R(j,i)*Q(:,j);  %| older steps projected onto previous norm
            end
            
            if norm(v_bar)< small*norm(vt(:,i))
                if isempty(V)
                   Vprev(:,i)=[];       %| eliminate singluar columns from past steps				
                   Wprev(:,i)=[];
                   Tprev(:,i)=[];
				   %warning('MATLAB:DataRemoved','Model removed data');
                   count(:,1) = count(:,1)+1;
                else
                    if i <= length(V(1,:))  
                        V(:,i)=[];           %| eliminate singluar columns from current data
                        W(:,i)=[];
                        T(:,i)=[];
                        %warning('MATLAB:DataRemoved','Model removed data'); 
                        count(:,1) = count(:,1)+1;
                    else                   
                        cols = length(V(1,:));
                        m = i - cols;
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
				R(i,i) = norm(v_bar);
				Q(:,i) = v_bar./R(i,i);                
            end
            
            if i == length(vt(1,:))
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