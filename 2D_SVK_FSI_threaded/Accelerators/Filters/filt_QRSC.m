function [Vprev, Wprev, Tprev, T2prev, Kprev, nRprev, nXprev, nVprev, nWprev, ...
    V, W, T, T2, K, nR, nX, nV, nW, Q, R, count, idx] = filt_QRSC(M, adaptv, small, count, k, r)
    
    [Vprev, Wprev, Tprev, T2prev, Kprev, nRprev, nXprev, nVprev, nWprev, ...
        V, W, T, T2, K, nR, nX, nV, nW] = local_vars(M);
    
%% filter interface data being collected   
    vt = [V Vprev];                   %| current and previous residual difference data   
    [idx, idxt] = sec_sel(adaptv, nR, K, nRprev, Kprev, r, k);    
    vt = vt(:,idx);
    
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
                if isempty(V)         %| eliminate singluar columns from past steps
                     
                     idx(i) = [];
                     idxt(i) = [];
                     
                     vt = [V Vprev];
                     vt = vt(:,idx);
                     
%                    ind1 = idx(i);
%                    ind2 = idxt(i)-2; 
                    
%                    Vprev(:,ind1)=[];       			
%                    Wprev(:,ind1)=[];
%                    Tprev(:,ind1)=[];
% 
%                    Kprev(:,ind2)  = [];       				
%                    T2prev(:,ind2) = [];                   
%                    nRprev(:,ind2) = [];
%                    nXprev(:,ind2) = [];
%                    nVprev(:,ind2) = [];                 
%                    nWprev(:,ind2) = [];   
%                      
%                    warning('MATLAB:DataRemoved','Model removed data');
%                    count(:,1) = count(:,1)+1;
                else 
                    ind1 = idx(i);
                    ind2 = idxt(i);
                    
                    if ind1 <= length(V(1,:))   %| eliminate singluar columns from current data
                        V(:,ind1)=[];          
                        W(:,ind1)=[];
                        T(:,ind1)=[];

                        T2(:,ind1) = [];
                        K(:,ind1)  = [];
                        nR(:,ind1) = [];
                        nX(:,ind1) = [];
                        nV(:,ind1) = [];
                        nW(:,ind1) = [];
                        
                        warning('MATLAB:DataRemoved','Model removed data'); 
                        count(:,1) = count(:,1)+1;
                    else                   
                        cols = length(V(1,:));
                        m1 = ind1 - cols;
                        m2 = ind2 - (cols + 1);

                        Vprev(:,m1) = [];      %| eliminate singluar columns from past steps
                        Wprev(:,m1) = [];
                        Tprev(:,m1) = [];

                        Kprev(:,m2)  = [];       %| eliminate singluar columns from past steps				
                        T2prev(:,m2) = [];                   
                        nRprev(:,m2) = [];
                        nXprev(:,m2) = [];
                        nVprev(:,m2) = [];                 
                        nWprev(:,m2) = [];  
                   
                        warning('MATLAB:DataRemoved','Model removed data');
                        count(:,1) = count(:,1)+1;
                    end
                
                    vt = [V Vprev]; 
                    [idx, idxt] = sec_sel(adaptv, nR, K, nRprev, Kprev, r, k);                
                    vt = vt(:,idx);
                end
                
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
                 
				V(:,end-2) = [];   
				W(:,end-2) = [];   
				T(:,end-2) = [];
				
				T2(:,end-3) = [];
				K(:,end-3)  = [];
				nR(:,end-3) = []; 
                nX(:,end-3) = [];
                nV(:,end-3) = [];
                nW(:,end-3) = [];
                        
				count(:,2) = count(:,2)+1;
             else 
                 
                old = find(Tprev == min(Tprev));   
                old2 = find(T2prev == min(T2prev));
                     
                 if (~isempty(old))                     
                    Vprev(:,old) = [];      
                    Wprev(:,old) = [];
                    Tprev(:,old) = [];

                    Kprev(:,old2) = [];                    
                    T2prev(:,old2) = [];                    
                    nRprev(:,old2) = [];
                    nXprev(:,old2) = [];
                    nVprev(:,old2) = [];                 
                    nWprev(:,old2) = [];                      
                 end                  
                       
                 count(:,2) = count(:,2)+1;
             end
             vt = [V Vprev];
        end
    end
end