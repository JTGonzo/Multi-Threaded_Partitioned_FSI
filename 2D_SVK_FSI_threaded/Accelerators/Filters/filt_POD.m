function [V, W, T, Vprev, Wprev, Tprev, vt, wt, Q, R, count, xmodes] = filt_POD(V, W, T, Vprev, Wprev, Tprev, small, count, nn)
%% filter interface data being collected   
    vt = [V Vprev];      %| current and previous residual difference data
    wt = [W Wprev];      %| current and previous solution difference data
    xmodes = 0;
    
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
             wt = [W Wprev];
        end 
    end

%% Perform POD Filter    
    if (~isempty(vt))           
         eta = length(vt(1,:));
         sig = (1/eta)*(vt'*vt);     %| form auto-correlation matrix
         [X,lam] = eig(sig);         %| evaluates its eigenvalues & vectors   

         lam_n = diag(sort(diag(lam),'descend'));  %| reorder eigenvalue matrix
         [p, ind] = sort(diag(lam),'descend');     %| store the indices of eigenvalue reordering
         X_n = X(:,ind);                           %| arrange the eigenvectors
         
         % find index where eigenvlaues become tiny
         for i = 1:eta
               if lam_n(i,i)< small
                    c = i-1;
                    %warning('MATLAB:DataRemoved','Model removed data');
                    xmodes = eta - c;
                    break
               else 
                   c = i;
               end
         end

         % project the difference matricies onto the eigenvector space
         V_pod = vt*X_n;  %re-verify this multiplication
         W_pod = wt*X_n;
        
         % only keep the most dominant modes
         vt = V_pod(:,1:c);
         wt = W_pod(:,1:c);
    end

%% Output filtered QR-decomposition 
[Q,R] = qr(vt,0); 
end
