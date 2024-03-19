function [V, W, T, Vprev, Wprev, Tprev, Q, R, count] = filt_QR1(V, W, T, Vprev, Wprev, Tprev, small, count)
%% filter interface data being collected    
    vt = [V Vprev];                   %| current and previous residual difference data
    
    singular = true;    
    % keep looping so long as singular condiiton is true                
    while (singular && (~isempty(vt)))    
       [Q,R] = qr(vt,0);                  
       [xi,i] = min(abs(diag(R)));     %| evaluate diag entries of decomposition   
		 if ((xi < small)||(xi < small*norm(R)))
             if isempty(V)
                 Vprev(:,i)=[];       %| eliminate singluar columns from past steps				
                 Wprev(:,i)=[];
				 Tprev(:,i)=[];
				 %warning('MATLAB:DataRemoved','Model removed data');
                 count(:,1) = count(:,1)+1;
             else
                  if i <= length(V(1,:))    %| eliminate singluar columns from current data
                     V(:,i)=[];              
                     W(:,i)=[];
                     T(:,i)=[];
                     %warning('MATLAB:DataRemoved','Model removed data');
                     count(:,1) = count(:,1)+1;
                  else
                     cols = length(V(1,:));
                     m = i - cols;
                     Vprev(:,m)=[];       %| eliminate singluar columns from past steps				
                     Wprev(:,m)=[];
                     Tprev(:,m)=[];
                     warning('MATLAB:DataRemoved','Model removed data');  
                     count(:,1) = count(:,1)+1;
                  end
             end
              vt = [V Vprev]; 
		  else
			 singular=false;
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
%% Output filtered QR-decomposition     
    [Q,R] = qr(vt,0); 
end