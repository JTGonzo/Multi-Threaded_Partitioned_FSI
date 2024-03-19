function [file2] = collector(V, Vprev, Q, R, n, k, file2)
%% operate on the difference data collected   
    vt = [V Vprev];       
    
    CND = cond(vt);
    
%% Evaluate the QR minimum values   
   %[Q,R] = qr(vt,0);                  
    minR = min(abs(diag(R)));        
    maxR = max(abs(diag(R)));
    Rnorm = norm(R);
    nRs = length(find(abs(diag(R)) < 1e-7*Rnorm));
    
%% Evaluate the V eigenvalues     
    sigE = vt'*vt;     
    [X,L] = eig(sigE);          
    minL = min(abs(diag(L))); 
    maxL = max(abs(diag(L))); 
    nLs = length(find(abs(diag(L)) < 1e-22));
    
    %% Evaluate the singular values and condition number  
     S = svd(sigE);  
     minS = min(abs(S)); 
     maxS = max(abs(S)); 
     nSs = length(find(abs(S) < 1e-22));
     
%% Print Data    
     if n >= 1  
         format2 = '%3i %3i %3i %e %e %e %e %3i %e %e %3i %e %e %3i \n';                 
         data2 = [n k size(vt,2) CND minR maxR Rnorm nRs minL maxL nLs minS maxS nSs]; 
         fprintf(file2,format2,data2);    
     end
end