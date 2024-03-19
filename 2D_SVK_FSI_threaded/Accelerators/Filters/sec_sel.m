%function [idx] = sec_sel(rang, V, W, T, T2, K, nR, Vprev, Wprev, Tprev, T2prev, Kprev, nRprev, r)
function [idx, idxt] = sec_sel(adaptv, nR, K, nRprev, Kprev, r, k)
    
    r_star = norm(r);
    
    rt = [nR nRprev];
    kt = [K Kprev];    
    
    if k == 1        
        idxt = find(kt == k);
        idx = zeros(length(idxt)-1,1);
        
        for i = 2:length(idxt)
            idx(i-1) = idxt(i)-i; 
        end
        
        idxt(1) = [];
    else
        idxt = [];
        idt = find(kt == 1);
        for i=length(rt):-1:1
           if (rt(i) >= adaptv.rangL*r_star) && (rt(i) <= (1/adaptv.rangU)*r_star)
                idxt = [i idxt];
           end       
        end
        
        rem = [];
        for j = length(idt):-1:1
            temp = find(idxt == idt(j));
            if ~isempty(temp)
                rem = [temp rem];
            end
        end
        
        if ~isempty(rem)
            idxt(rem) = [];
        end
        
        idx = zeros(length(idxt),1);
        for m = 1:length(idxt)           
            for n = 1:length(idt)
                if idxt(m) < idt(n)
                    idx(m) = idxt(m)-n+1;
                    break
                end
            end          
        end
    end
end