function [cSF, cFS] = rad_base(nodesS, nodesF, interfaceS_dofs, interfaceF_dofs, l_SF, l_FS, id_SF, id_FS, nf, k) 

shape_p = 50;

% dec = 1e-11;
% msh = .1;
% vrt = 5;
% shape_p = sqrt(-log(dec))/(msh*vrt)

%% Solid to Fluid Mapping
cSF = zeros(length(interfaceS_dofs),nf);
l_SF_norm = zeros(length(interfaceS_dofs),nf);
phi_SF = zeros(nf,length(interfaceS_dofs));
cond_SF   = zeros(length(interfaceS_dofs),1);
dref_SF = zeros(length(interfaceS_dofs),1);

for i = 1 : length(interfaceS_dofs)
    temp_SF = zeros(nf,nf);
    Phi_G_SF{i} = zeros(nf,nf);
end

for i = 1 : length(interfaceS_dofs)
    dref_SF(i) = l_SF(i,nf)*shape_p;
    for j = 1:nf
        l_SF_norm(i,j) = l_SF(i,j)/dref_SF(i);
%         phi_SF(j,i) = ((1 - l_SF_norm(i,j))^4 * (1 + 4 * l_SF_norm(i,j)))*heaviside(1-l_SF_norm(i,j));
        if l_SF_norm(i,j) >= 1    
            phi_SF(j,i) = 0;
        else
            phi_SF(j,i) = (1 - l_SF_norm(i,j))^4 * (1 + 4 * l_SF_norm(i,j));
        end
    end
end

for i = 1 : length(interfaceS_dofs)
    for j = 1: nf
        Pt = nodesF(:,interfaceF_dofs(id_SF(i,j)));
        for l = 1 : nf
            Pf =  nodesF(:,interfaceF_dofs(id_SF(i,l)));
            temp = norm(Pt - Pf);
            dbet = temp/dref_SF(i);
%             temp_SF(j,l) = ((1 - dbet)^4 * (1 + 4 * dbet))*heaviside(1-dbet);
            if dbet >= 1             
                temp_SF(j,l) = 0;
            else
                temp_SF(j,l) = (1 - dbet)^4 * (1 + 4 * dbet);
            end
        end
    end
    Phi_G_SF{i} = temp_SF;

    cond_SF(i) = cond(temp_SF);
    if cond_SF(i) > 1e13
		sfprintf('error') 
    end
    
end

for i = 1 : length(interfaceS_dofs)
    temp = Phi_G_SF{i};
   cSF(i,:) =  temp\phi_SF(:,i);
end

%% Fluid to Solid Mapping
cFS = zeros(length(interfaceF_dofs),nf);
l_FS_norm = zeros(length(interfaceF_dofs),nf);
phi_FS = zeros(nf,length(interfaceF_dofs));
cond_FS   = zeros(length(interfaceF_dofs),1);
dref_FS = zeros(length(interfaceF_dofs),1);

for i = 1 : length(interfaceF_dofs)
    temp_FS = zeros(nf,nf);
    Phi_G_FS{i} = zeros(nf,nf);
end

for i = 1 : length(interfaceF_dofs)
    dref_FS(i) = l_FS(i,nf)*shape_p;
    for j = 1:nf
        l_FS_norm(i,j) = l_FS(i,j)/dref_FS(i);
%         phi_FS(j,i) = ((1 - l_FS_norm(i,j))^4 * (1 + 4 * l_FS_norm(i,j)))*heaviside(1-l_FS_norm(i,j));
        if l_FS_norm(i,j) >= 1    
            phi_FS(j,i) = 0;
        else
            phi_FS(j,i) = (1 - l_FS_norm(i,j))^4 * (1 + 4 * l_FS_norm(i,j));
        end
    end
end

for i = 1 : length(interfaceF_dofs)
    for j = 1: nf
        Pt = nodesS(:,interfaceS_dofs(id_FS(i,j)));
        for l = 1 : nf
            Pf =  nodesS(:,interfaceS_dofs(id_FS(i,l)));
            temp = norm(Pt - Pf);
            dbet = temp/dref_FS(i);
%             temp_FS(j,l) = ((1 - dbet)^4 * (1 + 4 * dbet))*heaviside(1-dbet);
            if dbet >= 1             
                temp_FS(j,l) = 0;
            else
                temp_FS(j,l) = (1 - dbet)^4 * (1 + 4 * dbet);
            end
        end
    end
    Phi_G_FS{i} = temp_FS;

    cond_FS(i) = cond(temp_FS);
    if cond_FS(i) > 1e13
		sfprintf('error') 
    end
    
end

for i = 1 : length(interfaceF_dofs)
    temp = Phi_G_FS{i};
   cFS(i,:) =  temp\phi_FS(:,i);
end

end