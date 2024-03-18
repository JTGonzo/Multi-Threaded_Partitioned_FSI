function [MESH] = imap2(DATA, MESH, mapper)

%% Import solid & fluid boundary info 
dim            = MESH.dim;
nodesS         = MESH.Solid.nodes;
nodesF         = MESH.Fluid.nodes;
boundariesS    = MESH.Solid.boundaries;
boundariesF    = MESH.Fluid.boundaries;
dirichletS     = MESH.Solid.Dirichlet_dof_c;
dirichletF     = MESH.Fluid.Dirichlet_dof_c;
flag_interface = DATA.Fluid.flag_FSinterface;
flag_ALE       = DATA.Fluid.flag_ALE_fixed;
bcrow          = MESH.Fluid.bc_flag_row;
nbn            = MESH.Fluid.numBoundaryDof;

% Generate map between fluid and solid interface nodes
for k = 1 : dim
    %% identify all the solid domain interface nodes (in domain indexing)
    interfaceS_dofs{k} = [];
    for j = 1 : length(flag_interface{k})
        interfaceS_faces =  find(boundariesS(bcrow,:) == flag_interface{k}(j));
        tmp = boundariesS(1:nbn,interfaceS_faces);
        interfaceS_dofs{k}  = [interfaceS_dofs{k}; tmp(:)];
    end   
    interfaceS_dofs{k}  = unique(interfaceS_dofs{k}(:));  
    interfaceS_dofs{k}  = setdiff(interfaceS_dofs{k}, intersect(interfaceS_dofs{k},dirichletS{k}));
    
    % identify all the fluid domain interface nodes (in domain indexing)
    interfaceF_dofs{k} = [];
    for j = 1 : length(flag_interface{k})
        interfaceF_faces =  find(boundariesF(bcrow,:) == flag_interface{k}(j));
        tmp = boundariesF(1:nbn,interfaceF_faces);
        interfaceF_dofs{k}  = [interfaceF_dofs{k}; tmp(:)];
    end
    interfaceF_dofs{k}  = unique(interfaceF_dofs{k}(:));
    interfaceF_dofs{k}  = setdiff(interfaceF_dofs{k}, intersect(interfaceF_dofs{k},dirichletF{k}));
    
    % initialize mapping vector sizes
    Interface_SFmap{k} = zeros(length(interfaceS_dofs{k}),1);
    Interface_FSmap{k} = zeros(length(interfaceF_dofs{k}),1);
     
    %% New map without co-location assumption
    id_SF{k} = zeros(length(interfaceS_dofs{k}),length(interfaceF_dofs{k}));
    l_SF{k} = zeros(length(interfaceS_dofs{k}),length(interfaceF_dofs{k}));
    % to the fluid node from the solid nodes
    for i = 1 : length(interfaceS_dofs{k})
        iS_coord = nodesS(:,interfaceS_dofs{k}(i));
        for j = 1 : length(interfaceF_dofs{k})
            iF_coord = nodesF(:,interfaceF_dofs{k}(j));
            l_SF{k}(i,j) = norm(iF_coord-iS_coord);
            id_SF{k}(i,j) = j;
        end
    end
    
    % order the distance of the fluid nodes from the solid node
    for i = 1 : length(interfaceS_dofs{k})
        vecs = [l_SF{k}(i,:)', id_SF{k}(i,:)'];
        vec_or = sortrows(vecs,1);
        l_SF{k}(i,:) = vec_or(:,1)';
        id_SF{k}(i,:) = vec_or(:,2)';
    end

    id_FS{k} = zeros(length(interfaceF_dofs{k}),length(interfaceS_dofs{k}));
    l_FS{k} = zeros(length(interfaceF_dofs{k}),length(interfaceS_dofs{k}));   
    % to the solid node from the fluid node 
    for i = 1 : length(interfaceF_dofs{k})
        iF_coord = nodesF(:,interfaceF_dofs{k}(i));
        for j = 1 : length(interfaceS_dofs{k})
            iS_coord = nodesS(:,interfaceS_dofs{k}(j));
            l_FS{k}(i,j) = norm(iS_coord - iF_coord);
            id_FS{k}(i,j) = j;
        end
    end
    
    for i = 1 : length(interfaceF_dofs{k})
        vecs = [l_FS{k}(i,:)', id_FS{k}(i,:)'];
        vec_or = sortrows(vecs,1);
        l_FS{k}(i,:) = vec_or(:,1)';
        id_FS{k}(i,:) = vec_or(:,2)';
    end
    
    if mapper == 1          %| nearest-neighbour mapper       
        nf = 1;   
        cSF{k} = 1;
        cFS{k} = 1;
        
    elseif mapper == 3      %| linear mapper       
        nf = 2;
        [cSF{k}, cFS{k}] = lin_interp_3D(nodesS, nodesF, interfaceS_dofs, interfaceF_dofs, id_SF, id_FS, k);
        
    else
        nf = 12;
        [cSF{k}, cFS{k}] = rad_base(nodesS, nodesF, interfaceS_dofs{k}, interfaceF_dofs{k}, l_SF{k}, l_FS{k}, id_SF{k}, id_FS{k}, nf, k);      
    end
   
    Interface_SFmap{k} = id_SF{k}(:,1:nf);
    Interface_FSmap{k} = id_FS{k}(:,1:nf);
    
end

% Determines which of the mesh nodes are fixed 
for k = 1 : dim
      % Computes the Dirichlet sides of the domain
      nDir = length(flag_ALE{k});
      Dirichlet_side_tmp = [];
      for j = 1 : nDir         
            Dirichlet_side_j   = find(boundariesF(bcrow,:) == flag_ALE{k}(j));
            Dirichlet_side_tmp  = [Dirichlet_side_tmp, Dirichlet_side_j];
      end
      Dirichlet_side{k} = unique(Dirichlet_side_tmp);
      
      if ~isempty(Dirichlet_side{k})
            ALE_dirichlet{k} = boundariesF(1:nbn,Dirichlet_side{k});
            ALE_dirichlet{k} = unique(ALE_dirichlet{k}(:));
      else
            ALE_dirichlet{k} = [];
      end
end

%% Save the interface and boundary info to the global variable space
MESH.Fluid.dof_interface =  interfaceF_dofs;    % fluid interface node indicies
MESH.Solid.dof_interface =  interfaceS_dofs;    % solid interface node indicies
MESH.Interface_SFmap =  Interface_FSmap;        %| solid to fluid map
MESH.Interface_FSmap =  Interface_SFmap;        %| fluid to solid map
MESH.coeff_SFmap =  cSF;                        %| interp. coefficients for solid to fluid map
MESH.coeff_FSmap =  cFS;                        %| interp. coefficients for fluid to solid map
MESH.ALE_dirichlet   =  ALE_dirichlet;          % fluid domain nodes with fixed ALE boundary

return
