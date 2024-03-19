function [cSF, cFS] = lin_interp(nodesS, nodesF, interfaceS_dofs, interfaceF_dofs, id_SF, id_FS, k) 

    cSF = zeros(length(interfaceS_dofs{k}),2);
	for i = 1 : length(interfaceS_dofs{k})
		P0 = nodesS(:,interfaceS_dofs{k}(i));
		P1 = nodesF(:,interfaceF_dofs{k}(id_SF{k}(i,1)));
		P2 = nodesF(:,interfaceF_dofs{k}(id_SF{k}(i,2)));

		Vn = P1 - P0;
		Vt = (P2 - P1)/norm(P2 - P1);
		Pp = P0 + (Vn - Vt*dot(Vn,Vt));

		if dot(P1-Pp, P2-Pp) <= 0
			cSF(i,1) = norm(P2 - Pp)/norm(P2 - P1);
			cSF(i,2) = 1 - cSF(i,1);
		else
			cSF(i,1) = 1;
			cSF(i,2) = 1 - cSF(i,1);
		end        
	end

	cFS = zeros(length(interfaceF_dofs{k}),2);
	for i = 1 : length(interfaceF_dofs{k})
		P0 = nodesF(:,interfaceF_dofs{k}(i));
		P1 = nodesS(:,interfaceS_dofs{k}(id_FS{k}(i,1)));
		P2 = nodesS(:,interfaceS_dofs{k}(id_FS{k}(i,2)));

		Vn = P1 - P0;
		Vt = (P2 - P1)/norm(P2 - P1);
		Pp = P0 + (Vn - Vt*dot(Vn,Vt));

		if dot(P1-Pp, P2-Pp) <= 0
			cFS(i,1) = norm(P2 - Pp)/norm(P2 - P1);
			cFS(i,2) = 1 - cFS(i,1);
		else
			cFS(i,1) = 1;
			cFS(i,2) = 1 - cFS(i,1);
		end        
	end   
end