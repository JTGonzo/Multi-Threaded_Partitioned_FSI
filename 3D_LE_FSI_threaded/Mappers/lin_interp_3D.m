function [cSF, cFS] = lin_interp_3D(nodesS, nodesF, interfaceS_dofs, interfaceF_dofs, id_SF, id_FS, k) 

    cSF = zeros(length(interfaceS_dofs{k}),3);
	for i = 1 : length(interfaceS_dofs{k})
		P0 = nodesS(:,interfaceS_dofs{k}(i));
		P1 = nodesF(:,interfaceF_dofs{k}(id_SF{k}(i,1)));
		P2 = nodesF(:,interfaceF_dofs{k}(id_SF{k}(i,2)));
        P3 = nodesF(:,interfaceF_dofs{k}(id_SF{k}(i,3)));
        
        V12 = P2 - P1;
        V13 = P3 - P1;
        V23 = P3 - P2;

        area = norm(cross(V12,V13));
        long = max(norm(V12), norm(V13), norm(V23));
        AR = (long^2)/area;
        
        if ((area == 0)||(AR > 30))
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
        else    
               Vn = cross(P2-P1,P3-P1);
               Vn = Vn./norm(Vn);
               Pp = P0 + Vn*dot(P1-P0,Vn);
               
               a = cross(Pp-P1,P2-P1);
               b = cross(Pp-P2,P3-P2);
               c = cross(Pp-P3,P1-P3);
               
               if ((dot(a,b)>=0) && (dot(a,c)>=0))
                   a_ref = norm(cross(P2-P1,P3-P1)/2);
                   cSF(i,1) = norm(cross(P2-Pp,P3-Pp)/2)/a_ref;
                   cSF(i,2) = norm(cross(P1-Pp,P3-Pp)/2)/a_ref;
                   cSF(i,3) = 1 - cSF(i,1) - cSF(i,2);
               else
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
        end                              
	end

	cFS = zeros(length(interfaceF_dofs{k}),3);
	for i = 1 : length(interfaceF_dofs{k})
		P0 = nodesF(:,interfaceF_dofs{k}(i));
		P1 = nodesS(:,interfaceS_dofs{k}(id_FS{k}(i,1)));
		P2 = nodesS(:,interfaceS_dofs{k}(id_FS{k}(i,2)));
        P3 = nodesS(:,interfaceS_dofs{k}(id_FS{k}(i,3)));

        V12 = P2 - P1;
        V13 = P3 - P1;
        V23 = P3 - P2;

        area = norm(cross(V12,V13));
        long = max(norm(V12), norm(V13), norm(V23));
        AR = (long^2)/area;
        
        if ((area == 0)||(AR > 30))
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
        else    
		    Vn = cross(P2-P1,P3-P1);
		    Vn = Vn./norm(Vn);
		    Pp = P0 + Vn*dot(P1-P0,Vn);
		   
		    a = cross(Pp-P1,P2-P1);
		    b = cross(Pp-P2,P3-P2);
		    c = cross(Pp-P3,P1-P3);
		   
		    if ((dot(a,b)>=0) && (dot(a,c)>=0))
			   a_ref = norm(cross(P2-P1,P3-P1)/2);
			   cFS(i,1) = norm(cross(P2-Pp,P3-Pp)/2)/a_ref;
			   cFS(i,2) = norm(cross(P1-Pp,P3-Pp)/2)/a_ref;
			   cFS(i,3) = 1 - cFS(i,1) - cFS(i,2);
		    else
				Vn = P1 - P0;
				Vt = (P2 - P1)/norm(P2 - P1);
				Pp = P0 + (Vn - Vt*dot(Vn,Vt));

				if dot(P1-Pp,P2-Pp)<= 0
                    cFS(i,1) = norm(P2-Pp)/norm(P2-P1);
                    cFS(i,1) = 1 - cFS(i,1);
                else
                    cFS(i,1) = 1;
                    cFS(i,1) = 1 - cFS(i,1);
                end
		    end
        end                              
	end
end