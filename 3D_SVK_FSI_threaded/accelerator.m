function [uS, omega] = accelerator( MESH, uS, r, nLiter, extrapolator, model, Couple, omega, Acc_flag)   
%% Relaxation methods
%%(('explicit') || ('gaussseidel')|| ('constant')|| ('aitken')|| ('NIFC'))
if Acc_flag <= 4
	if nLiter == 1
         [uS_Ex] = extrapolator.predict();           % predict initial displacement
         uS(MESH.Solid.Gamma_global) = uS_Ex;       
	elseif nLiter == 2
		omega = sign(omega)*min(Couple.omegaMax,abs(omega));                % limit the max relaxation parameter
		if Acc_flag > 1                                                    
			if Acc_flag == 2        % constant relaxation
				r = omega * r;
			elseif Acc_flag == 3    % dynamic (Aitken) relaxation
			   r(MESH.Solid.Gamma_global) = omega * r(MESH.Solid.Gamma_global);
			else                    % NIFC update
				[~,q_it] = model.predict(r(MESH.Solid.Gamma_global));
				r(MESH.Solid.Gamma_global) = q_it;
			end
		 end
		uS = uS + r;                % Gauss-Seidel update if no relaxation applied
	else
		if Acc_flag > 1               
			if Acc_flag == 2
				r = omega * r;
			elseif Acc_flag == 3
				[dyn] = model.predict(r(MESH.Solid.Gamma_global));
				omega = -omega*dyn;
				r(MESH.Solid.Gamma_global) = omega * r(MESH.Solid.Gamma_global);
			else
				[~,q_it] = model.predict(r(MESH.Solid.Gamma_global));
				r(MESH.Solid.Gamma_global) = q_it;
			end               
		end
		uS = uS + r;
	end   
end   

%% Quasi-Newton methods   
%%(('andersonacceleration') || ('broyden2')|| ('genbroyden')|| ('mvls'))
if Acc_flag >= 4
	if nLiter == 1                                              % predict initial displacement
		 [uS_Ex] = extrapolator.predict();
		 uS(MESH.Solid.Gamma_global) = uS_Ex;
	elseif ((nLiter == 2) && (~model.ready()))
		r  = omega * r;                                         % relax update if info not available 
		uS = uS + r;
	else
        dx = zeros(size(r,1),1);

        temp = model.predict(-r(MESH.Solid.Gamma_global));      % output interface approx. Joacobian (mv proudct)
        dx(MESH.Solid.Gamma_global) = temp;  
        
        duS = dx + r;                                            % both components of residual information   
        uS = uS + duS;                                           % quasi-Newton update
	end       
end
