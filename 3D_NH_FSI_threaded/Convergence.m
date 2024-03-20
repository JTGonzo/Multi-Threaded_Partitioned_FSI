classdef Convergence < handle
%% Set in-function fixed (internal) and dynamic (interface) variable values 
    properties (SetAccess=private)
        reso = 0.0;   
        resk = 0.0;
		resD = 0.0;
        resV = 0.0;
        resP = 0.0;
		resT = 0.0;		
        k = 0;
        k_tot = 0;
        n = 0;
        file;
    end
    
    properties (SetAccess=immutable)
        tol;
        k_min;
        k_e;
        n_e;
        small;
        output;
        flag;
		type;
		problemString;
    end

%% Coupling iterations convergence checking 
    methods
        function C=Convergence(tol,k_min,k_e,n_e,small,output,flag,type,problemString)
        % Parameter definition/initialization
            C.tol = tol;                              %|  defined convergence tolerance
            C.k_min = k_min;                          %|  minimum number of coupling iterations
            C.k_e = k_e;                              %|  maximum number of coupling iterations
            C.n_e = n_e;                              %|  current time step 
            C.small = small;                          %|  smallest floating point number allowed
            C.output = output;                        %|  user defined output flag
            C.flag = flag;
            C.type = type; 
			
            if (C.output)
                filename4 = sprintf('Results/%s_residuals.txt',problemString);   %| iterations residual info data file
                C.file = fopen(filename4,'w');                  				 %| create residual output file if desired 
            end            
        end
 
 %% Add/Update info about current FP residual data      
        function add(C,r,rD,rV,rP,rT)   
            C.k = C.k + 1;                 %| update coupling iteration number             
            C.k_tot = C.k_tot + 1;         %| collect coupling iteration count
            C.resk = norm(r);              %| collect current coupling iteration residual norm
            C.resD = norm(rD);			   %| interface displacement resdiual norm
			C.resV = norm(rV); 			   %| interface velocity resdiual norm
            C.resP = norm(rP); 			   %| interface pressure resdiual norm
            C.resT = norm(rT); 			   %| interface traction resdiual norm
			
            if (C.k == 1)
                C.reso = C.resk;           %| initial FP iteration residual norm
            end
            
            if (C.output)
                    format = '%4i %4i %11.4e %11.4e %11.4e %11.4e %11.4e \n';    %| data format         
                    data = [C.n C.k C.resk C.resD C.resV C.resP C.resT];         %| time-step, coupling iteration, residual 
                    fprintf(C.file,format,data);     %| write results to file
                    fprintf(format,data);            %| print results to command line 
            end
        end

 %% Evaluate if the convergence criteria is achieved       
        function conv=is_satisfied(C)
			
			% check if solution converges per the set convergence requirement 
			if C.type == 1
				cvg_sat = ((C.k > C.k_min) && (C.resk < max(C.reso*C.tol,C.small))); 
			elseif C.type == 2
				cvg_sat = ((C.k > C.k_min) && (C.resk < C.small));
			elseif C.type == 3
				cvg_sat = ((C.k > C.k_min) && (C.resT < C.small));
			else
				cvg_sat = ((C.k == C.k_e) && (C.resk < C.reso*20));
			end

            if cvg_sat == 1
                conv=true;
            else
                conv=false;
            end
            
            if ((C.flag == 0) && (C.k >= 2))      % if explicit update is used (no FP iterations)
                 conv=true;
            end
            
            if ((C.k==C.k_e) && (~conv))          % display error message if convergence not achieved
                error(['Coupling not converged after maximal number' ...
                    ' of iterations']);
            end
        end
        
      function res = res_int(C)
           res = C.reso;
       end
        
 %% Updating time increment of handle     
        function increase_time(C)
            C.k = 0;
            C.n = C.n+1;          % increasing time step
        end
        
        % close output file once done
        function closefile(C)
            fclose(C.file);
        end
    end
end