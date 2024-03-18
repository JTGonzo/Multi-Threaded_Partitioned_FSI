classdef Convergence < handle
%% Set in-function fixed (internal) and dynamic (interface) variable values 
    properties (SetAccess=private)
        reso = 0.0;   
        resk = 0.0;
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
        problemString;
    end

%% Coupling iterations convergence checking 
    methods
        function C=Convergence(tol,k_min,k_e,n_e,small,output,flag,problemString)
        % Parameter definition/initialization
            C.tol = tol;                              %|  defined convergence tolerance
            C.k_min = k_min;                          %|  minimum number of coupling iterations
            C.k_e = k_e;                              %|  maximum number of coupling iterations
            C.n_e = n_e;                              %|  current time step 
            C.small = small;                          %|  smallest floating point number allowed
            C.output = output;                        %|  user defined output flag
            C.flag = flag;
            
            if (C.output)
                filename4 = sprintf('Results/%s_residuals.txt',problemString);   %| iterations residual info data file
                C.file = fopen(filename4,'w');                  %| create residual output file if desired 
            end            
        end
 
 %% Add/Update info about current FP residual data      
        function add(C,r,rT)   
            C.k = C.k + 1;                 %| update coupling iteration number             
            C.k_tot = C.k_tot + 1;         %| collect coupling iteration count
            C.resk = norm(r);              %| collect current coupling iteration residual norm
            C.resT = norm(rT); 
            
            if (C.k == 1)
                C.reso = C.resk;           %| initial FP iteration residual norm
            end
            
            if (C.output)
                    format = '%4i %4i %11.4e %11.4e \n';    %| data format         
                    data = [C.n C.k C.resk C.resT];         %| time-step, coupling iteration, residual 
                    fprintf(C.file,format,data);     %| write results to file
                    fprintf(format,data);            %| print results to command line 
            end
        end

 %% Evaluate if the convergence criteria is achieved       
        function conv=is_satisfied(C)
            
            if ((C.resk<max(C.reso*C.tol,C.small)) && (C.k>C.k_min))   % check if solution converges as desired
                conv=true;
            else
                conv=false;
            end
            
            if ((C.flag == 0) && (C.k >= 2))      % if explicit update is used (no FP iterations)
                 conv=true;
            end
            
            if ((C.k==C.k_e) && (~conv))                               % display error message if convergence not achieved
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
