classdef Extrapolator < handle
%% Set in-function fixed (internal) and dynamic (interface) variable values 
    properties (SetAccess=private)
        xp;
        x;
        xn;
        x2n;
        x3n;
        flag;
        
        n  = 0;
        n0 = 0;
        n1 = 1;
        n2 = 2;
        n3 = 3;
    end

%% New time step variable prediction
    methods
        function E = Extrapolator(x, flag)       
            E.initialize(x, flag);                        %| initialize the interface matrices              
        end 
        
        function [xe] = predict(E)
            
            if E.flag == 1 % Constant Extapolation
                xe = E.x;
            elseif E.flag == 2
                switch E.n                
                    case E.n0                                 %| No extrapolation for 0 time step 
                        error('Extrapolator cannot predict initial value');                    
                    case E.n1
                        xe = E.x;                             %| first time step prediction  (constant)                
                    otherwise    
                        xe = 2.0*E.x - E.xn;                  %| second time step prediction (linear)                    
                end    
            elseif E.flag == 3
                switch E.n                
                    case E.n0                                 %| No extrapolation for 0 time step 
                        error('Extrapolator cannot predict initial value');                    
                    case E.n1
                        xe = E.x;                             %| first time step prediction  (constant)                   
                    case E.n2
                        xe = 2.0*E.x - E.xn;                  %| second time step prediction (linear)                    
                    otherwise
                        xe = 2.5*E.x - 2.0*E.xn + 0.5*E.x2n;  %| second time step prediction (legacy) 
                end
            elseif E.flag == 4  
                switch E.n                
                    case E.n0                                 %| No extrapolation for 0 time step 
                        error('Extrapolator cannot predict initial value');                    
                    case E.n1
                        xe = E.x;                             %| first time step prediction  (constant)                   
                    case E.n2
                        xe = 2.0*E.x - E.xn;                  %| second time step prediction (linear)                    
                    otherwise
                        xe = 3.0*E.x - 3.0*E.xn + 1.0*E.x2n;  %| 3 subsequent time steps predictions (quadratic)
                end
            else   
               switch E.n                
                    case E.n0                                 %| No extrapolation for 0 time step 
                        error('Extrapolator cannot predict initial value');                    
                    case E.n1
                        xe = E.x;                             %| first time step prediction  (constant)                   
                    case E.n2
                        xe = 2.0*E.x - E.xn;                  %| second time step prediction (linear)                    
                    case E.n3
                        xe = 3.0*E.x - 3.0*E.xn + 1.0*E.x2n;  %| 3 subsequent time steps predictions (quadratic)
                   otherwise
                       xe = 4.0*E.x - 6.0*E.xn + 4.0*E.x2n - 1.0*E.x3n;  %| 4 subsequent time steps predictions (quadratic)
               end
            end
        end
        
        function add(E,x)          
            E.xp = x;      %| adding new time step information to the predictor array      
        end
 
%% Dynamically maintaining the data used for extrapolation
        function initialize(E,x,flag)   %| initializing the (dispacement) vector size
            E.xp  = x;
            E.x   = x;
            E.xn  = x;
            E.x2n = x;
            E.x3n = x;
            E.flag = flag;
        end
        
        function shift(E)          %| shifting what we define as the previous 3 time-step's data values            
            E.x3n = E.x2n;
            E.x2n = E.xn;         
            E.xn  = E.x;
            E.x   = E.xp;
        end
        
%% Updating time increment of handle         
        function increase_time(E)            
            E.n = E.n + 1;
            E.shift();            
        end
    end
end