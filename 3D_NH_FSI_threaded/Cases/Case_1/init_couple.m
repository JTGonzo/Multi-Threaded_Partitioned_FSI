%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  Initialize Coupling solver setting, export files %%%%%%%%%%
%%%%%% convergence/extrapolation schemes and accleration method  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fixed-Point iteration parameters
Couple.rtol = 1e-4;                 %| relative residual error
Couple.imin = 1;                    %| min number of fix-point iterations
Couple.imax = 50;                   %| max number of fix-point iterations
Couple.maxSteps = tf/dt;            %| number of time steps
Couple.probeOutFreq = 1;            %| physical variable output frequency
Couple.out = 1;                     %| output frequency for coupling data
Couple.probeOutFreq2 = 1;            %| physical variable output frequency

% Convergence acceleration parameters
algorithm = 'AndersonAcceleration';
Couple.omegaMax = 0.75; %0.1        %| relaxation factor
Couple.reuse = 10;                   %| time-steps to reuse
Couple.small = 1e-10;               
Couple.filter = 2;                  %| QR1 = 1; NM = 2; POD =3
Couple.count = zeros(1,3);          %| filtered/dropped columns
Couple.extra = 2;                   %| Const = 1; Linear = 2; Legacy =3; Quadratic =4; Cubic =5;

%% Initialize fixed-point iteration variables and data output files
r = zeros(length(uS),1);              % initial displacement residual 
Count = zeros(int32(tf/dt),2);              %| iteration counter
omega  = Couple.omegaMax;             %| initial relaxation factor
Wid = 33;                              %| watchpoint Id
WidF = 29;                            %| watchpoint Id

filename1 = sprintf('Results/%s_iters.oisd',problemString);   % iterations residual info data file
fileId1 = fopen(filename1,'w');
filename2 = sprintf('Results/%s_disp.othd',problemString);      % variable info data file
fileId2 = fopen(filename2,'w');
filename3 = sprintf('Results/%s_vel.othd',problemString);
fileId3 = fopen(filename3,'w');

%% Initalizing Vairables based on what coupling scheme is chosen
switch lower(algorithm)
    case 'explicit'
        Acc_flag = 0;
        model = Dummy();
    case 'gaussseidel'
        Acc_flag = 1;
        model = Dummy();
        omega  = 1;          
    case 'constant'  
        Acc_flag = 2; 
        model = Dummy();
        omega = .9;    
    case 'aitken'  
        Acc_flag = 3;
        model = Aitken(); 
    case 'nifc'   
        Acc_flag = 4;
        p_it = ones(length(r(MESH.Solid.Gamma_global)),1);
        q_it = zeros(length(r(MESH.Solid.Gamma_global)),1);
        model = NIFC(p_it, q_it);   
    case 'broydensecond'
        Acc_flag = 5;
        Couple.reuse = 0;                   
        model = Broyden(Couple.reuse);
    case 'genbroyden'
        Acc_flag = 6;
        Couple.reuse = 0;
        model = GB(Couple.small, Couple.reuse, Couple.count, problemString); %GB3
    case 'andersonacceleration'
        Acc_flag = 7;
        Couple.lindepen = 1e-4;     %| filtering threshold
        model = AA(Couple.lindepen, Couple.reuse, Couple.filter, Couple.count, problemString);
    case 'mvls'
        Acc_flag = 8;
        limit = 20;
        model = MVLS(Couple.small, limit, Couple.count, problemString);
    case 'mvlss'
        Acc_flag = 9;
        limit = 20;
        model = MVLSS(Couple.small, limit, Couple.count, problemString);
    otherwise
        error('Unknown algorithm');
end

%% Additional functional tools (step's inital variable approx. / convergence criteria)
extrapolator = Extrapolator(zeros(length(MESH.Solid.Gamma_global),1), Couple.extra);
convergence = Convergence(Couple.rtol,Couple.imin,Couple.imax,Couple.maxSteps,Couple.small,Couple.out, Acc_flag, problemString);
objects = {convergence,extrapolator,model};
