% -------------------------------------------------------------------------
% script to optimize controller terms for simulated flight
% -------------------------------------------------------------------------
% -------------------------------------------------
%% misc params
plotFlag = true ;
pertFlag = false ; % external perturbation or not
controlFlag = true ; % has to be true to test controller functions

solverType = 'lsqnonlin_unbnd' ;  % 'fmincon' | 'lsqnonlin_bnd' | 'lsqnonlin_unbnd'

runSolverFlag = true ;
parallelFlag = true ; % us parallel computing?

numStarts = 3  ; % number of starting points for MultiStart

% ---------------------------------------
%% delay type and time limits
% type of delay function to use:
% 'constant' = just deltaT (updates at each time point)
% 'average'  = average over previous wingbeat, also different at every time point
% 'function' = use "myControllerDelay.m" to get delay. should give same
%              value for multiple successive time points
delayFunType = 'function' ;  % 'constant' | 'average' | 'function'

% limit time of simulation run?
timeLimSec = 80 ; % nan ; %80 ;

% -------------------------------------------------
%% set initial guesses for controller terms
% (from JEB pitch paper)
% K_i0 = 0.0 ; %0.3 ; % unitless
% K_p0 = 0.020 ;%0.0157 ; %0.007 ; % seconds
% deltaT0 = 0.003 ; %0.006 ; % seconds

K_i0 = 0.437556 ;        % integral gain
K_p0 = 0.009728 ;        % prop gain  % 0.005728 ;
deltaT0 = 0.004296 ;   %time delay

% K_i0 = 0.633928 ;
% K_p0 = 0.006656 ;
% deltaT0 = 0.004491 ;

% K_i = 1.272641 ;
% K_p = 0.013782 ;
% deltaT = 0.004565 ;

% compile initial guesses into vector
x0 = [K_i0, K_p0, deltaT0] ;

% scale parameters to make them of similar maginitude
paramScale = [1e-2, 1, 1] ;

x0 = paramScale.*x0 ;

% level of noise which is added to x0 to generate MultiStart points
ms_noise_lvl = 0.0025 ;

% number of fit params
nvars = length(x0) ;

% -------------------------------------------------
%% time range, initial state, etc. for simulations
% general QS params
params = defineQuasiSteadyParams ; 
omega = params.omega ; 

% % start point
% t_int1 = -0.15 ;  % seconds
% t_int1 = t_int1 + (-1*rem(omega*t_int1, 2*pi) + pi/2)/omega; % ensure we start with sin(t_int1) = 1
% % end point
% t_int2 = 0.200 ;  % seconds

%t_int1 = 1.1e-3 ; %1.1e-3 ;
t_int1 = -150e-3 ;  
% t_int1 = t_int1 + (-1*rem(omega*t_int1, 2*pi) + pi/2)/omega; % ensure we start with sin(t_int1) = 1 
t_int1 = t_int1 + (-1*rem(omega*t_int1, 2*pi) - pi/2)/omega; % ensure we start with sin(t_int1) = 1 
t_int2 = 300e-3 ;

% interval(s)
t_int = [t_int1, t_int2] ;  % t_int = [-0.15, tf + 0.050] ;
thetaB0 = pi/4 ; % body pitch angle setpoint

s_0 = [0, 0, 0, 0, thetaB0 , 0] ; % initial state for sim

% ------------------------------------------------------------------------
% maximum time for local solvers in multistart (your mileage may vary)
sim_time = diff(t_int) ; 
sim_clock_time = 6000*sim_time ;
if parallelFlag
    ms_max_time = sim_clock_time*ceil(numStarts/5) ;
else
    ms_max_time = sim_clock_time*numStarts ;
end
% -------------------------------------------------------------------------
%% should we run solver or not? (if not, we're just manually testing params)
if runSolverFlag
    % -------------------------------------------------
    %% bounds on controller terms
    K_i_min = 0.0 ;  % integral term
    K_i_max = 0.0 ; % 2.0 ;
    
    K_p_min = 0.0 ;  % proportional term
    K_p_max = 0.05 ;
    
    deltaT_min = 0.002 ; % time delay
    deltaT_max = 0.020 ;
    
    lb = [K_i_min, K_p_min, deltaT_min] ;
    ub = [K_i_max, K_p_max, deltaT_max] ;
    
    % scale bounds
    lb = paramScale.*lb ;
    ub = paramScale.*ub ;
    % ------------------------------------------------
    %% set up solver
    switch solverType
        case 'fmincon'
            % options for fmincon
            options =  optimoptions(@fmincon, 'Display', 'iter', ...
                'FiniteDifferenceStepSize', 1e-5,...
                'FiniteDifferenceType','central', ...
                'UseParallel', parallelFlag) ;
            
            % objective function for fmincon
            fun = @(x) controllerStabilityCost(x,t_int, s_0, ...
                thetaB0, params, controlFlag, pertFlag, paramScale, ...
                delayFunType, timeLimSec) ;
            
            % set up optimization problem
            problem = createOptimProblem('fmincon', 'x0', x0,...
                'objective', fun, 'lb', lb, 'ub', ub, 'options', options) ;
            
        case 'lsqnonlin_unbnd'
            % options for UNBOUNDED nonlinear least squares
            options =  optimoptions(@lsqnonlin, 'Display', 'iter',...
                'Algorithm', 'levenberg-marquardt', ...
                'FiniteDifferenceStepSize', 1e-4,... % 1e-5
                'FiniteDifferenceType','central', ...
                'UseParallel', parallelFlag, ...
                'MaxIterations', 30, ...
                'MaxFunctionEvaluations', 300) ;
            
            % objective function for UNBOUNDED nonlinear least squares
            fun = @(x) controllerStabilityLSQ(x,t_int, s_0, ...
                thetaB0, params, controlFlag, pertFlag, paramScale, ...
                delayFunType, timeLimSec) ;
            
            % set up optimization problem
            problem = createOptimProblem('lsqnonlin', 'x0', x0,...
                'objective', fun, 'lb', [], 'ub', [], 'options', options) ;
        case 'lsqnonlin_bnd'
            % options for BOUNDED nonlinear least squares
            options =  optimoptions(@lsqnonlin, 'Display', 'iter',...
                'FiniteDifferenceStepSize', 1e-5,...
                'FiniteDifferenceType','central', ...
                'UseParallel', parallelFlag, ...
                'MaxIterations', 30, ...
                'MaxFunctionEvaluations', 300) ;
            
            % objective function for BOUNDED nonlinear least squares
            fun = @(x) controllerStabilityLSQ(x,t_int, s_0, ...
                thetaB0, params, controlFlag, pertFlag, paramScale, ...
                delayFunType, timeLimSec) ;
            
            % set up optimization problem
            problem = createOptimProblem('lsqnonlin', 'x0', x0,...
                'objective', fun, 'lb', lb, 'ub', ub, 'options', options) ;
    end
    
    % --------------------------------------------------
    %% start parallel pool?
    if parallelFlag  % && ~strcmp(solverType, 'simulatedannealing')
        if max(size(gcp)) == 0 % parallel pool needed
            parpool % create the parallel pool
        end
    end
    
    % ------------------------------------------------------------
    %% run MultiStart or solver
    if (numStarts > 1)
        % turn off solver display and local parallel
        problem.options.Display = 'off' ;
        problem.options.UseParallel = false ;
        
        % generate MultiStart object
        ms = MultiStart('FunctionTolerance',2e-5,...
            'XTolerance',5e-5, ...
            'UseParallel', parallelFlag, ...
            'Display', 'iter', ...
            'MaxTime', ms_max_time) ;
        
        % ---------------------------------------------------------------------
        % run MultiStart with custom start points. custom start points contains
        
        % pt_matrix, whose rows are the initial values for each start, contains
        % the user input and replicates of the default params
        pt_matrix = repmat(x0, numStarts,1) ;
        
        % add noise to the pt_matrix, depending on situation
        x0_noise = [ms_noise_lvl.*randn(numStarts-1, nvars-1), ...
            zeros(numStarts-1,1)] ;
        x0_noise = vertcat(zeros(1, nvars), x0_noise) ;
        
        % add the noise that we just generated
        pt_matrix = pt_matrix + x0_noise ;
        
        startpts = CustomStartPointSet(pt_matrix) ;
        [ms_x, ms_fval, ms_exitflag, ms_output, ms_solutions] = ...
            run(ms, problem, startpts) ;
        
        % ---------------------------------------------------------------------
        % after multistart, run regular solver to get extra output info
        % ---------------------------------------------------------------------
        fprintf('Running solver %s with output of MultiStart \n', solverType)
        % reset local solver settings
        problem.options.Display = 'iter' ;
        problem.options.UseParallel = parallelFlag ;
        
        % plug in output of MultiStart to problem structure
        problem.x0 = ms_x ;
    end
    
    % ------------------------------------------------------------
    %% run solver
    % NB: if we haven't used MultiStart, this takes the initial conditions as
    % the user input x0. if we have used MultiStart, this takes the initial
    % consitions as the best output of the MultiStart run
    switch solverType
        case 'fmincon'
            [x, fval, exitflag, output, lambda, grad, hessian] = ...
                fmincon(problem) ;
            
        case {'lsqnonlin_bnd', 'lsqnonlin_unbnd'}
            [x, fval, residual, exitflag, output, lambda, jacobian] = ...
                lsqnonlin(problem) ;
    end
    
    controllerTerms = (paramScale.^(-1)).*x ;
else
    controllerTerms = (paramScale.^(-1)).*x0 ;
    fval = 0.0 ;
end
% ------------------------------------------------
%% print results
fprintf('K_i = %f \n', controllerTerms(1))
fprintf('K_p = %f \n', controllerTerms(2))
fprintf('deltaT = %f \n', controllerTerms(3))
fprintf('RMSE = %f \n', fval)

% ------------------------------------------------
%% plot results?
if plotFlag
    % run solution with optimized results
    sol = simulatePitch_optim(controllerTerms, t_int, s_0, thetaB0, params, ...
        controlFlag, pertFlag, delayFunType, timeLimSec) ;
    
    % evaluate solution
    dt = 0.125e-3 ;
    t_eval = t_int(1) : dt : t_int(end) ;
    sint = deval(sol, t_eval) ;
    thetaB = (180/pi).*sint(5,:) ; %radians
    thetaB_dot = (180/pi).*sint(6,:) ;  %rad/s
    
    % filter output
    d_filt = designfilt('lowpassiir','FilterOrder',3,'SampleRate', (1/dt), ...
        'HalfPowerFrequency',params.omega/(4*pi),'DesignMethod','butter') ;
    figure ;
    
    thetaB_filt = filtfilt(d_filt, thetaB) ;
    thetaB_dot_filt = filtfilt(d_filt, thetaB_dot) ;
    
    % ----------------------------
    % pitch angle plot
    subplot(2,1,1)
    hold on
    
    plot(t_eval, thetaB)
    plot(t_eval, thetaB_filt, 'r-', 'LineWidth', 1.5)
    
    % axis properties
    xlabel('Time (s)')
    ylabel('Body Pitch Angle (deg)') ;
    axis tight
    
    % ----------------------------
    % pitch velocity plot
    subplot(2,1,2)
    hold on
    
    plot(t_eval, thetaB_dot)
    plot(t_eval, thetaB_dot_filt, 'r-', 'LineWidth', 1.5)
    
    % axis properties
    xlabel('Time (s)')
    ylabel('Body Pitch Vel (deg/s)') ;
    axis tight
end