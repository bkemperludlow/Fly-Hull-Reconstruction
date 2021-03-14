% -------------------------------------------------------------------------
% function to fit simulation controller parameters to measured data
% -------------------------------------------------------------------------
function sim_data_struct = fitSimToData_optim(data, x0, lb, ub, pulseDuration, ...
    paramScale, solverType, numStarts, parallelFlag, plotFlag, smoothFlag,...
    delayFunType, timeLimSec)
% ------------------------------
%% inputs and params
if ~exist('x0','var') || isempty(x0)
    % initial guess for fit params -- should try to get these from other
    % controller fit, but being lazy atm
    K_i = 0.437556 ;        % integral gain
    K_p = 0.005728 ;        % prop gain
    deltaT = 0.004296 ;   %time delay
    
    % NB: leaving initial pulseStrength guess as nan since we'll estimate
    % it later from data
    x0 = [K_i, K_p, deltaT, nan] ;
end
if ~exist('lb','var') || isempty(lb)
    lb = [0, 0, 0.002, -1e6] ;  % LOWER bounds on controller terms
end
if ~exist('ub','var') || isempty(ub)
    ub = [2.0, 0.5 , 0.03, 1e6] ;  % UPPER bounds on controller terms
end
if ~exist('pulseDuration','var') || isempty(pulseDuration)
    pulseDuration = 0.007 ;  % seconds
end
if ~exist('paramScale','var') || isempty(paramScale)
    paramScale = [1e-2, 1, 1, 1e-6] ;
end
if ~exist('numStarts','var') || isempty(numStarts)
    % solver type (either fmincon or lsqnonlin -- others too slow)
    numStarts = 1 ;  % if 1, solve normally. if >1, use MultiStart
end
if ~exist('parallelFlag','var') || isempty(parallelFlag)
    parallelFlag = false ;
end
if ~exist('parallelFlag','var') || isempty(parallelFlag)
    parallelFlag = false ;
end
if ~exist('plotFlag','var') || isempty(plotFlag)
    plotFlag = true ;
end
if ~exist('smoothFlag','var') || isempty(smoothFlag)
    smoothFlag = true ;
end
if ~exist('delayFunType','var') || isempty(delayFunType)
    delayFunType = 'constant' ;
end
if ~exist('timeLimSec','var') || isempty(timeLimSec)
    timeLimSec = 80 ;
end

% level of noise which is added to x0 to generate MultiStart points
ms_noise_lvl = 0.0025 ;

% number of fit params
nvars = length(x0) ;

% quasi-steady parameters
% NB: in future version, should try to get some morphological params from
% data
params = defineQuasiSteadyParams ;
omega = params.omega ; 
% ------------------------------------------------
%% read from data struct
% get smoothed body euler angles (in radians)
if smoothFlag
    [bodyPitch, bodyYaw, bodyRoll] = smoothBodyAngles(data) ;
else
    defineConstantsScript
    bodyPitch = data.anglesLabFrame(:,THETAB) ; 
    bodyYaw = data.anglesLabFrame(:, PHIB) ; 
    bodyRoll = data.anglesLabFrame(:, RHO) ; 
end
bodyYPR = (pi/180).*[bodyYaw, bodyPitch, bodyRoll] ;

% get time from movie data
dt = (1/data.params.fps) ;
t = dt.*(data.params.startTrackingTime : data.params.endTrackingTime) ;

if (t(1) > 0)
    fprintf('Invalid time range \n')
    keyboard
end

% define common time interval over which to compare simulation and data
ti = 0 ; % sec.
tf = min([t(end), 0.040]) ; % sec. either take 40 ms or time limit
t_ind = (t >= ti) & (t <= tf) ;

% use this to define DE solver range and to get data in correct format
% (left and right ends incorporates extra time to let simulated fly get to
% stable state)
pad_left = -0.15 ; % left and right interval time padding
pad_right = 0.015 ; 

% start point
t_int1 = (ti + pad_left) ;  
t_int1 = t_int1 + (-1*rem(omega*t_int1, 2*pi) + pi/2)/omega; % ensure we start with sin(t_int1) = 1 

% end point
t_int2 = tf + pad_right ; 

% interval(s)
t_int = [t_int1, t_int2] ;  % t_int = [-0.15, tf + 0.050] ; 
t_comp = [0, tf] ;
thetaB_data = bodyYPR(t_ind, 2) ;

% subtract off value of pitch angle at pert start
thetaB_data = thetaB_data - thetaB_data(1) ;

% % [shouldn't actually do] detrend data (same with simulated results)
% thetaB_data = detrend(thetaB_data) ;

% ------------------------------------------------------------------------
% for debugging
% plot(linspace(t_comp(1), t_comp(end), length(thetaB_data)), thetaB_data) 
% sim_data_struct = [] ; 
% return

% ------------------------------------------------------------------------
%% maximum time for local solvers in multistart (your mileage may vary)
sim_time = diff(t_int) ; 
sim_clock_time = 6000*sim_time ;
if parallelFlag
    ms_max_time = sim_clock_time*ceil(numStarts/5) ;
else
    ms_max_time = sim_clock_time*numStarts ;
end

% -----------------------------------------------------------------
%% set DE initial conditions
% get pre-perturbation pitch angle -- make that the set point
pre_pert_t_min = max([-0.01, t(1)]) ;
pre_pert_ind = (t < 0) & (t > pre_pert_t_min) ;
thetaB0 = mean(bodyYPR(pre_pert_ind,2)) ;

% % subtract off pre-pert mean from data
% thetaB_data = thetaB_data - thetaB0 ;

% since the fly is unstable to start, we can't match up the initial state
% vector between data and simulation (unless we add controllers for the
% other terms). so instead just start in standard initial conditions
s_0 = [0, 0, 0, 0, params.beta_0, 0] ;

% ----------------------------------------------------------------
%% estimate perturbation strength from data
params.pulseStart = 0 ;
params.pulseEnd = pulseDuration ;

if isnan(x0(4))
    % differentiate euler angles to get pulse strength
    pitchAccel = (1/dt)^2 .* [0; 0; diff(bodyYPR(:,2),2)];
    
    % pulse strength determined by peak of pitch acceleration
    pert_ind = (t >= 0) & (t <= pulseDuration) ;
    pitchAccelPert = pitchAccel(pert_ind) ;
    [~, locs] = findpeaks(pitchAccelPert) ;
    if ~isempty(locs)
        pulseStrength = pitchAccelPert(locs(1)) ;
    else
        [~, max_ind] = max(abs(pitchAccelPert)) ;
        pulseStrength = pitchAccelPert(max_ind) ;
    end
    
    x0(4) = pulseStrength ;
end

% -----------------------------------------------------------------
%% low pass filter to smooth simulation results (remove wb ringing)
d_filt = designfilt('lowpassiir','FilterOrder',3,'SampleRate', (1/dt), ...
    'HalfPowerFrequency',params.omega/(4*pi),'DesignMethod','butter') ;

% ----------------------------------------------------------------
%% scale initial guess (x0) and bounds (lb, ub)
% apply parameter scaling (if all params are same order of magnitude,
% should help solver? esp. since i'm setting finite difference step size?
x0 = paramScale.*x0 ;

% also apply to bounds
lb = paramScale.*lb ;
ub = paramScale.*ub ;

% ------------------------------------------------------------
%% define anonymous function for objective and solver options
switch solverType
    case 'fmincon'
        % options for fmincon
        options =  optimoptions(@fmincon, 'Display', 'iter', ...
            'FiniteDifferenceStepSize', 1e-5,...
            'FiniteDifferenceType','central', ...
            'UseParallel', parallelFlag) ;
        
        % objective function for fmincon
        fun = @(x) fitSimToDataCost_optim(x, t_int, t_comp, s_0, thetaB0, ...
            thetaB_data, params, d_filt, paramScale, delayFunType, timeLimSec) ;
        
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
        fun = @(x) fitSimToDataLSQ_optim(x, t_int, t_comp, s_0, thetaB0, ...
            thetaB_data, params, d_filt, paramScale, delayFunType, timeLimSec) ;
        
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
        fun = @(x) fitSimToDataLSQ_optim(x, t_int, t_comp, s_0, thetaB0, ...
            thetaB_data, params, d_filt, paramScale, delayFunType, timeLimSec) ;
        
        % set up optimization problem
        problem = createOptimProblem('lsqnonlin', 'x0', x0,...
            'objective', fun, 'lb', lb, 'ub', ub, 'options', options) ;
        
    otherwise
        fprintf('Invalid solver type: %s \n', solverType)
        sim_data_struct = nan ;
        return
end

% ------------------------------------------------
%% run MultiStart?
% if using MultiStart, do the following:
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
    % 1) the user provided initial guess
    % 2) the default sim values, as stored in the quasisteady params struct
    % 3) default sim values + gaussian noise
    % NB: all have same initial guess for pulse strength, as measured from
    % kinematics
    % ---------------------------------------------------------------------
    % first load in "default" params from qs params, and scale
    x0_default = [params.K_i_sim, params.K_p_sim, params.deltaT_sim, ...
        pulseStrength] ;
    x0_default = paramScale.*x0_default ;
    
    % pt_matrix, whose rows are the initial values for each start, contains
    % the user input and replicates of the default params
    pt_matrix = vertcat(x0, repmat(x0_default, numStarts-1,1)) ;
    
    % add noise to the pt_matrix, depending on situation
    if all(x0 == x0_default)
        % if x0 and x0_default are the same, add noise to all rows but 1st
        x0_noise = [ms_noise_lvl.*randn(numStarts-1, nvars-1), ...
            zeros(numStarts-1,1)] ;
        x0_noise = vertcat(zeros(1, nvars), x0_noise) ;
    elseif ~all(x0 == x0_default) && (numStarts == 2)
        % if x0 and x0 default are NOT the same but there are only 2 start,
        % don't add any noise
        x0_noise = zeros(numStarts, nvars) ;
    else
        % in this case, x0 and x0_default are different, and there are >2
        % starts. so we add noise rows 3-end of pt_matrix
        x0_noise = [ms_noise_lvl.*randn(numStarts-2, nvars-1), ...
            zeros(numStarts-2,1)] ;
        x0_noise = vertcat(zeros(2, nvars), x0_noise) ;
    end
    
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

% read out/rescale parameters back to original units
controllerTerms = (paramScale.^(-1)).*x ;

% print output
fprintf('K_i = %f \n', controllerTerms(1))
fprintf('K_p = %f \n', controllerTerms(2))
fprintf('deltaT = %f \n', controllerTerms(3))
fprintf('pulseStrength = %f \n', controllerTerms(4))

fprintf('\n SSE = %f \n', fval)
% ---------------------------------------------------
%% get kinematic data for solution
% re-run solution with optimized parameters
params.pulseStrength = controllerTerms(4) ;
sol = simulatePitch_optim(controllerTerms(1:3), t_int, s_0, thetaB0, ...
    params, true, true, delayFunType, nan) ;

% evaluate solution at time grid
dt = 0.125e-3 ;
t_eval = t_int(1) : dt : t_int(end) ;
sint = deval(sol, t_eval) ; %evaluate solution over full range

% read out body pitch
thetaB_sim = sint(5,:) ; %radians

% filter and detrend solver output (like subtracting pre-pert value)
thetaB_sim_filt = filtfilt(d_filt, thetaB_sim) ;
thetaB_sim_detrend = detrend(thetaB_sim_filt) ;

% restrict simulation data to just range of data
[~, ind1] = min(abs(t_eval - t_comp(1))) ;
[~, ind2] = min(abs(t_eval - t_comp(end))) ;
%     comp_ind = (t_eval >= t_comp(1)) & (t_eval <= t_comp(end)) ;
comp_ind = ind1:ind2 ;
thetaB_sim_final = thetaB_sim_detrend(comp_ind) ;
thetaB_sim_final = thetaB_sim_final - thetaB_sim_final(1) ;

% ---------------------------------------------------
%% compile info in structure for output
sim_data_struct = struct() ;

% controller coefficients
sim_data_struct.K_i = controllerTerms(1) ;
sim_data_struct.K_p = controllerTerms(2) ;
sim_data_struct.deltaT = controllerTerms(3) ;

% cost/SSE
sim_data_struct.fval = fval ;

% pulse data
sim_data_struct.pulseStrength = controllerTerms(4) ;
sim_data_struct.pulseStart = params.pulseStart ;
sim_data_struct.pulseEnd = params.pulseEnd ;

% kinematic data (both sim and data)
sim_data_struct.thetaB_sim = thetaB_sim_final ; % filteres/detrended pitch (sim)
sim_data_struct.thetaB_data = thetaB_data ; % filtered data pitch
sim_data_struct.t = t_eval(comp_ind) ;      % time range of data
sim_data_struct.t_eval = t_eval ;           % time range for full sim output
sim_data_struct.thetaB_sim_raw = thetaB_sim ;  % raw output of body angle from sim
sim_data_struct.sol = sol ;  % full simulation output struct
sim_data_struct.thetaB0 = thetaB0 ;  % controller angle set point

% initial guess, DE initial conditions, parameter bounds, and scaling
sim_data_struct.lb = (paramScale.^(-1)).*lb ;
sim_data_struct.ub = (paramScale.^(-1)).*ub ;
sim_data_struct.x0 = (paramScale.^(-1)).*x0 ;
sim_data_struct.s_0 = s_0 ;
sim_data_struct.paramScale = paramScale ;

% params and solver options
sim_data_struct.params = params ;
sim_data_struct.options = options ;
sim_data_struct.d_filt = d_filt ;   % filter used to get rid of wbf fluctuations
% sim_data_struct.solverType = solverType ;

% additional solver outputs, if they exist
solver_out_list = {'exitflag', 'output', 'lambda', 'grad', 'hessian', ...
    'residual', 'jacobian', 'ms_solutions', 'ms_x', 'ms_fval', ...
    'ms_exitflag', 'ms_output'} ;

% going to loop over possible solver outputs and assign them to structure
% fields if they exist. i think evaluating the strings to be a variable
% call is bad practice, but this is a pretty small application
for n = 1:length(solver_out_list)
    % current variable name
    solver_out = solver_out_list{n} ;
    
    % if current variable exists, add it to struct
    if exist(solver_out,'var')
        sim_data_struct.(solver_out) = eval(solver_out) ;
    end
end

% -----------------------------------------------------
%% plot results?
if plotFlag
    % make figure
    figure ;
    hold on
    % data
    plot(t_eval(comp_ind), thetaB_data)
    % sim results
    plot(t_eval(comp_ind), thetaB_sim_final, '--')
    
    % axis properties
    xlabel('Time (s)')
    ylabel('Body Pitch Angle (rad)')
    legend({'data','sim'})
    
end

end