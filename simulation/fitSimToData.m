% -------------------------------------------------------------------------
% function to fit simulation controller parameters to measured data
% -------------------------------------------------------------------------
function sim_data_struct = fitSimToData(data, x0, lb, ub, pulseDuration, ...
    paramScale, solverType, parallelFlag, plotFlag)
% ------------------------------
%% inputs and params
if ~exist('x0','var') || isempty(x0)
    % initial guess for fit params -- should try to get these from other
    % controller fit, but being lazy atm
    K_i = 0.437556 ;        % integral gain
    K_p = 0.005728 ;        % prop gain
    deltaT = 2*0.002148 ;   %time delay
    
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
if ~exist('solverType','var') || isempty(solverType)
    solverType = 'fmincon' ;
end
if ~exist('parallelFlag','var') || isempty(parallelFlag)
    parallelFlag = false ;
end
if ~exist('plotFlag','var') || isempty(plotFlag)
    plotFlag = true ;
end

% quasi-steady parameters
% NB: in future version, should try to get some morphological params from
% data
params = defineQuasiSteadyParams ;

% ------------------------------------------------
%% read from data struct
% get smoothed body euler angles (in radians)
[bodyPitch, bodyYaw, bodyRoll] = smoothBodyAngles(data) ;
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
t_int = [-0.15, tf + 0.050] ; % left and right ends incorporates extra time to let simulated fly get to stable state
t_comp = [0, tf] ;
thetaB_data = bodyYPR(t_ind, 2) ;

% subtract off value of pitch angle at pert start
thetaB_data = thetaB_data - thetaB_data(1) ;

% % detrend data (will do same with simulated results)
% thetaB_data = detrend(thetaB_data) ;
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

% ------------------------------------------------
%% run solver
% define anonymous function for objective
fun = @(x) fitSimToDataCost(x, t_int, t_comp, s_0, thetaB0, ...
    thetaB_data, params, d_filt, paramScale) ;

% set options and run solver (depending on selection)
switch solverType
    case 'fmincon'
        % fmincon options
        options =  optimoptions(@fmincon,'Display', 'iter',...
            'FiniteDifferenceStepSize', 1e-5,...
            'FiniteDifferenceType','central', ...
            'UseParallel', parallelFlag) ; % 'Algorithm', 'sqp',
        
        % fmincon solver
        [controllerTerms, fval] = fmincon(fun, x0, [], [], [], [], ...
            lb, ub, [], options) ;
    case 'patternsearch'
        % patternsearch options
        options =  optimoptions(@patternsearch,'Display', 'iter',...
            'Cache', 'on',...
            'UseParallel', parallelFlag) ;
        
        % patternsearch solver
        [controllerTerms, fval] = patternsearch(fun, x0,[], [], [], [], ...
            lb, ub, [], options) ;
    case 'simulatedannealing'
        % simulated annealing option
        options =  optimoptions(@simulannealbnd,'Display', 'iter') ;
        
        % simulated annealing solver
        [controllerTerms, fval] = simulannealbnd(fun, x0, lb, ub, options);
    case 'lsqnonlin_bnd'
        % if using least-squares, need to use different cost function
        fun = @(x) fitSimToDataLSQ(x, t_int, t_comp, s_0, thetaB0, ...
            thetaB_data, params, d_filt, paramScale) ;
        
        % options for BOUNDED nonlinear least squares
        options =  optimoptions(@lsqnonlin, 'Display', 'iter', ...
            'FiniteDifferenceStepSize', 1e-5,...
            'FiniteDifferenceType','central', ...
            'UseParallel', parallelFlag) ;
        
        % non linear least squares solver with bounds
        [controllerTerms, fval,residual, exitflag, output, lambda, ...
            jacobian] = lsqnonlin(fun, x0, lb, ub, options) ;
    case 'lsqnonlin_unbnd'
        % if using least-squares, need to use different cost function
        fun = @(x) fitSimToDataLSQ(x, t_int, t_comp, s_0, thetaB0, ...
            thetaB_data, params, d_filt, paramScale) ;
        
        % options for BOUNDED nonlinear least squares
        options =  optimoptions(@lsqnonlin, 'Display', 'iter',...
            'Algorithm', 'levenberg-marquardt', ...
            'FiniteDifferenceStepSize', 1e-5,...
            'FiniteDifferenceType','central', ...
            'UseParallel', parallelFlag) ;
        
        % non linear least squares solver with bounds
        [controllerTerms, fval, residual, exitflag, output, lambda, ...
            jacobian] = lsqnonlin(fun, x0, [], [], options) ;
    otherwise
        fprintf('Invalid solver type: %s \n', solverType)
        sim_data_struct = struct() ;
        return
end

% rescale parameters back to original values
controllerTerms = (paramScale.^(-1)).*controllerTerms ;

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
sol = simulatePitch(controllerTerms(1:3), t_int, s_0, thetaB0, params, ...
    true, true, false) ;

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
comp_ind = (t_eval >= t_comp(1)) & (t_eval <= t_comp(end)) ;
thetaB_sim_final = thetaB_sim_detrend(comp_ind) ;
thetaB_sim_final = thetaB_sim_final - thetaB_sim_final(1) ;

% ---------------------------------------------------
%% compile info in structure for output
sim_data_struct = struct() ;

% controller coefficients
sim_data_struct.K_i = controllerTerms(1) ;
sim_data_struct.K_p = controllerTerms(2) ;
sim_data_struct.deltaT = controllerTerms(3) ;

% pulse data
sim_data_struct.pulseStrength = controllerTerms(4) ;
sim_data_struct.pulseStart = params.pulseStart ;
sim_data_struct.pulseEnd = params.pulseEnd ;

% Sum of squared error
sim_data_struct.resnorm = fval ;

% kinematic data (both sim and data)
sim_data_struct.thetaB_sim = thetaB_sim_final ; % filteres/detrended pitch (sim)
sim_data_struct.thetaB_data = thetaB_data ; % filtered data pitch
sim_data_struct.t = t_eval(comp_ind) ;      % time range of data
sim_data_struct.thetaB_sim_raw = thetaB_sim ;  % raw output of body angle from sim
sim_data_struct.comp_ind = comp_ind ;  % indices over which raw sim output can be compared to real data
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
sim_data_struct.solverType = solverType ;

% additional solver outputs, if they exist
if exist('jacobian','var')
    sim_data_struct.residual = residual ;
    sim_data_struct.exitflag = exitflag ;
    sim_data_struct.output = output ;
    sim_data_struct.lambda = lambda ;
    sim_data_struct.jacobian = jacobian ;
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