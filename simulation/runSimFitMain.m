% -------------------------------------------------------------------------
% main script to run simulation fit on example data set
% -------------------------------------------------------------------------
% -------------------------
%% which data file to load
% can be altered, but default is 'fly_data.mat' in root subdir 'data'
rootPath = 'D:\Fly Data\VNC Motor Lines\94_11102019\' ; % 'D:\Fly Data\VNC Motor Lines\95_13102019\' ; %'D:\Fly Data\Janelia Flies\kir round 2\11_17082016\' ; %  'D:\Fly Data\VNC Motor Lines\54_21022019\' ; %
MovNum = 34 ; %2 ; % 43 ; %

normFlag = false ; % normalize data to perturbation amplitude?

% duration of magnetic pulse
pulseDuration = 7e-3 ;
% -------------------------
%% solver/output options
% first, choose solver method: UPDATE: DEFAULT TO LSQNONLIN_UNBND
% % 'fmincon' | 'lsqnonlin_bnd' | 'lsqnonlin_unbnd' 
% NB: 'patternsearch' and 'simulatedannealing' are too slow, so not
% included
solverType = 'lsqnonlin_unbnd' ; %'lsqnonlin_unbnd' ;  

% other solver options:
numStarts = 12 ;  % number of initializations to use. if >1, use multiStart 
optimFlag = true ;  % use optimized code?
parallelFlag = true ;  % use parallel computing?

% plot/save options:
plotFlag = true ;  % plot output results?
saveFlag = true ;  % save output results?

% % stop varying some parameters?
% varyFlags = [1, 1, 1, 1] ;

% get experiment number from folder name--too lazy to re-enter each time
pathStruct = generatePathStruct(rootPath) ;
pathSplit = strsplit(rootPath,'\') ;
folderSplit = strsplit(pathSplit{end-1},'_') ;
ExprNum = str2double(folderSplit{1}) ;

% smooth data?
smoothFlag = true ; 

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

% ---------------------------------------
%% import fly data structure
fprintf('Loading data for Expr %d, Mov %03d ... \n', ExprNum, MovNum)
% NB: data_fn will have weird suffix -- we're ignoring this
[data, ~, data_fn, ~, ~] = hierarchicalLoadData(pathStruct, MovNum) ;
[dataPath, ~, ~] = fileparts(data_fn) ;

% --------------------------
% normalize pitch data?
if normFlag
    % read out pitch time series
    defineConstantsScript ; 
    frames = (data.params.startTrackingTime : data.params.endTrackingTime) ;
    t = (1/data.params.fps).*frames ; 
    bodyPitch = data.anglesLabFrame(:, THETAB) ; 
    
    % subtract off initial value of pert
    [~, t0_ind] = min(abs(t)) ;  
    bodyPitch_initSub = bodyPitch - bodyPitch(t0_ind) ; 
    
    % find pert amp
    t_pert_ind = find((t > pulseDuration) & (t < 0.05)) ;
    [pert_amp, ~] = findpeaks(abs(bodyPitch_initSub(t_pert_ind)), ...
        t(t_pert_ind), 'Npeaks', 1) ;
    if isempty(pert_amp)
        [pert_amp, ~] = nanmax(abs(bodyPitch_initSub(t_pert_ind))) ;
    end
    
    % divide body pitch angle by pert amp
    bodyPitch_initSub_scale = (1/pert_amp).*bodyPitch_initSub ; 
    
    % add back initial value
    bodyPitch_scale = bodyPitch_initSub_scale + bodyPitch(t0_ind) ; 
    
    % insert back into data struct
    data.anglesLabFrame(:,THETAB) = bodyPitch_scale ; 
end
fprintf('Data loaded \n')
% --------------------------------------
%% params
fprintf('Guessing initial parameters ... \n')
% vector to scale each of 4 fit parameters by, so that they share similar
% order of magnitude. since i'm setting the finite difference step size,
% seems important to have consistent length scales for dimensions in param
% space. but maybe not?
paramScale = [1e-2, 1, 1, 1e-6] ;

% initial fit params guess (obtained from cure fit)
% controller_fit_struct = fitPitchController_dt(rootPath, MovNum, false, ...
%     false, false, false, false, false, false, false, pulseDuration) ;
% K_i = controller_fit_struct.K_i  ;        % integral gain
% K_p = controller_fit_struct.K_p ;        % prop gain
% deltaT = controller_fit_struct.deltaT ;   % time delay
K_i = 0.437556 ;
K_p = 0.005728 ;
deltaT = 0.00429 ;

% NB: leaving initial pulseStrength guess as nan since we'll estimate
% it later from data
x0 = [K_i, K_p, deltaT, nan] ;

fprintf('Initial parameters obtained \n')
% -------------------------------------
%% set bounds on fit parameters
% (controller gains and pulse amplitude)
K_i_min = 0.0 ;  % min/max for INTEGRAL gain coefficient. unitless
K_i_max = 2.0 ; %2.0 ;

K_p_min = 0 ; % min/max for PROPORTIONAL gain coefficient. seconds
K_p_max = 0.5 ;

deltaT_min = 2e-3 ; % min/max for time delay. seconds
deltaT_max = 30e-3 ;

pulseStrength_min = -1e6; % min/max for pulse amplitude. rad/s^2
pulseStrength_max =  1e6;

% vector of lower and upper bounds on fit params
lb = [K_i_min, K_p_min, deltaT_min, pulseStrength_min] ;  % LOWER bounds on controller terms
ub = [K_i_max, K_p_max, deltaT_max, pulseStrength_max] ;  % UPPER bounds on controller terms

% ------------------------------------------
%% start parallel pool?
if parallelFlag  % && ~strcmp(solverType, 'simulatedannealing')
    if max(size(gcp)) == 0 % parallel pool needed
        parpool % create the parallel pool
    end
end
% ---------------------------------------
%% run simulation fit

% time fit run
tic
if optimFlag
    fprintf('Running minimization with "optimized" functions ... \n')
    sim_data_struct = fitSimToData_optim(data, x0, lb, ub, pulseDuration, ...
        paramScale, solverType, numStarts, parallelFlag, plotFlag, ...
        smoothFlag, delayFunType, timeLimSec) ;
else
    fprintf('Running minimization ... \n')
    
    sim_data_struct = fitSimToData(data, x0, lb, ub, pulseDuration, ...
        paramScale, solverType, parallelFlag, plotFlag) ;
end
t_comp = toc ;
fprintf('Completed minmization \n')
fprintf('Elapsed time: %f \n', t_comp)

% add computation time to sim_data_struct
sim_data_struct.computation_time = t_comp ; % seconds
% ----------------------------------------
%% save results?
if saveFlag
    if normFlag
       suffixStr = 'norm_' ;
    else
        suffixStr = '' ; 
    end
    t_now = now ; 
    save(fullfile(dataPath, ...
        ['sim_data_struct_' suffixStr num2str(t_now) '.mat']), ...
        'sim_data_struct')
end