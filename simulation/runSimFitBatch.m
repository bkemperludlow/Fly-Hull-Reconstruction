% -------------------------------------------------------------------------
% script to run simulation fit to data on multiple perturbation movies
%
% look through aggregated data folder ("b1 paper data") and find all data
% folders with a controller fit in them. run sumulation fits on data in all
% such folders, and save output to those folders. a later script can
% compile the simulation results
% -------------------------------------------------------------------------
%% data path information
rootPath = 'D:\b1 paper data\' ;
mn_name = 'b1' ; % 'b1' | 'b2'
effector =  'new kir' ; %'UAS-GtACR1'  ; % 'new kir' | '5XUAS-DSCP-eGFPKir2.1'
pertType = 'Pitch' ; % 'Pitch', 'Roll', or 'Yaw'
KFlag = false ;
overWriteFlag = false ;

% folder where genotype data is stored
dataPath = fullfile(rootPath, mn_name, effector) ;

% duration of magnetic pulse (depends on what type of experiment we're
% looking at)
switch effector
    case 'UAS-GtACR1'
        pulseDuration = 15e-3 ;
    case {'new kir', 'weak kir'}
        pulseDuration = 7e-3 ; 
    otherwise
        fprintf('Invalid effector: %s \n', effector)
        keyboard
end

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

% vector to scale each of 4 fit parameters by, so that they share similar
% order of magnitude. since i'm setting the finite difference step size,
% seems important to have consistent length scales for dimensions in param
% space. but maybe not?
paramScale = [1e-2, 1, 1, 1e-6] ;

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

% ---------------------------------------------------------------------
%% simulation conditions for testing control parameters from curve fit
% (if initial guess doesn't converge, solver will fail)
t_int_test = [-0.15, 0.055] ;
s_0_test = [0, 0, 0, 0, pi/4, 0] ;
thetaB0_test = pi/4 ;
params_test = defineQuasiSteadyParams() ; 

% -------------------------------------------------------------
%% get directory structure for controller fit folders
% get search string for controller fit structs
baseStr = 'controller_fit_struct' ;
if KFlag && strcmp(pertType, 'Pitch')
    suffixStr = '_LM' ;
elseif ~KFlag && strcmp(pertType, 'Pitch') && ~strcmp(effector,'UAS-GtARC1')
    suffixStr = '_noK_new' ;
elseif ~KFlag && strcmp(pertType, 'Pitch') && strcmp(effector,'UAS-GtARC1')
    suffixStr = '_noK' ;
elseif strcmp(pertType,'Roll')
    suffixStr = '' ;
end
searchStr = [baseStr suffixStr '*'] ;

% exclude file paths that contain these strings
exclude_strs = {'old','combined','backup'} ;

% ---------------------------------------------------
% get all controller fits in folder
controllerDir = dir(fullfile(dataPath,'**',searchStr)) ;
pertIdx = arrayfun(@(x) contains(x.folder, pertType,'IgnoreCase',true),...
    controllerDir) ;
controllerDir = controllerDir(pertIdx) ;

% exclude controller fits in an "old" folder or that's already in a
% combined folder
exclude_idx = false(length(controllerDir),1) ;
for s_num = 1:length(exclude_strs)
    str_curr = exclude_strs{s_num} ;
    exclude_idx_curr = arrayfun(@(x) contains(x.folder, str_curr,...
        'IgnoreCase',true), controllerDir) ;
    exclude_idx = exclude_idx | exclude_idx_curr ;
end
controllerDir = controllerDir(~exclude_idx) ;

N_fits = length(controllerDir) ;

% --------------------------------------------------
%% start parallel pool?
if parallelFlag  % && ~strcmp(solverType, 'simulatedannealing')
    if max(size(gcp)) == 0 % parallel pool needed
        parpool % create the parallel pool
    end
end

% --------------------------------------------------------------
%% loop through controller fits and run simulation fit
% take controller curve fit guess as one of the starting points
% --------------------------------------------------------------
for ind = 1:N_fits
    % ---------------------------------------------------------
    %% load data
    % get path and filename for current controller fit
    folder = controllerDir(ind).folder ;
    fn = controllerDir(ind).name ;
    folder_split = strsplit(folder,'\') ;
    
    % check if simulation fit already exists
    if saveFlag
       saveNameFull = fullfile(folder, 'sim_data_struct.mat') ;
       if exist(saveNameFull,'file') && ~overWriteFlag
          fprintf('Already fit simulation to %s -- skipping \n', ...
              folder_split{end})
          continue
       end
    end
    fprintf('Loading data for %s ... \n', folder_split{end})
    
    % load current controller fit
    controller_fit_struct = importdata(fullfile(folder, fn)) ;
    
    % ------------------------------------------------------------
    % get controller parameters to use as intial guess for sim fit
    x0 = [controller_fit_struct.K_i ,...
        controller_fit_struct.K_p,...
        controller_fit_struct.deltaT,...
        nan] ;
  
    % -------------------------------------------------------------
    % test initial guess to see if it returns valid dde simulation
    sol_test = simulatePitch_optim(x0(1:3), t_int_test, s_0_test, ...
        thetaB0_test, params_test, true, true) ; 
    
    % if we don't get a valid solution, alter initial guess
    if max(sol_test.x) < t_int_test(end)
       x0 = [0.0, params_test.K_p_sim, params_test.deltaT_sim, nan] ;   
    end
    
    % ------------------------------------------------------
    % get other identifying info from controller fit struct
    MovNum = controller_fit_struct.MovNum ;
    ExprNum = controller_fit_struct.ExprNum ;
    if isfield(controller_fit_struct, 'effector') && ...
            isfield(controller_fit_struct, 'driver')
        effector = controller_fit_struct.effector ;
        driver = controller_fit_struct.driver ; 
    else
        genotype_str = folder_split{end-2} ;
        genotype_str_split = strsplit(genotype_str,'_') ; 
        effector = genotype_str_split{1} ; 
        driver = genotype_str_split{2} ; 
    end
    dataRootPath = controller_fit_struct.rootPath ;
    
    % ------------------------------------------------------------
    % load data struct for this movie (output of hullAnalysis)
    % NB: we can do this from either the "Fly Data" directory or the current
    % directory in the aggregated data folder since we just need to read out
    % time and body pitch angle -- these should be the same throughout
    
    pathStruct = generatePathStruct(dataRootPath) ;
    [data, ~, ~, ~, ~] = hierarchicalLoadData(pathStruct, MovNum) ;
    fprintf('Data loaded \n')
    
    % ----------------------------------------------------------
    %% run fit
    tic % time fit run
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
    
    % add fields to sim_data_struct
    sim_data_struct.computation_time = t_comp ; % computation time (seconds)
    sim_data_struct.MovNum = MovNum ;
    sim_data_struct.ExprNum = ExprNum ;
    
    sim_data_struct.effector = effector ;
    sim_data_struct.driver = driver ;
    sim_data_struct.rootPath = dataRootPath ;
    
    % -------------------------------------------------------------
    %% save results?
    if saveFlag
        saveNameFull = fullfile(folder, 'sim_data_struct.mat') ;
        save(saveNameFull, 'sim_data_struct')
        
        if plotFlag
            savefig(gcf, fullfile(folder, 'sim_fit_plot.fig')) 
            print(gcf, fullfile(folder, 'sim_fit_plot.png'),'-dpng','-r300') 
        end
    end
    
    % close any open figure windows
    close all
end