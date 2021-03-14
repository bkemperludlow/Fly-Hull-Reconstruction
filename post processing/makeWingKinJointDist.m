% -------------------------------------------------------------------------
% script to generate a joint distribution of wing kinematics. these joint
% distributions will be used to identify frames where the wing kinematics
% bear investigation
%
% Note: the method of collecting the kinematic data for forming this
% distribution is pretty specific to a given organization of data files,
% but the method should be generalizable. Specifically, here I'm going to
% collect data from pre-optogenetic-stimulus wingbeats across driver lines
% -------------------------------------------------------------------------
% ------------------------------------------------
%% path and data options
rootPath = 'D:\Fly Data\' ;
wingVarX = 'Phi' ;
wingVarY = 'Eta' ;
effector = 'Chrimson' ; % 'Chrimson' or 'GtACR1'
saveFlag = true ;
overWriteFlag = true ;
plotFlag = true ;
excludeFlag = false ; % exclude flagged wingbeats

% determine the binning to be used for each wing kinematic variable--should
% go a bit beyond typical range to include outliers (which should have low
% probability anyway)
phiMax = 220 ; % stroke (degrees)
phiMin = -10 ;
phiEdges = (phiMin - 0.5) : (phiMax + 0.5) ;

thetaMax = 75 ; % deviation (degrees)
thetaMin = -35 ;
thetaEdges = (thetaMin - 0.5) : (thetaMax + 0.5) ;

etaMax = 220 ; % wing pitch (degrees)
etaMin = -10 ;
etaEdges = (etaMin - 0.5) : (etaMax + 0.5) ;

% ------------------------------------------------------------------------
%% get binning edges for current variable pair
switch wingVarX
    case 'Phi'
        Xedges = phiEdges ;
    case 'Theta'
        Xedges = thetaEdges ;
    case 'Eta'
        Xedges = etaEdges ;
    otherwise
        fprintf('Invalid selection for wing var X: %s \n', wingVarX)
end
switch wingVarY
    case 'Phi'
        Yedges = phiEdges ;
    case 'Theta'
        Yedges = thetaEdges ;
    case 'Eta'
        Yedges = etaEdges ;
    otherwise
        fprintf('Invalid selection for wing var Y: %s \n', wingVarY)
end

% bin centers
Xbin_ctrs = (Xedges(1:end-1) + Xedges(2:end))./2 ;
Ybin_ctrs = (Yedges(1:end-1) + Yedges(2:end))./2 ;

% ------------------------------------------------------------------------
%% which data to pull from
% range of LED powers to calculate for
LED_pwr = 1 ; %[0.333, 0.5, 0.667, 0.8333, 1]  ; % 0.8333 ; % 1 ; % 0.667

% list of drivers that we care about and their corresponding MNs
driver_list = {'SS01062', 'ctrl' ; 'MB258C', 'b1' ; 'b2-gal4', 'b2' ; ...
    'SS41039', 'i1' ; 'SS37246', 'i2' ; 'SS01592', 'iii1' ; ...
    'SS48311', 'hg1' ; 'SS37253', 'hg2' ; 'SS47152', 'ps1' ; ...
    'SS41052', 'tp1' ; 'SS47120', 'tp2' ; 'SS51528', 'tpN' ; ...
    'SS31561', 'DLM' ; 'SS41068', 'DVM'} ;

% determine data directory based on effector type
if strcmp(effector, 'Chrimson')
    aggPath = fullfile(rootPath,'VNC MN Chrimson','aggregated') ;
elseif strcmp(effector, 'GtACR1')
    aggPath = fullfile(rootPath,'Opto Silencing','aggregated') ;
end


% -------------------------------------------------------------------------
%% initialize data storage then loop through drivers
counts_mat = zeros(length(Xedges)-1, length(Yedges)-1) ;
data_counter = 0 ;

for i = 1:size(driver_list,1)
    % get driver and MN name
    driver = driver_list{i,1} ;
    mn_name = driver_list{i,2} ;
    
    %------------------------------------------------------------
    % find data structure containing kinematics as a function of
    % phase
    dataPath = fullfile(aggPath, driver, [driver '_LED=' num2str(LED_pwr)...
        '_wingKinPhaseStruct.mat']) ;
    
    % check if data is there
    if ~exist(dataPath,'file')
        fprintf('no data folder for %s > %s \n', effector, driver)
        continue
    end
    
    % ----------------------------------------------------
    % load data
    try
        wingKinPhaseStruct = importdata(dataPath) ;
    catch
        fprintf('error loading data for %s > %s \n', effector, driver)
        continue
    end
    
    % ------------------------------------------
    % read out selected wing kinematic variables
    wb_num = [wingKinPhaseStruct.wb_num] ;
    prestim_idx = (wb_num < 0) ;
    
    % load in X and Y data
    Xdata = [[wingKinPhaseStruct(prestim_idx).([wingVarX 'R_interp'])] ;
        [wingKinPhaseStruct(prestim_idx).([wingVarX 'L_interp'])]] ;
    Ydata = [[wingKinPhaseStruct(prestim_idx).([wingVarY 'R_interp'])] ;
        [wingKinPhaseStruct(prestim_idx).([wingVarY 'L_interp'])]] ;
    
    % unwrap data
    Xdata = Xdata(:) ;
    Ydata = Ydata(:) ;
    
    % -------------------------------------------
    % bin data and add to running sum
    N = histcounts2(Xdata, Ydata, Xedges, Yedges,'normalization','count') ;
    
    counts_mat = counts_mat + N ;
    data_counter = data_counter + 1 ;
    
    clear wingKinPhaseStruct Xdata Ydata
    fprintf('Completed %d/%d drivers \n', i, size(driver_list,1))
end

% -------------------------------------------------------------------------
%% normalize so that each row (x value) sums to one
% also set rows to zero where there are so few counts that normalizing will
% artificially inflate the probability
joint_dist = counts_mat ;
% N_rows = size(joint_dist,1) ; 
% for j = 1:N_rows
%    rowSumCurr = sum(joint_dist(j,:)) ; 
%    if (rowSumCurr ~= 0)
%        joint_dist(j,:) = joint_dist(j,:)./rowSumCurr ; 
%    else
%        joint_dist(j,:) = zeros(size(joint_dist(j,:))) ; 
%    end
% end
row_sum = nansum(joint_dist,2) ; 
row_sum_z = (row_sum - median(row_sum))./ ...
    median(abs(row_sum - median(row_sum))) ; 
rows_to_norm = (row_sum_z > -3) & (row_sum ~= 0) ; 
joint_dist(rows_to_norm,:) = joint_dist(rows_to_norm,:)./ ...
    row_sum(rows_to_norm) ; 
joint_dist(~rows_to_norm,:) = 0 ; 


% ---------------------------------------------------------------
%% plot results?
if plotFlag
    h = figure ;
    imagesc(Xbin_ctrs, Ybin_ctrs, imgaussfilt(joint_dist')) ;
    
    set(gca,'ydir','normal') 
    xlabel(wingVarX)
    ylabel(wingVarY)
    colorbar
end

% ---------------------------------------------------------------
%% save results?
if saveFlag
    saveNameFull = fullfile(aggPath, [wingVarX '_' wingVarY ...
        '_jointDist.mat']) ;
    if ~exist(saveNameFull,'file') || overWriteFlag
        save(saveNameFull, 'joint_dist', 'Xbin_ctrs', 'Ybin_ctrs') ;
    end
end

%{
%load('D:\Fly Data\VNC MN Chrimson\04_10112017\Analysis\Unsorted\Expr_4_mov_003\Expr_4_mov_003_cleaned.mat')
defineConstantsScript
phiL = data.anglesBodyFrame(:,PHIL) ;
thetaL = data.anglesBodyFrame(:,THETAL) ; 

N_frames = length(phiL) ;
frames = 1:N_frames ; 
theta_joint_prob = zeros(N_frames,1) ; 
for ind = 1:N_frames 
	row_idx = (round(phiL(ind))==Xbin_ctrs) ;
    col_idx = (round(thetaL(ind)) == Ybin_ctrs) ; 
    if (sum(row_idx) == 0) || (sum(col_idx) == 0)
        theta_joint_prob(ind) = 0 ; 
    else
        theta_joint_prob(ind) = joint_dist(row_idx, col_idx) ; 
    end
end

% theta_joint_prob_z = (theta_joint_prob - mean(theta_joint_prob))./...
%     std(theta_joint_prob) ;
low_prob_thresh = 0.015 ; 
low_prob_idx = (theta_joint_prob < low_prob_thresh) ; 
figure ; 
subplot(2,1,1)
hold on
plot(frames, thetaL)
plot(frames(low_prob_idx), thetaL(low_prob_idx),'r.')

subplot(2,1,2)
hold on
plot(frames, theta_joint_prob,'r-')
plot([frames(1) frames(end)], low_prob_thresh*[1, 1], 'k--')
%}