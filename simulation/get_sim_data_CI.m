% -------------------------------------------------------------------------
% quick function to get CI for data fit to simulation (when fitting to
% normalized mean)
% -------------------------------------------------------------------------
function sim_data_struct_all = get_sim_data_CI(sim_data_struct_all, ...
    pertDataStructProc, pert_amp_mean, debugFlag)
% -------------------------
%% inputs and params
if ~exist('pert_amp_mean','var') || isempty(pert_amp_mean)
    % pert amp mean (across genotypes want to multiply average of normalized
    % pitch traces by consistent pert amplitude. this one is determined from
    % 'new kir' data
    pert_amp_mean = 0.2065 ;
end
if ~exist('debugFlag','var') || isempty(debugFlag)
    % make plot to check CI calculation?
    debugFlag = false ;
end

% number of bootstrap samples for CI calculation
ff_params = janelia_ff_params() ;
N_BOOT_SAMPLES = ff_params.N_BOOT_SAMPLES ;

% are we working with normalized data (probably yes, because we only need a
% CI for averaged data, and tough to average without normalizing first)
normFlag = pertDataStructProc.normFlag ; 

% -----------------------------------------------------
%% get CI for each genotype
% total number of genotypes
N_genotypes = length(sim_data_struct_all) ;
% check that this is likely simulation fit to averaged data
if (N_genotypes > 5) || (N_genotypes < 2)
    fprintf(['Warning: this might not be the right sim_data_struct for'...
        'CI calculation \n'])
    keyboard
end

% initialize figure for debugging?
if debugFlag
   h_debug = figure ;  
   ax_array = gobjects(N_genotypes,1) ; 
   plotColors = lines(3) ; % probably overkill but w/e
end

% get effector and driver info from pert data struct
drivers_all_data = [pertDataStructProc.drivers_all] ;
effectors_all_data = [pertDataStructProc.effectors_all] ;

% also get time range for pert data
t_data = pertDataStructProc.t ;

% loop over genotypes
for k = 1:N_genotypes
    % -------------------------------------------------
    %% calc CI for current genotype
    % current genotype of simulated data + time range
    genotype = sim_data_struct_all(k).genotype ;
    t_sim = sim_data_struct_all(k).t ;
    
    % find indices for matching genotype in pert data struct
    effector_idx = cellfun(@(y) strcmpi(y, genotype{1}), ...
        effectors_all_data) ;
    driver_idx = cellfun(@(y) strcmpi(y, genotype{2}), drivers_all_data) ;
    genotype_idx = (effector_idx & driver_idx) ;
    
    % read out angle data for current genotype
    if normFlag
        % if we fit simulation to mean of normalized data, will need to get
        % normalized data to calclulate CI. To start, need to use the same
        % data set used for average in pertDataStruct (in that code, we
        % exclude some vids)
        bad_idx = pertDataStructProc.bad_idx ; 
        small_pert_idx = pertDataStructProc.small_pert_idx ; 
        
        % use these to get all data that went into norm mean calculation
        idx = genotype_idx' & ~(bad_idx | small_pert_idx) ; 
        
        % read out raw data
        data_curr = pertDataStructProc.angle_aligned_mat(idx,:) ;
        
        % normalize raw data
        pert_amps = pertDataStructProc.pert_amps ; 
        data_curr = data_curr./abs(pert_amps(idx)) ; 
        
        % check if normalization led to weird results
        max_norm_angle = nanmax(abs(data_curr),[],2) ;
        bad_norm_idx = (max_norm_angle > 2) ;
        
        data_curr = data_curr(~bad_norm_idx,:) ;
        
        fprintf('N = %d movies for %s > %s \n', size(data_curr,1), ...
            genotype{1}, genotype{2})
        
    else
       % otherwise (if not normalized) just read data for current genotype
        data_curr = pertDataStructProc.angle_aligned_mat(genotype_idx,:) ;
        
    end
    
    % restrict data to just the time range of simulation
    t_match_idx = (t_data >= t_sim(1)) & (t_data <= t_sim(end)) ;
    t_data_match = t_data(t_match_idx) ; 
    data_curr = data_curr(:, t_match_idx) ;
     
    % calculate CI
    data_CI = bootci(N_BOOT_SAMPLES,@nanmean,data_curr) ;

    % ------------------------------------------------------------
    %% deal with propogating normalization scaling
    % since we likely normalized data to take mean, propogate this to CI
    if normFlag
        % get mean data as calculated from pertDataStruct -- to do this need
        % index of genotype for mean data mats
        mean_effector_idx = cellfun(@(y) strcmpi(y, genotype{1}), ...
            pertDataStructProc.genotypes(:,1)) ;
        mean_driver_idx = cellfun(@(y) strcmpi(y, genotype{2}), ...
            pertDataStructProc.genotypes(:,2)) ;
        mean_idx = mean_effector_idx & mean_driver_idx ;
        
        % check that we're getting the right mean data
        if sum(mean_idx) ~= 1
            fprintf('Error: could not find unique index for %s > %s\n',...
                genotype{1}, genotype{2})
            keyboard
        end
        
        % read out data mean and subtract off value at t=0
        data_mean = pertDataStructProc.angle_mean_mat(mean_idx,:) ;
        data_mean_t0 = data_mean(t_data == 0) ;
        data_mean = data_mean - data_mean_t0 ; 
        
        % get amplitude of data_mean -- use this to calculate scale factor
        % for data
        curr_pert_amp = findpeaks(-1.*data_mean,t_data,'SortStr','descend') ;
        curr_pert_amp = curr_pert_amp(1) ;
        scaleFactor = pert_amp_mean/curr_pert_amp ;
        
        % scale and translate CI to match mean data:
%         % subtract off initial value for CI
%         data_CI = data_CI - data_mean_t0 ; 
        
        % multiply CI by scale factor
%         data_CI_diff = data_CI - nanmean(data_curr) ; 
        data_CI_diff = data_CI - nanmean(data_curr) ; 
        data_CI_diff_scaled = scaleFactor.*data_CI_diff ; 
        data_CI = data_CI_diff_scaled + ...
            scaleFactor.*(nanmean(data_curr) - data_mean_t0) ; 
        
    end
    
%    % ----------------------------------------------------
%     %%  calculate CI
%     data_CI = bootci(N_BOOT_SAMPLES,@nanmean,data_curr) ;
    
    % --------------------------------------------------
    %% match up time ranges
    % seems like time from pert data struct and sim data struct have
    % different intervals, so want to make sure that arrays are all the
    % same length
    data_CI_interp = nan(size(data_CI,1), length(t_sim)) ; 
    for m = 1:size(data_CI,1)
       data_CI_interp(m,:) = interp1(t_data_match, data_CI(m,:), ...
           t_sim, 'spline');
    end
    
    % ------------------------------------------------
    %% add CI to sim_data_struct
    sim_data_struct_all(k).thetaB_data_CI = data_CI_interp ; 
    
    % ---------------------------------------
    %% make plot to debug?
    if debugFlag
        % initialize axis for current genotype
        ax_array(k) = subplot(N_genotypes,1,k, 'Parent', h_debug) ; 
        hold(ax_array(k),'on')
        
        % plot data CI 
        h_CI = fill(ax_array(k), [t_sim, fliplr(t_sim)],...
            [data_CI_interp(1,:),fliplr(data_CI_interp(2,:))],...
            plotColors(k,:),'linestyle','none','HandleVisibility','Off');
        h_CI.FaceAlpha = 0.3 ;
        
        % plot data mean
        data_from_struct = sim_data_struct_all(k).thetaB_data ; 
        plot(ax_array(k), t_sim, data_from_struct, '-', ...
            'Color', plotColors(k,:)) ; 
        
        % also plot simulation output
        sim_output = sim_data_struct_all(k).thetaB_sim ; 
        plot(ax_array(k), t_sim, sim_output, 'k:') ; 
        
        % clean up axes
        axis(ax_array(k),'tight')
        xlabel('time (s)')
        ylabel('angle (deg)')
        
    end
end


end