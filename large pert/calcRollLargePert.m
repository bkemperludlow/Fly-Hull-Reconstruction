%--------------------------------------------------------------------------
% if body roll is already sampled at appropriate times, just do the binning
% and smoothing part
%--------------------------------------------------------------------------
function [smoothed_rho, sp_rho, rho_t, rho_samp] = calcRollLargePert(data,t)
%% read in relevant data/params
rhoTimes = data.rhoTimes ;
rhoSamp = data.newRhoSamp ;

smoothingParams = setSmoothingParams ;
rollErr = smoothingParams.roll_est_err ;

% first make sure that we're doing the fit over the full t range, otherwise
% spaps will give zeros
if rhoTimes(end) < length(t)
    rhoTimes = [rhoTimes, length(t)] ; 
    rhoSamp = [rhoSamp, rhoSamp(end)] ; 
end

%% bin data
rhoTimes_diff = diff(rhoTimes) ;
chunk_ind = find(rhoTimes_diff > 2) + 1 ;
chunk_ind = [1, chunk_ind, (length(rhoTimes)+1)] ;

rhoTimes_old = rhoTimes ;
rhoSamp_old = rhoSamp ;

rhoTimes = nan(1,(length(chunk_ind)-1)) ;
rho_samp = nan(1,(length(chunk_ind)-1)) ;

for qq = 1:(length(chunk_ind)-1)
    idx1 = chunk_ind(qq) ;
    idx2 = chunk_ind(qq+1) - 1 ;
    rhoTimes(qq) = round(mean(rhoTimes_old(idx1:idx2))) ;
    rho_samp(qq) = nanmean(rhoSamp_old(idx1:idx2)) ;
end

%% get time information
rho_t    = rhoTimes + data.params.firstTrackableFrame - 1 ;
rho_t    = rho_t / data.params.fps  ;

%% fit spline

tol = rollErr ; % 4 ;
[sp_rho, ~] =  spaps(rho_t, rho_samp, tol) ;
rho0 = fnval(sp_rho, t) ;

smoothed_rho = rho0 ;
end