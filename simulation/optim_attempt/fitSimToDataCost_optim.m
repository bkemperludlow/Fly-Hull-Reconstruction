% -------------------------------------------------------------------------
% cost function for fitting real pitch data to simulation
%
% x is a state vector of controller parameters of the form
%   x = [K_i, K_p, deltaT, pulseStrength]
% -------------------------------------------------------------------------
function [cost, g, H] = fitSimToDataCost_optim(x, t_int, t_comp, s_0, ...
    thetaB0, thetaB_data, params, d_filt, paramScale, delayFunType, timeLimSec)
% --------------------------------------
%% rescale parameters in x
x = (paramScale.^(-1)).*x ;

% add pulse strength to params structure
params.pulseStrength = x(4) ;
% --------------------------------------
%% run simulation
% run pitch simulation
sol = simulatePitch_optim(x(1:3), t_int, s_0, thetaB0, params, true, ...
    true, delayFunType, timeLimSec) ;  % last 3 booleans: controlFlag=true, pertFlag=true, plotFlag=false

%fprintf('Evaluated DDE \n')
% ------------------------------------------------------------
%% evaluate solution after t = 0 (if dde23 ran without issue)
dt = 0.125e-3 ; % spacing of evaluation time (shouldn't matter much)
t_max = max(sol.x) ;
% t_eval = 0 : dt : t_max ;

% switch approach depending on solver output
if (t_max < t_int(end))
    % if DE couldn't be solved, set cost to NaN
    cost = nan ;
else
    % otherwise get mean SSE between simulated body pitch and DATA
    t_eval = t_int(1) : dt : t_int(end) ;
    sint = deval(sol, t_eval) ; %evaluate solution over full range
    
    % read out body pitch
    thetaB_sim = sint(5,:)' ; %radians
    
    % filter and detrend solver output (like subtracting pre-pert value)
    thetaB_sim_filt = filtfilt(d_filt, thetaB_sim) ;
    thetaB_sim_detrend = detrend(thetaB_sim_filt) ;
    
    % restrict simulation data to just range of data
    [~, ind1] = min(abs(t_eval - t_comp(1))) ;
    [~, ind2] = min(abs(t_eval - t_comp(end))) ;
    %     comp_ind = (t_eval >= t_comp(1)) & (t_eval <= t_comp(end)) ;
    comp_ind = ind1:ind2 ;
    thetaB_sim_final = thetaB_sim_detrend(comp_ind) ;
    
    % subtract off value at t == 0
    thetaB_sim_final = thetaB_sim_final - thetaB_sim_final(1) ;
    
    % evaluate cost (sum of squared errors)
    cost = sum((thetaB_sim_final - thetaB_data).^2) ;
    
    %     if (0)
    %        figure ;
    %        hold on
    %        plot(t_eval(comp_ind), thetaB_data)
    %        plot(t_eval(comp_ind), thetaB_sim_final)
    %        xlabel('Time (s)')
    %        ylabel('Body Pitch Angle (rad)')
    %        legend({'data','sim'})
    %     end
end

% ----------------------------------------------------------------
%% don't have a good sense for gradient or Hessian, so return 0
if (nargout > 1)
    g = zeros(size(x)) ;
    if (nargout > 2)
        H = zeros(length(x),length(x)) ;
    end
end

end