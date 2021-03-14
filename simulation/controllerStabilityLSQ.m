% -------------------------------------------------------------------------
% function to optimize controller coefficients for simulation of flight.
% optimization here just means levels out to a steady angle (set point?)
%
% x is a vector of controller coefficients: x = [K_i, K_p, deltaT] 
%
% NB: going to assume there's a long enough "burn in" time for the system
% to stabilize about a set point (initially there are some odd jumps) PRIOR
% TO T=0 (so we'll always evaluate after that
% -------------------------------------------------------------------------
function [cost, g, H] = controllerStabilityLSQ(x, t_int, s_0, ...
    thetaB0, params, controlFlag, pertFlag, paramScale, delayFunType,...
    timeLimSec)
% ----------------------------------------
%% run pitch simulation based on inputs
% sol = simulatePitch(x, t_int, s_0, thetaB0, params, controlFlag, ...
%     pertFlag, false) ; 
x = (paramScale.^(-1)).*x ; 

sol = simulatePitch_optim(x, t_int, s_0, thetaB0,  params, ...
    controlFlag, pertFlag, delayFunType, timeLimSec) ;
% ------------------------------------------------------------
%% evaluate solution after t = 0 (if dde23 ran without issue)
dt = 0.125e-3 ; % spacing of evaluation time (shouldn't matter much)
t_max = max(sol.x) ; 
t_eval = t_int(1) : dt : t_max ;
%t_eval = 0 : dt : t_max ;

% switch approach depending on solver output
if (t_max < t_int(end))
    % if DE couldn't be solved, set cost to NaN
   cost = nan ;  
else
    % otherwise get mean SSE between simulated body pitch and set point
    sint = deval(sol, t_eval) ; %evaluate solution at tspan
    thetaB = sint(5,:) ; %radians
    
    % thetaB_dot = sint(6,:) ; % radians/sec
    cost = (thetaB - thetaB0.*ones(size(thetaB))) ; 
%     cost = (thetaB - mean(thetaB)) ; 
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