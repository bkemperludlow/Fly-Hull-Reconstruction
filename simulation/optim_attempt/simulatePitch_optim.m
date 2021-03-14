% -------------------------------------------------------------------------
% updated function to run longitudinal flight simulations to test control
% parameters
%
% INPUTS:
%   - controllerTerms: 3-element vector specifying PI controller quantities
%       of the form [K_i, K_p, deltaT]
%   - t_int: interval over which to run de solver, specified as 2-element
%       vector. In seconds
%   - s_0: initial conditions state vector (also 'history' for dde). of the
%       form s_0 = [x, x_dot, z, z_dot, thetaB, thetaB_dot]
%   - thetaB0: controller angle set point
%   - params: quasi steady param structure, as defined by
%       defineQuasiSteadyParams.m
%   - controlFlag: boolean specifying whether or not to implement control
%   - pertFlag: boolean indicating whether or not to apply perturbation.
%       pert parameters are stored in params structure.
%   - plotFlag: boolean to plot results or not
%
% OUTPUT:
%   - sol: solution output structure from DE solver. can use deval.m to get
%       values of solution at arbitrary time points in range
%
% -------------------------------------------------------------------------
function sol = simulatePitch_optim(controllerTerms, t_int, s_0, thetaB0, ...
    params, controlFlag, pertFlag, delayFunType, timeLimSec)
% ----------------------------------------------
%% add controller terms to params struct
if controlFlag
    params.K_i = controllerTerms(1) ;
    params.K_p = controllerTerms(2) ;
    params.deltaT = controllerTerms(3) ;
else
    params.K_i = 0 ;
    params.K_p = 0 ;
    params.deltaT = 0 ;
end

% % ----------------------------------
% % type of delay function to use
% % 'constant' = just deltaT (updates at each time point)
% % 'average'  = average over previous wingbeat, also different at every time point
% % 'function' = use "myControllerDelay.m" to get delay. should give same
% %              value for multiple successive time points
% delayFunType = 'function' ;  % 'constant' | 'average' | 'function'
% 
% % limit time of simulation run?
% timeLimSec = 80 ; % nan ; %80 ;

% ----------------------------------------------
%% run solver
if controlFlag && (params.deltaT > 0)
    % values for timing solutions to dde
    timerVal = tic ;
    
    
    % set output function (timer)
    output_fun = @(t,y,flag) myTimerDDE(t, y, flag, timerVal, timeLimSec) ;
    options = ddeset('OutputFcn',output_fun) ;
    
    % define differential equation
    func = @(t,s,Z) longitudinalFlightODE_optim(t,s,Z,params,thetaB0,...
        pertFlag) ;
    
    % define delay function
    switch delayFunType
        case 'constant'
            delay_fun = params.deltaT ;
        case 'average'
            delay_fun = (linspace(params.deltaT, ...
                (params.deltaT + (2*pi)/params.omega), 30))' ;
        case 'function'
            delay_fun = @(t, s) myControllerDelay(t, s, params.deltaT, ...
                params.omega) ;
    end
    
    % run solver
    sol = ddesd(func, delay_fun, s_0, t_int, options) ; % t_int
    
    %     func = @(t,s,Z) longitudinalFlightODE_optim(t,s,Z,params,thetaB0,...
    %         pertFlag) ;
    %     sol = dde23(func, params.deltaT, s_0, t_int, options) ; % t_int
    
else
    % NB: want to use ode23 here (despite it being inferior to ode45 for
    % our purposes) for consistency with conditions for delayed vs
    % non-delayed cases
    
    func = @(t,s) longitudinalFlightODE_optim(t,s,[],params, thetaB0,...
        pertFlag) ;
    %     sol = ode45(@(t,s) longitudinalFlightODE(t,s,[],params), t_int, s_0) ;
    sol = ode23(func, t_int, s_0) ;
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function delta_phi_f = ...
%     pitch_controller_func(thetaB_lag, thetaB_dot_lag, thetaB0, K_i, K_p)
%
% % change in fwd stroke amplitude is determined by body pitch angle and vel
% delta_phi_f = K_i.*(thetaB_lag - thetaB0) + K_p.*(thetaB_dot_lag) ;
% delta_phi_f = nanmean(delta_phi_f) ;
%
% if abs(delta_phi_f) > 0.8727
%    delta_phi_f = sign(delta_phi_f)*0.8727 ;
% end
%
% end