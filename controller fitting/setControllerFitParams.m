% -------------------------------------------------------------------------
% function to load default parameters for controller fitting
% -------------------------------------------------------------------------
function fit_params = setControllerFitParams(pulseTiming)
% -----------------------------
%% magnetic pulse duration
% note that this can be modified as an input, in case folks change this
if ~exist('pulseTiming','var') || isempty(pulseTiming)
    pulseTiming = [0.0, 0.007] ; % in seconds
end
% -----------------------------
%% delta t spacing
deltaT_min = 0.002 ; % in seconds
deltaT_max = 0.013 ;
deltaT_N = 100 ;

% for yaw, delay time is much longer
deltaT_minYaw = 0.002 ;
deltaT_NYaw = 200 ;
deltaT_maxYaw = 0.020 ;  % for yaw, delay time is much longer
% ------------------------------------------
%% parameter guesses for fit initialization
% (based on wild-type data)

% ---------
% roll
% ---------
K_i_guessRoll = 0.7 ; % unitless. from Beatus, 2015
K_p_guessRoll = 0.006 ; % seconds. from Beatus, 2015
deltaTGuessRoll = 0.0046 ; % sec. from Beatus, 2015
K_guessRoll   = 0.0 ; % deg. trying to avoid using this
% ---------
% pitch
% ---------
K_i_guessPitch = 0.5 ; % unitless. from Whitehead, 2015
K_p_guessPitch = 0.008 ; % seconds. from Whitehead, 2015
K_guessPitch   = 0.0 ; % deg. trying to avoid using this
deltaTGuessPitch = 0.007 ; % sec. loosely from Whitehead, 2015 ( i think it was overestimated there)

% ---------
% yaw
% ---------
% NB: Leif's results were in different units (torque) so for yaw guess just
% taking roll values for gain and ~2.5 wingbeat deltaT until we have a
% better estimate
K_i_guessYaw = 0.7 ;  % unitless
K_p_guessYaw = 0.006 ; % seconds
deltaTGuessYaw = 0.012 ;
K_guessYaw   = 0.0 ; % deg. trying to avoid using this
% ------------------------------------------
%% time range for controller fit
% read out start time of magnetic pulse and convert to ms
pulseStartMS = 1000*pulseTiming(1) ;

%pulseEndMS = 1000*pulseTiming(2) ;
%pulseDurMS = pulseEndMS - pulseStartMS ;

% use this value to determine time range of fit
% pitchFitRangeMS = [pulseStartMS - 10, pulseStartMS + 35] ;
% rollFitRangeMS = [pulseStartMS - 10, pulseStartMS + 35] ; % +40
% yawFitRangeMS = [pulseStartMS - 10, pulseStartMS + 50]  ;

% NB: for some experiments with longer perturbations, will need wider fit
% windows


if pulseStartMS > 0
    pitchFitRangeMS = [pulseStartMS - 10, 50] ; % +45
    rollFitRangeMS = [pulseStartMS - 10, 50] ; % +45
    yawFitRangeMS = [pulseStartMS - 10, 50]  ;
else
    pitchFitRangeMS = [pulseStartMS - 10, pulseStartMS + 50] ; % +45
    rollFitRangeMS = [pulseStartMS - 10, pulseStartMS + 50] ; % +45
    yawFitRangeMS = [pulseStartMS - 10, pulseStartMS + 55]  ;
end

%{

if pulseDurMS < 10
    pitchFitRangeMS = [pulseStartMS - 10, pulseStartMS + 35] ;
    rollFitRangeMS = [pulseStartMS - 10, pulseStartMS + 35] ; % +40
    yawFitRangeMS = [pulseStartMS - 10, pulseStartMS + 50]  ;
else
    pitchFitRangeMS = [pulseStartMS - 10, pulseStartMS + 50] ; % +45
    rollFitRangeMS = [pulseStartMS - 10, pulseStartMS + 50] ; % +45
    yawFitRangeMS = [pulseStartMS - 10, pulseStartMS + 55]  ;
end

%}


% ----------------------------------------------
%% misc. params for different fit types
% ---------------------
% pitch fit params
% ---------------------
% fit for a constant term ("K") in the controller equation? (should
% probably never do this)
KFlag = false ;
% try to estimate steady-state (i.e. zero-torque) phi front value?
steadyPhiFrontFlag = true ;

% ----------------------------------------------
%% add to fit_params struct
fit_params = struct() ;

% magnetic pulse info
fit_params.pulseTiming = pulseTiming ;
fit_params.pulseDuration = diff(pulseTiming) ;

% time range for deltaT fitting pitch or roll
fit_params.deltaT_min = deltaT_min ;
fit_params.deltaT_max = deltaT_max ;
fit_params.deltaT_N = deltaT_N ;

% time range for deltaT fitting yaw
fit_params.deltaT_minYaw = deltaT_minYaw ;
fit_params.deltaT_maxYaw = deltaT_maxYaw ;
fit_params.deltaT_NYaw = deltaT_NYaw ;

% param guess for roll fits
fit_params.K_i_guessRoll = K_i_guessRoll ;
fit_params.K_p_guessRoll = K_p_guessRoll  ;
fit_params.K_guessRoll = K_guessRoll ;
fit_params.deltaTGuessRoll  = deltaTGuessRoll ;

% param guess for pitch fits
fit_params.K_i_guessPitch = K_i_guessPitch ;
fit_params.K_p_guessPitch = K_p_guessPitch ;
fit_params.K_guessPitch = K_guessPitch ;
fit_params.deltaTGuessPitch = deltaTGuessPitch ;

% param guess for yaw fits
fit_params.K_i_guessYaw = K_i_guessYaw ;
fit_params.K_p_guessYaw = K_p_guessYaw  ;
fit_params.K_guessYaw = K_guessYaw ;
fit_params.deltaTGuessYaw  = deltaTGuessYaw ;

% time range for fits
fit_params.pitchFitRangeMS  = pitchFitRangeMS ;
fit_params.rollFitRangeMS  = rollFitRangeMS ;
fit_params.yawFitRangeMS  = yawFitRangeMS ;

% misc params
fit_params.KFlag = KFlag ;
fit_params.steadyPhiFrontFlag = steadyPhiFrontFlag ;


end
