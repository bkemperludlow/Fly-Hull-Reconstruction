% -------------------------------------------------------------------------
% function to load default parameters for controller fitting
% -------------------------------------------------------------------------
function fit_params = setControllerFitParams(pulseDuration)
% -----------------------------
%% magnetic pulse duration
% note that this can be modified as an input, in case folks change this
if ~exist('pulseDuration','var') || isempty(pulseDuration)
    pulseDuration = 0.007 ; % in seconds
end
% -----------------------------
%% delta t spacing 
deltaT_min = 0.002 ; % in seconds
deltaT_max = 0.013 ; 
deltaT_N = 100 ; 

% ------------------------------------------
%% parameter guesses for fit initialization
% based on wild-type data
deltaTGuessPitch = 0.007 ; % sec. loosely from Whitehead, 2015 ( i think it was overestimated there)
deltaTGuessRoll = 0.0046 ; % sec. from Beatus, 2015

K_i_guessRoll = 0.7 ; % unitless. from Beatus, 2015
K_p_guessRoll = 0.006 ; % seconds. from Beatus, 2015

K_i_guessPitch = 0.5 ; % unitless. from Whitehead, 2015
K_p_guessPitch = 0.008 ; % seconds. from Whitehead, 2015
K_guessPitch   = 0.0 ; % deg. trying to avoid using this

% ------------------------------------------
%% time range for controller fit
pitchFitRangeMS = [-10, 35] ; 
rollFitRangeMS = [-10, 40] ; 

% ----------------------------------------------
%% add to fit_params struct 
fit_params = struct() ; 

fit_params.pulseDuration = pulseDuration ; 

fit_params.deltaT_min = deltaT_min ; 
fit_params.deltaT_max = deltaT_max ; 
fit_params.deltaT_N = deltaT_N ; 

fit_params.deltaTGuessPitch = deltaTGuessPitch ; 
fit_params.deltaTGuessRoll  = deltaTGuessRoll ; 

fit_params.K_i_guessRoll = K_i_guessRoll ; 
fit_params.K_p_guessRoll = K_p_guessRoll  ; 

fit_params.K_i_guessPitch = K_i_guessPitch ; 
fit_params.K_p_guessPitch = K_p_guessPitch ; 
fit_params.K_guessPitch = K_guessPitch ; 

fit_params.pitchFitRangeMS  = pitchFitRangeMS ; 
fit_params.rollFitRangeMS  = rollFitRangeMS ; 

end
