% ------------------------------------------------------------------------
% DDE output function to limit total integration time
% ------------------------------------------------------------------------
function stop = myTimerDDE(t, s, state, timerVal, timeLimSec ) %#ok<INUSL>
%MYOUTPUT stop intergration based on some criterion
stop = false;
%time_lim_sec = 80 ; 
switch state
    case []
        % several output values may be calculated in one major time
        % step
        time_elapsed = toc(timerVal) ; 
        
        if time_elapsed > timeLimSec
            % fprintf('Time limit of %d seconds exceeded \n', time_lim_sec)
            stop = true;
        end
    case 'init'
        % Setup for plots or guis
    case 'done'
        % Cleanup of plots, guis, or final plot
end
end