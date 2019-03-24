function [f, g, H] = controller_residuals_roll(x,phiAmpDiff,phiAmpTimes, sp_rho)
%---------------------------------------------------------
%Defines the fit residuals for a PI controller model
%
%Inputs:
%   x = input vector of the things I'm fitting for. Should
%       be of the form [K_i, K_p, deltaT] 
%   phiAmpDiff = difference in stroke amplitude
%   flipTimes = when phiAmpDiff is evaluated
%   sp_pitch = spline fit for the body roll angle
%
%Outputs:
%   f = chi squared for the fit parameters given in x
%   g = gradient of function (required by MATLAB's minimizer)
%   H = Hessian of function (required by MATLAB's minimizer)
%
%----------------------------------------------------------

K_i = x(1) ;
K_p = x(2) ;
deltaT = x(3) ;
%K = x(4) ;
%theta_0 = x(5) ;
%deltaT_2 = x(5) ; 

if length(phiAmpDiff) ~= length(phiAmpTimes)
    disp('phiAmpTimes and phiAmpDiff must have matching array lengths')
    return ;
end

bodyRoll = fnval(sp_rho, phiAmpTimes - deltaT) ;

%theta_0 = mean(fnval(sp_pitch, flipTimes(flipTimes < 0))) ;
rollVelocity = fnval( fnder(sp_rho,1), phiAmpTimes - deltaT) ;

controllerPrediction = K_i*bodyRoll + K_p*rollVelocity ;
%controllerPrediction = K_i*bodyPitch  + K_p*pitchVelocity + K ;
diff = phiAmpDiff - controllerPrediction ;

f = sum(diff.^2) ;

%need to define these later ; 
if nargout > 1
    g = zeros(size(x)) ;
    if nargout > 2
        H = zeros(length(x),length(x)) ;
    end
end