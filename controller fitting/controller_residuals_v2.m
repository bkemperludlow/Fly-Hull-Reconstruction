function [f, g, H] = controller_residuals_v2(x,phiFront,flipTimes, c_pitch)
%---------------------------------------------------------
%Defines the fit residuals for a PI controller model
%
%Inputs:
%   x = input vector of the things I'm fitting for. Should
%       be of the form [K_i, K_p, deltaT, K] 
%   phiFront = front stroke angle
%   flipTimes = wing forward flip times (when phiFront is evaluated)
%   sp_pitch = spline fit for the body pitch angle
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
K = x(4) ;
%theta_0 = x(5) ;
%deltaT_2 = x(5) ; 

if length(phiFront) ~= length(flipTimes)
    disp('Flip times and phiFront must have matching array lengths')
    return ;
end

bodyPitch = c_pitch(flipTimes - deltaT) ;
%theta_0 = fnval(sp_pitch, flipTimes - deltaT_2) ;
%theta_0 = mean(c_pitch(flipTimes(flipTimes < 0))) ;
theta_0 = c_pitch(0) ;
pitchVelocity = differentiate(c_pitch, flipTimes - deltaT) ;

controllerPrediction = K_i*(bodyPitch - theta_0) + K_p*pitchVelocity + K ;
%controllerPrediction = K_i*bodyPitch  + K_p*pitchVelocity + K ;
diff = phiFront - controllerPrediction' ;

f = sum(diff.^2) ;

%need to define these later ; 
if nargout > 1
    g = zeros(size(x)) ;
    if nargout > 2
        H = zeros(length(x),length(x)) ;
    end
end