% ----------------------------------------------------------------------
% calculate lift and drag coefficients
% based on (Dickinson, Lehmann, & Sane, 1999) and (Whitney & Wood, 2010)
% alpha is the wing angle of attack in RADIANS
% -----------------------------------------------------------------------
function [C_L, C_D] = getForceCoeffs(alpha, params)
% get empirically determined coefficient constants
if exist('params','var') && ~isempty(params)
    % read in param vals (if supplied)
    C_L_max = params.C_L_max ;
    C_D_0 = params.C_D_0 ;
    C_D_max = params.C_D_max ;
else
    % otherwise take hard-coded values
    C_L_max = 1.8 ;
    C_D_0 = 0.4 ;
    C_D_max = 3.4 ;
end

% get current lift and drag coefficients for current angle of attack
C_L = C_L_max.*sin(2*alpha) ; % lift coeff
C_D = ((C_D_max + C_D_0)/2)*(1 - cos(2*alpha)) ; % drag coeff

end
