% -----------------------------------------------------------------
% read out data from params to a wing_kin_params vector
% -----------------------------------------------------------------
function wing_kin_params = read_wing_kin_params(params)
if ~exist('params','var') || isempty(params)
    params = defineQuasiSteadyParams ; 
end
% -----------------------------------------
% fix delta phi front in this case, since we're explicitly varying phi_f
delta_phi_f = 0 ; 

% list of qs_params wing kin fields
field_list = {'omega', 'phi_f', 'phi_b', 'K', 'theta_0', 'theta_m', ...
    'del_theta','psi_0', 'psi_m', 'del_psi', 'C'} ; 

% factors to convert to radians, if necessary
scale_vec = ones(length(field_list),1) ;  
scale_vec([2,3,5,6,8,9]) = (pi/180).*scale_vec([2,3,5,6,8,9]) ;

% ------------------------------------------
% initialize wing kin params vector
wing_kin_params = zeros(length(field_list)+1,1) ; 

% loop through fields and fill vec
for k = 1:length(field_list)
    wing_kin_params(k) = scale_vec(k)*params.(field_list{k}) ; 
end

% -------------------------------
% add delta phi front
wing_kin_params(length(field_list)+1) = delta_phi_f ; 

end