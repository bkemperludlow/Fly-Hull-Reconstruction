% -------------------------------------------------------------------------
% function to determine parameter uncertainties from simulation fit to data
% -------------------------------------------------------------------------
function [param_CI, param_CI_struct] = ...
    estimate_sim_param_CI(sim_data_struct)
%% estimating uncertainty by using jacobian and "nlparci.m"
% read out optimized fit coefficients
beta = [sim_data_struct.K_i, sim_data_struct.K_p, sim_data_struct.deltaT, ...
    sim_data_struct.pulseStrength] ;

% need to scale data to units used for optimization
paramScale = sim_data_struct.paramScale ;
beta_scale = paramScale.*beta ;

% read out residual and jacobian
residual = sim_data_struct.residual ;
J = sim_data_struct.jacobian ;

% calculate CI and rescale
param_CI_scale = nlparci(beta_scale, residual, 'jacobian', J) ;
param_CI = repmat((paramScale.^(-1))', 1, 2).*param_CI_scale ;

% --------------------------------------------------
% arrange output in struture format if called for
if (nargout > 1)
    param_name_list = {'K_i', 'K_p', 'deltaT', 'pulseStrength'} ;
    param_CI_struct = struct() ;
    for k = 1:length(param_name_list)
        field_name = param_name_list{k} ;
        param_CI_struct.(field_name) = param_CI(k,:) ;
    end
end

end