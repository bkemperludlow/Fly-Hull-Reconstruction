% -----------------------------------------------------------------
% read out data from params to a body_params vector
% -----------------------------------------------------------------
function body_params = read_body_params(params)
if ~exist('params','var') || isempty(params)
    params = defineQuasiSteadyParams ; 
end
% list of qs_params body fields
field_list = {'span', 'r_hinge', 'thorax_width'} ; 

% initialize wing kin params vector
body_params = zeros(length(field_list),1) ; 

% loop through fields and fill vec
for k = 1:length(field_list)
    body_params(k) = params.(field_list{k}) ; 
end

end