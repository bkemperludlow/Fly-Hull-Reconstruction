% -------------------------------------------------------------------------
% function to test the predicted aerodynamic force based on simple body
% motions. the goal is to make sure i'm undertstanding how intuitive flight
% paths in the lab frame translate to body frame variables
% -------------------------------------------------------------------------
function [forcePred, torquePred] = testSyntheticData(fakeDataCaseList)
% ---------------------------
%% make sure input is cell
if ~iscell(fakeDataCaseList)
    fakeDataCaseList = {fakeDataCaseList} ; 
end
% ---------------------------
%% basic params 
params = defineQuasiSteadyParams() ; 
largePertFlag = false ; 
smoothAnglesFlag = false ; 
debugFlag = true ; 

body_length = 0.001/2 ; 
beta_0 = params.beta_0 ;
g = params.g;
body_mass = params.body_mass ;
% --------------------------------
%% fake data basics
% time of full "movie"
tms_start = 0 ; 
tms_end = 60 ; 
dtms = 0.125 ; 
tms = (tms_start:dtms:tms_end)' ; 
t = (1e-3).*tms ; 
N_pts = length(t) ; 

% wingbeat cutoff times
wbFreq = 220 ; % Hz
wbTimes = t(1) : (1/wbFreq) : t(end) ; 

% initialize arrays for body cm position and euler angles
pitchAdd = -pi/6 ; % +pi/8 ; %-pi/6 ; % 0 ; 
bodyCM = zeros(N_pts, 3) ; 
bodyYPR = [zeros(N_pts, 1), (beta_0 + pitchAdd).*ones(N_pts, 1), ...
    zeros(N_pts, 1)] ; 

% params determining basic motion (going to say y = p1*t^2 + p2*t)
CM_p1 = 1.4 ;
CM_p2 = 0.0 ; 

YPR_p1 = 50.0 ; %50 ; %20 ; 
YPR_p2 = 0.0 ; 

% structure translating dimension names to data indices
ind_struct = struct() ; 
ind_struct.x = 1 ; 
ind_struct.y = 2 ; 
ind_struct.z = 3 ; 
ind_struct.yaw = 1 ; 
ind_struct.pitch = 2 ; 
ind_struct.roll = 3 ; 

% also one for direction
sign_struct = struct() ; 
sign_struct.pos = +1 ; 
sign_struct.neg = -1 ; 
sign_struct.right = +1 ; 
sign_struct.left = -1 ; 
sign_struct.up = +1 ; 
sign_struct.down = -1 ; 
% --------------------------------------------------------------------
%% fill in translation/rotation based on option
% NB: the body CM positions are found in LAB FRAME
% first get sign and index of movement 
for k = 1:length(fakeDataCaseList)
    fakeDataCase = fakeDataCaseList{k} ;
    dataCaseSplit = strsplit(fakeDataCase,'_') ;
    var_sign = sign_struct.(dataCaseSplit{2}) ;
    var_ind = ind_struct.(dataCaseSplit{1}) ;
    
    % calculate motion and add to relevant array
    switch fakeDataCase
        case {'z_pos', 'z_neg', 'x_pos', 'x_neg', 'y_pos', 'y_neg'}
            var = CM_p1*t.^2 + CM_p2*t ;
            bodyCM(:, var_ind) = bodyCM(:, var_ind) +  var_sign.*var ;
        case {'pitch_up', 'pitch_down', 'roll_left', 'roll_right', ...
                'yaw_left', 'yaw_right'}
            var = YPR_p1*t.^2 + YPR_p2*t ;
            bodyYPR(:, var_ind) = bodyYPR(:, var_ind) + var_sign.*var ;
        otherwise
            fprintf('Invalid fake data selection: %s \n', fakeDataCase)
            keyboard
    end
end
% --------------------------------------------------------------------
%% plot fake trajectory?
if (0)
    % get body axis
    AHat = [ones(N_pts, 1), zeros(N_pts,2)] ; 
    for m = 1:N_pts
       rotM = eulerRotationMatrix(bodyYPR(m,1), bodyYPR(m,2), bodyYPR(m,3)) ; 
       AHat(m,:) = rotM'*AHat(m,:)' ; 
    end
    
    plot_ind = 1:20:N_pts ; 
    
    % body CM
   figure ; 
   hold on
   plot3(bodyCM(plot_ind,1), bodyCM(plot_ind,2), bodyCM(plot_ind,3),'o-')
   for k = plot_ind
      plot3(bodyCM(k,1) + body_length.*[0, AHat(k,1)],...
          bodyCM(k,2) + body_length.*[0, AHat(k,2)],...
          bodyCM(k,3) + body_length.*[0, AHat(k,3)],'r-')
      plot3(bodyCM(k,1) + body_length*AHat(k,1),...
          bodyCM(k,2) + body_length*AHat(k,2),...
          bodyCM(k,3) + body_length*AHat(k,3),...
          'ro','markerfacecolor','r')
   end
   axis equal
   box on 
   grid on
   xlabel('X')
   ylabel('Y') 
   zlabel('Z')
   
   % angles
   figure ; 
   plot(t, bodyYPR)
   
    
end
% --------------------------------------------------------------------
%% run force/torque prediction
[forcePred, torquePred, wb] = predictForceAndTorque(bodyCM, ...
    bodyYPR, wbTimes, t, [], [], params, largePertFlag, ...
    smoothAnglesFlag, debugFlag) ;

end