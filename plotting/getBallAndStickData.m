%--------------------------------------------------------------------------
% function to grab data necessary for a ball and stick plot and put that
% data in a combined structure. this way, a structure array can be the
% input into the plot function so we can plot multiple wingstrokes on top
% of each other.
% 
% INPUTS:
%       -data = movie analysis data structure
%       -wing_side = which wing to look at ('right', 'left', or 'average')
%       -t_val = a time value (in seconds) that occurs during the specific
%           wing stroke in question
%       -ind_spacing = how many time points to skip for the wingstroke
%       (i.e. ind_spacing=2 would mean take every other point) 
%--------------------------------------------------------------------------
function ballAndStickStruct = getBallAndStickData(data, wing_side, t_val, ...
    ind_spacing)
% --------------------
%% params and inputs
% --------------------
if ~exist('ind_spacing','var') || isempty(ind_spacing)
   ind_spacing = 1 ; 
end
% ---------------------
% morphology params
span = 2.5 ; %mm
chord = 0.5 ; %mm % 0.7/2 

% -------------------------------
% initialize data struct
ballAndStickStruct = struct() ; 

% ----------------------------------------------------------
%% pull out kinematics based on which side is being plotted
switch wing_side
    case 'right'
        fwdFlipTimes = data.fwdFlipTimesR ;
        backFlipTimes = data.backFlipTimesR ;
        [~, wingAngles, ~, ~, ~ ] =  smoothWingAngles(data, 'R') ;
        phi = wingAngles(1,:) ;
        theta = wingAngles(2,:) ;
        psi = wingAngles(3,:) ;
    case 'left'
        fwdFlipTimes = data.fwdFlipTimesL ;
        backFlipTimes = data.backFlipTimesL ;
        [~, wingAngles, ~, ~, ~ ] =  smoothWingAngles(data, 'L') ;
        phi = wingAngles(1,:) ;
        theta = wingAngles(2,:) ;
        psi = wingAngles(3,:) ;
    case 'average'
        [~, ~,fwdFlipTimes, backFlipTimes] = getPitchControllerData(data) ;
        [~, wingAnglesR, ~, ~, ~ ] =  smoothWingAngles(data, 'R') ;
        [~, wingAnglesL, ~, ~, ~ ] =  smoothWingAngles(data, 'L') ;
        phi = (wingAnglesR(1,:) + wingAnglesL(1,:))./2 ;
        theta = (wingAnglesR(2,:) + wingAnglesL(2,:))./2 ;
        psi = (wingAnglesR(3,:) + wingAnglesL(3,:))./2 ;
        
    otherwise
        keyboard
end

% time in seconds
t = (data.params.startTrackingTime : data.params.endTrackingTime ) / 8000 ;

% ----------------------------------------
%% define ball and stick coordinates
ball_x = span*cos((pi/180)*phi) ;
ball_y = span*sin((pi/180)*theta) ;

stick_x = chord*cos((pi/180)*psi) ;
stick_y = chord*sin((pi/180)*psi) ;

% ---------------------------------------
%% get wing stroke timing info
fwdFlipTimeStartInd = find(fwdFlipTimes < t_val, 1, 'last') ; 
fwdFlipTimeStart = fwdFlipTimes(fwdFlipTimeStartInd) ; 
fwdFlipTimeEndInd = find(fwdFlipTimes > t_val, 1, 'first') ;
fwdFlipTimeEnd = fwdFlipTimes(fwdFlipTimeEndInd) ; 
backFlipTimeMidInd = find(backFlipTimes > fwdFlipTimeStart, 1, 'first') ;
backFlipTimeMid = backFlipTimes(backFlipTimeMidInd) ; 

[~, ind1] = min(abs(t - fwdFlipTimeStart)) ;
[~, ind2] = min(abs(t - backFlipTimeMid)) ;
[~, ind3] = min(abs(t - fwdFlipTimeEnd)) ;

% ---------------------------------------
%% store info in struct
% ------------------
% back stroke info
ballAndStickStruct.ball_x_back = ball_x(ind1 : ind_spacing : ind2) ;
ballAndStickStruct.ball_y_back = ball_y(ind1 : ind_spacing : ind2) ;
ballAndStickStruct.stick_x_back = stick_x(ind1 : ind_spacing : ind2) ;
ballAndStickStruct.stick_y_back = stick_y(ind1 : ind_spacing : ind2) ;

% ------------------
% fwd stroke info
ballAndStickStruct.ball_x_fwd = ball_x(ind2 : ind_spacing : ind3) ;
ballAndStickStruct.ball_y_fwd = ball_y(ind2 : ind_spacing : ind3) ;
ballAndStickStruct.stick_x_fwd = stick_x(ind2 : ind_spacing : ind3) ;
ballAndStickStruct.stick_y_fwd = stick_y(ind2 : ind_spacing : ind3) ;

% ------------------
% timing info
ballAndStickStruct.fwdFlipTimeStart = fwdFlipTimeStart ; 
ballAndStickStruct.backFlipTimeMid = backFlipTimeMid ; 
ballAndStickStruct.fwdFlipTimeEnd = fwdFlipTimeEnd ; 

ballAndStickStruct.t_ind1 = ind1 ; 
ballAndStickStruct.t_ind2 = ind2 ; 
ballAndStickStruct.t_ind3 = ind3 ; 

ballAndStickStruct.t_val = t_val ; 
% ------------------
% other params 
ballAndStickStruct.span = wing_side ; 
ballAndStickStruct.span = span ; 
ballAndStickStruct.chord = chord ; 
ballAndStickStruct.ind_spacing = ind_spacing ; 

end