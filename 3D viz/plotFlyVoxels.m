%--------------------------------------------------------------------------
% function to generate a plot of the body and wing voxels for a specific
% frame
%--------------------------------------------------------------------------
function [h_vox, bodyCoords, rightWingCoords, leftWingCoords] = ...
    plotFlyVoxels(data, frame_num, h_vox, altChordFlag)
%--------------------------------------------------------------------------
%% params and inputs
if ~exist('h_vox','var') || isempty(h_vox)
    h_vox = figure('PaperPositionMode','auto') ;
else
    set(0,'CurrentFigure',h_vox)
end
if ~exist('altChordFlag','var') || isempty(altChordFlag)
    altChordFlag = true ; 
end
ax = gca ; 

scale = 4 ; % 2; % 4
Lstub = 3.0*scale ;
Lbar  = 10*scale ;

azview = -43 ;
elview =  26 ;

rightWingInd = data.rightWingInd ; 
leftWingInd = data.leftWingInd ; 
bodyInd = data.bodyInd ;
%--------------------------------------------------------------------------
%% get voxel data from struct
df = diff(data.res(:,1)) ;
frameStartInd = [1 ; find(df==1)+1] ;
frameEndInd   = [frameStartInd(2:end)-1 ; size(data.res,1)] ;
clear df ;

row1 = frameStartInd(frame_num) ;
row2 = frameEndInd(frame_num) ;

coords = data.res(row1:row2,2:4) ; % xyz positions of each voxel in frame=1st col value
IDX = data.RESIDX(row1:row2,:) ; % Classifying each voxel as L/R wing or body
bodyRows = (IDX(:,bodyInd)==1) ;
wingRows_R = (IDX(:,rightWingInd)==1) ;
wingRows_L = (IDX(:,leftWingInd)==1) ;

% voxel coordinates for fly
bodyCoords = coords(bodyRows,:) ; 
rightWingCoords = coords(wingRows_R, :) ; 
leftWingCoords = coords(wingRows_L, :) ; 
%--------------------------------------------------------------------------
%% get vector data from struct
cb = data.bodyCM(frame_num,:) ; % body cm
normRolls = data.rollVectors ; 
AHat = data.AHat(frame_num,:) ; 
psiCurr  = atan2(AHat(2), AHat(1)); % body angle with respect to x axis
psiHat  = [-sin(psiCurr)  cos(psiCurr)  0];

%------------------
% roll bar
checksum = sum(normRolls(frame_num,:));
if ((~isfinite(checksum) || norm(normRolls(frame_num,:))==0) && (frame_num>1))
    normRolls(frame_num,:) = normRolls(frame_num-1,:) ;
end

% if still zero, start with zero roll
if (norm(normRolls(frame_num,:))==0)
    normRolls(frame_num,:) = - cross( AHat, [0 0 1]) ;
    %disp('starting from rho=0') ;
end
normRoll = normRolls(frame_num,:) ;
% make sure the roll vector is perpendicular to AHat  and renormalize
normRoll = normRoll - AHat * dot(normRoll, AHat) ;
normRoll = normRoll / norm(normRoll) ;

%-------------------
% other body vectors
stubvec  = Lstub * RotatePoint(normRoll,[0 0 0],AHat,90) ;

%---------------------
% wing vectors
cr = data.rightWingCM(frame_num,:) ;
cl = data.leftWingCM(frame_num,:) ;

rightSpanHat = data.rightSpanHats(frame_num,:) ; 
leftSpanHat = data.leftSpanHats(frame_num,:) ;
rightChordHat = data.rightChordHats(frame_num,:) ; 
leftChordHat = data.leftChordHats(frame_num,:) ;
if altChordFlag
   rightChordAltHat = data.chord1AltHats(frame_num,:) ; 
   leftChordAltHat = data.chord2AltHats(frame_num,:) ; 
end
%--------------------------------------------------------------------------
%% make plot

% voxels
plot3(ax, bodyCoords(:,1), bodyCoords(:,2), bodyCoords(:,3),...
    'g.','MarkerSize',3) ;
hold(ax, 'on')
plot3(ax, rightWingCoords(:,1), rightWingCoords(:,2), ...
    rightWingCoords(:,3), 'r.','MarkerSize',3) ;
plot3(ax, leftWingCoords(:,1), leftWingCoords(:,2),...
    leftWingCoords(:,3), 'b.','MarkerSize',3) ;

% body vectors
line(ax,[cb(1), cb(1)+ stubvec(1)], [cb(2), cb(2)+ stubvec(2)], ...
    [cb(3), cb(3)+ stubvec(3)],'Color','r','LineWidth',8);
line(ax, [cb(1)+stubvec(1)-Lbar*normRoll(1), cb(1)+stubvec(1)+Lbar*normRoll(1)],...
    [cb(2)+stubvec(2)-Lbar*normRoll(2), cb(2)+stubvec(2)+Lbar*normRoll(2)], ...
    [cb(3)+stubvec(3)-Lbar*normRoll(3), cb(3)+stubvec(3)+Lbar*normRoll(3)],...
    'Color','k','LineWidth',8);
line(ax,[cb(1),cb(1)+scale*7*psiHat(1)], [cb(2),cb(2)+scale*7*psiHat(2)],...
    [cb(3),cb(3)+scale*7*psiHat(3)], 'Color','k','LineWidth',8);
line(ax, [cb(1), cb(1)+ scale*15*AHat(1)],...
    [cb(2), cb(2)+ scale*15*AHat(2)], ...
    [cb(3), cb(3)+ scale*15*AHat(3)], 'Color','r','LineWidth',8);

% wing vectors (span)
line(ax,[cr(1), cr(1)+ scale*10*rightSpanHat(1)], ...
    [cr(2), cr(2)+ scale*10*rightSpanHat(2)],...
    [cr(3), cr(3)+ scale*10*rightSpanHat(3)],...
    'Color','k','LineWidth',4);
line(ax, [cl(1), cl(1)+ scale*10*leftSpanHat(1)], ...
    [cl(2), cl(2)+ scale*10*leftSpanHat(2)], ...
    [cl(3), cl(3)+ scale*10*leftSpanHat(3)],...
    'Color','k','LineWidth',4);

% wing vectors (chord)
line(ax, [cr(1), cr(1)+ scale*8*rightChordHat(1)],...
    [cr(2), cr(2)+ scale*8*rightChordHat(2)],...
    [cr(3), cr(3)+ scale*8*rightChordHat(3)],...
    'Color','b','LineWidth',4);
line(ax,[cl(1), cl(1)+ scale*8*leftChordHat(1)],...
    [cl(2), cl(2)+ scale*8*leftChordHat(2)],...
    [cl(3), cl(3)+ scale*8*leftChordHat(3)],...
    'Color','r','LineWidth',4);

% wing vectors (alt chord)
if altChordFlag
    line(ax, [cr(1), cr(1)+ scale*8*rightChordAltHat(1)],...
        [cr(2), cr(2)+ scale*8*rightChordAltHat(2)],...
        [cr(3), cr(3)+ scale*8*rightChordAltHat(3)],...
        'Color','b','LineWidth',4,'LineStyle','--');
    line(ax,[cl(1), cl(1)+ scale*8*leftChordAltHat(1)],...
        [cl(2), cl(2)+ scale*8*leftChordAltHat(2)],...
        [cl(3), cl(3)+ scale*8*leftChordAltHat(3)],...
        'Color','r','LineWidth',4,'LineStyle','--');
end

hold(ax, 'off')
box(ax, 'on')
grid(ax, 'on')
axis(ax, 'equal')
rotate3d on 

%initialize view
view(azview,elview) ;


end