% -------------------------------------------------------------------------
% function to draw optogenetic LED + output beam for use in setup
% schematics
%{
fig = figure('Position',[2239, 346, 560, 420])  ;
hold on
ax = gca ; 
parent = hgtransform('Parent',ax);
grp = hgtransform('Parent',parent) ; 
%}
% -------------------------------------------------------------------------
function grp = draw3DoptoLED(ax, parent, ledLength, ledRadius, beamColorMap, ...
    beamLength, center, direction, lineStyleStr, objAlpha,...
    resolution) 
    
% -------------------
%% inputs and params
if ~exist('lineStyleStr', 'var') || isempty(lineStyleStr)
   lineStyleStr = 'none' ;  
end
if ~exist('objAlpha', 'var') || isempty(objAlpha)
   objAlpha = 1.0 ;  
end
if ~exist('resolution', 'var') || isempty(resolution)
   resolution = 100 ;  
end

% ---------------
% other params
%lineWidth = 0.5 * lineScale ; 
ledBodyColor = 0.1*[1, 1, 1] ; 
lensColor = [0.6196, 0.7922, 0.8824] ; % light blue 

beamRadius = 0.75*ledRadius ;
beamAlpha = 0.25 ; 

% ----------------------
% check colormap input
try
    valTest = brewermap([],beamColorMap) ; 
catch
    disp('invalid colormap selection')
    keyboard
end

% --------------------
% get axis info
axes(ax) ;
wasOnHold = ishold ;
hold on ;
% ------------------------------------------
%% generate LED body (just cylinder + circles)
% body is a cylinder that sits on the face of the camera body pointing
% along the positive x axis
[X, Y, Z ] = cylinder(ledRadius, resolution) ;
Z = Z * ledLength ; 
C = zeros([size(X) 3]) ;
C(:,:,1) = ledBodyColor(1) ;
C(:,:,2) = ledBodyColor(2) ;
C(:,:,3) = ledBodyColor(3) ;
hLED = surf(X, Y, Z) ;
set(hLED,'lineStyle','none') ;
set(hLED,'CData',C) ;
set(hLED,'FaceAlpha',objAlpha,'EdgeAlpha',objAlpha)

% lens is a disk that sits at the end of the lens tube
t = linspace(0,2*pi, resolution) ;
X = ledRadius * cos(t) ;
Y = ledRadius * sin(t) ;
Z = ledLength*ones(size(X))  ;
C = zeros([size(X) 3]) ;
C(:,:,1) = lensColor(1) ;
C(:,:,2) = lensColor(2) ;
C(:,:,3) = lensColor(3) ;
hLEDLens = patch(X,Y,Z,C,'lineStyle','none') ;
set(hLEDLens,'FaceAlpha',objAlpha,'EdgeAlpha',objAlpha)

% cap just closes out the bottom of the LED
t = linspace(0,2*pi, resolution) ;
X = ledRadius * cos(t) ;
Y = ledRadius * sin(t) ;
Z = zeros(size(X))  ;
C = zeros([size(X) 3]) ;
C(:,:,1) = ledBodyColor(1) ;
C(:,:,2) = ledBodyColor(2) ;
C(:,:,3) = ledBodyColor(3) ;
hLEDCap = patch(X,Y,Z,C,'lineStyle','none') ;
set(hLEDCap,'FaceAlpha',objAlpha,'EdgeAlpha',objAlpha)

% bundle body of LED into one group
grp = hgtransform('Parent',parent) ;
body_grp = hgtransform('Parent',grp);
set([hLED, hLEDLens, hLEDCap], 'Parent',body_grp) ;

% -----------------------------------------------------------
%% generate beam of light coming from LED
% this will extend upwards from the top of the lens
[X, Y, Z ] = cylinder(beamRadius, resolution) ;
Z = Z * (beamLength - ledLength) + ledLength ; 
%C = zeros([size(X) 3]) ;
%colorMat = flipud(brewermap(size(X,2),beamColorMap)) ; 
%C(:,:,1) = repmat(colorMat(:,1),1,2)' ;
%C(:,:,2) = repmat(colorMat(:,2),1,2)' ;
%C(:,:,3) = repmat(colorMat(:,3),1,2)' ;
hBeam = surf(X, Y, Z) ;
set(hBeam,'lineStyle','none') ;
%set(hBeam,'CData',C) ;
set(hBeam,'FaceColor','interp') ;
set(hBeam,'FaceAlpha',beamAlpha,'EdgeAlpha',beamAlpha)
colormap(flipud(brewermap([],beamColorMap)))

% give parentage to full group
beam_grp = hgtransform('Parent',grp) ;
set(hBeam, 'Parent',beam_grp) ; 

% -----------------------------------------------------------
%% perform group translation/rotation
% NB: cylindrical axis starts as z axis
direction_hat = direction./norm(direction) ; % make sure it's normalized
z_hat = [0, 0, 1] ; 
rot_ax = cross(z_hat, direction_hat) ; 
% for at least one camera, this will be parallel, giving zero cross product
if norm(rot_ax) <= 0 
    rot_ax = [1, 0, 0] ; 
else
    rot_ax = rot_ax./norm(rot_ax) ;
end
rot_ang = acos(dot(z_hat, direction_hat)) ; 

% make hg transforms
R = makehgtform('axisrotate',rot_ax, rot_ang) ; % rotation
T = makehgtform('translate', center) ;  % translation

set(grp, 'Matrix', T*R) ; 

% -----------------------------------------------------------
%% set material properties
material([hLED, hLEDLens, hLEDCap],'shiny')
set([hLED, hLEDLens, hLEDCap], 'ambientStrength',.9);  
set(hBeam, 'ambientStrength',1.0);  
% reset figure's hold properties
if (~wasOnHold)
    hold off ;
end

end