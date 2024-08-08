% -------------------------------------------------------------------------
% draw 3D schematic of an apparatus that could apply perturbative wind
% gusts (for MUIR grant figure). will settle for a tube that looks like 
% it could shoot air (so tapered kind of like a pipette)
%
%{
fig = figure ; 
hold on
ax = gca ;
parent = hgtransform('Parent',ax);

tubeLength = 8 ; 
tubeRadius = 1.0 ; 
taperLength = 6 ; 
center = [0,0,0] ; 
direction = [0,1,0] ; 

grp = draw3DwindTube(ax, parent, tubeLength, tubeRadius, ...
    taperLength, center, direction) ; 

%}
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function grp = draw3DwindTube(ax, parent, tubeLength, tubeRadius, ...
    taperLength, center, direction, tubeColor, objAlpha, resolution)

% -------------------
%% inputs and params
if ~exist('tubeColor', 'var') || isempty(tubeColor)
    tubeColor = 0.6*[1, 1, 1] ;
end
if ~exist('objAlpha', 'var') || isempty(objAlpha)
    objAlpha = 1.0 ;
end
if ~exist('resolution', 'var') || isempty(resolution)
    resolution = 100 ;
end


% -------------------------------------------------------------------------
% tapered portion will start with radius equal to tube and then decrease to
% a smaller radius:
taperRadius = 0.25*tubeRadius ; 

% --------------------
% get axis info
axes(ax) ;
wasOnHold = ishold ;
hold on ;

% ------------------------------------------
%% generate tbe body (just cylinder)
% body is a cylinder along the positive x axis
[X, Y, Z ] = cylinder(tubeRadius, resolution) ;
Z = Z * tubeLength ;
C = zeros([size(X) 3]) ;
C(:,:,1) = tubeColor(1) ;
C(:,:,2) = tubeColor(2) ;
C(:,:,3) = tubeColor(3) ;
hTube = surf(X, Y, Z) ;
set(hTube,'lineStyle','none') ;
set(hTube,'CData',C) ;
set(hTube,'FaceAlpha',objAlpha,'EdgeAlpha',objAlpha)

% bundle tube into a group
grp = hgtransform('Parent',parent) ;
body_grp = hgtransform('Parent',grp);
set(hTube, 'Parent',body_grp) ;

% ---------------------------------------------
%% generate tapered portion at the end of tube
% make tapered bit as cylinder with changing radius
t = linspace(tubeRadius,taperRadius,resolution) ;
[X, Y ,Z] = cylinder(t) ;
Z = Z * taperLength ;
C = zeros([size(X) 3]) ;
C(:,:,1) = tubeColor(1) ;
C(:,:,2) = tubeColor(2) ;
C(:,:,3) = tubeColor(3) ;
hTaper = surf(X, Y, Z) ;
set(hTaper,'lineStyle','none') ;
set(hTaper,'CData',C) ;

% add tapered bit to a group
taper_grp = hgtransform('Parent',grp); 
set(hTaper, 'Parent',taper_grp) ;

% translate tapered section to the end of the tube
T1 = makehgtform('translate', [0,0,tubeLength]) ;
set(taper_grp, 'Matrix', T1) ;

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
material([hTube, hTaper],'shiny')
set([hTube, hTaper], 'ambientStrength',0.9);

% reset figure's hold properties
if (~wasOnHold)
    hold off ;
end

end