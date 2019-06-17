function grp = myArrow(ax, parent, len, tipDiameter, tipHeight, rodDiam, ...
    basePoint, direction, rgb, resolution, lineStyleStr, arrowAlpha) %#ok<DEFNU>

if ~exist('arrowAlpha','var')
    arrowAlpha = 1.0 ; 
end

axes(ax) ;

wasOnHold = ishold ;

hold on ;

% plot cylinder
[X, Y, Z ] = cylinder(rodDiam/2, resolution) ;
Z = Z * len ;
C = zeros([size(X) 3]) ;
C(:,:,1) = rgb(1) ;
C(:,:,2) = rgb(2) ;
C(:,:,3) = rgb(3) ;
hr = surf(X, Y, Z) ;
set(hr,'lineStyle',lineStyleStr) ;
set(hr,'CData',C) ;
set(hr,'FaceAlpha',arrowAlpha,'EdgeAlpha',arrowAlpha)

% plot the cone
t = linspace(tipDiameter/2,0,resolution) ;
[X, Y ,Z] = cylinder(t) ;
Z = Z * tipHeight + len ;
C = zeros([size(X) 3]) ;
C(:,:,1) = rgb(1) ;
C(:,:,2) = rgb(2) ;
C(:,:,3) = rgb(3) ;
hc = surf(X, Y, Z) ;
set(hc,'lineStyle',lineStyleStr) ;
set(hc,'CData',C) ;
set(hc,'FaceAlpha',arrowAlpha,'EdgeAlpha',arrowAlpha)

% plot two rounded disks
% top disk
t = linspace(0,2*pi, resolution) ;
X = tipDiameter/2 * cos(t) ;
Y = tipDiameter/2 * sin(t) ;
Z = zeros(size(X)) + len ;
C = zeros([size(X) 3]) ;
C(:,:,1) = rgb(1) ;
C(:,:,2) = rgb(2) ;
C(:,:,3) = rgb(3) ;
htop = patch(X,Y,Z,C,'lineStyle',lineStyleStr) ;
set(htop,'FaceAlpha',arrowAlpha,'EdgeAlpha',arrowAlpha)

% bottom disk
t = linspace(0,2*pi, resolution) ;
X = rodDiam/2 * cos(t) ;
Y = rodDiam/2 * sin(t) ;
Z = zeros(size(X))  ;
C = zeros([size(X) 3]) ;
C(:,:,1) = rgb(1) ;
C(:,:,2) = rgb(2) ;
C(:,:,3) = rgb(3) ;
hbot = patch(X,Y,Z,C,'lineStyle',lineStyleStr) ;
set(hbot,'FaceAlpha',arrowAlpha,'EdgeAlpha',arrowAlpha)

grp = hgtransform('Parent',parent);
set([hc hr htop hbot], 'Parent',grp) ;

trans = makehgtform('translate',basePoint) ;
set(grp, 'Matrix',trans) ;

% set direction: direction(1) = yaw, direction(2) = pitch
% first make the arrow lie along the x axis
yawRad   = direction(1) * pi/180 ;
pitchRad = direction(2) * pi/180 ;

R1 = makehgtform('yrotate',pi/2) ;

% now pitch up
R2 = makehgtform('yrotate',-pitchRad) ;

% yaw
R3 = makehgtform('zrotate',yawRad) ;

set(grp, 'Matrix',trans*R3*R2*R1) ;

if (~wasOnHold)
    hold off ;
end

end
