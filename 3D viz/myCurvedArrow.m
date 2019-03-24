function grp = myCurvedArrow(ax, parent, arcRadius, startAngle, endAngle, ...
    tipDiameter, tipHeight, rodDiam, centerPoint, direction, rgb, resolution, lineStyleStr) %#ok<DEFNU>

axes(ax) ;
wasOnHold = ishold ;

hold on ;

grp     = hgtransform('Parent',parent);
coneGrp = hgtransform('Parent',grp);
rodGrp  = hgtransform('Parent',grp);

% correct endAngle to account for the length of the tip head
delta = tipHeight / arcRadius * 180 / pi ;
endAngle = endAngle - delta ;

% plot torus
a  = rodDiam/2 ;
c  = arcRadius ;
t1 = linspace(0*pi/180,(endAngle-startAngle)*pi/180, resolution) ;
t2 = linspace(0,2*pi, resolution) ;
[u,v]=meshgrid(t1,t2);
X = (c+a*cos(v)).*cos(u);
Y = (c+a*cos(v)).*sin(u);
Z = a*sin(v);

C = zeros([size(X) 3]) ;
C(:,:,1) = rgb(1) ;
C(:,:,2) = rgb(2) ;
C(:,:,3) = rgb(3) ;

hr = surf(X, Y, Z) ;
set(hr,'lineStyle',lineStyleStr) ;
set(hr,'CData',C) ;



% plot cone
t = linspace(tipDiameter/2,0,resolution) ;
[X, Y ,Z] = cylinder(t) ;
Z = Z * tipHeight;
C = zeros([size(X) 3]) ;
C(:,:,1) = rgb(1) ;
C(:,:,2) = rgb(2) ;
C(:,:,3) = rgb(3) ;
hc = surf(X, Y, Z) ;
set(hc,'lineStyle',lineStyleStr) ;
set(hc,'CData',C) ;

% disk under the cone
t = linspace(0,2*pi, resolution) ;
X = tipDiameter/2 * cos(t) ;
Y = tipDiameter/2 * sin(t) ;
Z = zeros(size(X))  ;
C = zeros([size(X) 3]) ;
C(:,:,1) = rgb(1) ;
C(:,:,2) = rgb(2) ;
C(:,:,3) = rgb(3) ;
htop = patch(X,Y,Z,rgb) ;

set([hc htop], 'Parent',coneGrp) ;

% rotate and move cone to the tip of the arrow
trans = makehgtform('translate',[arcRadius 0 0]) ;
R1    = makehgtform('xrotate',-pi/2) ;
R3    = makehgtform('zrotate',delta*pi/180) ;
R2    = makehgtform('zrotate',endAngle*pi/180) ;
set(coneGrp, 'Matrix',R2*trans*R3*R1) ;

% disk at the base of the rod
t = linspace(0,2*pi, resolution) ;
Z = rodDiam/2 * cos(t) ;
X = rodDiam/2 * sin(t) + arcRadius;
Y = zeros(size(X))  ;
C = zeros([size(X) 3]) ;
C(:,:,1) = rgb(1) ;
C(:,:,2) = rgb(2) ;
C(:,:,3) = rgb(3) ;
hbot = patch(X,Y,Z,rgb) ;
set(hbot,'lineStyle',lineStyleStr) ;
set([hbot hr],'parent',rodGrp) ;

% translate and rotate bottom disk to its right place
%trans = makehgtform('translate',[arcRadius 0 0]) ;
R1    = makehgtform('zrotate',startAngle*pi/180) ;
set(rodGrp, 'Matrix',R1) ;


trans = makehgtform('translate',centerPoint) ;

% set direction: direction(1) = yaw, direction(2) = pitch, direction(3) =
% roll
% first make the arrow lie along the x axis
yawRad = direction(1)   * pi / 180 ;
pitchRad = direction(2) * pi/180 ;
rollRad  = direction(3) * pi/180 ;

% roll
R1 = makehgtform('xrotate',rollRad) ;

% now pitch up
R2 = makehgtform('yrotate',-pitchRad) ;

% yaw
R3 = makehgtform('zrotate',yawRad) ;

set(grp, 'Matrix',trans*R3*R2*R1) ;

%set(grp, 'Matrix',trans) ;
if (~wasOnHold)
    hold off ;
end
end