% -------------------------------------------------------------------------
% function to draw 3D version of high speed camera for creating figure
% schecmatics
%
% INPUTS:
%   -ax: handle for axis to be drawn upon
%   -parent: object to which the output hgtransform (grp) belongs
%   -len, wid, dep: x, y, and z dimensions for camera body
%   -center: camera center coordinate, specified as a 1x3 vector
%   -direction: unit vector (1x3) specifying primary axis of camera
%   -lensLength: (optional) lenght of lens tube (scalar)
%   -lensDiam: (optional) diameter of lens tube (scalar)
%   -rgbCam, rgbLensTube, rgbLensCap: (optional) rgb vectors (1x3) 
%    descriping the color of the camera body, lens tube, and lens cap
%   -lineStyleStr: line style for camera body patch objects (tube and lens
%   cap have no lines)
% -------------------------------------------------------------------------
function grp = draw3Dcamera(ax, parent, len, wid, dep, center, direction, ...
    lineScale, lensLength, lensDiam, rgbCam, rgbLensTube, rgbLensCap, ...
    lineStyleStr, objAlpha, resolution)
% ---------------------
%% inputs and params (default preferences)
if ~exist('lensLength', 'var') || isempty(lensLength)
   lensLength = len/3 ;  
end
if ~exist('lensDiam', 'var') || isempty(lensDiam)
   lensDiam = wid/2 ;  
end
if ~exist('rgbCam', 'var') || isempty(rgbCam)
   rgbCam = 0.5*[1 1 1] ;  
end
if ~exist('rgbLensTube', 'var') || isempty(rgbLensTube)
   rgbLensTube = 0.25*[1 1 1] ;  
end
if ~exist('rgbLensCap', 'var') || isempty(rgbLensCap)
   rgbLensCap = [0.6196, 0.7922, 0.8824] ; % light blue 
end
if ~exist('lineStyleStr', 'var') || isempty(lineStyleStr)
   lineStyleStr = '-' ;  
end
if ~exist('objAlpha', 'var') || isempty(objAlpha)
   objAlpha = 1.0 ;  
end
if ~exist('resolution', 'var') || isempty(resolution)
   resolution = 100 ;  
end

% ---------------------
% other params
lightPos = [-0.928136 -1.27747 0.711783] ;
lineWidth = lineScale * 0.5 ; 

% axis controls (make sure we can add without overwriting)
axes(ax) ;
wasOnHold = ishold ;
hold on ;
%{
fig = figure('Position',[2239, 346, 560, 420])  ;
hold on
ax = gca ; 
parent = hgtransform('Parent',ax);
grp = hgtransform('Parent',parent) ; 
%}
% ------------------------------------------
%% generate camera body
% start with object centered (0,0,0) and facing along x axis
% object has dimensions len(gth) x wid(th) x dep(th)
Vertices = [0 0 0 ; len 0 0 ; len wid 0; 0 wid 0; 0 0 dep ; len 0 dep ; ...
    len wid dep ; 0 wid dep] - 0.5*[len wid dep] ; 
Faces = [1 2 6 5 ;2 3 7 6 ; 3 4 8 7 ; 4 1 5 8 ; 1 2 3 4 ; 5 6 7 8] ;
hCamBody =  patch('Vertices',Vertices,'Faces',Faces,...
    'FaceVertexCData',repmat(rgbCam,size(Faces,1), 1),'FaceColor','flat',...
    'lineStyle',lineStyleStr,'LineWidth',lineWidth) ;
set(hCamBody,'FaceAlpha',objAlpha,'EdgeAlpha',objAlpha)

% initialize the group that all camera objects will be part of 
grp = hgtransform('Parent',parent) ;
cam_grp = hgtransform('Parent',grp);
set(hCamBody, 'Parent',cam_grp) ;

% ------------------------------------------
%% generate lens body/cap
% body is a cylinder that sits on the face of the camera body pointing
% along the positive x axis
[X, Y, Z ] = cylinder(lensDiam/2, resolution) ;
Z = Z * lensLength ; 
C = zeros([size(X) 3]) ;
C(:,:,1) = rgbLensTube(1) ;
C(:,:,2) = rgbLensTube(2) ;
C(:,:,3) = rgbLensTube(3) ;
hLensTube = surf(X, Y, Z) ;
set(hLensTube,'lineStyle','none') ;
set(hLensTube,'CData',C) ;
set(hLensTube,'FaceAlpha',objAlpha,'EdgeAlpha',objAlpha)

% cap (glass bit) is a disk that sits at the end of the lens tube
t = linspace(0,2*pi, resolution) ;
X = lensDiam/2 * cos(t) ;
Y = lensDiam/2 * sin(t) ;
Z = lensLength*ones(size(X))  ;
C = zeros([size(X) 3]) ;
C(:,:,1) = rgbLensCap(1) ;
C(:,:,2) = rgbLensCap(2) ;
C(:,:,3) = rgbLensCap(3) ;
hLensCap = patch(X,Y,Z,C,'lineStyle','none') ;
set(hLensCap,'FaceAlpha',objAlpha,'EdgeAlpha',objAlpha)

% ----------------------------------------------------
%% translate and rotate lens parts
R1 = makehgtform('yrotate',pi/2) ;
T1 = makehgtform('translate',[len/2, 0, 0]) ;

% make lens part of the group and perform transformations
lens_grp = hgtransform('Parent',grp);
set([hLensTube, hLensCap], 'Parent',lens_grp) ;
set(lens_grp, 'Matrix',T1*R1) ;

% ----------------------------------------------------
%% translate/rotate full camera
% first get axis/angle for camera rotation so that it faces along
% 'direction'
direction_hat = direction./norm(direction) ; % in case it's not normalized
x_hat = [1, 0, 0] ; % initial camera direction
rot_ax = cross(x_hat, direction_hat) ; 
% for at least one camera, this will be parallel, giving zero cross product
if norm(rot_ax) <= 0 
    rot_ax = [0, 0, 1] ; 
else
    rot_ax = rot_ax./norm(rot_ax) ;
end
rot_ang = acos(dot(x_hat, direction_hat)) ; 

% make hg transforms
R2 = makehgtform('axisrotate',rot_ax, rot_ang) ; % rotation
T2 = makehgtform('translate', center) ;  % translation

set(grp, 'Matrix', T2*R2) ; 
% ----------------------------------------------------
%% set lighting and material properties
% don't have the firmest grasp on this, so worth re-checking
%hlight = light ;
lighting flat
%shading interp
%set(hlight,'position',lightPos) ;

material(hLensCap, 'shiny') ;
material([hCamBody, hLensTube],'metal')
%shading interp ;
set([hCamBody, hLensTube, hLensCap ], 'ambientStrength',.9);  

% ------------------------------------
%% turn off axis hold, if necessary
if (~wasOnHold)
    hold off ;
end
end