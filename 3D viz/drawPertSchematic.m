% -------------------------------------------------------------------------
% function to generate a figure showing the three camera views of the fly
% (with stills) as walls of a chamber with 3D-reconstructed flies in the 
% middle. this is an attempt to clean up "perturbationFigure.m" 
%{
rootPath = 'D:\Box Sync Old\Haltere Experiment\69_250720238\' ; 
MovNum = 1 ;
snapshotTimes = [-0.02 , 0.0 , 0.02 , 0.04] ; 
pertType = 'Pitch' ;

fig = drawPertSchematic(rootPath, MovNum, snapshotTimes, pertType) ;
%}
% -------------------------------------------------------------------------
function main_fig = drawPertSchematic(rootPath, MovNum, snapshotTimes, ...
    pertType, pix_pad_array, view_az_el)
% ----------------------------
%% params and inputs
if ~exist('pix_pad_array','var') || isempty(pix_pad_array)
    pix_pad_array = 20*[-1, 1, -1, 1, -1, 1] ; 
end
if ~exist('view_az_el','var') || isempty(view_az_el)
    view_az_el = [160, 20] ;
end
% scale for drawings
imScale = 5 ;
scale = 7*imScale ;
imageResolution = [1024 768] ; 
flyResolution =  100 ;

% view angles
az = view_az_el(1); %61 ; %54
el = view_az_el(2) ; %23 ; % 12

% pixel padding for when we trim images
pix_pad_array = imScale*pix_pad_array ;

% thickness of borders for image walls
wallLineWidth = 1.0 ; 

% lighting stuff
fly_material = 'shiny' ; 
fly_ambient_light = 0.9 ; 
wall_ambient_light = 0.7 ; 
light_pos = [-0.928136 -1.27747 0.711783] ; 

% path struct for directory
pathStruct = generatePathStruct(rootPath) ; 
ExprNum = pathStruct.ExprNum ; 


% --------------------------------------------------------------
%% get images (camera stills) to be used for walls of chamber
% this takes a little while :(
fprintf('Collecting fly images...\n')

rootPathSplit = strsplit(rootPath,'\') ;
cinRoot = strjoin(rootPathSplit(1:end-2),'\') ;
[~, ~, ~, im_xy, im_xz, im_yz] = getFlyVidStills(cinRoot, ExprNum, ...
    MovNum, snapshotTimes, imScale) ; 

% get axis limits for walls
axlim = get_wall_axlim(im_xy, im_xz, im_yz, pix_pad_array) ; 
% -----------------------------------------------------------------
%% initialize surfaces for the walls
[x1,y1] = meshgrid(axlim(1:2), axlim(3:4));
z1 = zeros(size(x1)) + axlim(5);

[y2, z2] = meshgrid(axlim(3:4), axlim(5:6));
x2 = zeros(size(y2)) + axlim(2); 

y2b = y2 ;
z2b = z2 ;
x2b = zeros(size(y2)) + axlim(1); 

[x3, z3] = meshgrid(axlim(1:2), axlim(5:6));
y3 = zeros(size(x3)) + axlim(3); 

x3b = x3 ;
z3b = z3 ;
y3b = zeros(size(x3)) + axlim(4); 

% --------------------------------------------------------------
%% load reconstruction data
% (we'll use if for the 3D kinematics for drawing the fly model)
fprintf('Loading data...\n')
analysisPath = findMovAnalysisPath(pathStruct, MovNum) ; 
data = hierarchicalLoad(analysisPath, ExprNum, MovNum) ;

% --------------------------------------------------------------------
%% smooth kinematic variables and get timing info
[~, smooth_anglesMat_R, ~, ~, ~ ] =  smoothWingAngles(data, 'R') ;
[~, smooth_anglesMat_L, ~, ~, ~ ] =  smoothWingAngles(data, 'L') ;
% fucking transpose
smooth_anglesMat_R = smooth_anglesMat_R' ; 
smooth_anglesMat_L = smooth_anglesMat_L' ; 

[bodyPitch, bodyYaw, bodyRoll] =  smoothBodyAngles(data) ;
bodyYPR = [180+bodyYaw, bodyPitch, bodyRoll] ; 

% [bodyCM_smooth, ~, ~] = smoothBodyCM(data.bodyCM,'lowpass') ;
% bodyCM_smooth = bodyCM_smooth ./ data.params.voxelSize ; 
bodyCM = getPixBodyCM(analysisPath, ExprNum, MovNum, imScale) ;

tms = (1/8)*(data.params.startTrackingTime : data.params.endTrackingTime) ;
tsec = tms./1000 ; 
snapshotIdx = zeros(size(snapshotTimes)) ; 
for i = 1:length(snapshotTimes)
   [~, snapshotIdx(i)] = min(abs(tsec - snapshotTimes(i))) ;  
end
% --------------------------------------------------------------------
%% initialize figure so we can start drawing on it
fprintf('Generating figure...\n')

main_fig = figure('position',[500    100   imageResolution],...
    'units','pixels') ; 
hold on
ax = gca ; % will be deleted and recreated later to make htimer on front
axpos = get(ax,'position') ;

set(gcf,'paperpositionMode','auto') ;
set(gcf,'color',[1 1 1 ]);% [16 16 16 ] / 255) ;
set(gcf,'invertHardcopy','off')
set(gcf,'renderer','opengl')

% set lighting properties
% camlight
% lighting gouraud
%shading interp
%set(hlight,'position',light_pos) ;
% --------------------------------------------------------------------
%% draw flies
fprintf('Drawing 3D flies...\n')
N_flies = length(snapshotTimes) ; 

% determine what kind of pin we want
if contains(pertType, 'pitch', 'IgnoreCase', true)
    pinType = 2 ; % 0==no pin, 1==roll pin, 2==pitch
elseif contains(pertType, 'roll', 'IgnoreCase', true)
    pinType = 1 ; % 0==no pin, 1==roll pin, 2==pitch
else
    pinType = 0 ; 
end
thetab0 = 45 * pi / 180 ;

% loop through and draw flies
flyGrpList = gobjects(N_flies) ; 
for ind = 1:N_flies 
    [flyGrp, ~, rightWingGrp, leftWingGrp, dL, ~, ~] = draw3Dfly(ax, ...
        scale, flyResolution, pinType, thetab0);
    setFlyDOF(flyGrp, rightWingGrp, leftWingGrp, ...
        bodyCM(snapshotIdx(ind),:), ...
        (pi/180)*bodyYPR(snapshotIdx(ind),:), ...
        (pi/180)*smooth_anglesMat_R(snapshotIdx(ind),:) , ...
        (pi/180)*smooth_anglesMat_L(snapshotIdx(ind),:), ...
        thetab0, dL ) ;
    
    material(flyGrp, fly_material)
    flyGrpList(ind) = flyGrp ; 

end

% ---------------------------------------------------------------------
%% set images to walls
fprintf('Setting wall images...\n')

% ---------------------------------------
% initialize surfaces
h1  = surface(x1,y1,z1) ; % z=0 plane
h2  = surface(x2,y2,z2) ; % x=0 plane
h2b = surface(x2b,y2b,z2b) ; % x=0 plane
h3  = surface(x3,y3,z3) ; % y=0 plane % old: surface(x,z+Lpix,y)
h3b = surface(x3b,y3b,z3b) ; 

material([h1 h2 h3 h2b h3b], 'shiny') ;
%shading interp ;
set([h1 h2 h3 h2b h3b], 'ambientStrength',wall_ambient_light);  

for indd = 1:N_flies
    c = findobj(flyGrpList(indd),'Type','surface');
    set(c,'ambientstrength',fly_ambient_light); % note that this line comes after "material shiny"
end
% --------------------------------------------
% crop images 
im1 = flipud(double(im_xy)) ;
im1 = im1(axlim(3):axlim(4), axlim(1):axlim(2)) ;

im2 = flipud(double(im_yz)) ;
im2 = im2(axlim(5):axlim(6), axlim(3):axlim(4)) ;

im3 = flipud(double(im_xz)) ;
im3 = im3(axlim(5):axlim(6), axlim(1):axlim(2)) ;

% --------------------------------------------
% update wall cdata
set(h1,'CData',(double(im1)),'FaceColor','texturemap','linestyle','-',...
    'LineWidth', wallLineWidth)
set(h2,'CData',(double(im2)),'FaceColor','texturemap','linestyle','-',...
    'LineWidth', wallLineWidth)
set(h2b,'CData',(double(im2)),'FaceColor','texturemap','linestyle','-',...
    'LineWidth', wallLineWidth)
set(h3,'CData',(double(im3)),'FaceColor','texturemap','linestyle','-',...
    'LineWidth', wallLineWidth)
set(h3b,'CData',(double(im3)),'FaceColor','texturemap','linestyle','-',...
    'LineWidth', wallLineWidth)

colormap gray ;

if (az<0)
    set(h2,'visible','on') ;
    set(h2b,'visible','off') ;
else
    set(h2,'visible','off') ;
    set(h2b,'visible','on') ;
end

if (az<=-90) || (az>=90)
    set(h3,'visible','on') ;
    set(h3b,'visible','off') ;
else
    set(h3,'visible','off') ;
    set(h3b,'visible','on') ;
end

fprintf('Finishing...\n')

axis equal
axis off
view(az, el)

fprintf('Complete\n')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
%% just move the stupid data filename stuff out of the main code
function data = hierarchicalLoad(analysisPath, ExprNum, MovNum)

dataFilename_corr = fullfile(analysisPath, ['Expr' num2str(ExprNum) 'mov' ...
    num2str(MovNum,'%03d') ...
    '_Data_manually_corrected.mat']) ;
dataFilename_cleaned = fullfile(analysisPath, ['Expr_' num2str(ExprNum) ...
    '_mov_'  num2str(MovNum,'%03d') '_cleaned.mat']) ;
dataFilename_test = fullfile(analysisPath, ['Expr_' num2str(ExprNum) ...
    '_mov_'  num2str(MovNum,'%03d') '_test.mat']) ;

% load data structure for current movie, trying first to get manually
% corrected data if it exists
if exist(dataFilename_corr, 'file')
    data = importdata(dataFilename_corr) ;
elseif exist(dataFilename_cleaned, 'file')
    data = importdata(dataFilename_cleaned) ;
elseif exist(dataFilename_test, 'file')
    data = importdata(dataFilename_test) ;
else
    fprintf('No file for Expr %s Mov %s \n', num2str(ExprNum), ...
        num2str(MovNum,'%03d') )
    data = [] ; 
    return
end

end

% -------------------------------------------------------------------------
%% function to get bounding box containing all objects in binary image
function bbox = totalBBox(BW)

s = regionprops(BW,'BoundingBox'); %get individual boundingboxes
s = cat(1, s.BoundingBox); %concatenate struct into matrix
%s has structure [left bottom width height]
t = [s(:,1) s(:,2) s(:,1)+s(:,3) s(:,2)+s(:,4)]; 
%t has structure [left bottom right top]
leftBottomCorner = min(t(:,1:2));
rightTopCorner = max(t(:,3:4));
%boundingbox, again in [left bottom width height] style
bbox = [leftBottomCorner rightTopCorner-leftBottomCorner];

end

% -------------------------------------------------------------------------
%% function to get bounding box containing all objects in binary image
function axlim = get_wall_axlim(im_xy, im_xz, im_yz, pix_pad_array)

% binarize images
bw_xy = imbinarize(imcomplement(im_xy)) ; 
bw_xz = imbinarize(imcomplement(im_xz)) ; 
bw_yz = imbinarize(imcomplement(im_yz)) ; 
imSize = size(bw_xy,1) ; 

% clean up any noise in images
bw_xy = bwareaopen(bw_xy, 5e3) ; 
bw_xz = bwareaopen(bw_xz, 5e3) ; 
bw_yz = bwareaopen(bw_yz, 5e3) ; 

% get bounding boxes that surround all objects
bbox_xy = totalBBox(bw_xy) ;
bbox_xz = totalBBox(bw_xz) ;
bbox_yz = totalBBox(bw_yz) ;

% get overall bounding coordinates
xmin = round(min([bbox_yz(2), bbox_xy(2)])) + pix_pad_array(1) ; 
xmin = max([xmin, 1]) ; 
xmax = round(max([bbox_yz(2)+bbox_yz(4), bbox_xy(2)+bbox_yz(4)])) + ...
    pix_pad_array(2)  ; 
xmax = min([xmax, imSize]) ;

ymin = round(min([bbox_xz(1), bbox_xy(1)])) + pix_pad_array(3) ; 
ymin = max([ymin, 1]) ; 
ymax = round(max([bbox_xz(1)+bbox_xz(3), bbox_xy(1)+bbox_xy(3)])) + ...
    + pix_pad_array(4); 
ymax = min([ymax, imSize]) ;

zmin = round(min([bbox_xz(2), bbox_yz(2)])) + pix_pad_array(5)   ; 
zmin = max([zmin, 1]) ; 
zmax = round(max([bbox_xz(2) + bbox_xz(4), bbox_yz(2)+bbox_yz(4)])) + ...
    + pix_pad_array(6); 
zmax = min([zmax, imSize]) ;

% combine into single array
%axlim = [xmin, xmax, ymin, ymax, zmin, zmax] ; 
axlim = [ymin, ymax, xmin, xmax, zmin, zmax] ; 

end

% -------------------------------------------------------------------------
%% get the 3D center of mass for the fly based on averages of pixel coords
% because we set up the walls to be perfectly orthogonal, the dlt stuff
% wouldn't look reasonable anyway
function bodyCM = getPixBodyCM(analysisPath, ExprNum, MovNum, imScale)

% load full analysis results
analysisOutput = load(fullfile(analysisPath, ['Expr_' num2str(ExprNum) ...
    '_mov_'  num2str(MovNum,'%03d') '_results.mat'])) ;

% average the CMs estimated from each camera
x_cm_avg = mean([analysisOutput.xcm_xy, analysisOutput.xcm_xz], 2) + 1 ;
y_cm_avg = mean([512-analysisOutput.ycm_xy, analysisOutput.xcm_yz], 2) + 1;
z_cm_avg = mean([512-analysisOutput.ycm_xz, 512-analysisOutput.ycm_yz], 2) + 1;

bodyCM = [x_cm_avg, y_cm_avg, z_cm_avg]*imScale ;


end