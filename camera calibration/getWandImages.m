%--------------------------------------------------------------------------
% function to pull out bw images of the wand spheres (separate images for
% large vs small sphere) in each view. much of this code taken from
% "getWandPoints"
%
% input parameters:
%    *  radius1, radius2 - radii (range?) of the spheres in pixels. if
%       inputs parameters are empty, set to default values of 17 and 45.
%
%    *  xyFilename, xzFilename, yzFilename - input file names.  can be
%       either cin or tif stacks.
%
%    *  outputName - name of the file to save image struct to
%
%  output parameters:
%       wandImageStruct - structure containing black and white images
%       (logical arrays) for each of the spheres in each of the views
%
% Example usage:
%   wandImageStruct = getWandImages('D:\Fly Data\VoxCleanUpTest\73_08062019\calibration\', 1)
%--------------------------------------------------------------------------
function wandImageStruct = getWandImages(calibrationPath, movieNum,...
    outputName, radius1, radius2)
% ---------------------------
%% inputs and params
if ~exist('outputName','var') || (isempty(outputName))
    outputName = ['wandImageStruct_' num2str(movieNum) '.mat'] ;
end
if ~exist('radius1','var') || (isempty(radius1))
    radius1 = 17 ;
end
if ~exist('radius2','var') || (isempty(radius2))
    radius2 = 45 ;
end

% --------------------
%% get cine info
xyFilename = fullfile(calibrationPath, ['xy_' num2str(movieNum) '.cin']) ;
yzFilename = fullfile(calibrationPath, ['yz_' num2str(movieNum) '.cin']) ;
xzFilename = fullfile(calibrationPath, ['xz_' num2str(movieNum) '.cin']) ;

md1 = getCinMetaData(xyFilename) ;
md2 = getCinMetaData(yzFilename) ;
md3 = getCinMetaData(xzFilename) ;

% open cin files
md1.cindata = myOpenCinFile(md1.filename) ;
md2.cindata = myOpenCinFile(md2.filename) ;
md3.cindata = myOpenCinFile(md3.filename) ;

t01 = md1.firstImage ;
t02 = md2.firstImage ;
t03 = md3.firstImage ;

imWidth = md1.width ;
imHeight = md1.height ;

info1 = [] ;
info2 = [] ;
info3 = [] ;

if (md1.firstImage~=md2.firstImage || md1.firstImage~=md3.firstImage)
    error('movies should have the same time axis. check this.') ;
end
Nim = md1.lastImage - md1.firstImage + 1 ;

% ----------------------------------------------
%% initialize data structure
wandImageStruct = struct() ;
%good = false(Nim,1) ; % keeps the indices where we had all circles in all frames

% ----------------------------------------------
%% loop through movie and grab images
parfor i=1:Nim
    disp(i) ;
    % read images
    t1 = t01 + i - 1 ;
    t2 = t02 + i - 1 ;
    t3 = t03 + i - 1 ;
    
    md1 = getCinMetaData(xyFilename) ;
    md2 = getCinMetaData(yzFilename) ;
    md3 = getCinMetaData(xzFilename) ;
    
    % open cin files
    md1.cindata = myOpenCinFile(md1.filename) ;
    md2.cindata = myOpenCinFile(md2.filename) ;
    md3.cindata = myOpenCinFile(md3.filename) ;
    
    %[im1, ~] = ReadCineFileImage(xyFilename, t1, false);
    %[im2, ~] = ReadCineFileImage(yzFilename, t2, false);
    %[im3, ~] = ReadCineFileImage(xzFilename, t3, false);
    
    im1 = myReadCinImage(md1.cindata, t1) ;
    im2 = myReadCinImage(md2.cindata, t2) ;
    im3 = myReadCinImage(md3.cindata, t3) ;
    
    % find circles with given radii range
    [centersxy, radiixy] = imfindcircles(im1,[radius1, radius2],'ObjectPolarity','dark');
    if ~(isequal(size(centersxy),[2,2]))
        continue ;
    end
    
    [centersyz, radiiyz] = imfindcircles(im2,[radius1, radius2],'ObjectPolarity','dark');
    if ~(isequal(size(centersyz),[2,2]))
        continue ;
    end
    
    [centersxz, radiixz] = imfindcircles(im3,[radius1, radius2],'ObjectPolarity','dark');
    if ~(isequal(size(centersxz),[2,2]))
        continue ;
    end
    
    % if we got to this line, it means that there are exactly two circles in each image
    
    % -------------------------------------------------------
    % ensure order of circle array goes small -> large
    if radiixy(1)>radiixy(2)
        centersxy = [centersxy(2,:); centersxy(1,:)];
        radiixy   = [radiixy(2),radiixy(1)]; %#ok<NASGU>
    end
    
    if radiiyz(1)>radiiyz(2)
        centersyz = [centersyz(2,:); centersyz(1,:)];
        radiiyz   = [radiiyz(2),radiiyz(1)]; %#ok<NASGU>
    end
    
    if radiixz(1)>radiixz(2)
        centersxz = [centersxz(2,:); centersxz(1,:)];
        radiixz   = [radiixz(2),radiixz(1)]; %#ok<NASGU>
        
        
    end
    
    % -----------------------------------------------------------
    % get bw image of circle
%     xy_im_small = drawBWcircle(centersxy(1,:), radiixy(1)) ;
%     xy_im_big = drawBWcircle(centersxy(2,:), radiixy(2)) ;
%     yz_im_small = drawBWcircle(centersyz(1,:), radiiyz(1)) ;
%     yz_im_big = drawBWcircle(centersyz(2,:), radiiyz(2)) ;
%     xz_im_small = drawBWcircle(centersxz(1,:), radiixz(1)) ;
%     xz_im_big = drawBWcircle(centersxz(2,:), radiixz(2)) ;
    
    if (0)
        figure(1) ; clf; %#ok<UNRCH>
        subplot(2,3,2) ; imshow(im1) ; hold on ; viscircles(centersxy(1,:),radiixy(1),'edgecolor','g') ;  viscircles(centersxy(2,:),radiixy(2),'edgecolor','r') ;title('xy')
        subplot(2,3,3) ; imshow(im2) ; hold on ; viscircles(centersyz(1,:),radiiyz(1),'edgecolor','g') ;  viscircles(centersyz(2,:),radiiyz(2),'edgecolor','r') ; title('yz');
        subplot(2,3,1) ; imshow(im3) ; hold on ; viscircles(centersxz(1,:),radiixz(1),'edgecolor','g') ;  viscircles(centersxz(2,:),radiixz(2),'edgecolor','r') ; title('xz') ;
        
        subplot(2,3,5) ; imshowpair(xy_im_big, xy_im_small) ;
        subplot(2,3,6) ; imshowpair(yz_im_big, yz_im_small) ;
        subplot(2,3,4) ; imshowpair(xz_im_big, xz_im_small) ;
        pause(0.01) ;
    end
    
    % ------------------------
    % add data to struct
    % ------------------------
    % xy 
    wandImageStruct(i).xy_radius_big = radiixy(2) ;
    wandImageStruct(i).xy_radius_small = radiixy(1) ;
    wandImageStruct(i).xy_center_small = centersxy(1,:) ;
    wandImageStruct(i).xy_center_big = centersxy(2,:) ;
    % yz 
    wandImageStruct(i).yz_radius_big = radiiyz(2) ;
    wandImageStruct(i).yz_radius_small = radiiyz(1) ;
    wandImageStruct(i).yz_center_small = centersyz(1,:) ;
    wandImageStruct(i).yz_center_big = centersyz(2,:) ;
    % xz
    wandImageStruct(i).xz_radius_big = radiixz(2) ;
    wandImageStruct(i).xz_radius_small = radiixz(1) ;
    wandImageStruct(i).xz_center_small = centersxz(1,:) ;
    wandImageStruct(i).xz_center_big = centersxz(2,:) ;
    
end

% ----------------------------------------
%%redefine cine data after parfor loop
md1 = getCinMetaData(xyFilename) ;
md2 = getCinMetaData(yzFilename) ;
md3 = getCinMetaData(xzFilename) ;

md1.cindata = myOpenCinFile(md1.filename) ;
md2.cindata = myOpenCinFile(md2.filename) ;
md3.cindata = myOpenCinFile(md3.filename) ;

% -----------------------------------
%% clear out empty points in struct 
empty_idx = arrayfun(@(x) isempty(x.xy_radius_big), wandImageStruct) ; 
wandImageStruct = wandImageStruct(~empty_idx) ; 

% ------------------------
%% close cine files
myCloseCinFile(md1.cindata) ;
myCloseCinFile(md2.cindata) ;
myCloseCinFile(md3.cindata) ;

% ------------------------------
%% save results
save(fullfile(calibrationPath, outputName),'wandImageStruct')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bw_circle = drawBWcircle(center, radius, imageSizeX, imageSizeY)

if ~exist('imageSizeX','var') || isempty(imageSizeX)
    imageSizeX = 512 ;
end
if ~exist('imageSizeY','var') || isempty(imageSizeY)
    imageSizeY = 512 ;
end

[columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
bw_circle = (rowsInImage - center(2)).^2 ...
    + (columnsInImage - center(1)).^2 <= radius.^2;
end