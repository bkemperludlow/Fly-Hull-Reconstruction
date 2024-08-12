function M = getWandPoints (radius1, radius2, xyFilename, xzFilename, yzFilename, outputCSVfile, saveFlag)
% find the centers of the two spheres that make up the "easy wand"
% calibrator.
%
% input parameters:
%    *  radius1, radius2 - radii (range?) of the spheres in pixels. if
%       inputs parameters are empty, set to default values of 17 and 45.
%
%    *  xyFilename, xzFilename, yzFilename - input file names.  can be
%       either cin or tif stacks.
%
%    *  outputCSVfile - name of the file to output points to as a CSV.
%       file extension is not needed. The has 12 columns
%
%       cam1pt1x cam1pt1y cam2pt1x cam2pt1y cam3pt1x cam3pt1y    cam1pt2x cam1pt2y cam2pt2x cam2pt2y cam3pt2x cam3pt2y
%
%
%  output parameters:
%       M - the coordinates of the points in the same format as the CSV file.
%
%
% if using cin files as inputs, the following two commands need to be
% called before running this function:
%
% 1. set phantom SDK path by executing the script:
% phantomSDK_setPath
% 2. Load the phantom SDK libraries using the function:
% LoadPhantomLibraries()
%
% M = getWandPoints (17, 45, 'xy.cin', 'xz.cin', 'yz.cin', 'wandPoints.csv') ;
%
tic
% spmd
%     LoadPhantomLibraries
% end

if (isempty(radius1))
    radius1 = 17 ;
end

if (isempty(radius2))
    radius2 = 45 ;
end

if ~exist('saveFlag', 'var') || isempty(saveFlag)
    saveFlag = true ;
end

% look at tracked wand points?
debugFlag = false ;

% get file type
[~, ~, xyFileExt] = fileparts(xyFilename) ;
isTif = strcmp(xyFileExt,'.tif') ;
isCin = strcmp(xyFileExt,'.cin') | strcmp(xyFileExt,'.cine');

if (isTif)
    info1=imfinfo(xyFilename);
    info2=imfinfo(yzFilename);
    info3=imfinfo(xzFilename);
    if (length(info1)~=length(info2) || length(info1)~=length(info3))
        error('movies should have the same time axis. check this.') ;
    end
    Nim = length(info1) ;
    t01 = [] ;
    t02 = [] ;
    t03 = [] ;
    md1 = [] ;
    md2 = [] ;
    md3 = [] ;
    imHeights = [info1.Height ; info2.Height ; info3.Height] ; 
    imWidths = [info1.Width ; info2.Width ; info3.Width] ; 
   
elseif (isCin) % read cin metadata
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
    
    info1 = [] ;
    info2 = [] ;
    info3 = [] ;
    
   % if (md1.firstImage~=md2.firstImage || md1.firstImage~=md3.firstImage)
   %     error('movies should have the same time axis. check this.') ;
   % end
   Nim = md1.lastImage - md1.firstImage + 1 ;
   
   imHeights = [md1.height ; md2.height ; md3.height] ; 
   imWidths = [md1.width ; md2.width ; md3.width] ; 
   % Nim = 720;
    %Nim = 1271 ;
else
    error ('Movie file type can be either TIF or CIN. Aborting.')
end

centers_xy = zeros(Nim,4);
centers_yz = zeros(Nim,4);
centers_xz = zeros(Nim,4);

good = false(Nim,1) ; % keeps the indices where we had all circles in all frames

parfor i=1:Nim %parfor i=1:Nim
    disp(i) ;
    % read images
    if (isTif)
        im1=imread(xyFilename,'Index', i, 'Info', info1);
        im2=imread(yzFilename,'Index', i, 'Info', info2);
        im3=imread(xzFilename,'Index', i, 'Info', info3);
    else
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
        
    end
    
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
    
    if debugFlag
        figure(1) ; clf; %#ok<UNRCH>
        subplot(1,3,2) ; imshow(im1) ; hold on ; viscircles(centersxy(1,:),radiixy(1),'edgecolor','g') ;  viscircles(centersxy(2,:),radiixy(2),'edgecolor','r') ;title('xy')
        subplot(1,3,3) ; imshow(im2) ; hold on ; viscircles(centersyz(1,:),radiiyz(1),'edgecolor','g') ;  viscircles(centersyz(2,:),radiiyz(2),'edgecolor','r') ; title('yz');
        subplot(1,3,1) ; imshow(im3) ; hold on ; viscircles(centersxz(1,:),radiixz(1),'edgecolor','g') ;  viscircles(centersxz(2,:),radiixz(2),'edgecolor','r') ; title('xz') ;
        pause(0.01) ;
    end
    good(i) = true ;
    
    centers_xy(i,:) = [centersxy(1,:), centersxy(2,:)] ;
    centers_yz(i,:)=  [centersyz(1,:), centersyz(2,:)] ;
    centers_xz(i,:) = [centersxz(1,:), centersxz(2,:)] ;
    
end

%redefine cine data after parfor loop
md1 = getCinMetaData(xyFilename) ;
md2 = getCinMetaData(yzFilename) ;
md3 = getCinMetaData(xzFilename) ;

md1.cindata = myOpenCinFile(md1.filename) ;
md2.cindata = myOpenCinFile(md2.filename) ;
md3.cindata = myOpenCinFile(md3.filename) ;

centers_xz = centers_xz(good,:) ;
centers_yz = centers_yz(good,:) ;
centers_xy = centers_xy(good,:) ;

count = size(centers_xy,1) ;

% arrrange point coordinates in the format that EasyWand likes
% the origin in every image in on the bottom left corner, so all y
% coordinates are transformed to (512-y+1)

% the format of the "points" file is:
% cam1pt1x cam1pt1y cam2pt1x cam2pt1y cam3pt1x cam3pt1y    cam1pt2x cam1pt2y cam2pt2x cam2pt2y cam3pt2x cam3pt2y

M      = zeros(count,12);

M(:,1)  = centers_xz(1:count,1);
M(:,3)  = centers_yz(1:count,1);
M(:,5)  = centers_xy(1:count,1);

M(:,7)  = centers_xz(1:count,3);
M(:,9)  = centers_yz(1:count,3);
M(:,11) = centers_xy(1:count,3);

M(:,2)  = imHeights(3) - centers_xz(1:count,2) + 1;
M(:,4)  = imHeights(2) - centers_yz(1:count,2) + 1;
M(:,6)  = imHeights(1) - centers_xy(1:count,2) + 1;

M(:,8)  = imHeights(3) - centers_xz(1:count,4) + 1;
M(:,10) = imHeights(2) - centers_yz(1:count,4) + 1;
M(:,12) = imHeights(1) - centers_xy(1:count,4) + 1;

% close cin files if needed
if (~isTif)
    myCloseCinFile(md1.cindata) ;
    myCloseCinFile(md2.cindata) ;
    myCloseCinFile(md3.cindata) ;
end

if saveFlag
    % save data to the same folder as the movies
    if (~exist('outputCSVfile', 'var'))
        outputCSVfile = [] ;
    end
    
    if (~isempty(outputCSVfile)) && (~strcmp(outputCSVfile(end-3:end),'.csv'))
        dlmwrite([outputCSVfile '.csv'],M) ;
    else
        dlmwrite(outputCSVfile,M) ;
    end
end
end