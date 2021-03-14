% ------------------------------------------------------------------------
% function to segment wings from any view -- Sam version (7/10/19)
% ------------------------------------------------------------------------
function [imWing1, imWing2, sameMasksFlag, noWingsFlag, wing1CM, wing2CM,...
    wingExtremaFlag] = segmentWings_Sam(imFly, imBody, extraCleanFlag, ...
    WING_AREA_THRESH, CONN, debugMode)
% ----------------------------------
%% params and inputs
if ~exist('extraCleanFlag', 'var') || isempty(extraCleanFlag)
    extraCleanFlag = false ;
end
if ~exist('WING_AREA_THRESH', 'var') || isempty(WING_AREA_THRESH)
    WING_AREA_THRESH = 20 ;
end
if ~exist('CONN', 'var') || isempty(CONN)
    CONN = 4 ;
end
if ~exist('debugMode', 'var') || isempty(debugMode)
    debugMode = false ;
end
WING_AREA_MIN = 35 ;
thicken_factor = 1 ; 
imMinLvl = 2 ; 
sameMasksFlag = false ;
noWingsFlag = false ;
wingExtremaFlag = false(2) ; 

% ----------------------------------------------------------
%% try to isolate wing only pixels
thickbody = bwmorph(imBody,'thicken',thicken_factor) ;
dif       = imFly - thickbody ;
dif       = bwmorph(dif,'clean') ; % remove isolated pixels
dif       = bwmorph(dif, 'fill') ;
%dif2       = bwmorph(dif, 'spur',Inf) ;
dif2       = bwareaopen(dif, WING_AREA_THRESH, CONN) ; %bwareaopen(dif2, WING_AREA_THRESH, CONN) ;
if extraCleanFlag
    dif2 = bwmorph(dif2, 'open') ;
    dif2 = bwmorph(dif2, 'close') ;
end

% use watershed segmentation to help in cases where wings are difficult to
% segment
D =  -bwdist(~dif2) ; 
mask = imextendedmin(D,imMinLvl);
D2 = imimposemin(D,mask);
L = watershed(D2) ;
L(~dif2) = 0 ;
dif2(L == 0) = 0 ;

%{
% some attempts at new ways to aid segmentation
if (0)
    figure ; imshow(dif2)
    % bridge background
   dif3 = ~bwmorph(~dif2, 'bridge') ;
   figure ; imshowpair(dif2,dif3) ;
   
   % watershed
   D =  -bwdist(~dif2) ; 
   mask = imextendedmin(D,2);
    figure ; imshowpair(dif2,mask,'blend')
    D2 = imimposemin(D,mask);
   L = watershed(D2) ; 
   L(~dif2) = 0 ; 
   rgb = label2rgb(L);
   figure ; imshowpair(rgb, dif2)
end
%}
% find connected componenets
CC  = bwconncomp(dif2,CONN);
Ncc = length(CC.PixelIdxList) ;

% sort connected components by size
svec = zeros(Ncc,1) ; % vector containing the size of each connected components
for j=1:Ncc
    svec(j) = length(CC.PixelIdxList{j}) ;
end

[~, idx] = sort(svec,'descend') ;
numObjects  = length(idx) ;
wingsIdx    = zeros(2,1) ;

% initialize images for each wing
imWing1 = false(size(imFly)) ;
imWing2 = false(size(imFly)) ;
wingAreas = zeros(1,2) ;

% ---------------------------------------------------------------------
%% determine whether we have at least 2 candidate wings, 1, or none
if (numObjects>=2) % if there are at least two objects, call biggest ones wings
    wingsIdx(1:2) = idx(1:2) ;
elseif (numObjects==1)
    wingsIdx = [ 1 ; 1 ] ;
    sameMasksFlag = true ;
else
    noWingsFlag = true ;
end

% get body bounding box
body_stats = regionprops(imBody,'BoundingBox') ;
bbox_body = body_stats(1).BoundingBox ; 

% --------------------------------------------------
%% fill in wing ONE image -- keep only large objects
if (wingsIdx(1)>0)
    % assign largest CC pixels to be wing1 (set them to 'true')
    imWing1(CC.PixelIdxList{wingsIdx(1)}) = true ;
    % get back some pixels we may have lost due to thickening/spur
    thickwing1 = bwmorph(imWing1,'thicken',thicken_factor) ; 
    imWing1 = imWing1 | (thickwing1 & thickbody) | (thickwing1 & dif) ; 
    % get wing pixel area
    wingAreas(1) = sum(imWing1(:)) ;
    % get wing centroid
    [wing1_idx1, wing1_idx2] = ind2sub(size(imWing1), find(imWing1)) ; 
    wing1CM = [mean(wing1_idx2), size(imWing1,1) - mean(wing1_idx1) + 1] ;
    
    % get wing 1 bounding box to see if wings occupy any extremal voxels
    wing1_stats = regionprops(imWing1,'BoundingBox') ;
    bbox_wing1 = wing1_stats(1).BoundingBox ; 
    
    for j = 1:2
        wingExtremaFlag(1,j) = ((bbox_wing1(j) <= bbox_body(j)) | ...
            ((bbox_wing1(j)+bbox_wing1(j+2)) >= ...
            (bbox_body(j)+bbox_body(j+2)))) & ...
            (wingAreas(1) > WING_AREA_MIN) ;
    end
else
    wing1CM = nan(1,2) ; 
    wingExtremaFlag(1,:) = false ; 
end

% --------------------------------------------------
%% fill in wing TWO image -- keep only large objects
if (wingsIdx(2)>0)
    % assign 2nd largest CC pixels to be wing2 (set them to 'true')
    imWing2(CC.PixelIdxList{wingsIdx(2)}) = true ;
    % get back some pixels we may have lost due to thickening/spur
    thickwing2 = bwmorph(imWing2,'thicken',thicken_factor) ; 
    imWing2 = imWing2 | (thickwing2 & thickbody) | (thickwing2 & dif) ; 
    % get wing pixel area
    wingAreas(2) = sum(imWing2(:)) ;
    
    % if we've killed off a wing2 blob
    if wingAreas(2) < 1
        imWing2 = imWing1 ;
        wingAreas(2) = wingAreas(1) ;
        wing2CM = wing1CM ;
        sameMasksFlag = true ;
    else
        [wing2_idx1, wing2_idx2] = ind2sub(size(imWing2), find(imWing2)) ;
        wing2CM = [mean(wing2_idx2), size(imWing2,1) - mean(wing2_idx1) + 1] ;
    end
    
    % get wing 2 bounding box to see if wings occupy any extremal voxels
    wing2_stats = regionprops(imWing2,'BoundingBox') ;
    bbox_wing2 = wing2_stats(1).BoundingBox ; 
    for k = 1:2
        wingExtremaFlag(2,k) = ((bbox_wing2(k) <= bbox_body(k)) | ...
            ((bbox_wing2(k)+bbox_wing2(k+2)) >= ...
            (bbox_body(k)+bbox_body(k+2)))) & ...
            (wingAreas(2) > WING_AREA_MIN) ;
    end
else
    wing2CM = nan(1,2) ; 
    wingExtremaFlag(2,:) = false ; 
end

% ---------------------------------------------
%% check results?
if  (debugMode) % debugMode
    L = labelmatrix(CC) ;
    figure ;
    % connected components from processed fly image
    subplot(1,3,1)
    imshow(label2rgb(L)) ;
    title('Connected components')
    
    % segmented wings
    subplot(1,3,2)
    imshowpair(imWing1, imWing2)
    title('Wing 1 and Wing 2')
    
    % comparison of wings to original image (make sure we didn't remove
    % features
    subplot(1,3,3)
    imshowpair(imFly, imWing1 | imWing2)
    title('Fly and Segmented Wings')
end

end