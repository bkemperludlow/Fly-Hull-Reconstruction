%--------------------------------------------------------------------------
% function to segment the wings in the XY (overhead) view of the fly. This
% function is also in 'hullReconstruction_mk*.m' 
%--------------------------------------------------------------------------
function [imWing1, imWing2, sameMasksFlag] = segementWings(imFly, imBody, WING_AREA_THRESH)

if ~exist('WING_AREA_THRESH', 'var')
    WING_AREA_THRESH = 12 ;
end

sameMasksFlag = false ;
oneWingEmptyFlag = false ;

thickbody = bwmorph(imBody,'thicken',1) ;
dif       = imFly - thickbody ;
dif       = bwmorph(dif,'clean') ; % remove isolated pixels
% find connected componenets and sort by size

CC  = bwconncomp(dif);
Ncc = length(CC.PixelIdxList) ;
    
svec = zeros(Ncc,1) ; % vector containing the size of each connected components
for j=1:Ncc
    svec(j) = length(CC.PixelIdxList{j}) ;
end

[svec, idx] = sort(svec,'descend') ;
numObjects  = length(idx) ;
wingsIdx    = zeros(2,1) ;

if (numObjects>=2)
    for j=1:2
        if (svec(j)>=WING_AREA_THRESH)
            wingsIdx(j) = idx(j) ;
        else
            wingsIdx(j) = -1 ;
        end
    end
elseif (numObjects==1)
    wingsIdx = [ 1 ; 1 ] ;
    sameMasksFlag = true ;
else
    disp('No objects... hmm...') ;
    keyboard ;
end

imWing1 = false(size(imFly)) ;
imWing2 = false(size(imFly)) ;

% keep only large objects
if (wingsIdx(1)>0)
    imWing1(CC.PixelIdxList{wingsIdx(1)}) = true ;
    %imWing1(objData(wingsIdx(1)).PixelIdxList) = true ; % old code
end

if (wingsIdx(2)>0)
    imWing2(CC.PixelIdxList{wingsIdx(2)}) = true ;
    %imWing2(objData(wingsIdx(2)).PixelIdxList) = true ; % old code
end

% if at least one wing is empty, use dif instead
% later on the two wings will be resolved using clustering in
% hullAnalysis_mk?
if (wingsIdx(1)==-1 || wingsIdx(2)==-1)
    
    imWing1 = false(size(dif)) ;
    
    [label, numObjects] = bwlabel(dif, 4);
    if (numObjects>1)
        objData = regionprops(label, 'area','PixelIdxList') ;
        [~, idx] = max([objData.Area]) ;
        imWing1(objData(idx).PixelIdxList) = true ;
    else
        imWing1 = dif ;
    end
    imWing2 = imWing1 ;
    sameMasksFlag = true ;
end

end