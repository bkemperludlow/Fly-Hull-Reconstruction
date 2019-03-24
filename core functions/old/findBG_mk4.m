function [bg, tin, tout, xcm, ycm] = findBG_mk4(cinFilename, thresh)
% find the background image for the given movie. THRESH is a rough binary
% threshold for. defaule value is 100 (works fine).
% bg is the image, 
% tin and tout are estimations for the time the fly goes in and out the
% frame.
% xcm and ycm are a ROUGH estimate for the center-of-mass of the fly

phantomSDK_setPath ;
LoadPhantomLibraries;

warning('off','MATLAB:gui:latexsup:UnableToInterpretTeXString')

if (~exist('thresh','var'))
    thresh = 100 ; %55 %100
end

DEBUG_FLAG = false;

metaData = getCinMetaData(cinFilename) ;
cindata  = myOpenCinFile(cinFilename) ;

%metaData.lastImage = -100  % specific for this movie

N    = double(metaData.lastImage) - double(metaData.firstImage) + 1 ;
%df   = zeros(N,1) ;
xcm  = zeros(N,1) ;
ycm  = zeros(N,1) ;
flyPixels_x = cell(N,1) ; %added by SW 7/4/16
flyPixels_y = cell(N,1) ; 
tvec = (double(metaData.firstImage)) : double(metaData.lastImage) ;
t0_ind = find(tvec == 0) ;

SE = strel('disk',4) ; % se size was 2
W = 60 ; %sets window size for making mask
% read first image
%[im1, ~] = ReadCineFileImage(cinFilename, metaData.firstImage, false);
im0 = myReadCinImage(cindata, 0) ;

bw0 = ( im0 < thresh ) ;
bw0 = imopen(bw0, SE) ;

centermask = false(size(bw0)) ;
x1 = metaData.width/2 - (W+15) ;
x2 = metaData.width/2 + (W+15) ;
y1 = metaData.height/2 - (W+15) ;
y2 = metaData.height/2 + (W+15) ;
centermask(y1:y2, x1:x2) = true ;

%test if center mask captures fly at t=0 (should, but sometimes things are
%misaligned)
bw0_test = bw0 & centermask ;
while sum(bw0_test(:)) < 10
    x1 = x1 - 1 ;
    x2 = x2 + 1 ;
    y1 = y1 - 1 ;
    y2 = y2 + 1  ;
    centermask(y1:y2, x1:x2) = true ;
    bw0_test = bw0 & centermask ;
end
    
bw0 = bw0 & centermask ;
bw0_permanent = bw0 ; 
clear im0 ;


if (~DEBUG_FLAG)
    hbar = waitbar(0,['findBG: Processing images in ' cinFilename(end-9:end)]) ;
end

%forward bit
mask = centermask ;
c = t0_ind ;

for it = 1 : metaData.lastImage
    c=c+1 ;
    %im2 = ReadCineFileImage(cinFilename, it, false);
    im2 = myReadCinImage(cindata, it) ;
    
    bw2 = ( im2 < thresh ) ;
    bw2 = imopen(bw2, SE) ;
    
    % remove objects smaller than something deleteObjSmallerThan
    
    bw2 = bw2 & mask ; % use the previous mask, since the fly did not move much
    
    M = abs( double(bw2) - double(bw0) ) ;
    
    % find the center-of-mass of points in which M is changing
    
    %[yind, xind] = find(bw2==true) ;
    [yind, xind] = find(M>0) ;
    
    %if (numel(yind)<Mthresh)
    %    yind = [] ;
    %    xind = [] ;
    %else
    if ~isempty(yind)
        ycm_init = round(mean(yind)) ; % initial guess
        xcm_init = round(mean(xind)) ;
    end
    % to find CM, mask bw2 in the neighborhood where M is changing
    
    mask = false(size(bw2)) ;
    W = 30 ;
    x1 = max([xcm_init-W, 1]) ;
    x2 = min([xcm_init+W, metaData.width]) ;
    y1 = max([ycm_init-W, 1]) ;
    y2 = min([ycm_init+W, metaData.height]) ;
    
    mask(y1:y2, x1:x2) = true ;
    %end
    bw3 = bw2 & mask ;
    
    [yind, xind] = find(bw3==true) ;
    ycm(c) = mean(yind) ; % initial guess
    xcm(c) = mean(xind) ;
    flyPixels_x{c} = xind ; 
    flyPixels_y{c} = yind ; 
    
    if (DEBUG_FLAG)
        figure(1) ;
        imshow(bw3) ; title(it) ;  hold on ;
        %imagesc(M) ; axis image ; caxis([-1 1]) ;
        hold on ;
        plot(xcm(c),ycm(c),'r+') ;
        hold off ;
        pause(0.01) ;
    end
    %df(c) = sum(M(:)) ;
    bw0 = bw3 ;
    if (~DEBUG_FLAG)
        waitbar(c/N, hbar) ;
    end
end

%backward bit
bw0 = bw0_permanent ; 
mask = centermask ;
c = t0_ind ;
backwardsLoopInd = metaData.firstImage : -1 ;
backwardsLoopInd = fliplr(backwardsLoopInd) ;

for it = backwardsLoopInd
    c=c-1 ;
    %im2 = ReadCineFileImage(cinFilename, it, false);
    im2 = myReadCinImage(cindata, it) ;
    
    bw2 = ( im2 < thresh ) ;
    bw2 = imopen(bw2, SE) ;
    
    % remove objects smaller than something deleteObjSmallerThan
    
    bw2 = bw2 & mask ; % use the previous mask, since the fly did not move much
    
    M = abs( double(bw2) - double(bw0) ) ;
    
    % find the center-of-mass of points in which M is changing
    
    %[yind, xind] = find(bw2==true) ;
    [yind, xind] = find(M>0) ;
    
    %if (numel(yind)<Mthresh)
    %    yind = [] ;
    %    xind = [] ;
    %else
    if ~isempty(yind)
        ycm_init = round(mean(yind)) ; % initial guess
        xcm_init = round(mean(xind)) ;
    end
    % to find CM, mask bw2 in the neighborhood where M is changing
    
    mask = false(size(bw2)) ;
    W = 30 ;
    x1 = max([xcm_init-W, 1]) ;
    x2 = min([xcm_init+W, metaData.width]) ;
    y1 = max([ycm_init-W, 1]) ;
    y2 = min([ycm_init+W, metaData.height]) ;
    
    mask(y1:y2, x1:x2) = true ;
    %end
    bw3 = bw2 & mask ;
    
    [yind, xind] = find(bw3==true) ;
    ycm(c) = mean(yind) ; % initial guess
    xcm(c) = mean(xind) ;
    flyPixels_x{c} = xind ; 
    flyPixels_y{c} = yind ; 
    
    if (DEBUG_FLAG)
        figure(1) ;
        imshow(bw3) ; title(it) ;  hold on ;
        %imagesc(M) ; axis image ; caxis([-1 1]) ;
        hold on ;
        plot(xcm(c),ycm(c),'r+') ;
        hold off ;
        pause(0.01) ;
    end
    %df(c) = sum(M(:)) ;
    bw0 = bw3 ;
    if (~DEBUG_FLAG)
        waitbar(c/N, hbar) ;
    end
end

if (~DEBUG_FLAG)
    close(hbar)
end

clear im0 im2 M
x_diff = abs(xcm(t0_ind-1) - xcm(t0_ind+1)) ;
y_diff = abs(ycm(t0_ind-1) - ycm(t0_ind+1)) ;
if (x_diff < 5) && (y_diff < 5)
    xcm(t0_ind) = mean([xcm(t0_ind-1) xcm(t0_ind+1)]) ;
    ycm(t0_ind) = mean([ycm(t0_ind-1) ycm(t0_ind+1)]) ;
else
    %keyboard ;
    
    splineWindow = 60 ;
    estErr = 2 ;
    t_spline = [(t0_ind-splineWindow):(t0_ind-1),  (t0_ind+1):(t0_ind+splineWindow)] ;
    
    c_x = mySplineSmooth(t_spline, xcm(t_spline),estErr) ;
    c_y = mySplineSmooth(t_spline, ycm(t_spline),estErr) ;
    
    xcm(t0_ind) = fnval(c_x, t0_ind) ; 
    ycm(t0_ind) = fnval(c_y, t0_ind) ;
    
    %this is just temporary. if this point is reached, it means that the
    %CoM estimations from the forward and backward loops are not
    %consistent. To automatically fix this, i could check which one is
    %closest to the center near t = 0, then interpolate between the "good
    %half" and the next point on the bad half that seems consistent
    %(assuming the "bad half" eventually gets to the write trajectory)
end

xcm(isnan(xcm)) = -1 ; 
ycm(isnan(ycm)) = -1 ; 

% the fly is outside the FOV if the CM has an abrupt change (center out).
% Want to use change point analysis here, but for now i'll just do a
% velocity thresh
xcm = hampel(xcm,4) ;
ycm = hampel(ycm,4) ; 
df = sqrt(diff(xcm).^2 + diff(ycm).^2) ;

%---------------------------------
%addedby SW 7/15/15

%out_candidates = abs(df)<= 4 ;
out_candidates = abs(df)<= .01 ; %0.3 ;
%{
%in_candidates = ~out_candidates ;

out_boxcar = zeros(length(out_candidates),1) ; 
winSize = 10 ;
for q = (1+winSize):(length(out_candidates)-winSize)  
    out_boxcar(q) = mean(out_candidates(q-winSize:q+winSize)) ;
end
out_boxcar(1:winSize) = out_boxcar(winSize+1)*ones(winSize,1) ; %not the most intelligent way to do this 
out_boxcar(end-winSize+1:end) = out_boxcar(end-winSize)*ones(winSize,1) ;

out = (out_boxcar > 0.05) ; %can set a real threshold (instead of zero)
%}

out = false(N-1,1) ;
%in  = false(N-1,1) ;
Wout   = 8 ; %4, 6
%Win    = 8 ;
for k=1:(N-1-Wout)
    if ( out_candidates(k:k+Wout-1) )
        out(k:k+Wout-1) = true ;
    end
end

% decide when the fly was in FOV
%tin = tvec(find(~out,1,'first')) ; %this could be changed so that it goes from zero outwards
%tout = tvec(find(~out,1,'last')) ;
%-------------------------------------

% if there is a time in which the fly is outside the FOV, take that as
% background
%out = medfilt1(double(out),10) ;
bgIndex = find(out,1,'first') ; % find the first "out" image and make it the BG image

if (~isempty(bgIndex))
    %bg = ReadCineFileImage(cinFilename, tvec(bgIndex), false);
    bg  = myReadCinImage(cindata, tvec(bgIndex)) ;
else
    % if the fly is always in FOV, fabricate a BG image (automatically, yay)
   
    % find the two most distant points along the trajectory
    D = squareform(pdist([xcm ycm])) ;
    [m1, ~]   = max(D) ;
    [~, ind2] = max(m1) ;
    [~, ind3] = max(D(:, ind2)) ;
    % max is in D(ind3, ind2) ;
    clear D
    
    t1 = ind3;
    t2 = ind2;
    
    if t1 > t2
        temp = t1 ; 
        t1 = t2 ; 
        t2 = temp ; 
    end
    
    xmid = round ( ( xcm(t1) + xcm(t2) ) /2 );
    ymid = round ( ( ycm(t1) + ycm(t2) ) /2 );
    
    %test added by SW 7/4/16
    flyPixels_x_t1 = [] ; 
    flyPixels_x_t2 = [] ; 
    flyPixels_y_t1 = [] ; 
    flyPixels_y_t2 = [] ; 
    for q = -5:5 
        i1 = max([1, t1 + q]) ;
        i2 = min([N, t2 + q]) ;
        flyPixels_x_t1 = [flyPixels_x_t1 ; flyPixels_x{i1} ] ; 
        flyPixels_x_t2 = [flyPixels_x_t2 ; flyPixels_x{i2} ] ; 
        flyPixels_y_t1 = [flyPixels_y_t1 ; flyPixels_y{i1} ] ; 
        flyPixels_y_t2 = [flyPixels_y_t2 ; flyPixels_y{i2} ] ; 
    end
    
    xdist = pdist2(flyPixels_x_t1, flyPixels_x_t2) ;
    ydist = pdist2(flyPixels_y_t1, flyPixels_y_t2) ;
    
    dx = min(xdist(:)) ; 
    dy = min(ydist(:)) ; 
    % find if the larger difference is along x or y
    %dx = abs(xcm(t1) - xcm(t2)) ;
    %dy = abs(ycm(t1) - ycm(t2)) ;
    
    % read the two images
    %imt1 = ReadCineFileImage(cinFilename, tvec(t1), false);
    %imt2 = ReadCineFileImage(cinFilename, tvec(t2), false);
    
    imt1 = myReadCinImage(cindata, tvec(t1)) ;
    imt2 = myReadCinImage(cindata, tvec(t2)) ;
    
    bg = imt1 ;
    
    if (dx>=dy)
        if (xcm(t1)<xmid) % if fly is left to the mid point at t1, take the left part from t2
            bg(:, 1:xmid) = imt2(:,1:xmid) ;
        else % if fly is right to the midpoint at t1, take the right part from t2
            bg(:, xmid:end) = imt2(:, xmid:end) ;
        end
    else
        if (ycm(t1)<ymid) % if fly is above to the mid point at t1, take the top part from t2
            bg(1:ymid,:) = imt2(1:ymid,:) ;
        else % if the fly is below the midpoint at t1, take the bottom part from t2
            bg(ymid:end,:) = imt2(ymid:end,:) ;
        end
    end
end
%------------------------------------------------------------------------
%added by SW 7/16/15
%{
ind = 0 ;
CM_temp = nan(N-1,2) ;
numPix = zeros(N-1,1) ;
buffer = 24 ; %amount inside image the CM has to be to say the fly is fully inside
for k = (metaData.firstImage+1) : metaData.lastImage
    ind = ind + 1 ;
    imtemp = myReadCinImage(cindata, k) ;
    imfly = bg - imtemp;
    bwtemp = ( imfly > thresh ) ;
    
    [yind_temp, xind_temp] = find(bwtemp > 0) ;
    if ~isempty(yind_temp)
        CM_temp(ind,1) = mean(yind_temp) ;
        CM_temp(ind,2) = mean(xind_temp) ;
    end
    numPix(ind) = sum(bwtemp(:)) ;
    if (0)
        figure(1) ;
        imshow(bwtemp) ; title(num2str(k)) ;
        pause(0.01) ;
    end
end

good_ind = find((1+buffer < CM_temp(:,1)) & (512-buffer > CM_temp(:,1)) & ...
    (1+buffer < CM_temp(:,2)) & (512-buffer > CM_temp(:,2)));
%}
tin_ind = find(out(1:t0_ind), 1, 'last') + 80;  
tout_ind = find(out(t0_ind:end), 1, 'first') + t0_ind - 1 - 80 ;

if isempty(tin_ind)
    tin_ind = 1 + 80 ;
end
if isempty(tout_ind)
    tout_ind = length(tvec) - 80 ;
end

tin = tvec(tin_ind) ;
tout = tvec(tout_ind) ;
%------------------------------------------------------------------------
if (DEBUG_FLAG)
    figure ;
    plot(tvec, df,'k.-') ;
    hold on ;
    plot(tvec, out*50,'r-','linewidth',2) ;
    %plot(tvec, in*50,'g-','linewidth',2) ;
    hold off ;
    
    figure ;
    plot(xcm, ycm,'bo-') ;
    axis equal ;
    
    %keyboard ;
end

myCloseCinFile(cindata) ;

warning('on','MATLAB:gui:latexsup:UnableToInterpretTeXString')

return
