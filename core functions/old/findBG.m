function [bg, tin, tout, xcm, ycm] = findBG(cinFilename, thresh)
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

DEBUG_FLAG = true;

metaData = getCinMetaData(cinFilename) ;
cindata  = myOpenCinFile(cinFilename) ;

%metaData.lastImage = -100  % specific for this movie

N    = double(metaData.lastImage) - double(metaData.firstImage) + 1 ;
df   = zeros(N-1,1) ;
xcm  = zeros(N-1,1) ;
ycm  = zeros(N-1,1) ;
tvec = (double(metaData.firstImage)+1) : double(metaData.lastImage) ;

SE = strel('disk',4) ; % se size was 2
% read first image
%[im1, ~] = ReadCineFileImage(cinFilename, metaData.firstImage, false);
im1 = myReadCinImage(cindata, metaData.firstImage) ;

bw1 = ( im1 < thresh ) ;
bw1 = imopen(bw1, SE) ;

clear im1 ;

mask = true(size(bw1)) ;
c=0;
if (~DEBUG_FLAG)
    hbar = waitbar(0,['findBG: Processing images in ' cinFilename(end-9:end)]) ;
end

for it = (metaData.firstImage+1) : metaData.lastImage
    c=c+1 ;
    %im2 = ReadCineFileImage(cinFilename, it, false);
    im2 = myReadCinImage(cindata, it) ;
    
    bw2 = ( im2 < thresh ) ;
    bw2 = imopen(bw2, SE) ;
    
    % remove objects smaller than something deleteObjSmallerThan
    
    bw2 = bw2 & mask ; % use the previous mask, since the fly did not move much
    
    M = abs( double(bw2) - double(bw1) ) ;
    
    % find the center-of-mass of points in which M is changing
    
    %[yind, xind] = find(bw2==true) ;
    [yind, xind] = find(M>0) ;
    
    %if (numel(yind)<Mthresh)
    %    yind = [] ;
    %    xind = [] ;
    %else
    ycm_init = round(mean(yind)) ; % initial guess
    xcm_init = round(mean(xind)) ;
    % to find CM, mask bw2 in the neighborhood where M is changing
    
    mask = false(size(bw2)) ;
    W = 60 ;
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
    
    if (DEBUG_FLAG)
        figure(1) ;
        imshow(bw3) ; title(it) ;  hold on ;
        %imagesc(M) ; axis image ; caxis([-1 1]) ;
        hold on ;
        plot(xcm(c),ycm(c),'r+') ;
        hold off ;
        pause(0.01) ;
    end
    df(c) = sum(M(:)) ;
    bw1 = bw3 ;
    if (~DEBUG_FLAG)
        waitbar(c/N, hbar) ;
    end
end

if (~DEBUG_FLAG)
    close(hbar)
end

clear im1 im2 M

% the fly is outside the FOV if there is a sequence of x frames in which
% df is small enough

%---------------------------------
%addedby SW 7/15/15

out_candidates = abs(df)<=4 ;
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
Wout   = 4 ;
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
out = medfilt1(double(out),10) ;
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
    
    xmid = round ( ( xcm(t1) + xcm(t2) ) /2 );
    ymid = round ( ( ycm(t1) + ycm(t2) ) /2 );
    
    % find if the larger difference is along x or y
    dx = abs(xcm(t1) - xcm(t2)) ;
    dy = abs(ycm(t1) - ycm(t2)) ;
    
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
tin = tvec(good_ind(1)) ;
tout = tvec(good_ind(end)) ;
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
