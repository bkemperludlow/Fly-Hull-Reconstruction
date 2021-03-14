function [all_fly_bw, body_only_bw, all_fly_thresholds, xcm_pass2, ...
     ycm_pass2, allAxlim, DELTA, with_legs_bw] = ...
    binaryThreshold(bg, cinFilename, tin, tout,twoflies,xcm_guess, ...
     ycm_guess, removeLegsFlag, stopWingsFlag)
% input parameters:
%    *  bg - the background image for this movie (UINT8)
%    *  cinFileName - a string of the full cine file name
%    *  tin, tout - start and end time of tracking. time is in frames such
%      that t=0 is the trigger time. negative time is, hence, possible.
%
% outputs:
%   * all_fly_bw - a sparse4D structure with the binary images of the entire-fly
%     (body + wings + legs). The dimension of all_fly_bw is:
%     [N, 1, height, width]
%   * body_only_bw = a sparse4D structure with the binary images of the
%     body-only. The first and last "DELTA" images are by-definition zeros.
%     The dimension of all_fly_bw is:  [N, 1, height, width]
%   * all_fly_thresholds - the binary thresholds used in each image for
%     binding the entire-fly bw image
%   * xcm_pass2, ycm_pass2 (double, not uint8)- center-of-mass position of the body-only.
%     These are two vectors of size (N-2*DELTA) that correspond to the
%     times between tin+DELTA ... tout-DELTA.
%     N is (tout-tin+1). The coordinates of CM are image coordinates, such
%     that x is the columns and y is rows. (0,0) is top left corner.
%   * allAxlim - Nx4 matrix with an estimated axis-limit for each frame.
%     When showing an image of frame k, you can use: axis(allAxLim(k,:)) to
%     center the view on the fly.
%   * DELTA - the number of frames taken before and after each frame.
%
%
% USAGE for example:
% ------------------
% load Luca_expr2_mov5_bg
% [all_fly_bw, body_only_bw, all_fly_thresholds, xcm_pass2, ycm_pass2, allAxlim, DELTA] = ...
%                 binaryThreshold(bgxy, 'C:\Users\Tsevi\Desktop\TEMP\xy_005.cin', -593, -438) ;
% figure ; hold on ;
% plot(xcm_pass2, ycm_pass2,'ks-','markerfacecolor','g') ;
% plot(xcm_pass2(DELTA+1), ycm_pass2(DELTA+1),'ks-','markerfacecolor','r') ;
% xlabel('x') ; ylabel('y') ; axis equal ; hold off ; grid on ;
%
%
%   note
% ---------
% this function does not cut the fly's legs off. if this is desired, then
% this function should be split into two. the first part should include
% finding the entire-fly bw images. then you should chop the legs in a new
% function or script, based on the existing leg-cutting function (or at
% least find the legs' coordinates). then, the rest of this function should
% become a separate function (pass1, pass2), that will process the
% entire-fly-no-legs images to find the center-of-mass positions
% accurately.

if ~exist('twoflies','var') || isempty(twoflies)
    twoflies =0;
end
if ~exist('removeLegsFlag','var') || isempty(removeLegsFlag)
    removeLegsFlag = true ; 
end
if ~exist('stopWingsFlag','var') || isempty(stopWingsFlag)
    stopWingsFlag = false ; 
end

loud = false ; % if true prints the progess to the screen too.

% get camera name 
%----------------------------------
% get camera name for current file
[~, fn_curr, ~] = fileparts(cinFilename) ; 
fn_curr = strrep(fn_curr, '_', '\_') ; 

warning('off','MATLAB:gui:latexsup:UnableToInterpretTeXString')

DEBUG_FLAG  = false ; % for pass 0
DEBUG_FLAG2 = false ; % for pass 1 false ;
DEBUG_FLAG3 = false ; % false % for pass 2;
DEBUG_FLAG4 = false ; % false % for wing stop fix;

FIT_SCOPE = 100 ;
DEGREE = 2  ; % polynomial degree used for fit

% --------------------------------------------------------
% params only used for "stopWingsFlag" case
SE = strel('disk', 4) ; % structural element used to open image to get body orientation
line_len_sub = 8 ; % when we create line strels, they'll have length = bodyLength - line_len_sub
orient_diff_thresh_low = 4.0 ; % max difference between body orientation in adjacent frames
orient_diff_thresh_high = 170.0 ; % max difference between body orientation in adjacent frames (to account for angle discontinuities)
thick_factor = 2 ; % amount to thicken image by
CC_size_thresh = 250 ; % used in pass 1 to find candidate connected components
centroid_dist_thresh = 2.0 ; % used in pass 1 to find candidate connected components
% --------------------------------------------------------

metaData = getCinMetaData(cinFilename) ;

if(~exist('tin','var'))
    tin = metaData.firstImage ;
    tout = metaData.lastImage ;
else
    if (isempty(tin))
        tin = metaData.firstImage ;
        tout = metaData.lastImage ;
    end
end

N = tout - tin + 1 ; % number of images

% set up output variables
% all_fly_bw   = init4D([N, 1, metaData.height, metaData.width]) ; % second dimension is 1 since this is only for 1 camera. combine later.
% body_only_bw = init4D([N, 1, metaData.height, metaData.width]) ; % second dimension is 1 since this is only for 1 camera. combine later.
all_fly_bw   = init4D([1, N, metaData.height, metaData.width]) ; % second dimension is 1 since this is only for 1 camera. combine later.
body_only_bw = init4D([1, N, metaData.height, metaData.width]) ; % second dimension is 1 since this is only for 1 camera. combine later.
all_fly_thresholds = zeros(N,1) ;

allAxlim = zeros(N,4) ;

W = 40 ;
maskWindow = 60 ; 

%xcm = round(xcm) ; ycm=round(ycm) ;

if (~DEBUG_FLAG)
    hbar = waitbar(0,['binaryThreshold: finding whole-body bw images for ' ...
        fn_curr ]) ;
end

% ------------------------------------------------------------------------
%% FIND THE ENTIRE FLY BW IMAGE IN EACH FRAME
% this part should be a sepapated function, after which the legs should be
% removed in a second function. the remaining part of this code (pass1,2)
% should be a third function.

cindata = myOpenCinFile(cinFilename) ;

c=0 ;

if twoflies>0
    xboxCell = cell(1,twoflies) ;
    yboxCell = cell(1,twoflies) ;
    for j = 1:twoflies
        im1=myReadCinImage(cindata, tin) ;
        figure;
        imshow(im1);
        [xbox,ybox]=ginput(4);
        xboxCell{1,j} = xbox ;
        yboxCell{1,j} = ybox ;
    end
end



for t=tin:tout
    
    c=c+1 ;
    
    if (loud)
        disp(['binaryThreshold: finding whole-body bw images for ' ...
            fn_curr ', frame ' num2str(c) ' / ' num2str(Nimages)]) ;
    end
    %im1 = ReadCineFileImage(cinFilename, t, false);
    im1 = myReadCinImage(cindata, t) ;
    
    im2 = imsubtract(bg, im1) ;
    
    %{
    % find bw-threshold level based on a sub image (same results, don't bother)
    x1 = max( [xcm(c)-W, 1]) ;
    x2 = min( [xcm(c)+W, metaData.width]) ;
    y1 = max( [ycm(c)-W, 1]) ;
    y2 = min( [ycm(c)+W, metaData.height]) ;
    
    sub_im2 = im2(y1:y2, x1:x2) ;
    %figure(3) ; imshow(sub_im2);
    level = graythresh(sub_im2)  ;
    c=c+1 ;
    %}
     
  
    if twoflies>0
        for j = 1:twoflies
            im2(int16(min(yboxCell{1,j})):int16(max(yboxCell{1,j})),int16(min(xboxCell{1,j})):int16(max(xboxCell{1,j})))...
                =zeros(int16(1-min(yboxCell{1,j})+max(yboxCell{1,j})),int16(-min(xboxCell{1,j})+max(xboxCell{1,j})+1));
        end
    end
    
    level = graythresh(im2) * 0.6 ; 
    all_fly_thresholds(c) = level ;
    
    bw2   = imbinarize(im2, level) ;
    
    ind_temp = t - metaData.firstImage ; 
    xcm_curr = xcm_guess(ind_temp) ; 
    ycm_curr = ycm_guess(ind_temp) ; 
    mask = false(size(bw2)) ;
    x1 = int16(max([xcm_curr - maskWindow, 1])) ;
    x2 = int16(min([xcm_curr + maskWindow, metaData.width])) ; 
    y1 = int16(max([ycm_curr - maskWindow, 1])) ;
    y2 = int16(min([ycm_curr + maskWindow, metaData.height])) ;
    mask(y1:y2, x1:x2) = true ;
    
    bw2 = bw2 & mask ; 
    % find center of mass of the largest cc of bw2
    % find largest CC
    CC  = bwconncomp(bw2);
    Ncc = length(CC.PixelIdxList) ;
    
   
    svec = zeros(Ncc,1) ; % vector containing the size of each connected components
    
    for j=1:Ncc
        svec(j) = length(CC.PixelIdxList{j}) ;
    end
    [~, idx] = sort(svec,'descend') ;
    
    %[idx1 idx2]  = ind2sub(size(bwtest2), CC.PixelIdxList{idx(1)}) ;
    bw3 = false(size(bw2)) ;
    %if t == 151
    %    disp('blerg') ;
    %end
    bw3(CC.PixelIdxList{idx(1)}) = true ;
    bw3 = imfill(bw3, 'holes') ;
    
    % ----
    % store binary images in the sparse 4d structure (see setImage4D)
    i1=1 ; i2=c ;
    ind1vec = all_fly_bw.dim(3)*(i1-1) +  (1:all_fly_bw.dim(3));
    ind2vec = all_fly_bw.dim(4)*(i2-1) +  (1:all_fly_bw.dim(4)) ;
    all_fly_bw.mat(ind1vec, ind2vec) = bw3 ;
    % stored
    % ---
    
    [idx1, idx2]  = ind2sub(size(bw3), CC.PixelIdxList{idx(1)}) ;
    
    %{
    % -- following part did not work. find body-only based on motion
    
    % apply the adaptive threshold for the fly-pixels only to find
    % body-only threshold
    pixvals     = im2(CC.PixelIdxList{idx(1)});
    level_body  = graythresh(pixvals) * 1.0;
    bw_body     = im2bw(im2, level_body) ;
    %}
    
    % find axis limits for display
    m1 = mean(idx1) ;
    m2 = mean(idx2) ;
    axlim = [m2, m2, m1, m1] + [-1 1 -1 1]*W;
    allAxlim(c,:) = axlim ;
    
    
    if (DEBUG_FLAG)
        % for display purposes
        BWoutline = bwperim(bw3);
        BWoutline_body = false & bwperim(bw_body);
        rbg = zeros([size(im2) 3],'uint8') ;
        im2temp = im2.*(uint8(~BWoutline)).*(uint8(~BWoutline_body)) ;
        rgb(:,:,1) = im2temp + uint8(BWoutline)*uint8(255) ;
        rgb(:,:,3) = im2temp + uint8(BWoutline_body)*uint8(255) ;
        rgb(:,:,2) = im2temp ;
        
        figure(1) ;
        subplot(1,4,1) ; imshow(rgb) ; title(t) ; axis(axlim) ;
        subplot(1,4,2) ; imshow(bw3) ; title(level) ; axis(axlim) ;
        subplot(1,4,3) ; imshow(bw_body) ; title(level_body) ; axis(axlim) ;
        subplot(1,4,4) ; imhist(im2) ; set(gca,'xlim',[0 255],'ylim',[0 30]) ;
        hold on ; plot(level*[255 255],[0 50],'r-','linewidth',2) ;
        pause(0.05);
    else
        waitbar(c/N, hbar) ;
    end
end

if (~DEBUG_FLAG)
    close(hbar);
else
    figure(1); clf ;
end


%% REMOVE LEGS (if not opening body for stopWings)
% ------------
if (removeLegsFlag) && (~stopWingsFlag)
    % if everything works, we should now work with the chopped legs pics
    with_legs_bw = all_fly_bw ; 
    [all_fly_bw, ~, ~, ~] = trackLegs_mk2(all_fly_bw, allAxlim, 1) ; 
end
%[legsRemovedChopLog, ~, ~, ~] = trackLegs_mk2(body_only_bw, allAxlim, 1) ; % XXX
%body_only_bw = legsRemovedChopLog ; % if everything works, we should now work with the chopped legs pics


%% FIND body-only images - PASS 1 (based on findBodyOnly_mk2.m)
DELTA = 18 ;

%xcm_pass1 = zeros(N-2*DELTA,1) ;
%ycm_pass1 = zeros(N-2*DELTA,1) ;
xcm_pass1 = zeros(N,1) + NaN ; % use N and not (N-2*DELTA) to be compatible with all_fly_bw and all_fly_thresholds
ycm_pass1 = zeros(N,1) + NaN ;
%N_bodyPix = nan(N,1) ; 

if (~DEBUG_FLAG2)
    hbar = waitbar(0,['binaryThreshold: finding body CM for ' fn_curr ...
        ' (pass 1)...']) ;
end

cc  = 0 ; % used only for hbar
Ndd = N - 2*DELTA ;

for t = tin+DELTA : tout-DELTA
    cc=cc+1 ; % for the hbar only
    
    if (loud)
        disp(['binaryThreshold: finding body CM for ' fn_curr ...
            ' (pass 1), frame ' num2str(cc) ' / ' num2str(Ndd)]) ;
    end
    combbw  = true(size(bw3)) ; %not actually used here
    comb2 = zeros(size(bw3),'uint8') ;
    %comb3 = zeros(size(bw3),'double') ;
    ind=0 ;
    
    % read images before and after time t
    for t2= (t-DELTA):(t+DELTA) % PASS 1
        k = t2 - tin + 1 ; %convert to frame number
        bw = getImage4D(all_fly_bw, 1,k) ;
        ind=ind+1 ;
        %stack(ind,:,:) = bw ;
        combbw = combbw & bw ;
        comb2 = comb2 + uint8(bw) ;
        
        if (t2==t && DEBUG_FLAG2)
            %im1 = ReadCineFileImage(cinFilename, t2, false);
            im1 = myReadCinImage(cindata, t2) ;
            keep = imsubtract(bg, im1) ;
            clear im1 ;
        end
    end
    %comb3 = comb3 / (2*DELTA+1) ;
    
    % check if there are "enough" pixels that appear in all frames
    % (t-DELTA ... t+DELTA). if not, reduce the threshold of (2*DELTA+1)
    % such that "enough" pixels are included.
    cont = true;
    thresh = 2*DELTA + 1 ;
    while (cont)
        comb2_bw = (comb2 >= thresh) ;
        % look for largest CC
        CC  = bwconncomp(comb2_bw);
        
        if (CC.NumObjects == 0)
            thresh = thresh - 1 ;
            continue ;
        end
        
        if (length(CC.PixelIdxList)>1)
            % take only largest CC
            Ncc = length(CC.PixelIdxList) ;
            svec = zeros(Ncc,1) ; % vector containing the size of each connected components
            for j=1:Ncc
                svec(j) = length(CC.PixelIdxList{j}) ;
            end
            [~, idx] = sort(svec,'descend') ;
        else
            idx = 1;
        end
        
        comb2_bw = false(size(comb2_bw)) ; % rest comb2_bw
        comb2_bw(CC.PixelIdxList{idx(1)}) = true ; % with the largest CC

        % if largest CC is still not large enough
        if (sum(comb2_bw(:)) < 400)
            thresh = thresh - 1 ;
            continue ;
        else
            cont = false ;
        end
    end
    
    % ----------------------------
    % get frame number 
    c2 = t - tin - DELTA + 1 ;
    c = t - tin + 1 ;
    
    % ----------------------------------------
    % try to deal with wing stopping issue
    
    if (stopWingsFlag) && (cc > (DELTA + 1)) && (t > 0) 
        % because we don't have body-only images to align yet, for the
        % first pass we'll just try to kill off any wing contributions with
        % a watershed segmentation
        comb2_bw = imopen(comb2_bw, SE) ; 
        
        D = bwdist(~comb2_bw) ;
        D = -D;
        D(~comb2_bw) = Inf;
        L = watershed(D);
        L(~comb2_bw) = 0;
        % figure ; imshow(label2rgb(L))
        
        if (numel(unique(L(:))) > 2)
            L_CC = bwconncomp(L) ;
            cc_sizes = cellfun(@(y) length(y), L_CC.PixelIdxList)' ;
            L_rprops = regionprops(L_CC, 'Centroid') ;
            cc_centroids = vertcat(L_rprops.Centroid) ;
            
            c = t - tin + 1 ;
            centroid_dist = myNorm(cc_centroids - ...
                [xcm_pass1(c-1), ycm_pass1(c-1)]) ;
            
            sizeCheck = (cc_sizes >= CC_size_thresh) ;
            centroidCheck = (centroid_dist < centroid_dist_thresh) ;
            
            if (sum(sizeCheck & centroidCheck) >= 1)
                comb3_bw = false(size(comb2_bw)) ;
                comb3_bw(L_CC.PixelIdxList{sizeCheck & ...
                    centroidCheck}) = true ;
                % figure ; imshowpair(comb2_bw, comb3_bw)
                comb2_bw = comb3_bw ;
            end
        end
    end
    
    % ----------------------------------------
    % store first-pass values of body cm
    [idx1, idx2]  = ind2sub(size(comb2_bw), CC.PixelIdxList{idx(1)}) ;
    
    xcm_pass1(c) = mean(double(idx2)) ; % previously used c2 instead of c.
    ycm_pass1(c) = mean(double(idx1)) ;
    %N_bodyPix(c) = length(idx1) ; 
    
    if (DEBUG_FLAG2)
        
        figure(1) ;
        axlim = allAxlim(c,:);
        imshow(comb2_bw); axis(axlim) ; title('after erosion') ;
        %subplot(1,3,1) ; imshow(combbw) ; axis(axlim) ;
        %subplot(1,3,2) ; imagesc(comb2); axis image ; axis(axlim) ; colorbar ;
        %subplot(1,3,3) ; imagesc(comb3) ;axis image ; axis(axlim) ; colorbar ;
        %colormap jet ;
        
        figure(2), imshow(keep) ; axis(axlim) ; colormap jet ; colorbar ; impixelinfo ; title('original') ;
        hold on;
        d1 = DELTA + 1 ;
        plot(xcm_pass1(d1:c), ycm_pass1(d1:c),'ko-','markerfacecolor','w','markersize',3) ;
        %figure(3), imagesc(comb3) ; axis image ; axis(axlim) ; colormap jet ; colorbar ; impixelinfo ; title('combined uint8') ;
        figure(4), imagesc(comb2) ; axis image ; axis(axlim) ; colormap jet ; colorbar ; impixelinfo ; title('combined binary') ;
        hold on ;
        plot(xcm_pass1(d1:c), ycm_pass1(d1:c),'ko-','markerfacecolor','w','markersize',3) ;
        pause(.01) ;
        hold off ;
    else
        waitbar(cc/Ndd, hbar) ;
    end
end


if (~DEBUG_FLAG2)
    close(hbar);
end

%figure ;
%plot(xcm_pass1, ycm_pass1,'bo-') ; axis equal ;


%% FIND BODY C.M. PASS-2 - shift the images to overlap, then do the same trick as in pass-1
% (see code in segmentWingSingleFrame.m)

if (~DEBUG_FLAG3) && (~DEBUG_FLAG4)
    hbar = waitbar(0,['binaryThreshold: finding body CM for ' fn_curr ...
        ' (pass 2)...']) ;
end

if (DEBUG_FLAG4) && (stopWingsFlag)
    h_fig = figure ;
    ax = gca ;
    hold on
    h_imshow = imshow(false(metaData.height, metaData.width), [],...
        'Parent', ax) ;
    h_ell = plot(ax, NaN, NaN, 'r-','LineWidth',1.0) ;
    h_dir = plot(ax, NaN, NaN, 'g-','LineWidth',1.0) ; 
    t_grid = linspace(0,2*pi,50);
    %h_dir_cap = plot(NaN, NaN, 'ko','MarkerFaceColor', 'g') ;
end
cc=0 ; % used only for hbar

xcm_pass2 = zeros(N,1) + NaN  ;
ycm_pass2 = xcm_pass2 ;


for t = tin+DELTA : tout-DELTA
    cc = cc + 1 ;
    if (loud)
        disp(['binaryThreshold: finding body CM for ' fn_curr ...
            ' (pass 2), frame ' num2str(cc) ' / ' num2str(Ndd)]) ;
    end
    t1 = max([t-FIT_SCOPE ; tin+DELTA]) ;
    t2 = min([t+FIT_SCOPE ; tout-DELTA]) ;
    
    tvec = (t1:t2)' ;
    tind = tvec - (tin+DELTA) + 1  + DELTA;
    
    xvec = xcm_pass1(tind) ;
    yvec = ycm_pass1(tind) ;
    
    px = polyfit(tvec, xvec, DEGREE) ;
    py = polyfit(tvec, yvec, DEGREE) ;
    
    %{
    figure ; hold on ;
    plot(xvec, yvec,'ko') ;
    plot( polyval(px,tvec), polyval(py, tvec),'r-') ;
    axis equal ; grid on ; box on ; title(t) ;
    %}
    
    tvec = (t-DELTA:t+DELTA)' ;
    
    % translate images - see segmentWingsSingleFrame lines 148 and on.
    cmvec  = [ polyval(px,tvec) polyval(py,tvec) ] ;
    cmcurr = cmvec(DELTA+1,:) ;
    drvec  = round(cmvec - repmat(cmcurr,2*DELTA+1,1)) ;
    
    combbw  = true(size(bw3)) ;
    comb2   = zeros(size(bw3),'uint8') ;
    %comb3 = zeros(size(bw3),'double') ;
    ind=0 ;
    
    for t2= (t-DELTA):(t+DELTA) % PASS 2
        ind=ind+1 ;
        k = t2 - tin + 1 ;
        bw = getImage4D(all_fly_bw, 1,k) ;
        bworig = bw ;
        %bw = translateImageBW(bworig, drvec(ind,:)) ;
        bw = translateImageBW2(bworig, drvec(ind,:)) ;
        
        %checksum = sum ( abs(double(bw(:)) - double(bb(:)))) ;
        %if(checksum~=0)
        %    disp('problem')
        %    keyboard ;
        %end
        combbw = combbw & bw ;
        comb2 = comb2 + uint8(bw) ;
        
        if (t2==t && DEBUG_FLAG3)
            im1 = ReadCineFileImage(cinFilename, t2, false);
            keep = imsubtract(bg, im1) ;
            clear im1 ;
        end
    end
    
    comb2_bw = (comb2 == 2*DELTA+1) ;
    
    % ------------------------------------------------------------------
    % if the wings are stopping mid flight, do some processing to remove
    % them from body pixels
    if (stopWingsFlag) && (cc > 1)
        comb2_bw_open = imopen(comb2_bw, se_line) ;
        comb2_bw = comb2_bw & bwmorph(comb2_bw_open,'thicken',thick_factor) ;
    end
    
    % ----------------------------------------------
    % get largest connected component
    CC  = bwconncomp(comb2_bw);
    if (length(CC.PixelIdxList)~=1)
        % take only largest CC
        Ncc = length(CC.PixelIdxList) ;
        svec = zeros(Ncc,1) ; % vector containing the size of each connected components
        for j=1:Ncc
            svec(j) = length(CC.PixelIdxList{j}) ;
        end
        [~, idx] = sort(svec,'descend') ;
        
        comb2_bw = false(size(comb2_bw)) ;
        
        comb2_bw(CC.PixelIdxList{idx(1)}) = true ;
        %disp(t)
    
    else
        idx = 1 ;
    end
    
    c = t - tin + 1 ;
    
    % --------------------------------------------------------------------
    % get orientation of body to update strel (only in stopWingsFlag case)
    % get orientation of body 
    if stopWingsFlag
        bw_open = imopen(comb2_bw,SE) ;
        % in case opening kills off image
        if (sum(bw_open(:)) < 1)
            strel_radius_temp = 4 ; 
            while sum(sum(bw_open(:)) < 1)
                strel_radius_temp = strel_radius_temp - 1; 
                se_temp = strel('disk',strel_radius_temp) ; 
                bw_open = imopen(comb2_bw,se_temp) ;
            end
        end
        % get region properties
        rprops = regionprops(bw_open, 'Orientation','Area','MajorAxisLength',...
            'MinorAxisLength', 'Centroid') ;
        [~, max_ind] = max([rprops.Area]) ;
        
        % update our estimate of the fly's current orientation
        if (cc > 1)
            orientation_diff = abs(rprops(max_ind).Orientation - ...
                orientation) ;
            if (orientation_diff < orient_diff_thresh_low) || ...
                    (orientation_diff > orient_diff_thresh_high)
                orientation = rprops(max_ind).Orientation ;
                line_len = round(rprops(max_ind).MajorAxisLength) - ...
                    line_len_sub ;
            else
                line_len = 14 ;
            end
        else
            orientation = rprops(max_ind).Orientation ;
            line_len = round(rprops(max_ind).MajorAxisLength) - ...
                line_len_sub  ;
        end
        
        % make line structuring element
        se_line = strel('line', line_len, orientation ) ;
        
    end
    % ----------------------------------------------------------------
    % SAVE comb2_bw (this is what the functin outputs)
    % store binary images in the sparse 4d structure (see setImage4D)
    i1=1 ; i2=c ;
    ind1vec = body_only_bw.dim(3)*(i1-1) +  (1:body_only_bw.dim(3));
    ind2vec = body_only_bw.dim(4)*(i2-1) +  (1:body_only_bw.dim(4)) ;
    body_only_bw.mat(ind1vec, ind2vec) = comb2_bw ;
    % stored.
    
    c2 = t - tin - DELTA + 1 ;
    [idx1, idx2]  = ind2sub(size(comb2_bw), CC.PixelIdxList{idx(1)}) ;
    xcm_pass2(c) = mean(double(idx2)) ; % previously used c2 intead of c.
    ycm_pass2(c) = mean(double(idx1)) ;
    
    if (DEBUG_FLAG3)
        
        axlim = allAxlim(c,:);
        %figure(1) ;
        %imshow(comb2_bw); axis(axlim) ; title('after erosion') ;
        %subplot(1,3,1) ; imshow(combbw) ; axis(axlim) ;
        %subplot(1,3,2) ; imagesc(comb2); axis image ; axis(axlim) ; colorbar ;
        %subplot(1,3,3) ; imagesc(comb3) ;axis image ; axis(axlim) ; colorbar ;
        %colormap jet ;
        
        figure(2), imshow(keep) ; axis(axlim) ; colormap jet ; colorbar ; impixelinfo ; title('original') ;
        hold on;
        d1 = DELTA + 1 ;
        plot(xcm_pass2(d1:c), ycm_pass2(d1:c),'ko-','markerfacecolor','w','markersize',3) ;
        %figure(3), imagesc(comb3) ; axis image ; axis(axlim) ; colormap jet ; colorbar ; impixelinfo ; title('combined uint8') ;
        figure(4), imagesc(comb2) ; axis image ; axis(axlim) ; colormap jet ; colorbar ; impixelinfo ; title('combined binary') ;
        hold on ;
        plot(xcm_pass2(d1:c), ycm_pass2(d1:c),'ko-','markerfacecolor','w','markersize',3) ;
        pause(.01) ;
        %keyboard  ;
    end
    if (DEBUG_FLAG4) && (stopWingsFlag)
        % update image
        set(h_imshow, 'CData', comb2_bw)
        
        % plot ellipse
        a = rprops(max_ind).MajorAxisLength/2;
        b = rprops(max_ind).MinorAxisLength/2;
        Xc = rprops(max_ind).Centroid(1);
        Yc = rprops(max_ind).Centroid(2);
        phi = deg2rad(-rprops(max_ind).Orientation);
        phi2 = deg2rad(-orientation) ; 
        % ellipse from region props 
        x = Xc + a*cos(t_grid)*cos(phi) - b*sin(t_grid)*sin(phi);
        y = Yc + a*cos(t_grid)*sin(phi) + b*sin(t_grid)*cos(phi);
        % ellipse using current orientation estimate
        x2 = Xc + a*cos(t_grid)*cos(phi2) - b*sin(t_grid)*sin(phi2);
        y2 = Yc + a*cos(t_grid)*sin(phi2) + b*sin(t_grid)*cos(phi2);
        
        % update ellipses
        set(h_ell, 'XData', x, 'YData', y) ; 
        set(h_dir, 'XData', x2, 'YData', y2) ; 
        title(ax, t)
        pause(0.05)
    end
    if ~(DEBUG_FLAG3) && ~(DEBUG_FLAG4)
        waitbar(cc/Ndd, hbar) ;
    end
    
    
end

if (~DEBUG_FLAG3) && ~(DEBUG_FLAG4)
    close(hbar)
end

% ------------------------------------------------------------------------
%% REMOVE LEGS (if the case of stopWings
if (removeLegsFlag) && (stopWingsFlag)
   with_legs_bw = all_fly_bw ;  
   [all_fly_bw, ~, ~, ~] = trackLegs_mk2(all_fly_bw, allAxlim, 1) ; 
end

if ~exist('with_legs_bw','var') && (nargout > 7)
    with_legs_bw = [] ; 
end
% -----------------------------------
%% Close cine and exit
myCloseCinFile(cindata) ;

warning('on','MATLAB:gui:latexsup:UnableToInterpretTeXString')

return

%% ---------------------------------------------------------------------


function B = translateImageBW2(A, dr)
%dr = (nextCM - currCM ); [dx dy] = [dcol drow]
drow = dr(2) ;
dcol = dr(1) ;
S = size(A) ;
%B = zeros(S,'uint8') ;
B = false(S); % zeros(S,'double') + NaN;

tmp = (1:S(1))+drow ;
rowind = tmp(tmp>=1 & tmp<=S(1));
tmp = (1:S(2))+dcol ;
colind = tmp(tmp>=1 & tmp<=S(2));

B(rowind-drow, colind-dcol) = A(rowind, colind) ;

% %disp('bla bla changed B from zeros to NaN''s in translateImage()') ;
%
% % can replace this loop by a more elegant matrix-style command
% for ix=1:S(2)
%     if (ix+dcol<=S(2) && ix+dcol>=1)
%         for iy=1:S(1)
%             if (iy+drow<=S(1) && iy+drow>=1)
%                 B(iy,ix) = double(A(iy+drow, ix+dcol)) ; % xxx
%             end
%         end
%     end
% end
return


function B = translateImageBW(A, dr)
%dr = (nextCM - currCM ); [dx dy] = [dcol drow]
drow = dr(2) ;
dcol = dr(1) ;
S = size(A) ;
%B = zeros(S,'uint8') ;
B = false(S); % zeros(S,'double') + NaN;

%disp('bla bla changed B from zeros to NaN''s in translateImage()') ;

% can replace this loop by a more elegant matrix-style command
for ix=1:S(2)
    if (ix+dcol<=S(2) && ix+dcol>=1)
        for iy=1:S(1)
            if (iy+drow<=S(1) && iy+drow>=1)
                B(iy,ix) = double(A(iy+drow, ix+dcol)) ; % xxx
            end
        end
    end
end


return

