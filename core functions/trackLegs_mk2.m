function [legsRemovedChopLog, res_t, startLines_t, endLines_t] = trackLegs_mk2(chopLog, allAxlim, cam)
% trackLegs_mk2(all_fly_bw, allAxlim_yz)

% handle only a single view.

Nframes     = chopLog.dim(2) ;
imageHeight = chopLog.dim(3) ;
imageWidth  = chopLog.dim(4) ;
allPoints   = zeros(Nframes*20,3) ; % over allocate
counter     = 0 ;

for f=1:Nframes 
    axlim = round(allAxlim(f,:)) ; 
    axlim(1) = max([axlim(1) 1]) ;
    axlim(2) = min([axlim(2) imageWidth]) ;
    axlim(3) = max([axlim(3) 1]) ;
    axlim(4) = min([axlim(4) imageHeight]) ;
    allAxlim(f,:) = axlim ;
    
    bw = getImage4D(chopLog, cam, f);
    
    
    [xendpoints, yendpoints, bwSkelForShow] = findEndPoints (bw, axlim) ;
    
    if (0)
        figure(1) ; clf ;
        set(gcf,'position',[309         273        1200         600]) ;
        subplot(1,2,1)  ; imshow(bw) ; axis(axlim) ; title(f) ;
        subplot(1,2,2)  ; imshow(bwSkelForShow) ; axis(axlim) ;hold on ;
        plot(xendpoints, yendpoints,'rx','markersize',10,'linewidth',2) ;
    end
    
    np = length(xendpoints) ;
    tvec = zeros(np,1) + f ;
    allPoints(counter+1:counter+np,:) = [  xendpoints, yendpoints, tvec ] ;
    counter = counter + np ;
end
%keyboard ;
allPoints = allPoints(1:counter,:) ;
par.mem = 8 ;  %4 %8
par.dim = 2 ; 
par.good = 30 ; %30
par.quiet= 0 ;
res = track(allPoints, 2.9, par) ;

% separate into particles "traj"
Nparticles = res(end,4) ;
traj = cell(Nparticles, 1) ;
startLines = find([ 1 ; diff(res(:,4)) ]) ;
endLines   = [startLines(2:end)-1 ; size(res,1) ] ;
for k=1:Nparticles
    traj{k} = res( startLines(k):endLines(k), 1:3) ;
end

res_t        = sortrows(res, [3 4]) ; % _t indicates "sorted by time"
startLines_t = find([ 1 ; diff(res_t(:,3)) ]) ;
endLines_t   = [startLines_t(2:end)-1 ; size(res_t,1) ] ;

legsRemovedChopLog = chopLog ;
legWidth = 4 ; %4

% remove the legs
for f=1:Nframes
    try
        coords = res_t( startLines_t(f):endLines_t(f), 1:2) ; % contains xy
    catch
        continue ;
    end
    currbw = getImage4D(chopLog, cam, f);
    axlim = allAxlim(f,:) ;
    for k=1:size(coords,1)
        legTip = coords(k,:) ;
        
        [~, removePixels] = removeOneLeg (currbw, legTip, legWidth, axlim);
        
        if (~isempty(removePixels))
            for q=1:size(removePixels,1)
                currbw(removePixels(q,1), removePixels(q,2)) = false ;
            end
        end
    end
    
    % find the largest CC in currbw
    CC = bwconncomp(currbw,4) ;
    Ncc = length(CC.PixelIdxList) ;
    svec = zeros(Ncc,1) ; % vector containing the size of each connected components
    for j=1:Ncc
        svec(j) = length(CC.PixelIdxList{j}) ;
    end
    [~, idx] = sort(svec,'descend') ;
    
    %[idx1 idx2]  = ind2sub(size(bwtest2), CC.PixelIdxList{idx(1)}) ;
    bw3 = false(size(currbw)) ;
    bw3(CC.PixelIdxList{idx(1)}) = true ;
 
    % save the currbw (legs cut) in legsRemovedChopLog
    
    % store binary images in the sparse 4d structure (see setImage4D)
    i1=cam ; i2=f ;
    ind1vec = legsRemovedChopLog.dim(3)*(i1-1) +  (1:legsRemovedChopLog.dim(3));
    ind2vec = legsRemovedChopLog.dim(4)*(i2-1) +  (1:legsRemovedChopLog.dim(4)) ;
    legsRemovedChopLog.mat(ind1vec, ind2vec) = bw3 ;
    
    if (0)
        figure(3); clf ;
        subplot(1,2,1) ; imshow( getImage4D(chopLog, cam, f) ) ; axis(axlim) ;
        hold on ; plot(coords(:,1), coords(:,2),'ro','linewidth',2) ;
        subplot(1,2,2) ; imshow( bw3) ; axis(axlim) ;
        hold on ; plot(coords(:,1), coords(:,2),'ro','linewidth',2) ;
        hold off ;
    end
    %keyboard ;
end


if (0)% show tracked points on movie
    figure(2)
    for f=1:Nframes
        bw = getImage4D(chopLog, cam, f);
        clf ; 
        subplot(1,2,1) ;
        imshow(bw) ; hold on ; title(f) ;
        for k=1:Nparticles
            % find if the current particle exists in frame f
            if (f>=traj{k}(1,3) && f<=traj{k}(end,3))
                % find f
                indf = find(traj{k}(:,3)==f) ;
                if (~isempty(indf))
                    plot( traj{k}(1:indf,1), traj{k}(1:indf,2),'ro-') ;
                end
            end
        end
        axis(allAxlim(f,:)) ;
        
        subplot(1,2,2) ;
        imshow( getImage4D(legsRemovedChopLog, cam, f)) ;
        axis(allAxlim(f,:)) ;
        
        pause(0.01);
    end
end



return


% ========================================================================

function [xendpoints, yendpoints, bwSkelForShow] = ...
    findEndPoints (bw, axlim)

bwcropped = bw( axlim(3):axlim(4), axlim(1):axlim(2)) ;
% pad with 1 row/column of zeros from all sides,
% to handle the case where the leg is at the boundary of the frame
s   = size(bwcropped) ;
bw2 = false(s+2) ;
bw2(2:s(1)+1, 2:s(2)+1) = bwcropped ;

bwskel = bwmorph(bw2, 'skel', inf) ;
bwspur = bwmorph(bwskel, 'spur', 1) ;
ed     = imsubtract(bwskel, bwspur) ;
[yendpoints, xendpoints] = find(ed==true) ;

% shift back
xendpoints    = xendpoints + axlim(1) - 2;
yendpoints    = yendpoints + axlim(3) - 2;

bwSkelForShow = false(size(bw)) ;
bwSkelForShow(axlim(3):axlim(4), axlim(1):axlim(2)) = ...
    bwskel(2:s(1)+1, 2:s(2)+1) ;

return

% ========================================================================

function axlim = findAxisLimits (bw1) %#ok<DEFNU>

[label, numObjects] = bwlabel(bw1, 4);
objData = regionprops(label, 'area','BoundingBox') ;
if (numObjects~=1)
    objInd = find([objData.Area]==max([objData.Area])) ;
else
    objInd = 1 ;
end
bbox = objData(objInd).BoundingBox ;
axlim = [bbox(1)-2 bbox(1)+bbox(3)+2 bbox(2)-2 bbox(2)+bbox(4)+2] ;
axlim = round(axlim) ;

H = size(bw1,1) ;
W = size(bw1,2) ;
axlim(1) = max([axlim(1) 1]) ;
axlim(2) = min([axlim(2) W]) ;
axlim(3) = max([axlim(3) 1]) ;
axlim(4) = min([axlim(4) H]) ;


return
