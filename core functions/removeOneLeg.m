function [ bw, removePixels] = removeOneLeg (bwcurr, legTip, legWidth, axlim)

Nb = 30 ; %30
plotFlag = false ;

% this is the length the algorithm travels along the boundary. it should be
% longer than a leg.

if (~exist('axlim','var'))
    s = size(bwcurr) ;
    axlim = [1 s(2) 1 s(1)] ;
end

W = 4; % floor(legWidth/2)*2 ;

bwcropped = bwcurr( axlim(3):axlim(4), axlim(1):axlim(2)) ;
% pad with W row/column of zeros from all sides,
% to handle the case where the leg is at the boundary of the frame
s   = size(bwcropped) ;
bw2 = false(s+W*2) ;
bw2( (W+1):s(1)+W, (W+1):s(2)+W) = bwcropped ;

bwthick = bwmorph(bw2,'thicken') ; 
bwthick = bwmorph(bwthick,'diag',1) ;

%bwskel = bwmorph(bwthick, 'skel',Inf) ;
bwskel = bwmorph(bwthick, 'thin',Inf) ;
bwspur = bwmorph(bwskel, 'spur',1) ;

bwedges = bwskel - bwspur ;

[ir, ic] = find(bwedges==true) ;

% find the on-pixel in bwedges which is closest to legTip
r1 = legTip(2) - axlim(3) + (W+1) ;
c1 = legTip(1) - axlim(1) + (W+1) ;

dst2 = (ir-r1).^2 + (ic-c1).^2 ;
[~, ind] = min(dst2) ;

legTip = [ic(ind(1)) ir(ind(1))] ;

% -----------------------------------------
% FIND PIXELS ALONG THE SKELETON OF THE LEG
% -----------------------------------------

% go over the neighbors of the tip pixel until we get to a pixel that has
% more than two neighbors

rowList = zeros(100,1) ; % allocate 100, then trim later
colList = zeros(100,1) ;

c=0 ;
cont = true ;
currRow = legTip(2) ;
currCol = legTip(1) ;

% transform to cropped frame of reference
%currRow = currRow - axlim(3) + (W+1) ;
%currCol = currCol - axlim(1) + (W+1) ;


% the leg is given in the non-thickened image, and hence the tip is in the
% wrong place in the thickened image

while (cont)
    % need to use the zero-padded version
    neigh = bwskel(currRow-1:currRow+1, currCol-1:currCol+1) ;
    neigh(2,2) = false ; % remove current pixel (will be used later to find neighbor)
    s = sum(neigh(:))  ; %
    
    if (s>=2 || s==0)
        cont = false ;
    elseif (s==1)
        c=c+1 ;
        rowList(c) = currRow ;
        colList(c) = currCol ;
        bwskel(currRow, currCol) = false ; % delete current pixel on bwskel
        % move to the neighbour
        [ir ic] = find(neigh==true) ;
        currRow = currRow + ir(1) - 2 ;
        currCol = currCol + ic(1) - 2 ;
    end
end

%{
% shift back
xendpoints    = xendpoints + axlim(1) - 2;
yendpoints    = yendpoints + axlim(3) - 2;
%}

% if there are two few pixels in skeletonized leg, abort.
if (c<=2)
    disp('removeOneLeg_mk2: leg it too short. Aborting.')
    bw = [] ;
    removePixels = [] ;
    return
end


rowList = rowList(1:c) ;
colList = colList(1:c) ;

% ROWLIST AND COLLIST CONTAIN THE LEG-SKELETON PIXELS

% -----------------------
% FIND THE IMAGE BOUNDARY
% -----------------------

bwboundary = bwmorph(bwthick, 'remove') ;

% check that the first pixel in (rowList,colList) is TRUE in bwboundary.
% if not add a pixel from the boundary

r1 = rowList(1) ;
c1 = colList(1) ;
if (bwboundary(r1,c1) == false )
    % calculate slope to the next pixel and go the opposite way
    r2 = rowList(2) ;
    c2 = colList(2) ;
    dr = r2 - r1 ;
    dc = c2 - c1 ;
    newr = r1 - dr ;
    newc = c1 - dc ;
    
    % if this does not work, find the closest on-pixel in the neighborhood
    % of (newr, newc)
    if ( bwboundary(newr, newc) == false )        
        neigh = bwboundary(newr-1:newr+1 , newc-1:newc+1) ;        
        [rows cols] = find(neigh) ;
        if (isempty(rows))
            %disp('error in removeOneLeg_mk2 - need to take larger neighborhood?') ;
            disp('removeOneLeg_mk2 failed. Aborting without removing leg.') ;
            bw = [] ;
            removePixels = [] ;
            return 
        end
        rows = rows + newr - 2 ;
        cols = cols + newc - 2 ;        
        dst2 = (rows-newr).^2 + (cols-newc).^2 ;
        [~, ind] = min(dst2) ;
        newr = rows(ind(1)) ;
        newc = cols(ind(1)) ;
    end
    %{
    if ( bwboundary(newr, newc) == false )
        newr = r1-dr ;
        newc = c1 ;
    end
    if ( bwboundary(newr, newc) == false )
        newr = r1 ;
        newc = c1-dr;
    end    
    %}
    if ( bwboundary(newr, newc) == false )
        disp('error')
        keyboard ;
    end
    rowList = [newr ; rowList] ;
    colList = [newc ; colList] ;
end

if (plotFlag) ;
    hfig=figure(8);
    imshow(bwboundary,'initialmagnification',400) ; title('boundary');
    hold on ;
    plot(colList, rowList,'r+') ;
    hold off ;
end

%{
start from the tip pixel, such that it is on the thickened boundary
this pixel should have two closest neighbors. find them, then for each one do the
same (two at a time), until the whole leg boundary is found.

for each pair of bounary pixels, calculate their distance.

When the distance goes above some threshold OR when the distance derivative
grows consistently, declare we found the base of the leg.

remove the leg (somehow, it should not be so hard).
%}

list1 = zeros(Nb,2) + NaN ; % NaN will indicate we stopped the search along the boundary
list2 = zeros(Nb,2) + NaN ;

r1 = rowList(1) ;
c1 = colList(1) ;

% initial with tip pixel
list1(1,:) = [ r1 c1 ] ;
list2(1,:) = [ r1 c1 ] ;

% find the tip pixel's two closest neighbors on the boundary

bwboundary(r1,c1)  = false ;  % remove current pixel
neigh = bwboundary(r1-1:r1+1 , c1-1:c1+1) ;

[rows cols] = find(neigh) ;

rows = rows + r1 - 2 ;
cols = cols + c1 - 2 ;

dst2 = (rows-r1).^2 + (cols-c1).^2 ;
[dst2 ind] = sort(dst2) ;

% take one closest pixels + offest
list1(2,:) = [rows(ind(1)) cols(ind(1))] ; % + [ r1 c1] - 2;

% find the second closest pixel that is also farthest from list1(2,:)
dst2 = dst2(2:end) ;
ind  = ind(2:end) ;

if (isempty(dst2))
    disp('removeOneLeg_mk2 failed. Aborting without removing leg.') ;
    bw = [] ;
    removePixels = [] ;
    return
end

d = find(dst2==dst2(1)) ;
if (numel(d)==1)
   % if these is only one "second closest", take it
    list2(2,:) = [rows(ind(1)) cols(ind(1))] ; % + [ r1 c1] - 2;
else
    % d contains indices to "dst2" and "ind" with the closest pixels 
    dst3 = ( rows(ind(d)) - list1(2,1) ).^2 + (cols(ind(d))- list1(2,2) ).^2 ;
    % among these, find the farthest from list1(2,:)
    [~, ind3] = max(dst3) ;
    ind3 = ind3(1) ;
    list2(2,:) = [ rows(ind(d(ind3))) cols(ind(d(ind3))) ] ;
end

% list1 and list2 contain the tip pixel in their first row
% each second row has a pixel of each side along the boundary

% combine lists to make loop more elegant
combList = [ list1 list2 ] ;



%distVec = zeros(1,Nb) ; % not needed anymore

% walk on the boundary 
for k=3:Nb
    
    % remove curret TWO pixels from bwboundary  
    for listIndex = [1 2]
        j = listIndex*2 - 1 ;
        r1 = combList(k-1, j) ;
        c1 = combList(k-1, j+1) ;
        if (isnan(r1))
            continue ;
        end
        bwboundary(r1,c1) = false ;
    end
    
    for listIndex = [1 2]
              
        j = listIndex*2 - 1 ;        
        r1 = combList(k-1, j) ;
        c1 = combList(k-1, j+1) ;
        
        if (isnan(r1))
            continue ;
        end
        
        %bwboundary(r1,c1) = false ; % done before in a separate loop
        
        neigh = bwboundary(r1-1:r1+1 , c1-1:c1+1) ;
        [rows cols] = find(neigh) ;

        rows = rows + r1 - 2 ;
        cols = cols + c1 - 2 ;

        if (numel(rows)==1)
           combList(k, j) = rows ;
           combList(k, j+1) = cols ;
        else
            % if more than one neighbor, find the closest
            dst2 = (rows-r1).^2 + (cols-c1).^2 ;
            if (isempty(dst2)) 
                disp('cont...' ) 
                %keyboard ;
                continue ;
            end
            
            ind = find(dst2==min(dst2));
            
            if (numel(ind)==1)
                combList(k, j) = rows(ind(1)) ;
                combList(k, j+1) = cols(ind(1)) ;
            else
                ne = length(ind) ;
                % if two are closest, find the one which is in the general
                % direction of the previous lines
                % previous slope is (dr,dc)
                % this part assumes k>2 which is ok (see loop above)
                if (k>4)
                    dk = 4 ;
                else
                    dk = k-1 ;
                end
                    
                dr = ( combList(k-1,j) - combList(k-dk,j) ) / (dk-1) ;
                dc = ( combList(k-1,j+1) - combList(k-dk,j+1)) / (dk-1) ;
                % find slopes for the pixels found in "ind"
                drind = rows(ind) - combList(k-1,j) ;
                dcind = cols(ind) - combList(k-1,j+1) ;
                dotvec = zeros(ne,1) ;
                v=[dc dr] ;
                for q=1:ne
                    dotvec(q) = dot ( v, [dcind(q) drind(q)]) ;
                end
                selectedInd = find(dotvec==max(dotvec)) ;
                
                if (numel(selectedInd)==1)
                    combList(k, j) = rows(ind(selectedInd)) ;
                    combList(k, j+1) = cols(ind(selectedInd)) ;
                else
                    disp(['k=' num2str(k) ' - Stop traveling along boundary.']);
                    %eyboard ;
                end                
                
                %keyboard ;
            end
            
        end
        
    end
    
    if (plotFlag)
        figure(hfig) ; clf ;
        imshow(bwboundary,'initialmagnification',400) ;
        hold on ;
        plot(colList, rowList,'r+') ;
        plot( combList(1:k, [2 4]) , combList(1:k,[1 3]),'c+');
        hold off ;
        %pause(0.1) ;
    end
    
end

% for each point on the midline (rowList colList)
% find the closest two points along the two parts of combList
Nc = size(rowList,1) ;
pairsList = zeros(Nc,4) ; % [r1 c1 r2 c2]
pairsDist = zeros(Nc,1) ;

% the following part inheretnly ignores any NaN entries in combList.

for k=1:Nc
    r1 = rowList(k) ;
    c1 = colList(k) ;
    for listIndex = [1 2]
        j = listIndex*2 - 1 ;
        %r1 = combList(k-1, j) ;
        %c1 = combList(k-1, j+1) ;
        dst2 = (combList(:,j)-r1).^2 + (combList(:,j+1)-c1).^2 ;
        [~, ind] = min(dst2) ;
        ind = ind(1) ; % think if it's better to use ind(end)
        pairsList(k,j:j+1) = combList(ind,j:j+1) ;
        if (plotFlag)
            hold on ;
            plot([ c1 pairsList(k,j+1)], [r1 pairsList(k,j)],'yo-') ;
        end
    end
    % calculate distance between each pair of points
    pairsDist(k) = norm( pairsList(k,1:2) - pairsList(k,3:4)) ;
    if (plotFlag)
        text(c1,r1,num2str( round(pairsDist(k)*10)/10),'BackgroundColor',[.7 .9 .7]);
    end
    %pause(0.1) ;
end

% find last pair with distance <= legWidth

ind = find(pairsDist<=legWidth,1,'last'); 

r1  = pairsList(ind,1) ;
c1  = pairsList(ind,2) ;
r2  = pairsList(ind,3) ;
c2  = pairsList(ind,4) ;

% remove the line connecting (r1,c1) to (r2,c2) on bw
bw = bwcurr ;
n  = round(legWidth*5) ;
rows = round(linspace(r1,r2,n)) ;
cols = round(linspace(c1,c2,n)) ;

% transform to the 'full image' frame of reference
rows = rows + axlim(3) - (W+1) ;
cols = cols + axlim(1) - (W+1) ;

if (plotFlag)
    hold on ;
    plot(cols, rows,'gx') ;
end

% re-arrange pixels otfor removal

removePixels = unique( [rows'  cols'],'rows') ;
goodInd = find(removePixels(:,1) <= 512 & removePixels(:,1) >= 1 & removePixels(:,2) <= 512 & removePixels(:,2) >= 1) ;
removePixels = removePixels(goodInd,:) ; 
n = size(removePixels,1) ;

% remove
for k=1:n
    bw(removePixels(k,1), removePixels(k,2)) = false ;
end
        
% plot stuff
if (plotFlag)
    figure(11) ;
    imshow(bw) ; title('result') ; axis(axlim) ;
end

return ;
