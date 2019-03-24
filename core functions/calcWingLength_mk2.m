function [newWingRes, meanWingLength] =  calcWingLength_mk2(wingRes) 
% calculate the average wing length and also return a new wingRes that's
% entirely a connected component
% [newWingRes, meanWingLength] =  calcWingLength_mk2(wing1Res) ;
plotFlag = false ;

Nimages = max(wingRes(:,1)) - min(wingRes(:,1)) + 1 ;
wingLength = zeros(Nimages,1) ;
newWingRes = [] ;
wingCM = zeros(Nimages+1,1) ; %PADDED BY ONE!!!

df = diff(wingRes(:,1)) ; % [t x y z ]
frameStartIndWing = [1 ; find(df==1)+1] ;
frameEndIndWing   = [frameStartIndWing(2:end)-1 ; size(wingRes,1)] ;

fraction1 = 0.03 ; % how many of wing voxels to consider for extremal ends
fraction2 = 0.05 ; % how many wing length values to ignore
delta = 8 ; % difference allowed for extremal voxels from pca value

gmreplicates = 10 ;
AICtol = -.02 ;

if plotFlag
    hfig = figure ;
end

for j = 1:Nimages
    i1 = frameStartIndWing(j) ;
    i2 = frameEndIndWing(j) ;
    
    %{
    minX = min(wingRes(i1:i2,2));
    maxX = max(wingRes(i1:i2,2));
    sizeX = maxX - minX + 1; 
    minY = min(wingRes(i1:i2,3));
    maxY = max(wingRes(i1:i2,3));
    sizeY = maxY - minY + 1; 
    minZ = min(wingRes(i1:i2,4));
    maxZ = max(wingRes(i1:i2,4));
    sizeZ = maxZ - minZ + 1; 
    
    tempSpace = false(sizeX,sizeY,sizeZ) ;
    onInd = sub2ind(size(tempSpace), wingRes(i1:i2,2) - minX + 1, wingRes(i1:i2,3) - minY + 1, wingRes(i1:i2,4) - minZ + 1) ;
    tempSpace(onInd) = true ; 
    
    CC = bwconncomp(tempSpace) ; 
    Ncc = length(CC.PixelIdxList) ;
    svec = zeros(Ncc,1) ; % vector containing the size of each connected components
    for k=1:Ncc
        svec(k) = length(CC.PixelIdxList{k}) ;
    end
    [~, idx] = sort(svec,'descend') ;
    %}
    
    [largestCC, ~, ~, ~] = findLargestHullCC (wingRes(i1:i2,2:4), 26) ;
    largestCC = double(largestCC) ;
    numPix = size(largestCC,1) ;
    
    if numPix > 800 %arbitrarily check to see if largest CC is large enough 
        %{
        tempSpace2 = false(size(tempSpace)) ;
        tempSpace2(CC.PixelIdxList{idx(1)}) = true ;
        D = -bwdist(~tempSpace2) ;
        %D(~tempSpace2) = -Inf ;
        %L = watershed(D) ;
        %if plotFlag
        %    figure ; 
        %    PATCH_3Darray(L == 0)
        %end
        %CC = bwconncomp(L == 0) ;
        %Ncc = length(CC.PixelIdxList) ;
        %svec = zeros(Ncc,1) ; % vector containing the size of each connected components
        %for k=1:Ncc
        %    svec(k) = length(CC.PixelIdxList{k}) ;
        %end
        %[~, idx] = sort(svec,'descend') ;
        %tempSpace3 = false(size(tempSpace2)) ;
        %tempSpace3(CC.PixelIdxList{idx(1)}) = true ;
        %tempSpace3 = imfill(tempSpace3,'holes') ;
        %newInd = find(tempSpace2 == true) ;
        
        [newX, newY, newZ] = ind2sub(size(tempSpace2),CC.PixelIdxList{idx(1)}) ;
        newX = newX + double(minX) + 1 ; newY = newY + double(minY) + 1 ; newZ = newZ + double(minZ) + 1 ; 
        newHull = [newX, newY, newZ] ;
        %}
        
        AIC = zeros(1,4);
        GMModels = cell(1,4);
        clust = cell(1,4) ;
        options = statset('MaxIter',500);
        for k = 1:4
            GMModels{k} = fitgmdist(largestCC,k,'Options',options);
            AIC(k)= GMModels{k}.AIC;
        end
        
        AICslope = diff(AIC)./AIC(2:end) ;
        numComponents = 1 + find(AICslope < AICtol, 1, 'first')  ;
        if isempty(numComponents)
            numComponents = 1 ;
        end
        idx = cluster(GMModels{numComponents}, largestCC) ;
            
        %want to check which cluster is closest to previous wingCM and
        %then see if that needs clustering as well
        
        figure ;
        plotType = {'bo', 'r+', 'gx', 'ksq', 'c.'} ;
        hold on
        for q = 1:numComponents
            currClust = largestCC(idx == q,:) ;
            plot3(currClust(:,1),currClust(:,2),currClust(:,3),plotType{q});
        end
        axis equal
        
    disp('blah')    
    end   
    %{   
        currWingRes = [double(wingRes(i1,1))*ones(numPix,1), currLargestCC] ;
        currWingRes = int16(currWingRes) ;
    else
        currWingRes = wingRes(i1:i2,:) ;
    end
    
    Nvox = size(currWingRes,1) ;
    pcaCoeffs = pca(double(currWingRes(:,2:4))) ;
    Rcm = round(mean(currWingRes(:,2:4))) ;
   
    mat1 = double(currWingRes(:,2:4)) - repmat(Rcm, Nvox,1)  ;
    distFromRcm = ( sum (mat1.^2,2) ) .^ 0.5 ;
    clear mat1
   
    [~, sortedInd] = sort(distFromRcm,'descend') ;
    Nvox1          = ceil(Nvox*fraction1) ;
    selectedInd    = sortedInd(1:Nvox1) ;
  
    mat1 = double(currWingRes(selectedInd,2:4)) - repmat(Rcm, Nvox1, 1) ;
    mat2 = repmat (pcaCoeffs(:,1)', Nvox1, 1) ;
   
    dotprod = sum( mat1 .* mat2 , 2) ;
   
    % find the voxels with distance around the maximum/minimum
    maxdist = max(dotprod) ;
    mindist = min(dotprod) ;
   
    grp1Flag = (dotprod>maxdist-delta) ;
    grp2Flag = (dotprod<mindist+delta) ;
   
    clear dotprod
  
    voxelGrp1 = currWingRes(selectedInd(grp1Flag),2:4) ;
    voxelGrp2 = currWingRes(selectedInd(grp2Flag),2:4) ;
   
    grp1CM = mean (voxelGrp1,1) ;
    grp2CM = mean (voxelGrp2,1) ;
    
    if plotFlag
        figure(hfig) ; clf; 
        hold on
        plot3(currWingRes(:,2),currWingRes(:,3),currWingRes(:,4),'bo')
        axis equal
        vecScale = 20 ;
        plot3([Rcm(1) Rcm(1)+vecScale *pcaCoeffs(1,1)], [Rcm(2) Rcm(2)+vecScale *pcaCoeffs(2,1)],...
            [Rcm(3) Rcm(3)+vecScale *pcaCoeffs(3,1)], 'r','LineWidth',5)
        plot3(grp1CM(1), grp1CM(2), grp1CM(3), 'gx', 'LineWidth',5)
        plot3(grp2CM(1), grp2CM(2), grp2CM(3), 'gx', 'LineWidth',5)
    end
    newWingRes = [newWingRes ; currWingRes] ;
    wingLength(j) = myNorm(grp2CM - grp1CM) ;
        %}
end    
    
% get rid of the 5% smallest and 5% largest values to estimate mean

v = sort(wingLength) ;
n1 = round(fraction2*Nimages) ;
n2 = round( (1-fraction2)*Nimages) ;

meanWingLength = mean(v(n1:n2))  ; % mean(wingLengthVecR) ;


return
