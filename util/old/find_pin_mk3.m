function [ tip ] = find_pin_mk3( data )


%--------------------------------------------------------------------------
%Calculates the vector corresponding to the tip of the pin glued to the fly
%
%Developed by Luca and Sam
%
%Input
%-----
%data - the output of hullAnalysis_mk3. Notably contains the coordinates of
%the body voxels of the fly, as well as the body axis
%
%Output
%------
%pin - array containing the location of the pin tip
%
%--------------------------------------------------------------------------

%% First initialize array and load relevant aspects of 'data'
Nframes = data.Nimages ;
bodyCM=data.bodyCM;
AHat=data.AHat;
tip = zeros(Nframes,3);          %the tip of the pin

try 
    body_axis=data.body_axis;
catch
    body_axis = 50*data.AHat;
end

%Find the first and last frame of the movie
df = diff(data.res(:,1)) ;
frameStartInd = [1 ; find(df==1)+1] ;
frameEndInd   = [frameStartInd(2:end)-1 ; size(data.res,1)] ;
clear df ;


%% Establish reference frame information

FrameRef=input('Enter the reference frame number frame \n');

row1       = frameStartInd(FrameRef) ;
row2       = frameEndInd(FrameRef) ;
coords     = data.res(row1:row2,2:4) ; % xyz
IDX        = data.RESIDX(row1:row2,:) ;
bodyRows   = (IDX(:,data.bodyInd)==1) ;
bodyVoxels = coords(bodyRows,:) ;
Nvox       = sum(bodyRows) ;

figure;
plot3(bodyVoxels(:,1),bodyVoxels(:,2),bodyVoxels(:,3),'bo');
axis equal;
hold on;
plot3([bodyCM(FrameRef,1)-body_axis(FrameRef,1) bodyCM(FrameRef,1)+body_axis(FrameRef,1)],[bodyCM(FrameRef,2)-body_axis(FrameRef,2) body_axis(FrameRef,2)+bodyCM(FrameRef,2)],[bodyCM(FrameRef,3)-body_axis(FrameRef,3) bodyCM(FrameRef,3)+body_axis(FrameRef,3)],'g','LineWidth',7);
tip(FrameRef,:) = input('Enter coordinates of pin tip in the reference frame \n');
plot3(tip(FrameRef,1),tip(FrameRef,2),tip(FrameRef,3),'r+','MarkerSize',15);

pinDist = norm(tip(FrameRef,:)-bodyCM(FrameRef,:));
tol = 10;
%% Track pin using pin location from the reference frame
for i=FrameRef+1:Nframes
    
    row1       = frameStartInd(i) ;
    row2       = frameEndInd(i) ;
    coords     = data.res(row1:row2,2:4) ; % xyz
    IDX        = data.RESIDX(row1:row2,:) ;
    bodyRows   = (IDX(:,data.bodyInd)==1) ;
    bodyVoxels = coords(bodyRows,:) ;
    Nvox       = sum(bodyRows) ;
    
    old_tip = tip(i-1,:);
    
    diffFromOldTip = double(bodyVoxels)-repmat(old_tip,size(bodyVoxels,1),1);
    diffSum = sqrt(sum(diffFromOldTip .* diffFromOldTip,2));
    %     closePointsIdx = find(diffSum<5);
    %     closePoints = bodyVoxels(closePointsIdx,:);
    %
    %     diffFromCM = double(closePoints) - repmat(bodyCM(i,:),size(closePoints,1),1);
    %     diffFromCMNorm = norm(diffFromCM,1);
    %     maxDistFlag = find(diffFromCMNorm == max(diffFromCMNorm));
    Idx=find(diffSum == min(diffSum));
    %Idx=find(diffSum <= pinDist*(1+tol) | diffSum >= pinDist*(1-tol) );
    
    %Idx=find(diffSum < 5);
    
    %diffFromCM = double(bodyVoxels)-repmat(bodyCM(i,:),size(bodyVoxels,1),1);
    %diffSum = sqrt(sum(diffFromCM .* diffFromCM,2));
    %Idx2 =find(diffSum <= pinDist+tol & diffSum >= pinDist-tol );
    
    %Idx3 = intersect(Idx,Idx2);
    
    tip(i,:) = mean(bodyVoxels(Idx,:),1);
    %     tip(i,:) = closePoints(maxDistFlag,:);
    
    
    plot3(bodyVoxels(:,1),bodyVoxels(:,2),bodyVoxels(:,3),'bo');
    axis equal;
    hold on;
    plot3([bodyCM(i,1)-body_axis(i,1) bodyCM(i,1)+body_axis(i,1)],[bodyCM(i,2)-body_axis(i,2) body_axis(i,2)+bodyCM(i,2)],[bodyCM(i,3)-body_axis(i,3) bodyCM(i,3)+body_axis(i,3)],'g','LineWidth',7);
    plot3(tip(i,1),tip(i,2),tip(i,3),'r+','MarkerSize',15);
    
    clf
end


for i=1:FrameRef-1
    
    j=FrameRef-i;
    
    row1       = frameStartInd(j) ;
    row2       = frameEndInd(j) ;
    coords     = data.res(row1:row2,2:4) ; % xyz
    IDX        = data.RESIDX(row1:row2,:) ;
    bodyRows   = (IDX(:,data.bodyInd)==1) ;
    bodyVoxels = coords(bodyRows,:) ;
    Nvox       = sum(bodyRows) ;
    
    old_tip = tip(j+1,:);
    
    diffFromOldTip = double(bodyVoxels)-repmat(old_tip,size(bodyVoxels,1),1);
    diffSum = sqrt(sum(diffFromOldTip .* diffFromOldTip,2));
    %     closePointsIdx = find(diffSum<5);
    %     closePoints = bodyVoxels(closePointsIdx,:);
    %
    %     diffFromCM = double(closePoints) - repmat(bodyCM(i,:),size(closePoints,1),1);
    %     diffFromCMNorm = norm(diffFromCM,1);
    %     maxDistFlag = find(diffFromCMNorm == max(diffFromCMNorm));
    Idx=find(diffSum == min(diffSum));
    %Idx=find(diffSum < 5);
    
    %diffFromCM = double(bodyVoxels)-repmat(bodyCM(j,:),size(bodyVoxels,1),1);
    %diffSum = sqrt(sum(diffFromCM .* diffFromCM,2));
    %Idx2 =find(diffSum <= pinDist+tol & diffSum >= pinDist-tol );
    
    %Idx3 = intersect(Idx,Idx2);
    
    tip(j,:) = mean(bodyVoxels(Idx,:),1);
    %     tip(i,:) = closePoints(maxDistFlag,:);
    
    clf
    plot3(bodyVoxels(:,1),bodyVoxels(:,2),bodyVoxels(:,3),'bo');
    axis equal;
    hold on;
    plot3([bodyCM(j,1)-body_axis(j,1) bodyCM(j,1)+body_axis(j,1)],[bodyCM(j,2)-body_axis(j,2) body_axis(j,2)+bodyCM(j,2)],[bodyCM(j,3)-body_axis(j,3) bodyCM(j,3)+body_axis(j,3)],'g','LineWidth',7);
    plot3(tip(j,1),tip(j,2),tip(j,3),'r+','MarkerSize',15);
    
    
end

figure, plot3(tip(:,1), tip(:,2),tip(:,3),'b.-') 
hold on;
plot3(bodyCM(:,1), bodyCM(:,2),bodyCM(:,3),'g.-') 
axis equal;

end

