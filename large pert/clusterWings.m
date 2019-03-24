%--------------------------------------------------------------------------
% Function to deal with cases when the voxels that should be assigned to
% multiple wings are only assigned to one. This happens often in extreme
% perturbations
%
% NB: a version of this is performed in hullAnalysis, so I'll try to build
% off of that
%
% Things to try:
%   - image opening
%   -cluster gap metric
%--------------------------------------------------------------------------
function [wingVox, wingRowsOut, idx, centroids, badClusterFlag] = ...
    clusterWings(data, frameNum, wing_str, debugFlag)
%--------------------------------------------------------------------------
%% params and inputs
if ~exist('debugFlag','var')
    debugFlag = false ;
end
numReplicates = 20 ;
downSampleStep = 2 ; % 1 would mean no down sampling, 2 reduces by 50%, ...
k_list = [1:3] ; % how many clusters to test
%numTestClusts = 4 ;
%clustCriterion = 'DaviesBouldin' ; % 'CalinskiHarabasz', 'DaviesBouldin', 'gap', 'silhouette'

wingLength = data.wingLength ;
LL = 0.55 ;
lengthThresh = wingLength * LL ;

if strcmp(wing_str,'right')
    wing_str_alt = 'left' ;
    plt_color = 'r' ;
    plt_color_alt = 'b' ;
elseif strcmp(wing_str,'left')
    wing_str_alt = 'right' ;
    plt_color = 'b' ;
    plt_color_alt = 'r' ;
else
    disp('Invalid wing side selection')
    keyboard
end
%--------------------------------------------------------------------------
%% find appropriate voxels for frame
% get the indices for frame voxels
df = diff(data.res(:,1)) ;
frameStartInd = [1 ; find(df==1)+1] ;
frameEndInd   = [frameStartInd(2:end)-1 ; size(data.res,1)] ;
clear df ;

row1 = frameStartInd(frameNum) ;
row2 = frameEndInd(frameNum) ;

% coordinates that correspond to current frame
coords = data.res(row1:row2,2:4) ; % xyz positions of each voxel in frame=1st col value
IDX = data.RESIDX(row1:row2,:) ; % Classifying each voxel as L/R wing or body
wingRows = (IDX(:,data.([ wing_str 'WingInd']))==1) ;
wingRowsAlt = (IDX(:,data.([ wing_str_alt 'WingInd']))==1) ;
wingVox = double([coords(wingRows,1), coords(wingRows,2), ...
    coords(wingRows,3)]);
wingVoxAlt = double([coords(wingRowsAlt,1), coords(wingRowsAlt,2), ...
    coords(wingRowsAlt,3)]);
bodyCM = data.bodyCM(frameNum, :) ;

%--------------------------------------------------------------------------
%% show initial voxel configuration
if debugFlag
    h_debug = figure ;
    hold on
    plot3(wingVox(1:downSampleStep:end,1), ...
        wingVox(1:downSampleStep:end,2), ...
        wingVox(1:downSampleStep:end,3), '.', ...
        'Color',plt_color)
    plot3(wingVoxAlt(1:downSampleStep:end,1), ...
        wingVoxAlt(1:downSampleStep:end,2), ...
        wingVoxAlt(1:downSampleStep:end,3), 'o', ...
        'Color',plt_color_alt)
    axis equal
    grid on
    box on
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    
end
%--------------------------------------------------------------------------
%% first test to see if clustering is appropriate
%{
%tic
% use a faster gmm cluster evaluation to check if we're likely to need 3
% clusters
eva_gmm = evalclusters(wingVox,'gmdistribution','DaviesBouldin',...
    'KList',k_list_init) ; 
if eva_gmm.OptimalK > 2
    k_list = 1:3 ; 
else
    k_list = 1:2 ; 
end
% based on gmm evaluation, try different kmeans clustering
eva = evalclusters(wingVox(1:downSampleStep:end,:),'kmeans','gap',...
    'KList',k_list) ;
%disp(k_list(end)) 
%toc
%}
eva = evalclusters(wingVox(1:downSampleStep:end,:),'kmeans','CalinskiHarabasz',...
    'KList',k_list) ;
if eva.OptimalK < 2
    badClusterFlag = true ;
    idx = [] ;
    centroids = [] ;
    wingRowsOut = [] ; 
    return
else
    badClusterFlag = false ;
end
%--------------------------------------------------------------------------
%% perform clustering
% NB: now we're potentially clustering with more than two components, so
% need to find the most "wing like" if there are >2 clusters
numClusters = eva.OptimalK ;
fprintf('Optimal cluster number: %d \n',numClusters) 
idx = kmeans(wingVox,numClusters,'Replicates',numReplicates) ;

if debugFlag
    plot_mrkr_cell = {'gx', 'cx', 'yx', 'kx' } ;
    set(0, 'CurrentFigure',h_debug)
    for q = 1:numClusters
        plot3(wingVox(idx==q,1),wingVox(idx==q,2),wingVox(idx==q,3),...
            plot_mrkr_cell{q})
    end
end

%--------------------------------------------------------------------------
%% check to see if this gives rise to reasonable wings (and calc centroids)
wingClust_cell = cell(numClusters,1) ; 
wingLargestCC_cell = cell(numClusters,1) ; 
centroids_all = nan(numClusters,3) ; 
farPoints_all = nan(numClusters,3) ; 
farPointDist_all = nan(numClusters,1) ; 

for j = 1:numClusters
    wingClust = wingVox(idx==j,:) ;
    [wingLargestCC, farPoint, ~, farPointDist] = ...
        refineWingVoxels(wingClust, bodyCM , lengthThresh) ;
    newCentroid = ...
        calcCentroidFromTopView(wingLargestCC, farPoint, wingLength/3) ;
    
    wingClust_cell{j} = wingClust ;
    wingLargestCC_cell{j} = wingLargestCC ;
    centroids_all(j,:) = newCentroid ;
    farPoints_all(j,:) = farPoint ;
    farPointDist_all(j) = farPointDist ; 
end

% re-shuffle the cluster indices to make sure we only output 2 clusters
[~, sort_ind] = sort(farPointDist_all,'descend') ; 
good_cluster_inds = sort_ind(1:2) ; 
%idx((idx ~= good_cluster_inds(1)) & (idx ~= good_cluster_inds(2))) = 0 ; 
idx_temp = zeros(size(idx)) ; 
idx_temp(idx == good_cluster_inds(1)) = 1 ; 
idx_temp(idx == good_cluster_inds(2)) = 2 ; 
idx = idx_temp ; 

centroids = centroids_all(good_cluster_inds, :) ; 

% need to deal with voxels and the issues of nested logical indexing :(
wingRowsOut = false(size(wingRows,1),2) ; 
wingRowsOut((wingRows ==1),1) = (idx == 1) ; 
wingRowsOut((wingRows ==1),2) = (idx == 2) ; 

end

% %--------------------------------------------------------------------------
% %% try image opening? (would need to move this up obvi)
% SE = strel('sphere',1); %Structuring element used for image morphologies
% target_size = [120,120,120]; %target fixed im size
% buffer = [60,60,60]; %padding when converting from voxels to im to make all voxel coord>0 %[50 50 50] ;
% 
% wingVox_meanSub = wingVox-repmat(mean(wingVox),size(wingVox,1),1) ; 
% wing_im = voxel_to_3d_array(wingVox_meanSub,buffer,target_size);
% [opened_im,~] = open_image3D(wing_im,SE) ;
% %eroded_wing_vox =  im2voxel(eroded_im);
% opened_wing_vox =  im2voxel(opened_im);
% %eroded_wing_vox =  eroded_wing_vox - meshgrid(buffer,ones(size(eroded_wing_vox,1),1));
% opened_wing_vox =  opened_wing_vox - meshgrid(buffer,ones(size(opened_wing_vox,1),1));
% %COM = find_COM_vox(opened_wing_vox);
% opened_wing_vox =  opened_wing_vox -repmat(COM,size(opened_wing_vox,1),1);
