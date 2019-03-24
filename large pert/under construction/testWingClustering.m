dataPath = ['I:\Fly Data\VNC Sensory Lines\08_09042018\Analysis\' ...
    'Unsorted\Expr_8_mov_008\' ] ;

dataFilename = 'Expr_8_mov_008_test.mat' ;
dataFilename_LPC = 'Expr_8_mov_008_LPC_03.mat' ;
flagStructName = 'flagStruct_03.mat' ;

data = importdata(fullfile(dataPath, dataFilename)) ;
data_LPC = importdata(fullfile(dataPath, dataFilename_LPC)) ;
flagStruct = importdata(fullfile(dataPath, flagStructName)) ;

k_max = 3 ;
clust_criterion_cell = {'CalinskiHarabasz', 'DaviesBouldin', 'gap',...
    'silhouette'} ;
error_wing_side_cell = {'right','left','left','left'} ;
%--------------------------------------------------------------------------
clusterRightFlag = flagStruct.clusterRightFlag ;
clusterLeftFlag = flagStruct.clusterLeftFlag ;
clusterErrorFlag = flagStruct.clusterErrorFlag ;

clusterRightInd = find(clusterRightFlag) ;
clusterLeftInd = find(clusterLeftFlag) ;
clusterErrorInd = find(clusterErrorFlag) ;

clusterInd_cell = {clusterRightInd, clusterLeftInd, clusterErrorInd} ;

%clust_eval_array = nan(length(clusterInd_cell), length(clust_criterion_cell)) ;
cluster_eval_array = nan(data.Nimages, 3) ;
clust_optimal_k = nan(data.Nimages,1) ; 
%--------------------------------------------------------------------------
df = diff(data.res(:,1)) ;
frameStartInd = [1 ; find(df==1)+1] ;
frameEndInd   = [frameStartInd(2:end)-1 ; size(data.res,1)] ;
clear df ;

wing_str = 'right' ;
for i = 1:length(frameStartInd)
    tic
    row1 = frameStartInd(i) ;
    row2 = frameEndInd(i) ;
    
    % coordinates that correspond to current frame
    coords = data.res(row1:row2,2:4) ; % xyz positions of each voxel in frame=1st col value
    IDX = data.RESIDX(row1:row2,:) ; % Classifying each voxel as L/R wing or body
    wingRows = (IDX(:,data.([ wing_str 'WingInd']))==1) ;
    %wingRowsAlt = (IDX(:,data.([ wing_str_alt 'WingInd']))==1) ;
    wingVox = double([coords(wingRows,1), coords(wingRows,2), ...
        coords(wingRows,3)]);
    %wingVoxAlt = double([coords(wingRowsAlt,1), coords(wingRowsAlt,2), ...
    %    coords(wingRowsAlt,3)]);
    
    try
        eva = evalclusters(wingVox,'kmeans',clust_criterion_cell{1},...
            'KList',1:k_max) ;
    catch
        continue
    end
    cluster_eval_array(i,:) = eva.CriterionValues ;
    clust_optimal_k(i) = eva.OptimalK ; 
    toc
    disp(i)
end
%--------------------------------------------------------------------------

% for i = 1:length(clusterInd_cell)
%    clusterIndCurr = clusterInd_cell{i} ;
%    clust_eval_curr = nan(length(clusterIndCurr), k_max) ;
%    for j = clusterIndCurr
%
%    end
%
% end