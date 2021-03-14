%--------------------------------------------------------------------------
% function to switch the hull reconstruction/analysis variables that are
% specific to the right vs left wing. operates on the data structure that
% is output by the analysis
%
% NB: This does not handle angle variables (e.g. data.anglesBodyFrame), so
% angles should be re-calculated after performing swaps.
%
% Q: do the roll vectors need to be switched?
%--------------------------------------------------------------------------
function data_out = swapWingLeftRight(data_in,swapFrames)
% make data_out, a copy of data_in
data_out = data_in ; 

% define frame boundaries in res structure
df = diff(data_in.res(:,1)) ;
frameStartInd = [1 ; find(df==1)+1] ;
frameEndInd   = [frameStartInd(2:end)-1 ; size(data_in.res,1)] ;
clear df ;

% get indices for the two wings 
rightWingInd = data_in.rightWingInd ; 
leftWingInd = data_in.leftWingInd ; 

for i = 1:length(swapFrames)
    swapFrame = swapFrames(i) ;
    row1 = frameStartInd(swapFrame) ;
    row2 = frameEndInd(swapFrame) ;
    
    % swap left and right voxel indices
    data_out.RESIDX(row1:row2, rightWingInd) = ...
        data_in.RESIDX(row1:row2, leftWingInd) ;
    data_out.RESIDX(row1:row2, leftWingInd) = ...
        data_in.RESIDX(row1:row2, rightWingInd) ;
    
    % swap span and chord vectors
    data_out.rightChordHats(swapFrame,:)= data_in.leftChordHats(swapFrame,:) ;
    data_out.rightSpanHats(swapFrame,:)= data_in.leftSpanHats(swapFrame,:) ;
    data_out.leftChordHats(swapFrame,:) = data_in.rightChordHats(swapFrame,:) ;
    data_out.leftSpanHats(swapFrame,:) = data_in.rightSpanHats(swapFrame,:) ;
    
    % swap wing centers of mass and wing tips
    data_out.rightWingCM(swapFrame,:) = data_in.leftWingCM(swapFrame,:) ;
    data_out.leftWingCM(swapFrame,:) = data_in.rightWingCM(swapFrame,:) ;
    data_out.rightWingTips(swapFrame,:) = data_in.leftWingTips(swapFrame,:) ;
    data_out.leftWingTips(swapFrame,:)  = data_in.rightWingTips(swapFrame,:) ;
    
    data_out.rightChordTopProjections(swapFrame) = ...
        data_in.leftChordTopProjections(swapFrame) ;
    data_out.leftChordTopProjections(swapFrame)  = ...
        data_in.rightChordTopProjections(swapFrame) ;
    data_out.diag11Right(swapFrame) = data_in.diag21Left(swapFrame) ;
    data_out.diag12Right(swapFrame) = data_in.diag22Left(swapFrame) ;
    data_out.diag21Left(swapFrame)  = data_in.diag11Right(swapFrame) ;
    data_out.diag22Left(swapFrame)  = data_in.diag12Right(swapFrame) ;
end

end