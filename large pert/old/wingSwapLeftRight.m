function data = wingSwapLeftRight(data,swapFrames)
% define frame boundaries in res structure
df = diff(data.res(:,1)) ;
frameStartInd = [1 ; find(df==1)+1] ;
frameEndInd   = [frameStartInd(2:end)-1 ; size(data.res,1)] ;
clear df ;

% swap left and right
for j = swapFrames
    row1 = frameStartInd(j) ;
    row2 = frameEndInd(j) ;
    tempResidx = data.RESIDX(row1:row2, 2) ;
    data.RESIDX(row1:row2, 2) = data.RESIDX(row1:row2, 3) ;
    data.RESIDX(row1:row2, 3) = tempResidx ;
    
    tempChord = data.rightChordHats(j,:) ;
    tempSpan = data.rightSpanHats(j,:) ;
    data.rightChordHats(j,:)= data.leftChordHats(j,:) ;
    data.rightSpanHats(j,:)= data.leftSpanHats(j,:) ;
    data.leftChordHats(j,:) = tempChord ;
    data.leftSpanHats(j,:) = tempSpan ;
    
    tempWingCM = data.rightWingCM(j,:) ;
    data.rightWingCM(j,:) = data.leftWingCM(j,:) ;
    data.leftWingCM(j,:) = tempWingCM ;
    
    data.rightChordTopProjections(j) = data.leftChordTopProjections(j) ;
    data.leftChordTopProjections(j)  = data.rightChordTopProjections(j) ;
    data.diag11Right(j) = diag21Left(j) ;
    data.diag12Right(j) = diag22Left(j) ;
    data.diag21Left(j)  = diag11Right(j) ;
    data.diag22Left(j)  = diag12Right(j) ;
    data.rightWingTips(j,:) = data.leftWingTips(j,:) ;
    data.leftWingTips(j,:)  = data.rightWingTips(j,:) ;
end
end