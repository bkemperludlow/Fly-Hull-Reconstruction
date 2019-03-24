function fakeFly = setFakeFlyDOF(fakeFly, CM, bodyYPR, rightYPR, leftYPR, plotFlag) 

%This should take a set of fake fly voxels generated by makeFakeFlyData and
%set the position and orientation of those voxels

if nargin < 6
    plotFlag = false ;
end

deg2rad = (pi/180) ; 

%CM = [10 10 10] ; %voxels. temporary, for testing
bodyYPR = bodyYPR * deg2rad ; %radians
rightYPR = rightYPR * deg2rad ; % [stroke, wing pitch, elevation] 
leftYPR  = leftYPR * deg2rad ;

% Rotate right wing
N_R = size(fakeFly.wingRVoxels,1) ; 
hingeCenteredWingR = fakeFly.wingRVoxels - repmat(fakeFly.hingeR, N_R,1) ;
R1 = eulerRotationMatrix(-rightYPR(1), -rightYPR(2), rightYPR(3)) ;
hingeCenteredWingR = R1*hingeCenteredWingR' ;
fakeFly.wingRVoxels = hingeCenteredWingR' + repmat(fakeFly.hingeR, N_R,1) ;

%Rotate left wing

N_L = size(fakeFly.wingLVoxels,1) ; 
hingeCenteredWingL = fakeFly.wingLVoxels - repmat(fakeFly.hingeL, N_L,1) ;
R2 = eulerRotationMatrix(leftYPR(1), -leftYPR(2), -leftYPR(3)) ;
hingeCenteredWingL = R2*hingeCenteredWingL' ;
fakeFly.wingLVoxels = hingeCenteredWingL' + repmat(fakeFly.hingeL, N_L,1) ;


%Rotate whole body

R3 = eulerRotationMatrix(bodyYPR(1), bodyYPR(2), bodyYPR(3)) ;
fakeFly.bodyVoxels = (R3'*fakeFly.bodyVoxels')' ;
fakeFly.wingRVoxels = (R3'*fakeFly.wingRVoxels')' ;
fakeFly.wingLVoxels = (R3'*fakeFly.wingLVoxels')' ;
fakeFly.hingeR = (R3'*fakeFly.hingeR')' ;
fakeFly.hingeL = (R3'*fakeFly.hingeL')' ;

%Translate fly

fakeFly.bodyVoxels = fakeFly.bodyVoxels + repmat(CM, size(fakeFly.bodyVoxels,1),1) ; 
fakeFly.wingRVoxels = fakeFly.wingRVoxels + repmat(CM, size(fakeFly.wingRVoxels,1),1) ; 
fakeFly.wingLVoxels = fakeFly.wingLVoxels + repmat(CM, size(fakeFly.wingLVoxels,1),1) ; 
fakeFly.hingeR = fakeFly.hingeR + CM ;
fakeFly.hingeL = fakeFly.hingeL + CM ; 

if plotFlag
    h = plotFakeFly(fakeFly) ;
end


end










