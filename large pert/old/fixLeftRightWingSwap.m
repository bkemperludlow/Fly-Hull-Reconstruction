% try to fix issues of left and right wing flipping when roll is extreme.
% Should be run on data structure with angles already estimated

dataPath = 'H:\Reconstruction Problems\10_20072017\Analysis\Expr_10_mov_039' ; 
filename = 'Expr_10_mov_039_test.mat' ; 

defineConstantsScript 
data = importdata(fullfile(dataPath, filename)) ; 

vel_thresh = 15 ; %should figure out a way to pull from data

leftWingCM = data.leftWingCM ; 
rightWingCM = data.rightWingCM ; 
rho = data.anglesLabFrame(:,RHO) ; 
frames = 1:size(rightWingCM,1) ; 

nanIndL_all = false(size(frames))' ; 
nanIndR_all = false(size(frames))' ; 

for dim =1:3 
    nanIndL = isnan(leftWingCM(:,dim)) ;
    nanIndR = isnan(rightWingCM(:,dim)) ;
    leftWingCM(nanIndL,dim) = interp1(frames(~nanIndL),leftWingCM(~nanIndL),frames(nanIndL)) ; 
    rightWingCM(nanIndR,dim) = interp1(frames(~nanIndR),rightWingCM(~nanIndR),frames(nanIndR)) ; 
    
    nanIndL_all = (nanIndL_all | nanIndL) ;
    nanIndR_all = (nanIndR_all | nanIndR) ;
end

leftWingCM_velNorm = [0 ; myNorm(diff(leftWingCM,1)) ] ; 
rightWingCM_velNorm = [0 ; myNorm(diff(rightWingCM,1)) ] ;

% figure ;
% hold on 
% plot(frames, leftWingCM_velNorm ,'b-') ; 
% plot(frames, rightWingCM_velNorm,'r-') ; 

highVelInd = (leftWingCM_velNorm > vel_thresh) & (rightWingCM_velNorm > vel_thresh) ;
highVel_or_nan = highVelInd | nanIndR_all | nanIndL_all ; 
%goodRoll = (abs(rho) < 45) ;
badInd = find(highVel_or_nan) ; 

ignoreFrames = [] ; 
potentialSwapFrames = [] ;
for i = badInd'
    rcm = rightWingCM(i,:) ; 
    lcm = leftWingCM(i,:) ; 
    rcm_prev = rightWingCM(i-1,:) ; 
    lcm_prev = leftWingCM(i-1,:) ; 
    
    if (norm(rcm-lcm_prev) < vel_thresh) && (norm(lcm-rcm_prev) < vel_thresh) 
        potentialSwapFrames = [potentialSwapFrames , i] ;
    else
        ignoreFrames = [ignoreFrames , i ] ; 
    end

end

swapFrameStart = potentialSwapFrames(1:2:end) ; 
swapFrameEnd = potentialSwapFrames(2:2:end) ; 
% for j = 1:2:length(potentialSwapFrames)
%     swapFrames = [swapFrames, potentialSwapFrames(j):(potentialSwapFrames(j+1)-1)] ; 
% end

data.ignoreFrames = [data.ignoreFrames , ignoreFrames ] ; 
data = wingSwapLeftRight(data,swapFrames) ;

%rollBar = leftWingCM - rightWingCM ; 
%rollBar = rollBar./myNorm(rollBar) ; 
%LR_check = rollBar


