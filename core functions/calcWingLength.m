function [meanWingLengthR meanWingLengthL wingLengthVecR wingLengthVecL ] =  calcWingLength (data) 
% calculate the length of each wing as the distance from the tip to the
% hinge. return the mean lengths as well as the length calculated for each
% frame

fraction = 0.05 ;

d = data.rightWingTips - data.rightHinges ;
wingLengthVecR = myNorm(d) ;
% get rid of the 5% smallest and 5% largest values to estimate mean

v = sort(wingLengthVecR) ;
N = length(v) ;
n1 = round(fraction*N) ;
n2 = round( (1-fraction)*N) ;

meanWingLengthR = mean(v(n1:n2))  ; % mean(wingLengthVecR) ;

d = data.leftWingTips - data.leftHinges ;
wingLengthVecL = myNorm(d) ;
v = sort(wingLengthVecL) ;
meanWingLengthL = mean(v(n1:n2))  ; % mean(wingLengthVecL) ;

clear d v

return
