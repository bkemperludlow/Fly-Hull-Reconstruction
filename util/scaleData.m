function [scaledX, meanX, stdX] = scaleData(X) 
%-------------------------------------------------------------------------

%Rescales data so that x' = x - mean(x) / std(x) . This should make the
%data look normally distributed and standardize the units. 
%
%Assumes rows correspond to oberservations and columns correspond to
%different features (variables, measurement types, etc.)
%
%same as scikit.learn.preprocessing.scale
%-------------------------------------------------------------------------
[M, ~] = size(X) ;

meanX = mean(X, 1) ;
stdX = std(X,1) ;

scaledX = (X - repmat(meanX, M, 1))./ repmat(stdX, M, 1) ;

end