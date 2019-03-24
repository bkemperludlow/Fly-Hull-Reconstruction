function tip = findWingTip(coords, spanVec, wingCM)

% find distances of coords from wingCM
Nvox = size(coords,1) ;

fraction = 0.05 ;

% 1st screening
% -------------
% consider only voxels in a cone of 30 degrees around the line starting
% from wingCM along the direction of spanVec
mat1 = coords - repmat(wingCM, Nvox, 1) ;
mat1 = mat1 ./ repmat( myNorm(mat1), 1,3) ;
mat2 = repmat(spanVec, Nvox, 1) ;
dotprod = dot(mat1, mat2, 2) ;

ind1 = dotprod > 0.5 ; % cos(30 deg) = 0.5

% 2nd screening
% -------------
% take the farthest "fraction" of the voxels and calculate their mean position

coords = coords(ind1,:) ;
Nvox   = size(coords,1) ;

dst  = myNorm ( coords - repmat(wingCM,Nvox,1) ) ;

[~, sortedInd] = sort(dst,'descend') ;
Nvox1 = ceil(Nvox*fraction) ;
selectedInd    = sortedInd(1:Nvox1) ;
if (isempty(selectedInd))
    tip = [NaN NaN NaN] ;
else
    if (numel(selectedInd)==1)
        tip = coords(selectedInd,:) ;
    else
        tip = mean(coords(selectedInd,:)) ;
    end
end

end