%finds set intersection between two vectors A and B, where elements are
%considered the same if their euclidean difference is lower than some
%tolerance level (tol)
%
%Note: returns the elements from the A. These are obviously not the same as
%the corresponding ones from B, but are closer than tol.
function [C, indA, indB] = intersectTol(A,B,tol)

if size(A,1) < size(A,2)
    A = A' ;
end
if size(B,1) < size(B,2)
    B = B' ;
end

if (size(A,2) ~= 1) || (size(B,2) ~= 1)
    disp('Error: inputs must be in the form of vectors')
    return ;
end

D = pdist2(A,B) ; 
[indA, indB] = find(D < tol) ; 

C = A(indA) ; 

end