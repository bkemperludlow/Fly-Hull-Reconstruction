function m = init4D(dim, nzmax)
% returns m which is a structure with two fields
% m.dim
% m.mat - a 2d sparse matrix that emulates a 4D matrix
%
% the input parameter nzmax is optional. it indicates how much memory will
% be allocated to the sparse matrix

% check size of dim
if (size(dim,1)~=1 || size(dim,2)~=4)
    error('init4D : input parameter ''dim'' should have size of [1 4]') ;
end

m.dim = dim ;

s1 = dim(1) * dim(3) ;
s2 = dim(2) * dim(4) ;

m.dim = dim ;
if (exist('nzmax','var'))
    m.mat = sparse([],[],false,s1, s2, nzmax) ;
else
    m.mat = sparse([],[],false,s1, s2) ; 
end


return
