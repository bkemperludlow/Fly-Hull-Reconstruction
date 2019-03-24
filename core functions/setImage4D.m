function m = setImage4D(m, i1, i2, bw)

% data checkups:
% i1 should be between 1 and m.dim(1)
% i2 should be between 1 and m.dim(2)
% bw shoule be a logical matrix of size m.dim(3:4)


ind1vec = m.dim(3)*(i1-1) +  (1:m.dim(3));
ind2vec = m.dim(4)*(i2-1) +  (1:m.dim(4)) ;

m.mat(ind1vec, ind2vec) = bw ;

return

