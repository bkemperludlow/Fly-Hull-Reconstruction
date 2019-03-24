function m = convert4D2sparse(choplogmatrix)

% find what dim is
dim = size(choplogmatrix) ;
nzmax = sum(choplogmatrix(:)) ;

m = init4D(dim, nzmax) ;

for it=1:dim(1)
    for cam=1:dim(2)
        bw = squeeze(choplogmatrix(it, cam, :, :));
        
        ind1vec = m.dim(3)*(it-1) +  (1:m.dim(3));
        ind2vec = m.dim(4)*(cam-1) +  (1:m.dim(4)) ;

        m.mat(ind1vec, ind2vec) = bw ;
    end
end

return

