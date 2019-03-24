function ellipsoidRes = myEllipsoid(xc, yc, zc, xr, yr, zr)
    %use integer inputs so this can output voxels
    size = ceil((4/3)*pi*xr*yr*zr) ; %pad the volume a little bit
    ellipsoidRes = nan(size,3) ;
    counter = 1 ;
    
    for i = (xc-xr):(xc+xr)
        for j = (yc-yr):(yc+yr)
            for k = (zc-zr):(zc+zr)
                if ((i-xc)/xr)^2 + ((j-yc)/yr)^2 + ((k-zc)/zr)^2 <= 1
                    ellipsoidRes(counter, :) = [i, j, k] ; 
                    counter = counter + 1 ;
                end    
            end
        end
    end

    ind = find(~isnan(ellipsoidRes(:,1))) ;
    ellipsoidRes = ellipsoidRes(ind, :) ;
    
end