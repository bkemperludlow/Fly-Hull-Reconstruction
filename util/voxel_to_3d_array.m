function A = voxel_to_3d_array(V,buffer,target_size)
%Takes a n*3 array of 3D coordinates(V) and convert into a finemesh 3D array
%with values representing intensity.
%Buffer is a list of positive translations along each dimension to make
%sure all coordinate values are positive.
%The array is then padded with 0s at the end to reach the target size
%--------------------------------------------------------------------------
%tic
V = double(V);
buffer_grid = meshgrid(buffer,ones(size(V,1),1));
V_t = round(V)+ buffer_grid;
V = V + buffer_grid;
V_floor = floor(V) ;
V_ceil = ceil(V) ;

%round_diff = V_t - V ;
% round_diff_pos_idx = (round_diff > 0.45) ;
%round_diff_neg_idx = (round_diff < -0.45) ;
round_diff_pos_idx = ((V-V_t) > 0.45)  ;
round_diff_neg_idx = ((V_t-V) > 0.45)  ;
V_pos = round_diff_pos_idx.*V_ceil + ~round_diff_pos_idx.*V_t ;
V_neg = round_diff_neg_idx.*V_floor + ~round_diff_neg_idx.*V_t ;

A_max = max(V_ceil) ;
A=zeros(A_max(1),A_max(2),A_max(3));

try
    ind_coords = sub2ind(size(A), V_t(:,1), V_t(:,2), V_t(:,3)) ;
catch
    keyboard
end
ind_coords_pos = sub2ind(size(A), V_pos(:,1), V_pos(:,2), V_pos(:,3)) ;
ind_coords_neg = sub2ind(size(A), V_neg(:,1), V_neg(:,2), V_neg(:,3)) ;

A(ind_coords) = 1 ;
A(ind_coords_pos) = 1 ;
A(ind_coords_neg) = 1 ;
%A = imfill(A,6,'holes') ;
A = padarray(A,target_size-size(A),0,'post');
A=logical(A);
%toc
%--------------------------------------------------------------------------

end

%--------------------------------------------------------------------------
% tic
% V = double(V);
% buffer_grid = meshgrid(buffer,ones(size(V,1),1));
% V_t = round(V)+ buffer_grid;
% V = V + buffer_grid;
% A_max = max(V_t) ;
% A=zeros(A_max(1),A_max(2),A_max(3));
% for i=1:size(V,1)
%     coords=uint8(V_t(i,:)); %rounds nicely
%     A(coords(1),coords(2),coords(3))=1;
%     %If voxel near middle, fill both surrounding cell. Prevents fake
%     %missing cells in the middle of image due to rounding
%     low = uint8(floor(V(i,:)));
%     if V(i,1)-double(low(1))>0.45 && V(i,1)-double(low(1))<0.55
%         A(low(1):low(1)+1,coords(2),coords(3)) = 1;
%     end
%     if V(i,2)-double(low(2))>0.45 && V(i,2)-double(low(2))<0.55
%         A(coords(1),low(2):low(2)+1,coords(3)) = 1;
%     end
%     if V(i,3)-double(low(3))>0.45 && V(i,3)-double(low(3))<0.55
%         A(coords(1),coords(2),low(3):low(3)+1) = 1;
%     end
% end
% A = padarray(A,target_size-size(A),0,'post');
% A=logical(A);
% toc
%--------------------------------------------------------------------------

