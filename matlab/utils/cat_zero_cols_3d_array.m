function A_out = cat_zero_cols_3d_array(A_in, num_added_cols)
%CAT_ZERO_COLS_3D_ARRAY A_out is a deep copy of A_in input 2D matrix, with an extra 
% 0s row at the bottom
% One row is added by default i.e., if only one argument is passed
if (nargin < 2)
    A_out = zeros(size(A_in, 1), size(A_in, 2)+1, size(A_in, 3));
    for ii = 1:size(A_in, 3)
        A_out(:,:,ii) = cat_zero_col(A_in(:,:,ii));
    end
else
    A_out = zeros(size(A_in, 1), size(A_in, 2)+num_added_cols, size(A_in, 3));
    for ii = 1:size(A_in, 3)
        A_out(:,:,ii) = cat_zero_col(A_in(:,:,ii), num_added_cols);
    end
end

end %function

