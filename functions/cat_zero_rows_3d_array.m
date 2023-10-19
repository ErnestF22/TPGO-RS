function A_out = cat_zero_rows_3d_array(A_in, num_added_rows)
%CAT_ZERO_ROWS_3D_ARRAY A_out is a deep copy of A_in input 2D matrix, with an extra 
% 0s row at the bottom
% One row is added by default i.e., if only one argument is passed
if (nargin < 2)
    A_out = zeros(size(A_in, 1)+1, size(A_in, 2), size(A_in, 3));
    for ii = 1:size(A_in, 3)
        A_out(:,:,ii) = cat_zero_row(A_in(:,:,ii));
    end
else
    A_out = zeros(size(A_in, 1)+num_added_rows, size(A_in, 2), size(A_in, 3));
    for ii = 1:size(A_in, 3)
        A_out(:,:,ii) = cat_zero_row(A_in(:,:,ii), num_added_rows);
    end
end

end %function

